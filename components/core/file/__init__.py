''' Provides helpers for working with files

Newly created files should be upserted with the upsert_file function,
this will ensure the file is accessible to any worker that needs to
access it.

Files are referenced by url, most commonly `drs://` urls to this (or another) DRS server,
which we resolve to their access url (e.g. a presigned S3 url) and read with HTTP Range
requests: only the parts of a file which are actually read are ever transferred, files
are never copied in full.
'''
import io
import os
import typing
import tempfile
import requests
import contextlib
import urllib.parse

class File(typing.TypedDict):
  url: str
  filename: str
  description: typing.Optional[str]
  size: int
  sha256: typing.Optional[str]

class TemporaryFile(dict):
  def __init__(self, suffix='', description=None):
    super().__init__(
      url='',
      filename='',
      description=description,
      size=0,
    )
    self.file = tempfile.mktemp(suffix=suffix)
    self['filename'] = self.file

  def __repr__(self) -> str:
    return f"TemporaryFile({repr(self)})"

  def unlink(self):
    os.unlink(self.file)

@contextlib.contextmanager
def upsert_file(suffix='', description=None):
  ''' Usage:
  from components.core.file import upsert_file
  with upsert_file('.csv') as f:
    # or, for example: df.to_csv(f.file)
  # f.url contains the url to the file uploaded
  print(f.url)
  '''
  # NOTE: this context manager will remove the temporary file when done
  tmp = TemporaryFile(suffix, description=description)
  try:
    # we give you a temporary file to write to
    yield tmp
    # we upload the file to the main server to get a persistent url
    req = requests.post(
      f"{os.environ['PUBLIC_URL']}/api/file/upload",
      headers={'Authorization': f"Token {os.environ.get('WORKER_AUTH_TOKEN') or os.environ['BETTER_AUTH_SECRET']}"},
      files=[('file', open(tmp.file, 'rb'))],
    )
    assert req.status_code >= 200 and req.status_code < 300, f"Error ({req.status_code}): {req.json()}"
    res = req.json()
    # we return register url to the file
    tmp.update(res['file'][0])
  finally:
    # we remove the temporary file
    tmp.unlink()

class _Resolved(typing.TypedDict, total=False):
  # a local file
  path: str
  # or a url (to read with Range requests)
  url: str
  headers: typing.Dict[str, str]
  size: typing.Optional[int]
  # whether the url is already encoded (presigned urls must not be re-encoded)
  encoded: bool

def _drs_object_url(netloc: str, object_id: str):
  # like ufs: a host with an explicit port (i.e. a development server) is served over http
  scheme = 'http' if ':' in netloc else 'https'
  return f"{scheme}://{netloc}/ga4gh/drs/v1/objects/{object_id}"

def _drs_headers(headers) -> typing.Dict[str, str]:
  ''' DRS access headers are a list of "Name: value" strings '''
  if not headers: return {}
  if isinstance(headers, dict): return headers
  return dict(header.split(': ', 1) for header in headers)

def _local_storage_path(key: str):
  ''' storage:// urls are only resolvable by python with the (default) local filesystem storage '''
  if os.environ.get('S3_BUCKET') or os.environ.get('UFS_STORAGE') or os.environ.get('CAVATICA_STORAGE_PROJECT'):
    raise ValueError("storage:// urls are internal to the app, reference files by their drs:// url")
  root = os.environ.get('STORAGE_PATH') or os.path.join(os.environ.get('PYTHON_ROOT') or os.getcwd(), 'data', 'ufs')
  return os.path.join(root, key)

def _resolve(url: str) -> _Resolved:
  parsed = urllib.parse.urlparse(url)
  if parsed.scheme == 'file':
    return dict(path=urllib.parse.unquote(parsed.path))
  elif parsed.scheme == 'storage':
    return dict(path=_local_storage_path(parsed.netloc))
  elif parsed.scheme in ('http', 'https'):
    return dict(url=url, headers={}, size=None, encoded=False)
  elif parsed.scheme == 'drs':
    object_url = _drs_object_url(parsed.netloc, parsed.path.lstrip('/'))
    req = requests.get(object_url)
    assert req.ok, f"DRS object {url} could not be resolved ({req.status_code})"
    obj = req.json()
    size = obj.get('size')
    access_methods = [m for m in obj.get('access_methods', []) if m.get('type') in ('https', 'http')]
    # prefer an access id: the access url it resolves to is typically direct (i.e. presigned)
    for method in access_methods:
      if method.get('access_id'):
        req = requests.get(f"{object_url}/access/{method['access_id']}")
        if req.ok:
          access = req.json()
          return dict(url=access['url'], headers=_drs_headers(access.get('headers')), size=size, encoded=True)
    for method in access_methods:
      access_url = method.get('access_url')
      if access_url:
        if isinstance(access_url, str):
          return dict(url=access_url, headers={}, size=size, encoded=True)
        return dict(url=access_url['url'], headers=_drs_headers(access_url.get('headers')), size=size, encoded=True)
    raise ValueError(f"DRS object {url} has no supported access method")
  else:
    raise NotImplementedError(f"Unsupported file url {url}")

def _open_url(resolved: _Resolved, block_size: int, cache_type: str):
  import fsspec
  fs = fsspec.filesystem('http', encoded=resolved.get('encoded', False), headers=resolved.get('headers') or {})
  # with a known size, the file is seekable without a HEAD request (presigned GET urls can't be HEAD)
  return fs.open(resolved['url'], 'rb', block_size=block_size, cache_type=cache_type, size=resolved.get('size'))

@contextlib.contextmanager
def file_as_stream(file: File, mode: str = 'rb', encoding: str = 'utf-8'):
  ''' Read a file sequentially (mode 'rb', or 'r' for text) '''
  resolved = _resolve(file['url'])
  if 'path' in resolved:
    with open(resolved['path'], mode, **({} if 'b' in mode else dict(encoding=encoding))) as fr:
      yield fr
  else:
    with _open_url(resolved, block_size=8 << 20, cache_type='readahead') as fr:
      yield fr if 'b' in mode else io.TextIOWrapper(fr, encoding=encoding)

@contextlib.contextmanager
def file_as_seekable(file: File, block_size: int = 1 << 20):
  ''' Random access to a (binary) file, e.g. for HDF5 (h5py/anndata) or excel files.
  Only the blocks which are read are transferred (and cached for re-reads). '''
  resolved = _resolve(file['url'])
  if 'path' in resolved:
    with open(resolved['path'], 'rb') as fr:
      yield fr
  else:
    with _open_url(resolved, block_size=block_size, cache_type='blockcache') as fr:
      yield fr

def file_read_stream(file: File, chunk_size=8192):
  with file_as_stream(file, 'rb') as fr:
    while True:
      buf = fr.read(chunk_size)
      if not buf: break
      yield buf

def one_or_none(it):
  try:
    one = next(it)
  except StopIteration:
    return None
  try:
    next(it)
  except StopIteration:
    return one
  else:
    return None

def file_stats_from_url(url: str, chunk_size=8192):
  import hashlib
  h = hashlib.sha256()
  s = 0
  with file_as_stream(dict(url=url), 'rb') as fr:
    while True:
      buf = fr.read(chunk_size)
      if not buf:
        break
      h.update(buf)
      s += len(buf)
  return dict(
    size=s,
    sha256=h.hexdigest(),
  )

def file_from_url(url: str):
  parsed = urllib.parse.urlparse(url)
  filename = urllib.parse.unquote(parsed.path).rpartition('/')[-1]
  if parsed.scheme == 'drs':
    req = requests.get(_drs_object_url(parsed.netloc, parsed.path.lstrip('/')))
    assert req.ok
    file_info = req.json()
    filename = file_info.get('name') or filename
    sha256 = one_or_none(checksum['checksum'] for checksum in file_info.get('checksums', []) if checksum.get('type') == 'sha-256')
    if sha256:
      return dict(url=url, filename=filename, size=file_info['size'], sha256=sha256)
  return dict(url=url, filename=filename, **file_stats_from_url(url))
