''' Storage in a CAVATICA project with the Seven Bridges API client (sevenbridges-python).

This backs `CavaticaStorage` (components/core/file/server/storage.ts), used by CAVATICA workers
(`CAVATICA_STORAGE_PROJECT`): files are kept in the user's project, named by their sha256.

Environment:
- CAVATICA_API_TOKEN: the user's CAVATICA API token
- CAVATICA_STORAGE_PROJECT: the project (`{owner}/{project}`) to store files in
- CAVATICA_API_URL: the API (default https://cavatica-api.sbgenomics.com/v2)
'''
import os
import functools

@functools.cache
def _api():
  import sevenbridges as sb
  return sb.Api(
    url=os.environ.get('CAVATICA_API_URL') or 'https://cavatica-api.sbgenomics.com/v2',
    token=os.environ['CAVATICA_API_TOKEN'],
  )

def _find(key: str):
  for file in _api().files.query(project=os.environ['CAVATICA_STORAGE_PROJECT'], names=[key], limit=1):
    if file.type != 'folder':
      return file
  return None

def stat(key: str):
  file = _find(key)
  return dict(size=file.size) if file else None

def download_url(key: str):
  ''' A (presigned) url to GET the file directly, supporting Range requests '''
  file = _find(key)
  return file.download_info().url if file else None

def put(key: str, path: str):
  ''' Move a local file into the project (content addressed: existing files aren't re-uploaded) '''
  if _find(key) is None:
    _api().files.upload(path, project=os.environ['CAVATICA_STORAGE_PROJECT'], file_name=key, overwrite=False, wait=True)
  os.unlink(path)

def delete(key: str):
  file = _find(key)
  if file: file.delete()

def copy(src_key: str, dst_key: str):
  ''' Copy a file within the project (the platform copies it, nothing is transferred) '''
  if _find(dst_key) is not None: return
  file = _find(src_key)
  if file is None: raise FileNotFoundError(src_key)
  file.copy(os.environ['CAVATICA_STORAGE_PROJECT'], name=dst_key)

# Multipart uploads (parts are PUT to the platform's own presigned urls, e.g. directly by the browser)

def create_multipart(key: str, size: int, part_size: int):
  res = _api().post('/upload/multipart', data=dict(
    project=os.environ['CAVATICA_STORAGE_PROJECT'],
    name=key,
    size=size,
    part_size=part_size,
  )).json()
  return dict(upload_id=res['upload_id'], part_size=res.get('part_size') or part_size)

def part_url(upload_id: str, part_number: int):
  res = _api().get(f'/upload/multipart/{upload_id}/part/{part_number}').json()
  return dict(url=res['url'], headers=res.get('headers') or {})

def complete_multipart(upload_id: str, parts: list):
  for part in parts:
    _api().post(f'/upload/multipart/{upload_id}/part', data=dict(
      part_number=part['part_number'],
      response=dict(headers=dict(ETag=part['etag'])),
    ))
  _api().post(f'/upload/multipart/{upload_id}/complete')

def abort_multipart(upload_id: str):
  _api().delete(f'/upload/multipart/{upload_id}')
