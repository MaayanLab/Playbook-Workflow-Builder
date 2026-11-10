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
