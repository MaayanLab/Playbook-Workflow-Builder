''' Provides helpers for working with files

Newly created files should be upserted with the upsert_file function,
this will ensure the file is accessible to any worker that needs to
access it.
'''
import os
import tempfile
import requests
import contextlib

class TemporaryFile:
  def __init__(self, suffix='') -> None:
    self.file = tempfile.mktemp(suffix=suffix)
    self.url = None

  def __repr__(self) -> str:
    return f"TemporaryFile({repr(self)})"
  
  def unlink(self):
    os.unlink(self.file)

@contextlib.contextmanager
def upsert_file(suffix=''):
  ''' Usage:
  from components.core.file import upsert_file
  with upsert_file('.csv') as f:
    # or, for example: df.to_csv(f.file)
  # f.url contains the url to the file uploaded
  print(f.url)
  '''
  # NOTE: this context manager will remove the temporary file when done
  tmp = TemporaryFile(suffix)
  try:
    # we give you a temporary file to write to
    yield tmp
    # we upload the file to the main server to get a persistent url
    req = requests.post(
      f"{os.environ['PUBLIC_URL']}/api/components/core/file/upload",
      headers={'Authorization': f"Token {os.environ['NEXTAUTH_SECRET']}"},
      files=[('file', open(tmp.file, 'rb'))],
    )
    assert req.status_code >= 200 and req.status_code < 300, f"Error ({req.status_code}): {req.json()}"
    res = req.json()
    # we return register url to the file
    tmp.url = res['file'][0]
  finally:
    # we remove the temporary file
    tmp.unlink()
