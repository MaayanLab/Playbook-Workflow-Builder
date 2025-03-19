from components.core.file import file_as_stream

def read_dmt(fr):
  return {
    term: dict(description=description, set=drugs)
    for line in fr
    for line_split in (line.decode().strip().split('\t'),)
    if len(line_split) >= 3
    for term, description, *drugs in (line_split,)
  }

def load_drug_matrix_transpose(file):
  if file['filename'].endswith('.dmt') or file['filename'].endswith('.gmt') or file['filename'].endswith('.txt'):
    with file_as_stream(file, 'r') as fr:
      return read_dmt(fr)
  elif file['filename'].endswith('.dmt.gz') or file['filename'].endswith('.gmt.gz') or file['filename'].endswith('.txt.gz'):
    import gzip
    with file_as_stream(file, 'rb') as fr_compressed:
      with gzip.open(fr_compressed, 'rt') as fr:
        return read_dmt(fr)
  else:
    raise NotImplementedError
