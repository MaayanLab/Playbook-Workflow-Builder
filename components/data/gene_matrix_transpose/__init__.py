from components.core.file import file_as_stream

def read_gmt(fr):
  return {
    term: dict(description=description, set=genes)
    for line in fr
    for line_split in (line.strip().split('\t'),)
    if len(line_split) >= 3
    for term, description, *genes in (line_split,)
  }

def load_gene_matrix_transpose(file):
  if file['filename'].endswith('.gmt') or file['filename'].endswith('.txt'):
    with file_as_stream(file, 'r') as fr:
      return read_gmt(fr)
  elif file['filename'].endswith('.gmt.gz') or file['filename'].endswith('.txt.gz'):
    import gzip
    with file_as_stream(file, 'rb') as fr_compressed:
      with gzip.open(fr_compressed, 'rt') as fr:
        return read_gmt(fr)
  else:
    raise NotImplementedError
