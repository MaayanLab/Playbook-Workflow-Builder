from components.core.file import file_as_stream

def load_gene_matrix_transpose(file):
  with file_as_stream(file['url']) as fr:
    return {
      term: dict(description=description, set=genes)
      for line in fr
      for line_split in (line.decode().strip().split('\t'),)
      if len(line_split) >= 3
      for term, description, *genes in (line_split,)
    }
