from components.core.file import file_as_stream

def load_drug_matrix_transpose(file):
  with file_as_stream(file) as fr:
    return {
      term: dict(description=description, set=drugs)
      for line in fr
      for line_split in (line.decode().strip().split('\t'),)
      if len(line_split) >= 3
      for term, description, *drugs in (line_split,)
    }
