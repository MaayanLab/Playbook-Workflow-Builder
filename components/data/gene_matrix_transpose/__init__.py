import fsspec
import contextlib

@contextlib.contextmanager
def fsspec_open_as_iterator(url, *args, **kwargs) -> str:
  with fsspec.open(url, *args, **kwargs) as fr:
    yield fr

def load_gene_matrix_transpose(file):
  with fsspec_open_as_iterator(file['url']) as fr:
    return {
      term: dict(description=description, set=genes)
      for line in fr
      for line_split in (line.decode().strip().split('\t'),)
      if len(line_split) >= 3
      for term, description, *genes in (line_split,)
    }
