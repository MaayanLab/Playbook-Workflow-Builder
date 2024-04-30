import tempfile
from components.core.file import File
from components.data.gene_count_matrix import anndata_from_file

def csv_read_stream(file: File, chunk_size=8192):
  ''' Stream the AnnData file as a csv
  '''
  d = anndata_from_file(file).to_df()
  with tempfile.TemporaryFile('wb+') as fh:
    d.to_csv(fh)
    fh.flush()
    fh.seek(0)
    while True:
      buf = fh.read(chunk_size)
      if not buf: break
      yield buf
