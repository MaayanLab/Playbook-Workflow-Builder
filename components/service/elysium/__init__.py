import os
import uuid
import time
import requests
import pandas as pd
from core.file import upsert_file
from data.gene_count_matrix import gene_count_matrix
from collections import Counter

def process_single_end(upload_uid, filenames, organism='human', polling_interval=10):
  alignment_uid = uuid.uuid4()
  # start jobs with uploaded files
  jobs = []
  for filename in filenames:
    cloudalignmentCreatejobReq = requests.get('https://maayanlab.cloud/cloudalignment/createjob', params=dict(
      username=os.environ['ELYSIUM_USERNAME'],
      password=os.environ['ELYSIUM_PASSWORD'],
      organism=organism, # or mouse
      file1=filename,
      outname=f"{alignment_uid}-{filename.rsplit('.', 2)[0]}-{organism}",
    ))
    cloudalignmentCreatejobRes = cloudalignmentCreatejobReq.json()
    jobs.append(cloudalignmentCreatejobRes)
  # wait for alignment to complete
  while True:
    time.sleep(polling_interval)
    cloudalignmentProgressReq = requests.get('https://maayanlab.cloud/cloudalignment/progress', params=dict(
      username=os.environ['ELYSIUM_USERNAME'],
      password=os.environ['ELYSIUM_PASSWORD'],
      prefix=alignment_uid,
    ))
    cloudalignmentProgressRes = cloudalignmentProgressReq.json()
    # failed, submitted, waiting
    # TODO: log status
    print(cloudalignmentProgressRes)
    status_count = Counter([alignment['progress'] for alignment in cloudalignmentProgressRes])
    if 'failed' in status_count: break
    elif 'submitted' in status_count: continue
    elif 'waiting' in status_count: continue
    else: break
  # identify aligned files
  charonFilesReq = requests.get('https://maayanlab.cloud/charon/files', params=dict(
    username=os.environ['ELYSIUM_USERNAME'],
    password=os.environ['ELYSIUM_PASSWORD'],
    prefix=alignment_uid,
  ))
  charonFilesRes = charonFilesReq.json()
  alignedFiles = [
    filename 
    for filename in charonFilesRes['filenames']
    if filename.startswith(upload_uid) and filename.endswith('_gene.tsv')
  ]
  # download and merge aligned files
  df = pd.concat({
    filename: pd.DataFrame(f"https://s3.amazonaws.com/biodos/c095573dc866f2db2cd39862ad89f074/${encodeURIComponent(filename)}", sep='\t')
    for f in alignedFiles
  }, axis=1)
  # add file to the platform and return as gene count matrix
  with upsert_file('.tsv', description='Aligned gene count matrix') as f:
    df.to_csv(f.file, sep='\t')
  return gene_count_matrix(f)
