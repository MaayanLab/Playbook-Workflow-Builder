import os
import re
import sys
import uuid
import time
import requests
import pandas as pd
from urllib.parse import quote
from components.core.file import upsert_file
from components.data.gene_count_matrix import gene_count_matrix
from collections import Counter

def process_alignment(upload_uid, filenames, paired=False, organism='human', polling_interval=10):
  alignment_uid = str(uuid.uuid4())
  # start jobs with uploaded files
  jobs = []
  if paired:
    pairs = {}
    for filename in filenames:
      m = re.match(r'^(.+)_(1|2)\.fastq\.gz$', filename)
      if not m: continue
      prefix = m.group(1)
      suffix = m.group(2)
      if prefix not in pairs: pairs[prefix] = { 'file1': None, 'file2': None }
      if suffix == '1': pairs[prefix]['file1'] = f"{upload_uid}-{filename}"
      elif suffix == '2': pairs[prefix]['file2'] = f"{upload_uid}-{filename}"
    #
    for suffix, files in pairs.items():
      print(f"Creating job for {prefix}...", file=sys.stderr)
      cloudalignmentCreatejobReq = requests.get('https://maayanlab.cloud/cloudalignment/createjob', params=dict(
        username=os.environ['ELYSIUM_USERNAME'],
        password=os.environ['ELYSIUM_PASSWORD'],
        organism=organism,
        **files,
        outname=f"{alignment_uid}-{suffix}-{organism.replace('human', 'hs').replace('mouse', 'mm')}",
      ))
      cloudalignmentCreatejobRes = cloudalignmentCreatejobReq.json()
      jobs.append(cloudalignmentCreatejobRes)
  else:
    for filename in filenames:
      print(f"Creating job for {filename}...", file=sys.stderr)
      cloudalignmentCreatejobReq = requests.get('https://maayanlab.cloud/cloudalignment/createjob', params=dict(
        username=os.environ['ELYSIUM_USERNAME'],
        password=os.environ['ELYSIUM_PASSWORD'],
        organism=organism,
        file1=f"{upload_uid}-{filename}",
        outname=f"{alignment_uid}-{filename.rsplit('.', 2)[0]}-{organism.replace('human', 'hs').replace('mouse', 'mm')}",
      ))
      cloudalignmentCreatejobRes = cloudalignmentCreatejobReq.json()
      jobs.append(cloudalignmentCreatejobRes)
  # wait for alignment to complete
  status = -2
  while True:
    time.sleep(polling_interval)
    cloudalignmentProgressReq = requests.get('https://maayanlab.cloud/cloudalignment/progress', params=dict(
      username=os.environ['ELYSIUM_USERNAME'],
      password=os.environ['ELYSIUM_PASSWORD'],
      prefix=alignment_uid,
    ))
    cloudalignmentProgressRes = cloudalignmentProgressReq.json()
    status_count = Counter([alignment['status'] for alignment in cloudalignmentProgressRes.values()])
    if 'failed' in status_count:
      print(f"\nJob failed", file=sys.stderr, flush=True)
      break
    elif 'submitted' in status_count:
      if status < 0:
        status = 0
        print(f"\nJob started", end='', file=sys.stderr, flush=True)
      else:
        while status < status_count.get('completed', 0):
          print('!', end='', file=sys.stderr)
          status += 1
        print('.', end='', file=sys.stderr, flush=True)
      continue
    elif 'waiting' in status_count:
      if status != -1:
        status = -1
        print(f"\nWaiting for job to start", end='', file=sys.stderr, flush=True)
      else:
        print('.', end='', file=sys.stderr, flush=True)
      continue
    else:
      print('\nDone.', file=sys.stderr, flush=True)
      break
  print('Assembling gene count matrix...', file=sys.stderr)
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
    if filename.startswith(alignment_uid) and filename.endswith('_gene.tsv')
  ]
  # download and merge aligned files
  df = pd.concat([
    pd.read_csv(
      f"https://s3.amazonaws.com/biodos/bda6d89b6a9bfc05956b33dbf3eeff6b/{quote(filename)}",
      sep='\t',
      index_col=0,
      header=None,
      names=[filename.split('-', 5)[-1][:-len('-hs_gene.tsv')]],
    ).iloc[:, 0]
    for filename in alignedFiles
  ], axis=1)
  # add file to the platform and return as gene count matrix
  with upsert_file('.tsv', description='Aligned gene count matrix') as f:
    df.to_csv(f.file, sep='\t')
  return gene_count_matrix(f)
