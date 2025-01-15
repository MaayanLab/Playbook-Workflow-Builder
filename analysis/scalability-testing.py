#%%
import time
import random
import requests
import pandas as pd
import seaborn as sns

ENDPOINT = 'http://192.168.1.2'

#%%
def test_depth(depth: int):
  start = time.time()
  res = requests.post(ENDPOINT + '/api/db/fpl', json={
    "data":{
      "input":{"type":"TestView","value":f"{random.random()}"}
    },
    "workflow":[
      {"id":"0","type":"InitialTestProcess","data":{"id":"input"}},
    ] + [
      {"id":f"{i+1}","type":"TestProcess", "inputs":{"input":{"id":f"{i}"}}}
      for i in range(depth)
    ]
  })
  playbook_id = res.json()
  submit = time.time()
  res = requests.get(ENDPOINT + '/api/db/fpl/' + playbook_id + '/output/export')
  while res.status_code == 504:
    res = requests.get(ENDPOINT + '/api/db/fpl/' + playbook_id + '/output/export')
  res.raise_for_status()
  playbook = res.json()
  resolve = time.time()
  return {
    'submission': submit - start,
    'resolution': resolve - submit,
  }

def test_parallel(users: int, depth: int):
  from concurrent.futures import ThreadPoolExecutor
  with ThreadPoolExecutor(max_workers=users) as pool:
    futures = [pool.submit(test_depth, depth) for _ in range(users)]
    results = [f.result() for f in futures]
  return results

#%%
results = []

for u in range(1, 50+1):
  for d in range(1, 10+1):
    results += [dict(result, n_users=u, depth=d) for result in test_parallel(u, d)]

df = pd.DataFrame(results)
df.to_csv('results-2x5.tsv', sep='\t')

#%%
(
  sns.FacetGrid(df[df['depth'].isin([5,10])&df['parallel users'].isin([1,10,20,30,40,50])], col='depth', hue='parallel users')
    .map(sns.stripplot, 'processes', 'resolution')
    .add_legend()
)

#%%
(
  sns.FacetGrid(df[df['depth'].isin([5,10])&df['parallel users'].isin([1,10,20,30,40,50])], col='depth', hue='parallel users')
    .map(sns.stripplot, 'processes', 'submission')
    .add_legend()
)

#%%
# sns.heatmap(df.pivot_table(columns='depth', index='n_users', values='resolution', aggfunc='mean'))
