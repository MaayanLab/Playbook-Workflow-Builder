#%%
import os
import re
import time
import openai
import pandas as pd
from tqdm import tqdm
from textwrap import dedent
openai.api_key = os.getenv("OPENAI_API_KEY")

def chunked(it, chunk_size=100):
  chunk = []
  for el in it:
    chunk.append(el)
    if len(chunk) >= chunk_size:
      yield chunk
      chunk = []
  if chunk:
    yield chunk

def chat_completion_with_retries(*args, **kwargs):
  for i in range(10):
    try:
      return openai.ChatCompletion.create(*args, **kwargs)
    except KeyboardInterrupt: raise
    except Exception as e: exc = e
    time.sleep(5)
    print(f"retry {i+1}")
  raise exc

#%%
# random-paths.tsv was generated from the playbook partnership codebase, it consists of
#  stories associated with workflows, which are simply paths of process nodes
stories = []
workflows = []
with open('data/random-paths.tsv', 'r') as fr:
  for line in fr:
    line_split = line.strip().split('\t')
    if len(line_split) < 3: continue
    story, _, *workflow = line_split
    stories.append(story)
    workflows.append(workflow)

#%%
# We want to build instructions a user would have provided which
#  would have resulted in the random path. We do this with chat completions.
gpt_instruction_story = []
gpt_instructions = []
for story_index, story in enumerate(tqdm(stories)):
  chat_completion = chat_completion_with_retries(
    model="gpt-3.5-turbo",
    messages=[
      {
        "role": "user",
        "content": dedent(f"""
          Consider a biomedical workflow with the following description:
          {story}

          Please provide several possible user instructions that would have prompted us to create the above workflow.
        """)
      },
    ]
  )
  instructions = [
    m.group(3).strip('"')
    for line in map(str.strip, chat_completion['choices'][0]['message']['content'].splitlines())
    if line
    for m in (re.match(r'^((\d+)\.\s*)?(.+)$', line),)
  ]
  gpt_instructions.extend(instructions)
  gpt_instruction_story.extend([story_index]*len(instructions))

#%%
# Here we fetch embeddings of the user instructions GPT generated
gpt_instructions_embedded = [
  result['embedding']
  for chunk in tqdm(chunked(gpt_instructions, 100))
  for result in openai.Embedding.create(
    model='text-embedding-ada-002',
    input=chunk,
  )['data']
]

#%%
# we save these things for use machine learning models to power the chat feature
store = pd.HDFStore('data/store.h5', 'w')
store['stories'] = pd.Series(stories)
store['workflows'] = pd.Series(workflows)
store['gpt_instruction_story'] = pd.Series(gpt_instruction_story)
store['gpt_instructions'] = pd.Series(gpt_instructions)
store['gpt_instructions_embedded'] = pd.DataFrame(gpt_instructions_embedded)
store.close()
