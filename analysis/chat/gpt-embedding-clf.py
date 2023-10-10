#%%
import json
import pandas as pd
import tensorflow as tf
from sklearn.feature_extraction.text import CountVectorizer

#%%
store = pd.HDFStore('store.h5', 'r')
gpt_instructions = store['gpt_instructions'].values
gpt_instruction_story = store['gpt_instruction_story'].values
gpt_instructions_embedded = store['gpt_instructions_embedded'].values
stories = store['stories'].values
workflows = store['workflows'].values
store.close()
gpt_instructions_stories = [stories[i] for i in gpt_instruction_story]
gpt_instructions_workflows = [workflows[i] for i in gpt_instruction_story]

#%%
vectorizer = CountVectorizer(ngram_range=(1, 2), tokenizer=lambda v: v, preprocessor=lambda v: v)
vectorizer.fit(workflows)
_, vectorizer_vocabulary = zip(*sorted((i, v) for v, i in vectorizer.vocabulary_.items()))
workflows_vectorized = (vectorizer.fit_transform(gpt_instructions_workflows)>0).astype(int)
workflows_vectorized

#%%
model = tf.keras.Sequential([
  tf.keras.layers.Input((gpt_instructions_embedded.shape[1],)),
  tf.keras.layers.Dense(100, activation='relu'),
  tf.keras.layers.Dense(100, activation='relu'),
  tf.keras.layers.Dense(workflows_vectorized.shape[1], activation='sigmoid'),
])
model.compile(loss='binary_crossentropy', optimizer='adam')
model.fit(gpt_instructions_embedded, workflows_vectorized.toarray(), epochs=10)

#%%
model.save('data/gpt-embedding-clf')

#%%
with open('data/gpt-embedding-clf.json', 'w') as fw:
  json.dump(dict(
    url='http://tensorflow:8501/v1/models/gpt-embedding-clf',
    vocab=list(vectorizer_vocabulary)
  ), fw)
