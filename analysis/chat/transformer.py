# based on https://keras.io/examples/nlp/neural_machine_translation_with_transformer/

#%%
import json
import string
import random
import pandas as pd
import tensorflow as tf

#%%
store = pd.HDFStore('store.h5', 'r')
gpt_instructions = store['gpt_instructions'].values
gpt_instruction_story = store['gpt_instruction_story'].values
stories = store['stories'].values
workflows = store['workflows'].values
store.close()
gpt_instructions_stories = [stories[i] for i in gpt_instruction_story]
gpt_instructions_workflows = [workflows[i] for i in gpt_instruction_story]

#%%
workflow_steps = pd.Series(workflows).explode().unique()

#%%
instruct_pairs = []
for i, j in enumerate(gpt_instruction_story):
  instruct_pairs.append((
    '[start] ' + gpt_instructions[i] + ' [end]',
    '\t'.join(['[start]', *workflows[j], '[end]'])
  ))

for _ in range(5):
    print(random.choice(instruct_pairs))

random.shuffle(instruct_pairs)
num_val_samples = int(0.10 * len(instruct_pairs))
num_train_samples = len(instruct_pairs) - 2 * num_val_samples
train_pairs = instruct_pairs[:num_train_samples]
val_pairs = instruct_pairs[num_train_samples : num_train_samples + num_val_samples]
test_pairs = instruct_pairs[num_train_samples + num_val_samples :]

print(f"{len(instruct_pairs)} total pairs")
print(f"{len(train_pairs)} training pairs")
print(f"{len(val_pairs)} validation pairs")
print(f"{len(test_pairs)} test pairs")

#%%
batch_size = 32
text_vocab_size = 10000
text_sequence_length = 64
workflow_vocab_size = workflow_steps.size
workflow_sequence_length = 16

strip_chars = string.punctuation
strip_chars = strip_chars.replace("[", "")
strip_chars = strip_chars.replace("]", "")
def custom_standardization(input_string):
    import re
    lowercase = tf.strings.lower(input_string)
    return tf.strings.regex_replace(lowercase, "[%s]" % re.escape(strip_chars), "")

eng_vectorization = tf.keras.layers.TextVectorization(
    max_tokens=text_vocab_size,
    output_mode="int",
    output_sequence_length=text_sequence_length + 1,
    standardize=custom_standardization,
)
workflow_vectorization = tf.keras.layers.TextVectorization(
    max_tokens=workflow_vocab_size,
    output_mode="int",
    output_sequence_length=workflow_sequence_length + 1,
    standardize=None,
    split=lambda tokens: tf.strings.split(tokens, sep='\t'),
)
train_text, train_workflows = map(list, zip(*train_pairs))
eng_vectorization.adapt(train_text)
workflow_vectorization.adapt(train_workflows)

def format_dataset(text, workflow):
    text = eng_vectorization(text)
    workflow = workflow_vectorization(workflow)
    return (
        {
            "encoder_inputs": text,
            "decoder_inputs": workflow[:, :-1],
        },
        workflow[:, 1:],
    )
def make_dataset(pairs):
    texts, workflows = zip(*pairs)
    texts = list(texts)
    workflows = list(workflows)
    dataset = tf.data.Dataset.from_tensor_slices((texts, workflows))
    dataset = dataset.batch(batch_size)
    dataset = dataset.map(format_dataset)
    return dataset.shuffle(2048).prefetch(16).cache()

train_ds = make_dataset(train_pairs)
val_ds = make_dataset(val_pairs)
#%%
def TransformerEncoder(embed_dim, dense_dim, num_heads, inputs, mask=None):
    attention = tf.keras.layers.MultiHeadAttention(
        num_heads=num_heads, key_dim=embed_dim
    )
    dense_proj = tf.keras.Sequential(
        [
            tf.keras.layers.Dense(dense_dim, activation="relu"),
            tf.keras.layers.Dense(embed_dim),
        ]
    )
    layernorm_1 = tf.keras.layers.LayerNormalization()
    layernorm_2 = tf.keras.layers.LayerNormalization()
    supports_masking = True

    # def call(inputs, mask=None):
    attention_output = attention(query=inputs, value=inputs, key=inputs)
    proj_input = layernorm_1(inputs + attention_output)
    proj_output = dense_proj(proj_input)
    return layernorm_2(proj_input + proj_output)
    # return call

def PositionalEmbedding(sequence_length, vocab_size, embed_dim, inputs):
    token_embeddings = tf.keras.layers.Embedding(
        input_dim=vocab_size, output_dim=embed_dim
    )
    position_embeddings = tf.keras.layers.Embedding(
        input_dim=sequence_length, output_dim=embed_dim
    )

    # def call(inputs):
    length = tf.shape(inputs)[-1]
    positions = tf.range(start=0, limit=length, delta=1)
    embedded_tokens = token_embeddings(inputs)
    embedded_positions = position_embeddings(positions)
    return embedded_tokens + embedded_positions

    # return call

def TransformerDecoder(embed_dim, latent_dim, num_heads, inputs, encoder_outputs, mask=None):
    attention_1 = tf.keras.layers.MultiHeadAttention(
        num_heads=num_heads, key_dim=embed_dim
    )
    attention_2 = tf.keras.layers.MultiHeadAttention(
        num_heads=num_heads, key_dim=embed_dim
    )
    dense_proj = tf.keras.Sequential(
        [
            tf.keras.layers.Dense(latent_dim, activation="relu"),
            tf.keras.layers.Dense(embed_dim),
        ]
    )
    layernorm_1 = tf.keras.layers.LayerNormalization()
    layernorm_2 = tf.keras.layers.LayerNormalization()
    layernorm_3 = tf.keras.layers.LayerNormalization()
    add = tf.keras.layers.Add()  # instead of `+` to preserve mask
    supports_masking = True

    # def call(inputs, encoder_outputs, mask=None):
    attention_output_1 = attention_1(
        query=inputs, value=inputs, key=inputs, use_causal_mask=True
    )
    out_1 = layernorm_1(add([inputs, attention_output_1]))

    attention_output_2 = attention_2(
        query=out_1,
        value=encoder_outputs,
        key=encoder_outputs,
    )
    out_2 = layernorm_2(add([out_1, attention_output_2]))

    proj_output = dense_proj(out_2)
    return layernorm_3(add([out_2, proj_output]))

    # return call

#%%
embed_dim = 64
latent_dim = 1024
num_heads = 16

encoder_inputs = tf.keras.Input(shape=(None,), dtype="int64", name="encoder_inputs")
x = PositionalEmbedding(text_sequence_length + 1, text_vocab_size, embed_dim, encoder_inputs)
encoder_outputs = TransformerEncoder(embed_dim, latent_dim, num_heads, x)
encoder = tf.keras.Model(encoder_inputs, encoder_outputs)

decoder_inputs = tf.keras.Input(shape=(None,), dtype="int64", name="decoder_inputs")
encoded_seq_inputs = tf.keras.Input(shape=(None, embed_dim), name="decoder_state_inputs")
x = PositionalEmbedding(workflow_sequence_length + 1, workflow_vocab_size, embed_dim, decoder_inputs)
x = TransformerDecoder(embed_dim, latent_dim, num_heads, x, encoded_seq_inputs)
x = tf.keras.layers.Dropout(0.5)(x)
decoder_outputs = tf.keras.layers.Dense(workflow_vocab_size, activation="softmax")(x)
decoder = tf.keras.Model([decoder_inputs, encoded_seq_inputs], decoder_outputs)

decoder_outputs = decoder([decoder_inputs, encoder_outputs])
transformer = tf.keras.Model(
    [encoder_inputs, decoder_inputs], decoder_outputs, name="transformer"
)

#%%
epochs = 100

transformer.summary()
transformer.compile(
    "adam", loss="sparse_categorical_crossentropy", metrics=["accuracy"]
)
transformer.fit(train_ds, epochs=epochs, validation_data=val_ds)

#%%
text_input = tf.keras.Input(shape=tuple(), dtype="string", name="text")
workflow_input = tf.keras.Input(shape=tuple(), dtype="string", name="workflow")
workflow_output = transformer([eng_vectorization(text_input), workflow_vectorization(workflow_input)])

model = tf.keras.models.Model(
    inputs=[text_input, workflow_input],
    outputs=workflow_output,
)

tf.keras.models.save_model(
    model,
    'data/transformer',
    include_optimizer=False,
)

#%%
eng_vocab = eng_vectorization.get_vocabulary()
workflow_vocab = workflow_vectorization.get_vocabulary()
with open('data/transformer.json', 'w') as fw:
  json.dump(dict(
    url='http://tensorflow:8501/v1/models/transformer',
    vocab=dict(eng=eng_vocab, workflow=workflow_vocab),
  ), fw)
