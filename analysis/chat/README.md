# Analysis: Chat

These scripts were used to build machine learning models to power the Playbook Partnership chat feature.

1. `gpt-data-generation.py`: Use random walks through the playbook partnership and the stories generated to prompt OpenAI to build instructions.
2. `gpt-embedding-clf.py`: Given instruction embeddings (from openai), predict the most likely workflow steps in the original workflow.
3. `transformer.py`: Given instructions in english, predict the most likely workflow in the original workflow. i.e. "translate" from english to the workflow DSL.
