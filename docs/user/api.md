# API User Guide

The Playbook Partnership API allows you to construct workflows via a REST API. Workflows are specified using a JSON serialized DSL, they are then registered into the database and either the output can be resolved or the resulting playbook can be visited.

```bash
# Construct a workflow Gene => GTEx Tissue => Barplot with the input: ACE2
curl -H'Content-Type: application/json' -XPOST https://playbook-workflow-builder.cloud/api/db/fpl -d @- << EOF
[
  {"id":"1","type":"Input[Gene]","data":{"type":"Term[Gene]","value":"ACE2"}},
  {"id":"2","type":"GTExTissueExpressionFromGene","inputs":{"gene":{"id":"1"}}},
  {"id":"3","type":"BarplotFrom[Scored[Tissue]]","inputs":{"terms":{"id":"2"}}}
]
EOF

# Construct a workflow FileInput => AnnDataFromFile => Limma-Voom with the example.h5ad file
curl -H'Content-Type: application/json' -XPOST https://playbook-workflow-builder.cloud/api/db/fpl -d @- << EOF
[
  {"id": "0", "type": "FileInput", "data": { "type": "FileURL", "value": { "url": "https://playbook-workflow-builder.cloud/api/v1/components/core/file/example.h5ad", "filename": "example.h5ad" } } },
  {"id": "1", "type": "AnnDataFromFile", "inputs": { "file": { "id": "0" } } },
  {"id": "2", "type": "Limma-Voom", "inputs": { "anndata": { "id": "1" } } }
]
EOF
```

Both of these examples will return an `id` for the workflow. This id can be used as follows:
- `https://playbook-workflow-builder.cloud/graph/{id}` => The graph mode of the workflow
- `https://playbook-workflow-builder.cloud/report/{id}` => The report mode of the workflow
- `https://playbook-workflow-builder.cloud/api/db/fpl/{id}` => See the resulting workflow as registered in the db
- `https://playbook-workflow-builder.cloud/api/db/fpl/{id}/output` => See the resulting workflow along with its output (for all nodes)

Each process (step of the workflow) can be accessed independently, those ids are at `[{"id": "fpl_id", "process": { "id": "{process_id}", ... }, ... }]` and can be accessed like:
- `https://playbook-workflow-builder.cloud/api/db/process/{id}` => See the resulting process as registered in the db
- `https://playbook-workflow-builder.cloud/api/db/process/{id}/output` => See the resulting process along with its output (for all nodes)
