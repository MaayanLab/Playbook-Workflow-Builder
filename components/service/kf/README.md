# Knowledge Finder Gene Expression in Tumors Service

This directory contains the Knowledge Finder gene expression in tumors service, which is designed to retrieve and process gene expression data for a given Ensembl ID. The service uses the OpenPedCan API to obtain the gene expression data and returns the data as a dictionary in the format:

```python
[
  {
    "TPM_mean": <float>,
    "TPM_sd": <float>,
    "TPM_median": <float>,
    "Disease": <str>,
    "Gene_symbol": <str>,
    "Gene_Ensembl_ID": <str>,
    "Dataset": <str>
  },
  ...
]
```

# Dependencies
This service has the following dependencies:

* Python 3.6+
* pandas
* requests

You can install these dependencies by running:
```shell
pip install -r requirements.txt
```

# Usage
To use this service, simply call the main function in kf/__init__.py, passing in an Ensembl ID as a string:

This will return a list of dictionaries, with each dictionary representing gene expression data for a specific cancer type.