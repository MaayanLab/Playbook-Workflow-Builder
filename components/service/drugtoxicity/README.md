# Cytotoxicity, BBB Permeability, and Ranked List Drugs

This directory contains Cytotoxicity, BBB Permeability, and Ranked List Drugs Nodes.

## Cytotoxicity CHEMBL API

Queries the ChEMBL database for cytotoxicity assay data for a given list of drugs. Then filters the data with a PCHEMBL_VALUE >=1 (Note that pCHEMBL = -log10(IC50, XC50, AC50, Ki, Kd, potency). Greater than or equal to 5 (10um) is a typical minimum rule of thumb for binding activity between a compound and a protein target.)

```python
[
    {
        "activity_id": <str>,
        "assay_description": <str>,
        "assay_type": <str>,
        "drug_name": <str>,
        "standard_type": <str>,
        "standard_units": <str>,
        "standard_value": <int>,
        "target_pref_name":<str>,
        "pchembl_value": <float>,
    }
]
```

## BBB Permeability

Queries the B3DB (Blood-Brain Barrier Database) for BBB permeability data and merges it with drug scores.

```python
[
    {
        "drug_name": <str>,
        "logBB": <float>,
        "bbb_permeable": <str>,
    }
]
```

## Ranked List of Drugs

Produces a ranked list of drug candidates based on the combination of cytotoxicity, BBB and confidence z-scores.

```python
[
    {
        "drug_name": <str>,
        "confidence_zscore": <str>,
        "logBB": <float>,
        "bbb_permeable": <str>,
        "cytotoxicity_activity_id": <int>,
        "cytotoxicity_assay_description": <str>,
        "cytotoxicity_assay_type": <str>,
        "cytotoxicity_standard_type": <str>,
        "cytotoxicity_standard_units": <str>,
        "cytotoxicity_standard_value": <int>,
        "cytotoxicity_target_pref_name": <str>,
        "cytotoxicity_pchembl_value": <float>,
    }
]
```
