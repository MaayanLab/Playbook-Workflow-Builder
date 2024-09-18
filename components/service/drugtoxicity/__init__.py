"""
Module: Drug Cytotoxicity, Blood-Brain Barrier (BBB) Permeability, and Candidate Ranking

This module provides a set of functions to query external databases like ChEMBL and B3DB to retrieve
drug-related data such as cytotoxicity assay information and blood-brain barrier (BBB) permeability.
Additionally, the module includes functionality to rank drug candidates based on a combination of
cytotoxicity data and user-provided confidence scores.

Key Functions:
--------------
1. query_drug_cytotoxicty_from_chembl(drugs: list[dict]) -> list[dict]:
   - Queries the ChEMBL database for cytotoxicity assay data for a given list of drugs.

2. batch_generator(large_drug_list: list[str], batch_size: int) -> list[str]:
   - A helper generator function to split a large list of drug names into batches.

3. query_chembl_api(drugs_to_query: list[str]) -> list[dict]:
   - Queries the ChEMBL API for cytotoxicity assay data for a specific batch of drugs.

4. produce_ranked_drug_candidates(drug_scores: list[dict], drug_cytotoxicity_chembl: list[dict]) -> list[dict]:
   - Produces a ranked list of drug candidates based on the combination of cytotoxicity and confidence z-scores.

5. query_drug_bbb_from_b3db(drug_scores: list[dict]) -> list[dict]:
   - Queries the B3DB (Blood-Brain Barrier Database) for BBB permeability data and merges it with drug scores.

Dependencies:
-------------
- pandas: For DataFrame operations and data manipulation.
- requests: For making API calls to external databases like ChEMBL.

Usage Example:
--------------
1. Provide a list of drugs with confidence scores and cytotoxicity data from ChEMBL.
2. Use `produce_ranked_drug_candidates()` to rank drugs based on their cytotoxicity and confidence z-scores.
3. Optionally, query B3DB using `query_drug_bbb_from_b3db()` to retrieve additional data about blood-brain barrier
   permeability.

"""

import pandas as pd
import requests


def query_drug_cytotoxicty_from_chembl(drugs: list[dict]):
    """
    Query the ChEMBL database for cytotoxicity assay data for a given list of drugs.

    Args:
        drugs (list[dict]): List of drugs with a "term" key, each containing a drug name.

    Returns:
        list[dict]: A list of dictionaries containing cytotoxicity data sorted by drug name.
    """
    # Extract the drug terms from the input list of dictionaries
    large_drug_list = [drug["term"] for drug in drugs]

    batch_size = 100
    results = None
    # Query ChEMBL API in batches for large lists of drugs
    for batch in batch_generator(large_drug_list, batch_size):
        results = query_chembl_api(batch)

    # Specify columns to keep from the API response
    cols = [
        "activity_id",
        "assay_description",
        "assay_type",
        "molecule_pref_name",
        "standard_type",
        "standard_units",
        "standard_value",
        "target_pref_name",
        "pchembl_value",
    ]

    drug_cytotoxicity_df = pd.DataFrame(columns=cols)

    if results:
        # Create a DataFrame from the API response
        drug_cytotoxicity_df = pd.DataFrame(results)[cols]
        # Convert specific columns to numeric types for further analysis
        drug_cytotoxicity_df[["pchembl_value", "standard_value"]] = (
            drug_cytotoxicity_df[["pchembl_value", "standard_value"]].apply(
                pd.to_numeric
            )
        )

    # Sort by molecule name and return the results as a list of dictionaries
    return drug_cytotoxicity_df.sort_values("molecule_pref_name").to_dict(
        orient="records"
    )


def batch_generator(large_drug_list: list[str], batch_size: int):
    """
    Generator function to yield batches from a large list.

    Args:
        large_drug_list (list[str]): List of drug names to query.
        batch_size (int): Number of items per batch.

    Yields:
        list[str]: A batch of drug names.
    """
    for i in range(0, len(large_drug_list), batch_size):
        yield large_drug_list[i : i + batch_size]


def query_chembl_api(drugs_to_query: list[str]):
    """
    Query the ChEMBL API for cytotoxicity assay data for a given list of drugs.

    Args:
        drugs_to_query (list[str]): List of drug names to query the ChEMBL API.

    Returns:
        list[dict]: A list of dictionaries containing cytotoxicity assay data from ChEMBL.
    """
    # Join drug names into a single comma-separated string
    drugs_str = ",".join(drug.upper() for drug in drugs_to_query)

    # Set threshold for pChEMBL activity value and limit to Toxicity Assays
    pchembl_value = 1  # Specify a minimum threshold of the pChEMBL activity value. Note that pCHEMBL = -log10(IC50, XC50, AC50, Ki, Kd, potency). Greater than or equal to 5 (10um) is a typical minimum rule of thumb for binding activity between a compound and a protein target.
    assay_type = "T"  # Only look for Toxicity Assays (Cytotoxicity)

    # Define the base URL for the ChEMBL API
    host = "https://www.ebi.ac.uk"  # This is the stem of the url
    url = f"{host}/chembl/api/data/activity?molecule_pref_name__in={drugs_str}&pchembl_value__gte={pchembl_value}&assay_type={assay_type}&format=json&limit=1000"

    # Query the API and convert the response to JSON
    response = requests.get(
        url
    ).json()  # This calls the information back from the API using the 'requests' module, and converts it to json format

    # Extract the list of activities from the response
    results = response["activities"]

    # If there are additional pages, fetch them
    while response["page_meta"]["next"]:
        response = requests.get(host + response["page_meta"]["next"]).json()
        results = results + response["activities"]
    return results


def produce_ranked_drug_candidates(
    drug_scores: list[dict],
    drug_cytotoxicity_chembl: list[dict],
    drug_bbb: list[dict],
):
    """
    Produce a ranked list of drug candidates based on confidence z-scores and cytotoxicity.

    Args:
        drug_scores (list[dict]): List of drugs with z-scores indicating confidence levels.
        drug_cytotoxicity_chembl (list[dict]): Cytotoxicity data from the ChEMBL database.

    Returns:
        list[dict]: Ranked list of drug candidates based on cytotoxicity and confidence z-scores.
    """
    # Convert drug scores and cytotoxicity data into DataFrames
    drug_scores_df = pd.DataFrame.from_dict(drug_scores).rename(
        columns={"zscore": "confidence_zscore", "term": "drug_name"}
    )
    drug_cytotoxicity_df = (
        pd.DataFrame.from_dict(drug_cytotoxicity_chembl)
        .add_prefix("cytotoxicity_")
        .rename(columns={"cytotoxicity_molecule_pref_name": "drug_name"})
    )

    bbb_df = pd.DataFrame.from_dict(drug_bbb)

    # Convert drug names to lowercase for consistency
    drug_cytotoxicity_df["drug_name"] = drug_cytotoxicity_df[
        "drug_name"
    ].str.lower()

    # Merge cytotoxicity data with drug scores based on drug names
    ranked_drug_candidates_df = pd.merge(
        drug_scores_df, drug_cytotoxicity_df, how="left", on="drug_name"
    )

    ranked_drug_candidates_df = pd.merge(
        ranked_drug_candidates_df, bbb_df, how="left", on="drug_name"
    )

    ranked_drug_candidates_df = ranked_drug_candidates_df.fillna("")

    # Specify columns to keep
    # cytotoxicity_cols
    drug_cytotoxicity_cols = [
        "cytotoxicity_activity_id",
        "cytotoxicity_assay_description",
        "cytotoxicity_assay_type",
        "cytotoxicity_standard_type",
        "cytotoxicity_standard_units",
        "cytotoxicity_standard_value",
        "cytotoxicity_target_pref_name",
        "cytotoxicity_pchembl_value",
    ]

    drug_score_cols = [
        "drug_name",
        "confidence_zscore",
    ]

    bbb_cols = [
        "logBB",
        "bbb_permeable",
    ]

    ranked_drug_candidates_df = ranked_drug_candidates_df[
        drug_score_cols + drug_cytotoxicity_cols + bbb_cols
    ]

    # Return the ranked list as a list of dictionaries
    return ranked_drug_candidates_df.sort_values(
        "confidence_zscore", ascending=False
    ).to_dict(orient="records")


def query_drug_bbb_from_b3db(drug_scores: list[dict]):
    """
    Query the B3DB database for blood-brain barrier (BBB) permeability data and merge it with drug scores.

    Args:
        drug_scores (list[dict]): List of drugs with their respective scores, where each drug has a "term" key.

    Returns:
        list[dict]: A list of dictionaries containing drug names, logBB values, and BBB permeability status.
    """
    # Convert the list of drug scores into a DataFrame and rename "term" column to "drug_name"
    drug_scores_df = pd.DataFrame.from_dict(drug_scores).rename(
        columns={"term": "drug_name"}
    )

    # URL to the B3DB classification data (Blood-Brain Barrier permeability information)
    b3db_url = "https://raw.githubusercontent.com/theochem/B3DB/refs/heads/main/B3DB/B3DB_classification.tsv"

    # Read the B3DB data from the URL as a DataFrame, specifying tab as the delimiter
    b3db_df = pd.read_csv(b3db_url, sep="\t")

    # Perform a left join to merge the drug scores with the B3DB data on the drug names
    # Retain only the relevant columns: drug name, logBB value (permeability coefficient), and BBB status
    df = pd.merge(
        drug_scores_df, b3db_df, left_on="drug_name", right_on="compound_name"
    )[["drug_name", "logBB", "BBB+/BBB-"]]

    # Replace any NaN values with empty strings for cleaner output
    df = df.fillna("")
    # Remove any duplicate rows that may have been introduced during the merge
    df.drop_duplicates()

    # Filter values out that don't have both logBB and BBB+/BBB- fields
    df = df[(df["logBB"] != "") & (df["BBB+/BBB-"] != "")]

    # Rename the "BBB+/BBB-" column to a more readable format, "bbb_permeable"
    df = df.rename(columns={"BBB+/BBB-": "bbb_permeable"})
    # Sort the DataFrame by drug name for better readability and convert the result to a list of dictionaries
    return df.sort_values("drug_name").to_dict(orient="records")
