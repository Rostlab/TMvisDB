from dataclasses import dataclass
import os
import logging

import pymongo
import pandas as pd
import numpy as np

FIELDS = {
    "_id": "UniProt ID",
    "sequence": "Sequence",
    "predictions.transmembrane": "Prediction",
    "annotations.tm_categorical": "Alpha, Beta, Signal",
    "seq_length": "Sequence length",
    "organism_name": "Organism name",
    "organism_id": "Organism ID",
    "uptaxonomy.Lineage_all": "Lineage",
    "uptaxonomy.Domain": "Domain",
    "uptaxonomy.Kingdom": "Kingdom",
}

DATA_FORM = {elem: 1 for elem in FIELDS}


@dataclass
class MembraneAnnotation:
    tmbed_annotation: list[str] = None
    topdb_annotation: list[str] = None
    membdb_annotation: list[str] = None

    @property
    def has_tmbed(self):
        return self.tmbed_annotion is not None

    @property
    def has_topdb(self):
        return self.tobdb_annotion is not None

    @property
    def has_membranome(self):
        return self.membdb_annotion is not None


def construct_query(
    selected_organismid,
    selected_domain,
    selected_kingdom,
    selected_type,
    show_signal_peptide,
    selected_length,
):
    sp = int(show_signal_peptide)

    query_form = dict()
    # sequence length
    if selected_length != (16, 5500):
        query_form["seq_length"] = {
            "$gt": selected_length[0],
            "$lt": selected_length[1],
        }

    # add topology filter if selected_type not "All"
    if "Both" in selected_type:
        query_form["annotations.tm_categorical"] = [1, 1, sp]
    elif "Alpha-helix" in selected_type:
        query_form["annotations.tm_categorical"] = [1, 0, sp]
    elif "Beta-strand" in selected_type:
        query_form["annotations.tm_categorical"] = [0, 1, sp]

    if (
        selected_organismid != ""
        and selected_organismid != "0"
        and selected_organismid.isnumeric()
    ):
        query_form["organism_id"] = int(selected_organismid)

    # add filter for domain and kingdom
    if "All" not in selected_domain and selected_organismid == "0":
        query_form["uptaxonomy.Domain"] = selected_domain

    if "All" not in selected_kingdom and selected_organismid == "0":
        query_form["uptaxonomy.Kingdom"] = selected_kingdom

    return query_form


def init_connection():
    host = os.environ.get("TMVIS_MONGO_HOST", "localhost")
    port = int(os.environ.get("TMVIS_MONGO_PORT", 27017))
    username = os.environ.get("TMVIS_MONGO_USERNAME", "")
    password = os.environ.get("TMVIS_MONGO_PASSWORD", "")
    auth_source = os.environ.get("TMVIS_MONGO_DB", "admin")

    # Log the connection parameters (excluding password for security reasons)
    logging.debug(
        f"Connecting to MongoDB with host: {host}, port: {port}, username: {username}, authSource: {auth_source}"
    )

    return pymongo.MongoClient(
        host=host,
        port=port,
        username=username,
        password=password,
        authSource=auth_source,
        appname="TMvis Frontend",
        serverSelectionTimeoutMS=5000,  # 5 seconds timeout for server selection
        connectTimeoutMS=5000,  # 5 seconds timeout for connection
    )


def get_random_data(db_conn, sample_size: int):
    items = db_conn.tmvis.aggregate(
        [
            {"$sample": {"size": sample_size}},
            {"$project": DATA_FORM},
        ]
    )

    logging.debug(items)

    df = pd.json_normalize(items)
    df = df.rename(columns=FIELDS)
    return df


def get_filtered_data(db_conn, query: dict[str, str | int], sample_size: int):
    items = db_conn.find(query, DATA_FORM).limit(sample_size)
    df = pd.json_normalize(items)

    if df.empty:  # No items matching criteria
        return None

    df = df.rename(columns=FIELDS)

    # Check for non-existing columns and add them
    for col_id, col_name in FIELDS.items():
        if col_name not in df.columns:
            df[col_name] = None

    # Ensure the order of the columns is consistent with FIELDS
    df = df[[col_name for col_name in FIELDS.values()]]

    return df


def construct_membranome_annotation(item: dict[str, str | dict[str, str]]):
    """
    Extracts Membranome data from the given item.
    """
    if "membranomedb" not in item:
        return None

    seq_length = int(item["seq_length"]) - 1
    pos_start = int(item["membranomedb"]["tm_seq_start"]) - 1
    pos_end = int(item["membranomedb"]["tm_seq_end"]) - 1

    membdb_annotation = ["*"] * seq_length
    membdb_annotation[pos_start:pos_end] = ["AH"] * (pos_end - pos_start + 1)

    return membdb_annotation


# TODO Error Message:
# st.warning(
#     "We are having trouble finding the predicted transmembrane topology of your protein in TMvisDB. "
#     "This could mean, e.g., (1) your protein is outside the length restrictions of TMvisDB (see FAQ), (2) your protein is not predicted as a transmembrane protein, or (3) the UniProt ID is misspelled. "
#     "If an AlphaFold structure is displayed below, it is without transmembrane topology annotation.",
#     icon="🚨",
# )


def get_membrane_annotation_for_id(db, selected_id):
    query = {"_id": selected_id}
    membranome_form = {
        "seq_length": 1,
        "topdb.TopDB_Entry": 1,
        "membranomedb.tm_seq_start": 1,
        "membranomedb.tm_seq_end": 1,
    }

    tmbed_annotation: list[str] = None
    topdb_annotation: list[str] = None
    membdb_annotation: list[str] = None

    item = list(db.find(query, membranome_form))

    if not item:
        raise ValueError("Item not found in database.")

    # TMbed prediction
    tmbed_annotation = list(item[0]["predictions"]["transmembrane"])

    # TopDB annotation
    topdb_annotation = list(item[0].get("topdb", {}).get("TopDB_Entry", None))

    # Membranome annotation
    membdb_annotation = construct_membranome_annotation(item[0])

    return MembraneAnnotation(tmbed_annotation, topdb_annotation, membdb_annotation)
