from dataclasses import dataclass
from enum import Enum
import os
import logging

import pymongo
import pandas as pd
import numpy as np

from .protein_info import MembraneAnnotation


# Options: Filter
class Topology(Enum):
    ALL = "All"
    BOTH = "Both"
    ALPHA_HELIX = "Alpha-helix"
    BETA_STRAND = "Beta-strand"


class TaxaSelectionCriterion(Enum):
    ORGANISM = "Organism ID"
    DOMAIN = "Domain/Kingdom"


class Domain(Enum):
    ALL = "All"
    BACTERIA = "Bacteria"
    EUKARYOTA = "Eukaryota"
    ARCHAEA = "Archaea"
    UNCLASSIFIED = "unclassified sequences"


class Kingdom(Enum):
    pass


class Archaea(Kingdom):
    ALL = "All Archaea"
    ASGARD_GROUP = "Asgard group"
    HYDROTHERMARCHAEOTA = "Candidatus Hydrothermarchaeota"
    THERMOPLASMATOTA = "Candidatus Thermoplasmatota"
    DPANN_GROUP = "DPANN group"
    EURYARCHAEOTA = "Euryarchaeota"
    TACK_GROUP = "TACK group"
    ARCHA_INCERTAE_SEDIS = "Archaea incertae sedis"
    UNCLASSIFIED_ARCHAEA = "unclassified Archaea"
    ENVIRONMENTAL_SAMPLES = "environmental samples"


class Eukaryota(Kingdom):
    ALL = "All Eukaryota"
    AMOEBOZOA = "Amoebozoa"
    ANCYROMONADIDA = "Ancyromonadida"
    APUSOZOA = "Apusozoa"
    BREVITEA = "Breviatea"
    CRUMS = "CRuMs"
    CRYPTOPHYCEAE = "Cryptophyceae (cryptomonads)"
    DISCOBA = "Discoba"
    GLAUCOCYSTOPHYCEAE = "Glaucocystophyceae"
    HAPTISTA = "Haptista"
    HEMIMASTIGOPHORA = "Hemimastigophora"
    MALAWIMONADIDA = "Malawimonadida"
    METAMONADA = "Metamonada"
    OPISTHOKONTA = "Opisthokonta"
    RHODELPHEA = "Rhodelphea"
    RHODOPHYTA = "Rhodophyta (red algae)"
    SAR = "Sar"
    VIRIDIPLANTAE = "Viridiplantae"
    EUKARYOTA_INCERTAE_SEDIS = "Eukaryota incertae sedis"
    UNCLASSIFIED_EUKARYOTES = "unclassified eukaryotes"
    ENVIRONMENTAL_SAMPLES = "environmental samples"


class Bacteria(Kingdom):
    ALL = "All Bacteria"
    ACIDOBACTERIA = "Acidobacteria"
    AQUIFICAE = "Aquificae"
    ATRIBACTEROA = "Atribacterota"
    CALDISERICA_CRYOSERICOTA_GROUP = "Caldiserica/Cryosericota group"
    CALDITRICHAEOTA = "Calditrichaeota"
    CANDIDATUS_KRUMHOLZIBACTERIOTA = "Candidatus Krumholzibacteriota"
    CANDIDATUS_THARPELLOTA = "Candidatus Tharpellota"
    CHRYSIOGENETES = "Chrysiogenetes"
    COLEOSPERMUM = "Coleospermum"
    COPROTHERMOBACTEROTA = "Coprothermobacterota"
    DEFERRIBACTERES = "Deferribacteres"
    DESULFOBACTEROTA = "Desulfobacterota"
    DICTYOGLOMI = "Dictyoglomi"
    ELUSIMICROBIA = "Elusimicrobia"
    FCB_GROUP = "FCB group"
    FUSOBACTERIA = "Fusobacteria"
    MYXOCCOCOTA = "Myxococcota"
    NITROSFINAE_TECTOMICROBIA_GROUP = "Nitrospinae/Tectomicrobia group"
    NITROSPIRAE = "Nitrospirae"
    PROTEOBACTERIA = "Proteobacteria"
    PVC_GROUP = "PVC group"
    SPIROCHAETES = "Spirochaetes"
    SYNERGISTETES = "Synergistetes"
    TERRABACTERIA_GROUP = "Terrabacteria group"
    THERMODESULFOBACTERIA = "Thermodesulfobacteria"
    THERMOTOGAE = "Thermotogae"
    BACTERIA_INCERTAE_SEDIS = "Bacteria incertae sedis"
    UNCLASSIFIED_BACTERIA = "unclassified Bacteria"
    ENVIRONMENTAL_SAMPLES = "environmental samples"


AllKingdoms = Kingdom(
    "AllKingdoms",
    {
        "ALL": "All",
        **{item.name: item.value for item in Archaea if item.name != "ALL"},
        **{item.name: item.value for item in Eukaryota if item.name != "ALL"},
        **{item.name: item.value for item in Bacteria if item.name != "ALL"},
    },
)

DOMAIN_MAP = {
    Domain.ARCHAEA: Archaea,
    Domain.BACTERIA: Bacteria,
    Domain.EUKARYOTA: Eukaryota,
    Domain.ALL: AllKingdoms,
    Domain.UNCLASSIFIED: AllKingdoms,
}


def get_kingdom_for_domain(domain: Domain):
    if domain in Domain:
        kingdom_type = DOMAIN_MAP[domain]
    else:
        kingdom_type = AllKingdoms
    return kingdom_type


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
class DBFilter:
    taxonomy_selection: TaxaSelectionCriterion = TaxaSelectionCriterion.ORGANISM
    organism_id: int | None = 9606
    domain: Domain = Domain.ALL
    kingdom: Kingdom = AllKingdoms.ALL
    topology: Topology = Topology.ALL
    signal_peptide: bool = False
    sequence_lengths: tuple[int, int] = (16, 5500)
    num_sequences: int = 100
    random_selection: bool = True

    def construct_query(self):
        sp = int(self.signal_peptide)

        query_form = dict()
        # sequence length
        if self.sequence_lengths != (16, 5500):
            query_form["seq_length"] = {
                "$gt": self.sequence_lengths[0],
                "$lt": self.sequence_lengths[1],
            }

        if self.topology == Topology.BOTH:
            query_form["annotations.tm_categorical"] = [1, 1, sp]
        elif self.topology == Topology.ALPHA_HELIX:
            query_form["annotations.tm_categorical"] = [1, 0, sp]
        elif self.topology == Topology.BETA_STRAND:
            query_form["annotations.tm_categorical"] = [0, 1, sp]
        elif self.topology == Topology.ALL:
            query_form["annotations.tm_categorical"] = [1, 1, sp]

        if (
            self.taxonomy_selection == TaxaSelectionCriterion.ORGANISM
            and self.organism_id is not None
        ):
            query_form["organism_id"] = int(self.organism_id)

        else:
            if self.domain != Domain.ALL:
                query_form["uptaxonomy.Domain"] = self.domain.value

            kingdom_type = get_kingdom_for_domain(self.domain)
            if self.kingdom != kingdom_type.ALL:
                query_form["uptaxonomy.Kingdom"] = self.kingdom.value

        logging.debug(query_form)
        return query_form

    def __str__(self):
        parts = []

        # Topology
        if self.topology is not None:
            parts.append(
                f"topology: {self.topology.value} ({' not' if (not self.signal_peptide and self.topology is not Topology.ALL) else ''} including signal peptides )"  # noqa: E501
            )

        # Taxonomy
        taxonomy_parts = []
        if self.organism_id is not None:
            taxonomy_parts.append(f"Organism ID: {self.organism_id}")
        if self.domain is not None:
            taxonomy_parts.append(f"Domain: {self.domain.value}")
        if self.kingdom is not None:
            taxonomy_parts.append(f"Kingdom: {self.kingdom.value}")

        if taxonomy_parts:
            parts.append(f"taxonomy ({', '.join(taxonomy_parts)})")

        # Sequence Lengths
        if self.sequence_lengths != (16, 5500):
            parts.append(
                f"length range: {self.sequence_lengths[0]}-{self.sequence_lengths[1]}"
            )
        else:
            parts.append("all lengths")

        # Final String
        return ", ".join(parts)


def init_connection():
    host = os.environ.get("TMVIS_MONGO_HOST", "localhost")
    port = int(os.environ.get("TMVIS_MONGO_PORT", 27017))
    username = os.environ.get("TMVIS_MONGO_USERNAME", "")
    password = os.environ.get("TMVIS_MONGO_PASSWORD", "")
    auth_source = os.environ.get("TMVIS_MONGO_DB", "admin")

    # Log the connection parameters (excluding password for security reasons)
    logging.debug(
        f"Connecting to MongoDB with host: {host}, port: {port}, username: {username}, authSource: {auth_source}"  # noqa: E501
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

    df = pd.json_normalize(items)
    df = df.rename(columns=FIELDS)
    return df


def get_filtered_data(db_conn, query: dict[str, str | int], sample_size: int):
    items = db_conn.tmvis.find(query, DATA_FORM).limit(sample_size)
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


def get_data_for_id(db_conn, selected_id: str):
    query = {"_id": selected_id}
    return get_filtered_data(db_conn, query, sample_size=1)


def get_membrane_annotation_for_id(db_conn, selected_id: str):
    query = {"_id": selected_id}
    membranome_form = {
        "seq_length": 1,
        "predictions": 1,
        "topdb.TopDB_Entry": 1,
        "membranomedb.tm_seq_start": 1,
        "membranomedb.tm_seq_end": 1,
    }

    tmbed_annotation: list[str] = None
    topdb_annotation: list[str] = None
    membdb_annotation: list[str] = None

    item = list(db_conn.tmvis.find(query, membranome_form))

    if not item:
        raise ValueError("Item not found in database.")

    # TMbed prediction
    tmbed_annotation = list(item[0]["predictions"]["transmembrane"])

    # TopDB annotation
    topdb_annotation = construct_topdb_annotation(item[0])

    # Membranome annotation
    membdb_annotation = construct_membranome_annotation(item[0])

    return MembraneAnnotation(tmbed_annotation, topdb_annotation, membdb_annotation)


def construct_topdb_annotation(item: None | dict[str, str | dict[str, str]]):
    topdb_value = item.get("topdb", {}).get("TopDB_Entry", None)
    return None if not topdb_value else list(topdb_value)


def construct_membranome_annotation(item: None | dict[str, str | dict[str, str]]):
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
