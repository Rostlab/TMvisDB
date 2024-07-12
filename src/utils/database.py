from dataclasses import dataclass
import os
from functools import reduce
from operator import and_

import pandas as pd
from peewee import (
    SqliteDatabase,
    Model,
    CharField,
    ForeignKeyField,
    IntegerField,
    AutoField,
    TextField,
    DateField,
    FloatField,
    BooleanField,
    fn,
)

from .protein_info import MembraneAnnotation
import utils.lineage_definitions as lineage_definitions
from .lineage_definitions import (
    TaxaSelectionCriterion,
    Domain,
    Kingdom,
    Topology,
    AllKingdoms,
)


# Read the DATABASE_URL environment variable
DATABASE_URL = os.getenv("DATABASE_URL", "sqlite:///data/tmvis.db")

# Database connection
DATABASE = SqliteDatabase(
    DATABASE_URL.split("///")[-1],
    pragmas={"journal_mode": "off", "cache_size": -1024 * 64},
)


class BaseModel(Model):
    class Meta:
        database = DATABASE


class Organism(BaseModel):
    id = AutoField(primary_key=True)
    taxon_id = CharField(index=True, unique=True)
    name = CharField()
    super_kingdom = CharField(index=True)
    clade = CharField(index=True, null=True)


class Sequence(BaseModel):
    id = AutoField(primary_key=True)
    uniprot_id = CharField()
    uniprot_accession = CharField(index=True, unique=True)
    organism = ForeignKeyField(Organism, backref="sequences")
    sequence = TextField()
    seq_length = IntegerField()


class TMInfo(BaseModel):
    id = AutoField(primary_key=True)
    sequence = ForeignKeyField(Sequence, backref="tm_info")
    tm_helix_count = IntegerField()
    tm_helix_percent = FloatField()
    tm_strand_count = IntegerField()
    tm_strand_percent = FloatField()
    signal_count = IntegerField()
    signal_percent = FloatField()
    generated_at = DateField()
    has_alpha_helix = BooleanField()
    has_beta_strand = BooleanField()
    has_signal = BooleanField()


class Annotation(BaseModel):
    id = AutoField(primary_key=True)
    sequence = ForeignKeyField(Sequence, backref="annotations")
    start = IntegerField()
    end = IntegerField()
    label = CharField(
        max_length=100,
    )
    date_added = DateField()
    source_db = CharField(choices=["topdb", "membranome", "tmvis"])
    source_db_ref = CharField(null=True)
    source_db_url = CharField(max_length=400, null=True)


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
        query = Sequence.select()

        if not self.random_selection:
            filters = (
                self.sequence_length_filter()
                + self.topology_filter()
                + self.taxonomy_filter()
                + self.signal_peptide_filter()
            )

            if filters:
                combined_filters = reduce(and_, filters)
                query = query.where(combined_filters)

        query = query.join(TMInfo).join(Organism)

        if self.random_selection:
            query = query.order_by(fn.Random())

        query = query.limit(self.num_sequences)
        return query

    def sequence_length_filter(self):
        filters = []
        if self.sequence_lengths != (16, 5500):
            filters.append(
                Sequence.seq_length.between(
                    self.sequence_lengths[0], self.sequence_lengths[1]
                )
            )
        return filters

    def topology_filter(self):
        filters = []
        if self.topology != Topology.ALL:
            if self.topology == Topology.BOTH:
                filters.append(
                    (TMInfo.has_alpha_helix == True)  # noqa: E712
                    & (TMInfo.has_beta_strand == True)  # noqa: E712
                )
            elif self.topology == Topology.ALPHA_HELIX:
                filters.append(TMInfo.has_alpha_helix == True)  # noqa: E712
            elif self.topology == Topology.BETA_STRAND:
                filters.append(TMInfo.has_beta_strand == True)  # noqa: E712
                filters.append(TMInfo.has_signal == self.signal_peptide)
        return filters

    def taxonomy_filter(self):
        filters = []
        if (
            self.taxonomy_selection == TaxaSelectionCriterion.ORGANISM
            and self.organism_id is not None
        ):
            filters.append(Organism.taxon_id == str(self.organism_id))
        else:
            if self.domain != Domain.ALL:
                filters.append(Organism.super_kingdom == self.domain.value)
            kingdom_type = lineage_definitions.get_kingdom_for_domain(self.domain)
            if self.kingdom != kingdom_type.ALL:
                filters.append(Organism.clade == self.kingdom.value)
        return filters


def get_sequence_data(db_conn: SqliteDatabase, db_filter: DBFilter):
    query = db_filter.construct_query()
    items = db_conn.execute(query)
    dict_items = items.dicts()
    if dict_items is not None:
        return pd.DataFrame(dict_items)
    return None


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
