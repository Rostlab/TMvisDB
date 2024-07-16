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
DATABASE = SqliteDatabase(DATABASE_URL.split("///")[-1])


def set_pragma_settings(db):
    db.execute_sql(f"PRAGMA cache_size = {-1024 * 64}")  # Set cache size to 64MB
    db.execute_sql("PRAGMA journal_mode = off")  # Set journal mode to WAL
    db.execute_sql("PRAGMA synchronous = NORMAL")  # Set synchronous mode to NORMAL
    db.execute_sql("PRAGMA temp_store = MEMORY")  # Store temporary tables in memory
    db.execute_sql(
        "PRAGMA mmap_size = 268435456"
    )  # Use memory-mapped I/O for performance
    db.execute_sql("PRAGMA read_uncommitted = true")  # Allow read uncommitted
    db.execute_sql("PRAGMA optimize")  # Run the optimize pragma


def initialize_database_connection():
    DATABASE.connect()
    set_pragma_settings(DATABASE)
    return DATABASE


class BaseModel(Model):
    class Meta:
        database = DATABASE


class Organism(BaseModel):
    id = AutoField(primary_key=True)
    taxon_id = CharField(index=True, unique=True)
    name = CharField()
    super_kingdom = CharField(index=True)
    clade = CharField(index=True, null=True)

    class Meta:
        indexes = (
            (
                ("super_kingdom", "clade"),
                False,
            ),  # Composite index on super_kingdom and clade
        )


class Sequence(BaseModel):
    id = AutoField(primary_key=True)
    uniprot_id = CharField()
    uniprot_accession = CharField(index=True, unique=True)
    organism = ForeignKeyField(Organism, backref="sequences", index=True)
    sequence = TextField()
    seq_length = IntegerField(index=True)


class TMInfo(BaseModel):
    id = AutoField(primary_key=True)
    sequence = ForeignKeyField(Sequence, backref="tm_info", index=True)
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

    class Meta:
        indexes = (
            (
                ("has_alpha_helix", "has_beta_strand", "has_signal"),
                False,
            ),  # Composite index
        )


class Annotation(BaseModel):
    id = AutoField(primary_key=True)
    sequence = ForeignKeyField(Sequence, backref="annotations", index=True)
    start = IntegerField()
    end = IntegerField()
    label = CharField(
        max_length=100,
    )
    date_added = DateField()
    source_db = CharField(choices=["topdb", "membranome", "tmvis"])
    source_db_ref = CharField(null=True)
    source_db_url = CharField(max_length=400, null=True)

    class Meta:
        indexes = (
            (
                ("sequence", "start", "end"),
                False,
            ),  # Composite index on sequence, start, and end
        )


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
        query = Sequence.select(
            Sequence.uniprot_id,
            Sequence.uniprot_accession,
            Sequence.seq_length,
            Organism.name,
            Organism.super_kingdom,
            Organism.clade,
            TMInfo.tm_helix_count,
            TMInfo.tm_strand_count,
            TMInfo.signal_count,
        )

        if not self.random_selection:
            filters = (
                self.sequence_length_filter()
                + self.topology_filter()
                + self.taxonomy_filter()
            )

            if filters:
                combined_filters = reduce(and_, filters)
                query = query.where(combined_filters)

        query = query.join(TMInfo).switch(Sequence).join(Organism)

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


def get_sequence_data(db_filter: DBFilter):
    query = db_filter.construct_query()
    return query


def get_sequence_data_for_id(selected_id: str):
    return Sequence.get_or_none(Sequence.uniprot_accession == selected_id)


def get_membrane_annotation_for_id(selected_id: str):
    sequence = Sequence.get_or_none(Sequence.uniprot_accession == selected_id)

    annotations = Annotation.select().where(Annotation.sequence == sequence)

    if len(annotations) == 0:
        raise ValueError("Could not find any annotations for the given sequence.")

    return annotations
