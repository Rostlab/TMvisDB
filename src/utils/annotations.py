# Copyright 2024 Tobias Olenyi.
# SPDX-License-Identifier: Apache-2.0

from dataclasses import dataclass, field
from collections import defaultdict
from enum import Enum

import pandas as pd

from .database import Annotation


class AnnotationSource(Enum):
    MEMBRANOME = "membranome"
    TOPDB = "topdb"
    TMBED = "tmbed"
    UNIPROT = "uniprot"
    TMALPHAFOLD = "tmalphafold"


DISPLAY_NAMES = {
    AnnotationSource.MEMBRANOME: "Membranome Annotation",
    AnnotationSource.TOPDB: "TopDB Annotation",
    AnnotationSource.TMBED: "TMbed Prediction",
    AnnotationSource.UNIPROT: "UniProt Annotation",
    AnnotationSource.TMALPHAFOLD: "TmAlphaFold Annotation",
}


def annotations_from_db(annotations: Annotation):
    parsed_annotations = defaultdict(list)
    reference_urls = dict()

    for annotation in annotations:
        parsed_annotations[annotation.source].append(
            ResidueAnnotation(annotation.start, annotation.end, annotation.label)
        )
        reference_urls[annotation.source] = (
            annotation.source_db_url
        )  # FIXME This is a bit uggly and should be refactored; Constant overwrites.

    return parsed_annotations, reference_urls


@dataclass
class ResidueAnnotation:
    start: int
    end: int
    label: str


@dataclass
class MembraneAnnotation:
    annotations: dict[AnnotationSource, list[ResidueAnnotation]] = field(
        default_factory=dict
    )
    reference_urls: dict[AnnotationSource, str] = field(default_factory=dict)


def construct_df_from_annotation(annotation: MembraneAnnotation, sequence: list[str]):
    return pd.DataFrame(
        zip(
            sequence,
            *[
                ANNOTATION_HANDLER_MAP[source](
                    sequence, annotation.annotations.get(source)
                )
                for source in annotation.annotations
            ],
        ),
        columns=["Sequence"]
        + [DISPLAY_NAMES[source] for source in annotation.annotations],
    )


# TODO Fix this


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


def construct_tmembd_annotation(item: None | dict[str, str | dict[str, str]]):
    pass


ANNOTATION_HANDLER_MAP = {
    "membranome": construct_membranome_annotation,
    "topdb": construct_topdb_annotation,
    "tembed": construct_tmembd_annotation,
}
