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


def construct_df_from_annotation(annotation: MembraneAnnotation, sequence: str):
    return pd.DataFrame(
        zip(
            list(sequence),
            *[
                construct_annotation_vector(
                    annotation.annotations[source], len(sequence)
                )
                for source in annotation.annotations
                if annotation.annotations[source] is not None
                and len(annotation.annotations[source]) > 0
            ],
        ),
        columns=["Sequence"]
        + [
            DISPLAY_NAMES[source]
            for source in annotation.annotations
            if annotation.annotations[source] is not None
            and len(annotation.annotations[source]) > 0
        ],
    )


def construct_annotation_vector(annotations: list[ResidueAnnotation], seuence_len: int):
    annotation_vector = ["*"] * seuence_len
    for annotation in annotations:
        for i in range(annotation.start - 1, annotation.end):
            annotation_vector[i] = annotation.label
    return annotation_vector
