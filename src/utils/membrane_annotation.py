# Copyright 2024 Tobias Olenyi.
# SPDX-License-Identifier: Apache-2.0

import logging
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


def annotations_from_db(annotations: list[Annotation]):
    parsed_annotations = defaultdict(list)
    reference_urls = {}

    for annotation in annotations:
        try:
            source = AnnotationSource(annotation.source_db)
        except ValueError:
            # Log a warning if the source is not recognized
            logging.warning(f"Unrecognized annotation source: {annotation.source_db}")
            continue

        parsed_annotations[source].append(
            ResidueAnnotation(annotation.start, annotation.end, annotation.label)
        )

        # Only set the reference URL if it's not already set for this source
        if source not in reference_urls and annotation.source_db_url:
            reference_urls[source] = annotation.source_db_url

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

    def add_annotation(
        self, source: AnnotationSource, annotation: list[ResidueAnnotation]
    ):
        self.annotations[source] = annotation

    def add_reference_url(self, source: AnnotationSource, url: str):
        self.reference_urls[source] = url

    def update_annotations(
        self, new_annotations: dict[AnnotationSource, list[ResidueAnnotation]]
    ):
        self.annotations.update(new_annotations)

    def update_reference_urls(self, new_urls: dict[AnnotationSource, str]):
        self.reference_urls.update(new_urls)

    @property
    def has_annotations(self):
        return any(len(self.annotations[source]) > 0 for source in self.annotations)


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
