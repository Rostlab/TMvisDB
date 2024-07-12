# Copyright 2023 RostLab.
# SPDX-License-Identifier: 	AGPL-3.0-or-later


from dataclasses import dataclass, field, fields
import logging

import pandas as pd

from utils import db, api
from utils.api import UniprotACCType


@dataclass
class MembraneAnnotation:
    predicted: list[str] | None = field(
        default=None, metadata={"col_name": "TMbed Prediction"}
    )  # noqa: E501
    topdb: list[str] | None = field(
        default=None, metadata={"col_name": "TopDB Annotation"}
    )  # noqa: E501
    membdb: list[str] | None = field(
        default=None, metadata={"col_name": "Membranome Annotation"}
    )  # noqa: E501
    uniprot: list[str] | None = field(
        default=None, metadata={"col_name": "UniProt Annotation"}
    )  # noqa: E501
    alphafold: list[str] | None = field(
        default=None, metadata={"col_name": "TmAlphaFold Annotation"}
    )  # noqa: E501

    @property
    def has_an_predicted(self):
        return self.predicted is not None

    @property
    def has_an_topdb(self):
        return self.topdb is not None

    @property
    def has_an_membdb(self):
        return self.membdb is not None

    @property
    def has_an_uniprot(self):
        return self.uniprot is not None

    @property
    def has_an_alphafold(self):
        return self.alphafold is not None

    @property
    def __annotation_fields(self):
        return [
            annotation
            for annotation in fields(self)
            if "col_name" in annotation.metadata
        ]

    @property
    def has_annotations(self):
        return any(
            getattr(self, f"has_an_{annotation.name}")
            for annotation in self.__annotation_fields
        )

    @property
    def available_annotations(self):
        return [
            annotation.metadata.get("col_name")
            for annotation in self.__annotation_fields
            if getattr(self, f"has_an_{annotation.name}")
        ]

    def construct_annotation_table(self, sequence: list[str]):
        used_annotation_fields = [
            annotation
            for annotation in self.__annotation_fields
            if getattr(self, annotation.name) is not None
        ]

        return pd.DataFrame(
            zip(
                sequence,
                *[
                    getattr(self, annotation.name)
                    for annotation in used_annotation_fields
                ],
            ),
            columns=["Sequence"]
            + [
                annotation.metadata.get("col_name", "")
                for annotation in used_annotation_fields
            ],
        )


@dataclass
class ProteinInfo:
    supplied_accession: str
    uniprot_accession: str
    uniprot_name: str
    sequence: str
    structure: bytes
    annotations: MembraneAnnotation
    info_df: pd.DataFrame

    @property
    def has_annotations(self):
        return (
            self.annotations.has_an_predicted
            or self.annotations.has_an_topdb
            or self.annotations.has_an_membdb
            or self.annotations.has_an_uniprot
            or self.annotations.has_an_alphafold
        )

    @staticmethod
    def collect_for_id(db_conn, selected_id: str, uniprot_acc_type: UniprotACCType):
        # get uniprot annotation
        (
            uniprot_accession,
            uniprot_name,
            uniprot_annotation,
            uniprot_seq_length,
        ) = api.get_uniprot_tmvec(selected_id, uniprot_acc_type)

        db_annotation = db.get_membrane_annotation_for_id(db_conn, selected_id)

        alphafold_annotation = api.get_tmalphafold_annotation(
            uniprot_name if uniprot_name is not None else selected_id,
            uniprot_seq_length,
        )

        db_annotation.uniprot = uniprot_annotation
        db_annotation.alphafold = alphafold_annotation

        sequence_info_df = db.get_data_for_id(db_conn, selected_id)
        sequence, structure = api.get_af_structure(
            uniprot_accession if uniprot_accession is not None else selected_id
        )

        return ProteinInfo(
            supplied_accession=selected_id,
            uniprot_accession=uniprot_accession,
            uniprot_name=uniprot_name,
            sequence=sequence,
            structure=structure,
            annotations=db_annotation,
            info_df=sequence_info_df,
        )
