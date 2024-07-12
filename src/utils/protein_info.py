# Copyright 2023 RostLab.
# SPDX-License-Identifier: 	AGPL-3.0-or-later
from dataclasses import dataclass, field, fields
import logging

import pandas as pd

from utils import database, api
from utils.api import UniprotACCType
from utils import annotations
from utils.annotations import MembraneAnnotation, AnnotationSource


def collect_annotations_for_id(selected_id: str, uniprot_acc_type: UniprotACCType):
    annotation = MembraneAnnotation()

    (
        uniprot_accession,
        uniprot_name,
        uniprot_annotation,
        uniprot_seq_length,
    ) = api.get_uniprot_tmvec(selected_id, uniprot_acc_type)

    if uniprot_annotation is not None:
        annotation[AnnotationSource.UNIPROT] = uniprot_annotation
        annotation.reference_urls[AnnotationSource.UNIPROT] = (
            f"https://www.uniprot.org/uniprotkb/{uniprot_accession}/entry"
        )

    alphafold_annotation, url = api.get_tmalphafold_annotation(
        uniprot_name if uniprot_name is not None else selected_id,
        uniprot_seq_length,
    )

    if alphafold_annotation is not None:
        annotation[AnnotationSource.ALPHAFOLD] = alphafold_annotation
        annotation.reference_urls[AnnotationSource.ALPHAFOLD] = url

    db_annotations = database.get_membrane_annotation_for_id(selected_id)
    parsed_db_annoations, parsed_db_refs = annotations.annotations_from_db(
        db_annotations
    )

    annotation.annotations |= parsed_db_annoations
    annotation.reference_urls |= parsed_db_refs

    return annotation, uniprot_accession, uniprot_name


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
    def collect_for_id(
        db_conn,
        selected_id: str,
    ):
        # get uniprot annotation

        db_annotation, uniprot_name, uniprot_accession = (
            database.get_membrane_annotation_for_id(db_conn, selected_id)
        )

        sequence_info_df = database.get_sequence_data_for_id(db_conn, selected_id)
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
