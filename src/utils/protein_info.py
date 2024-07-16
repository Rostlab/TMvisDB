# Copyright 2023 RostLab.
# SPDX-License-Identifier: 	AGPL-3.0-or-later
from dataclasses import dataclass

import pandas as pd

from utils import database, api
from utils import membrane_annotation
from utils.membrane_annotation import MembraneAnnotation, AnnotationSource


FIELDS = {
    "uniprot_accession": "UniProt Accession",
    "uniprot_id": "UniProt ID",
    "seq_length": "Sequence length",
    "name": "Organism name",
    "taxon_id": "Organism ID",
    "super_kingdom": "Domain",
    "clade": "Kingdom",
    "tm_has_alpha_helix": "Has alpha helix",
    "tm_has_beta_strand": "Has beta strand",
    "tm_has_signal": "Has signal peptide",
    "tm_helix_count": "Transmembrane helix residues",
    "tm_strand_count": "Transmembrane strand residues",
    "signal_count": "Signal peptide residues",
}


def db_to_df(query_result: database.Sequence):
    data = query_result.dicts()
    df = pd.DataFrame(data)
    df.rename(columns=FIELDS, inplace=True)
    return df


def fetch_membrane_annotations(selected_id: str):
    annotation = MembraneAnnotation()
    uniprot_response = api.uniprot_fetch_annotation(selected_id)

    if uniprot_response is not None:
        annotation.add_annotation(
            AnnotationSource.UNIPROT, uniprot_response.membrane_annotations
        )
        annotation.add_reference_url(
            AnnotationSource.UNIPROT,
            f"https://www.uniprot.org/uniprotkb/{uniprot_response.accession}/entry",
        )

    tmalphafold_annotation = api.tmalphafold_fetch_annotation(
        uniprot_response.name
        if uniprot_response is not None and uniprot_response.name is not None
        else selected_id
    )

    if tmalphafold_annotation is not None:
        annotation.add_annotation(AnnotationSource.TMALPHAFOLD, tmalphafold_annotation)
        annotation.add_reference_url(
            AnnotationSource.TMALPHAFOLD,
            f"https://tmalphafold.ttk.hu/entry/{uniprot_response.accession if uniprot_response is not None and uniprot_response.accession is not None else selected_id}",
        )

    db_annotations = database.get_membrane_annotation_for_id(selected_id)
    parsed_db_annotations, parsed_db_refs = membrane_annotation.annotations_from_db(
        db_annotations
    )

    annotation.update_annotations(parsed_db_annotations)
    annotation.update_reference_urls(parsed_db_refs)

    return annotation, uniprot_response


def fetch_sequence_data(selected_id: str):
    sequence_info_df = db_to_df(database.get_sequence_data_for_id(selected_id))
    return sequence_info_df


@dataclass
class ProteinInfo:
    supplied_accession: str
    uniprot_accession: str
    uniprot_name: str
    sequence: str
    structure: bytes
    annotation: MembraneAnnotation
    info_df: pd.DataFrame

    @property
    def has_annotations(self):
        return len(self.annotation.annotations) > 0

    @staticmethod
    def collect_for_id(
        selected_id: str,
    ):
        annotation, uniprot_info = fetch_membrane_annotations(selected_id)

        sequence_info_df = fetch_sequence_data(selected_id)

        sequence, structure = api.alphafolddb_fetch_structure(
            uniprot_info.accession
            if uniprot_info.accession is not None
            else selected_id
        )

        return ProteinInfo(
            supplied_accession=selected_id,
            uniprot_accession=uniprot_info.accession,
            uniprot_name=uniprot_info.name,
            sequence=sequence,
            structure=structure,
            annotation=annotation,
            info_df=sequence_info_df,
        )
