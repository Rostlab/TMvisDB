import logging

import py3Dmol
from stmol import showmol
import streamlit as st
import pandas as pd
from st_aggrid import AgGrid


from utils.protein_visualization import (
    ALPHAFOLD_LEGEND_DF,
    ANNOTATION_LEGEND_DF,
)
from utils import protein_visualization

from utils.protein_info import ProteinInfo
from utils.protein_visualization import Style, ColorScheme
from utils import annotations


def display_legend(color_scheme: ColorScheme, has_no_pred):
    """Displays color code explanation based on the provided color protocol."""
    st.write("Color code")
    if color_scheme == ColorScheme.ALPHAFOLD_PLDDT_SCORE or has_no_pred:
        st.write(
            ALPHAFOLD_LEGEND_DF.style.applymap(
                protein_visualization.alphafold_legend_coloring, subset=["pLDDT score"]
            )
        )
    else:
        st.write(
            ANNOTATION_LEGEND_DF.style.applymap(
                protein_visualization.annotation_legend_coloring, subset=["Color"]
            )
        )
        st.caption(
            "Inside/outside annotations of TMbed are not optimized and must be interpreted with caution."  # noqa: E501
        )


def display_links(protein_info: ProteinInfo):
    # TOD refactor based on links
    """Displays links to various resources for further evaluation."""
    links = [
        f"- UniProt entry: [{protein_info.uniprot_accession}](https://www.uniprot.org/uniprotkb/{protein_info.uniprot_accession}/entry)",
        f"- Evaluate protein-specific phenotype predictions: [LambdaPP](https://lambda.predictprotein.org/o/{protein_info.uniprot_accession})",
        "- Generate structural alignments: [Foldseek](https://search.foldseek.com/search)",
        "- Experimentally derived topology information: [Topology Data Bank of Transmembrane Proteins](https://topdb.unitmp.org//)",  # noqa: E501
        "- Membranome database for single-helix transmembrane proteins: [Membranome](https://membranome.org/)",
        f"- Alpha-helical transmembrane proteins: [TmAlphaFold database](https://tmalphafold.ttk.hu/entry/{selected_id})",
    ]

    st.markdown("Resources to evaluate your selection further:")
    st.markdown("\n".join(links))


def format_available_annotations(annotations):

    # TODO refactor according to sources

    def prefix_annotation(annotation):
        """Prefixes the annotation with 'a' or 'an'."""
        first_word = annotation.split()[0].lower()
        prefix = "an" if first_word[0] in ["a", "e", "i", "o", "u"] else "a"
        return f"{prefix} {annotation.title()}"

    if not annotations:
        return "Found no available annotations."
    else:
        formatted_annotations = [prefix_annotation(ann) for ann in annotations]
        if len(formatted_annotations) == 1:
            connected = formatted_annotations[0]
        elif len(formatted_annotations) == 2:
            connected = f"{formatted_annotations[0]} and {formatted_annotations[1]}"
        else:
            connected = f"{', '.join(formatted_annotations[:-1])}, and {formatted_annotations[-1]}"
        return f"Found {connected}."


def display_membrane_annotation(protein_info: ProteinInfo):
    """Visualizes membrane annotations using the provided data."""

    st.Markdown("## Membrane Annotations")
    if not protein_info.annotation.has_annotations:
        st.caption("Could not find any annotations for this protein.")
        return

    pred_table = annotations.construct_df_from_annotation(
        protein_info.annotation, protein_info.sequence
    )
    styled_table = pred_table.T.style.apply(
        protein_visualization.color_prediction, axis=0
    )

    st.write(
        format_available_annotations(protein_info.annotation.available_annotations)
    )

    st.write(styled_table)

    st.caption(
        "If entries in a row are '*', there are no annotations for these residues."  # noqa: E501
    )


def display_other_annotations(protein_info: ProteinInfo):
    # Further sequence annotation
    st.write("Protein Annotation")
    AgGrid(
        protein_info.info_df.drop(columns=["Sequence", "Prediction"]),
        height=100,
        fit_columns_on_grid_load=True,
    )


def display_protein_structure(protein_info: ProteinInfo, style: Style):
    structure = protein_info.structure
    view = py3Dmol.view(js="https://3dmol.org/build/3Dmol.js")
    view.addModelsAsFrames(structure)
    view.setBackgroundColor("#262730")

    # add color
    if (
        style.color_scheme == ColorScheme.ALPHAFOLD_PLDDT_SCORE
        or not protein_info.annotation.has_an_predicted
    ):
        view.setStyle(
            {"model": -1},
            {
                style.visualization_style.value: {
                    "colorscheme": {
                        "prop": "b",
                        "gradient": "roygb",
                        "min": 50,
                        "max": 90,
                    }
                }
            },
        )
    else:
        tm_color = protein_visualization.map_annotation_to_color(
            protein_info.annotation.predicted
        )
        view.setStyle(
            {"model": -1},
            {
                style.visualization_style.value: {
                    "colorscheme": {"prop": "resi", "map": tm_color}
                }
            },
        )

    view.spin(style.spin)

    view.zoomTo()
    showmol(view, height=500, width=800)


def create_visualization_for_id(
    protein_info: ProteinInfo,
    style: Style,
):
    display_protein_structure(protein_info, style)

    st.markdown("---")

    display_membrane_annotation(protein_info)

    display_other_annotations(protein_info)

    # Explain colors used in the visualization
    display_legend(
        color_scheme=style.color_scheme,
        has_no_pred=not protein_info.annotation.has_an_predicted,
    )

    st.markdown("---")

    # Display links to other resources
    display_links(protein_info)
