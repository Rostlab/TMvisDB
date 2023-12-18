from dataclasses import dataclass
from enum import Enum
import re
import logging

import streamlit as st

from utils.db import TaxaSelectionCriterion, Domain, Topology, DBFilter
from utils import db
from utils import api
from utils.coloring import ProteinStyle, ColorScheme, Style

sb = st.sidebar


def display_sidebar_header():
    sb.markdown(
        "<style>div.block-container{padding-top:1rem;}</style>", unsafe_allow_html=True
    )
    sb.markdown("---")
    sb.subheader("Search TMvisDB")
    sb.caption("Please open the 'Database' tab to see results.")


def create_random_form():
    emp = sb.empty()
    emp.button(
        "Show random selection",
        help="Click here to show 100 random proteins of TMvisDB.",
        on_click=handle_random_selection,
    )


def create_filter_form():
    with sb.expander("Click here to access filters for TMvisDB."):
        taxonomy_selection = st.radio(
            "Select Taxonomy via",
            TaxaSelectionCriterion,
            format_func=lambda x: x.value,
            key="taxonomy_selection",
        )

        disable_domain = True
        if taxonomy_selection == TaxaSelectionCriterion.DOMAIN:
            disable_domain = False

        st.number_input(
            "Enter Organism ID",
            min_value=0,
            help="Type in UniProt Organism ID.",
            placeholder=9606,
            disabled=not disable_domain,
            value=None,
            key="organism_id",
        )

        # Select Taxonomy: Domain/Kingdom
        domain = st.selectbox(
            "Select Domain",
            Domain,
            format_func=lambda x: x.value,
            help="Type domain or select from list.",
            disabled=disable_domain,
            key="domain",
        )

        kingdom_type = db.get_kingdom_for_domain(domain)

        st.selectbox(
            "Select Kingdom",
            kingdom_type,
            format_func=lambda x: x.value,
            help="Type kingdom or select from list.",
            disabled=disable_domain,
            key="kingdom",
        )

        topology = st.selectbox(
            "Filter by Transmembrane Topology ",
            Topology,
            format_func=lambda x: x.value,
            help="TMbed predicts per-residue transmembrane topology as either alpha-helical or beta-stand.",  # noqa: E501
            key="topology",
        )
        st.checkbox(
            "Show sequences with signal peptides",
            value=False,
            help="TMbed predicts whether a sequence contains signal peptides. Change 'Filter by Transmembrane Topology' to access.",  # noqa: E501
            disabled=(topology == Topology.ALL),
            key="signal_peptide",
        )

        st.slider(
            "Select sequence length",
            16,
            5500,
            (16, 5500),
            help="Select a minimum and maximum value for sequence length.",
            key="sequence_lengths",
        )

        st.number_input(
            "Select limit of shown sequences",
            1,
            10000,
            value=100,
            help="As TMvisDB is a large database, you may want to set a limit for your table.",
            key="num_sequences",
        )

        # Submit results
        st.button(
            "Apply filters",
            help="Click here to show your selection.",
            on_click=handle_db_filter,
        )


def create_vis_form():
    st.sidebar.subheader("Visualize predicted transmembrane proteins")
    sb.caption("Please open the 'Visualization' tab to see results.")

    with sb.expander("Click here to access 3D visualization for single proteins."):
        st.text_input("Insert Uniprot ID", placeholder="Q9NVH1")

        st.selectbox(
            "Style",
            ProteinStyle,
            key="style",
            format_func=lambda x: x.value,
        )
        # select color
        st.selectbox(
            "Color Scheme",
            ColorScheme,
            format_func=lambda x: x.value,
            key="color_scheme",
        )
        # select spin
        st.checkbox("Spin", value=False, key="spin")

        # Submit results
        st.button(
            "Visualize Protein",
            help="Click here to visualize your selection.",
            on_click=handle_vis_changes,
        )


def end():
    sb.markdown("---")
    st.sidebar.write(
        "Authors: [Céline Marquet](https://github.com/C-Marquet), [Tobias Olenyi](https://github.com/t03i)"
    )
    st.sidebar.write("Source: [Github](https://github.com/rostlab/TMvisDB)")


def handle_db_filter():
    attributes = [
        "taxonomy_selection",
        "sequence_lengths",
        "topology",
        "organism_id",
        "domain",
        "kingdom",
        "signal_peptide",
        "num_sequences",
    ]
    filter_kwargs = {
        attr: getattr(st.session_state, attr)
        for attr in attributes
        if hasattr(st.session_state, attr)
    }
    st.session_state.filter = DBFilter(**filter_kwargs, random_selection=False)


def handle_random_selection():
    st.session_state.filter = DBFilter()


def handle_vis_changes():
    selected_id = getattr(st.session_state, "selected_id", "Q9NVH1").strip().upper()
    input_format = api.check_input_format(selected_id)
    # TODO fix this
    protein_info = api.get_protein_info(selected_id, input_format)
    attributes = ["style", "color_scheme", "spin"]
    style_kwargs = {
        attr: getattr(st.session_state, attr)
        for attr in attributes
        if hasattr(st.session_state, attr)
    }
    style = Style(**style_kwargs)
    st.session_state.visualization_protein = protein_info
    st.session_state.visualization_style = style


def display_sidebar():
    display_sidebar_header()
    create_random_form()
    sb.markdown("---")
    create_filter_form()
    sb.markdown("---")
    create_vis_form()
    end()
