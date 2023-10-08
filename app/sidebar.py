from dataclasses import dataclass
from enum import Enum

import streamlit as st
import re

from utils.db import TaxaSelectionCriterion, Domain, Topology, DBFilter
from utils import db

sb = st.sidebar


####################################################################


def display_sidebar_header():
    sb.markdown(
        "<style>div.block-container{padding-top:1rem;}</style>", unsafe_allow_html=True
    )
    sb.markdown("---")
    sb.subheader("Search TMvisDB")
    sb.caption("Please open the 'Database' tab to see results.")


def handle_random_selection():
    emp = sb.empty()
    rand_selected = emp.button(
        "Show random selection",
        help="Click here to show 100 random proteins of TMvisDB.",
        disabled=st.session_state.rndm,
        key="1",
    )
    return rand_selected


def handle_taxonomy_selection():
    taxonomy_selection = st.radio("Select Taxonomy via", TaxaSelectionCriterion)

    disable_domain = True
    if taxonomy_selection == TaxaSelectionCriterion.DOMAIN:
        disable_domain = False

    # Select Taxonomy: Organism ID
    organism_id = st.text_input(
        "Enter Organism ID",
        help="Type in UniProt Organism ID.",
        placeholder="9606",
        disabled=not disable_domain,
        value="",
    )

    # Select Taxonomy: Domain/Kingdom
    domain = st.selectbox(
        "Select Domain",
        Domain,
        help="Type domain or select from list.",
        disabled=disable_domain,
    )

    kingdom_type = db.get_kingdom_for_domain(domain)

    kingdom = st.selectbox(
        "Select Kingdom",
        kingdom_type,
        help="Type kingdom or select from list.",
        disabled=disable_domain,
    )

    return taxonomy_selection, organism_id, domain, kingdom


def filters():
    display_sidebar_header()

    random_selected = handle_random_selection()
    if random_selected:
        st.session_state.rndm = True
        return DBFilter()
    else:
        with sb.expander("Click here to access filters for TMvisDB."):
            (
                taxonomy_selection,
                organism_id,
                domain,
                kingdom,
            ) = handle_taxonomy_selection()

            # Select TMP type
            topology = st.selectbox(
                "Filter by Transmembrane Topology ",
                Topology,
                help="TMbed predicts per-residue transmembrane topology as either alpha-helical or beta-stand.",
            )
            signal_peptide = st.checkbox(
                "Show sequences with signal peptides",
                value=False,
                help="TMbed predicts whether a sequence contains signal peptides. Change 'Filter by Transmembrane Topology' to access.",
                disabled=(topology != Topology.ALL),
            )

            # Sequence length range
            sequence_length = st.slider(
                "Select sequence length",
                16,
                5500,
                (16, 5500),
                help="Select a minimum and maximum value for sequence length.",
            )

            # Number of shown sequences
            num_sequences = st.number_input(
                "Select limit of shown sequences",
                1,
                10000,
                value=100,
                help="As TMvisDB is a large database, you may want to set a limit for your table.",
            )

            # Submit results
            subm = st.button(
                "Submit filters",
                help="Click here to show your selection.",
                disabled=st.session_state.filt,
            )
            if subm:
                st.session_state.filt = True

        return DBFilter(
            taxa_selection=taxonomy_selection,
            organism_id=organism_id,
            domain=domain,
            kingdom=kingdom,
            topology=topology,
            signal_peptide=signal_peptide,
            sequence_lengths=sequence_length,
            num_sequences=num_sequences,
            random_selection=False,
        )


def vis():
    sb.markdown("---")
    st.sidebar.subheader("Visualize predicted transmembrane proteins")
    sb.caption("Please open the 'Visualization' tab to see results.")

    with sb.expander("Click here to access 3D visualization for single proteins."):
        # select ID
        selected_id = st.text_input("Insert Uniprot ID", placeholder="Q9NVH1")
        selected_id = selected_id.strip().upper()
        # selected_id = re.sub(r'[^a-zA-Z0-9]','', selected_id).upper()
        # select style
        style = st.selectbox(
            "Style", ["Cartoon", "Line", "Cross", "Stick", "Sphere"]
        ).lower()
        # select color
        color_prot = st.selectbox(
            "Color Scheme", ["Transmembrane Prediction", "Alphafold pLDDT score"]
        )
        # select spin
        spin = st.checkbox("Spin", value=False)
        if selected_id == "":
            selected_id = "Q9NVH1"
    return selected_id, style, color_prot, spin


def end():
    sb.markdown("---")
    st.sidebar.write(
        "Authors: [Céline Marquet](https://github.com/C-Marquet), [Tobias Olenyi](https://github.com/t03i)"
    )
    st.sidebar.write("Source: [Github](https://github.com/rostlab/TMvisDB)")
