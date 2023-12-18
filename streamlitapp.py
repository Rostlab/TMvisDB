import os
import streamlit as st
import pymongo
import logging
from pymongo.errors import ConnectionFailure
import pandas as pd
import logging

from app import database, faq, overview, visualization, about, sidebar, header
from utils import db, api
from utils.db import DBFilter
from utils.api import UniprotACCType
from utils.protein_info import ProteinInfo
from utils.coloring import Style


def db_error():
    st.warning(
        "We are having trouble finding the predicted transmembrane topology of your protein in TMvisDB. "  # noqa: E501
        "This could mean, e.g., (1) your protein is outside the length restrictions of TMvisDB (see FAQ), (2) your protein is not predicted as a transmembrane protein, or (3) the UniProt ID is misspelled. "  # noqa: E501
        "If an AlphaFold structure is displayed below, it is without transmembrane topology annotation.",  # noqa: E501
        icon="🚨",
    )


def collect_and_display_protein_info(db_conn, selected_id):
    uniprot_acc_type = api.check_input_format(selected_id)

    if uniprot_acc_type == UniprotACCType.UNKNOWN:
        st.error(
            "The input format of your selected ID ** "
            + selected_id
            + " ** is not correct.",
            icon="🚨",
        )

        return

    protein_info: ProteinInfo = ProteinInfo.collect_for_id(
        db_conn, selected_id, uniprot_acc_type
    )
    if not protein_info.has_annotations:
        st.warning(
            "We found no transmembrane annotation; neither predicted nor in UniProt or TmAlphaFold.",  # noqa: E501
            icon="🚨",
        )

    if protein_info.has_annotations and protein_info.uniprot_accession is None:
        st.warning(
            f"We could not find a protein in UniProtKB/TrEMBL matching your ID: {selected_id}."  # noqa: E501
        )

    # get and protein structure
    if protein_info.sequence is not None:
        if protein_info.uniprot_accession != protein_info.uniprot_name:
            st.write(
                "Displaying protein with UniProt accession number: ",
                protein_info.uniprot_accession,
                " and UniProt entry name:",
                protein_info.uniprot_name,
            )
        else:
            st.write(
                "Displaying protein with ID: ",
                protein_info.uniprot_accession
                if protein_info.uniprot_accession is not None
                else selected_id,
            )

    return protein_info


def acknowledge_statistics_warning():
    """
    Display a warning about usage statistics and ask for user acknowledgement.
    """
    if "user_acknowledged_stats" not in st.session_state:
        st.session_state.user_acknowledged_stats = False

    if not st.session_state.user_acknowledged_stats:
        stats_warning_message = st.warning(
            "Welcome to TMvisDB. The authors of TMvisDB opted out of gathering any usage summary statistics.  \n"  # noqa: E501
            "However, this web application is implemented with Streamlit. "
            "Please familiarize yourself with [Streamlit's privacy policy](https://streamlit.io/privacy-policy) before proceeding. "  # noqa: E501
            "The authors of TMvisDB have no insight into or control over Streamlit's data collection process and, thus, cannot accept any liability for said process.",  # noqa: E501
            icon="🚨",
        )
        continue_button_placeholder = st.empty()
        user_clicked_continue = continue_button_placeholder.button(
            "Click here to continue to TMvisDB."
        )
        if user_clicked_continue:
            st.session_state.user_acknowledged_stats = True
            continue_button_placeholder.empty()
            stats_warning_message.empty()


@st.cache_resource
def initialize_database_connection():
    """
    Initialize the MongoDB connection.
    Returns the database connection if successful or None if unsuccessful.
    """
    client = db.init_connection()

    try:
        client.admin.command("ismaster")
        return client.microscope
    except (
        pymongo.errors.ServerSelectionTimeoutError,
        pymongo.errors.ConnectionFailure,
    ):
        logging.error("Failed to connect to MongoDB.")
        st.error(
            "Error establishing a connection to TMvisDB! Please try again later, and/or contact us here: tmvisdb@rostlab.org",  # noqa: E501
            icon="🚨",
        )
        return None
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        return None


def initialize_session_state():
    """
    Initialize the Streamlit session state variables.
    """
    default_state = {
        "data": pd.DataFrame(),
        "user_display": "",
        "filter": DBFilter(),
        "visualization_style": Style(),
        "visualization_protein_info": None,
    }
    for key, value in default_state.items():
        if key not in st.session_state:
            st.session_state[key] = value


@st.cache_data(ttl=600, show_spinner=False)
def get_initial_protein_info(_db_conn):
    selected_id = getattr(st.session_state, "selected_id", "Q9NVH1").strip().upper()
    input_format = api.check_input_format(selected_id)
    protein_info = ProteinInfo.collect_for_id(_db_conn, selected_id, input_format)
    return protein_info


def handle_database_tab(db_conn):
    """
    Handle the logic for the 'Database' tab.
    """
    if db_conn is None:
        db_error()
        return

    db_filter = st.session_state.filter

    if db_filter.random_selection:
        database.display_random_data(db_conn, db_filter)
    else:
        database.display_filtered_data(db_conn, db_filter)

    if not st.session_state.data.empty:
        database.show_table(st.session_state.data)
        st.download_button(
            "Download selection",
            database.convert_df(st.session_state.data),
            "file.csv",
            "text/csv",
            key="download-csv",
        )


def show_3d_visualization(db_conn, protein_info: ProteinInfo, style: Style):
    """
    Display 3D visualization of the protein.
    """
    if db_conn is None:
        db_error()
        return

    # try:
    visualization.create_visualization_for_id(protein_info, style)
    # except Exception as e:
    #     logging.error(f"An error occurred: {e}")
    #     st.error(
    #         "Sorry, we could not visualize your selected protein. Please contact us, so we can help you with your search.",  # noqa: E501
    #         icon="🚨",
    #     )


def main():
    st.set_page_config(page_title="TMvisDB", page_icon="⚛️", layout="wide")
    logging.basicConfig(level=logging.DEBUG)

    acknowledge_statistics_warning()

    if st.session_state.user_acknowledged_stats:
        db_conn = initialize_database_connection()

        initialize_session_state()

        if st.session_state.visualization_protein_info is None:
            st.session_state.visualization_protein_info = get_initial_protein_info(
                db_conn
            )

        # Sidebar
        sidebar.display_sidebar()

        # Header
        header.title()
        tab_overview, tab_database, tab_visualization, tab_faq, tab_about = st.tabs(
            ["Overview", "Database", "Visualization", "FAQ", "About"]
        )

        # Tabs handling
        with tab_overview:
            overview.intro()

        with tab_database:
            handle_database_tab(db_conn)

            st.markdown("---")
            if not st.session_state.data.empty:
                local_id = st.selectbox(
                    "Choose an ID to visualize predicted transmembrane topology below",
                    st.session_state.data["UniProt ID"],
                    0,
                )

                with st.spinner("Loading Protein Data"):
                    protein_info = collect_and_display_protein_info(db_conn, local_id)
                    if protein_info is not None:
                        show_3d_visualization(
                            db_conn, protein_info, st.session_state.visualization_style
                        )
            st.markdown("---")

        with tab_visualization:
            show_3d_visualization(
                db_conn,
                st.session_state.visualization_protein_info,
                st.session_state.visualization_style,
            )
            st.markdown("---")

        with tab_faq:
            faq.quest()

        with tab_about:
            about.handle_about()


if __name__ == "__main__":
    main()
