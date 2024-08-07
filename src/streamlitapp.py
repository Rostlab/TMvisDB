import logging
import os

import pandas as pd
import streamlit as st
from peewee import OperationalError

from views import (
    faq,
    overview,
    protein_detail,
    protein_list,
    about,
    sidebar,
    header,
)
from utils import database, api
from utils.api import UniprotACCType
from utils.protein_visualization import ColorScheme, VizFilter, Style
from utils.database import DBFilter
from utils.protein_info import ProteinInfo


def db_error():
    st.warning(
        "We are having trouble finding the predicted transmembrane topology of your protein in TMvisDB. "  # noqa: E501
        "This could mean, e.g., (1) your protein is outside the length restrictions of TMvisDB (see FAQ), (2) your protein is not predicted as a transmembrane protein, or (3) the UniProt ID is misspelled. "  # noqa: E501
        "If an AlphaFold structure is displayed below, it is without transmembrane topology annotation.",  # noqa: E501
        icon="🚨",
    )


def collect_and_display_protein_info(db_conn, selected_id):
    uniprot_acc_type = api.uniprot_get_input_type(selected_id)

    if uniprot_acc_type == UniprotACCType.UNKNOWN:
        st.error(
            f"The input format of your selected ID ** {selected_id} ** is not correct.",
            icon="🚨",
        )
        return None

    protein_info: ProteinInfo = ProteinInfo.collect_for_id(selected_id)
    if not protein_info.has_annotations:
        st.warning(
            "We found no transmembrane annotation; neither predicted nor in UniProt or TmAlphaFold.",  # noqa: E501
            icon="🚨",
        )

    if protein_info.has_annotations and not protein_info.uniprot_accession:
        st.warning(
            f"We could not find a protein in UniProtKB/TrEMBL matching your ID: {selected_id}."  # noqa: E501
        )

    # get and protein structure
    if protein_info.sequence is not None:
        if protein_info.uniprot_accession != protein_info.uniprot_name:
            st.write(
                f"Displaying protein with UniProt accession number: {protein_info.uniprot_accession} and UniProt entry name: {protein_info.uniprot_name}"  # noqa: E501
            )
        else:
            st.write(
                f"Displaying protein with ID: {protein_info.uniprot_accession if protein_info.uniprot_accession is not None else selected_id}"
            )

    return protein_info


def display_gdpr_banner():
    """
    Display a removable GDPR compliance banner using Streamlit components.
    """
    # Initialize session state for banner visibility
    if "show_banner" not in st.session_state:
        st.session_state.show_banner = True

    # Display the banner if it hasn't been dismissed
    if st.session_state.show_banner:
        with st.container():
            cols = st.columns([8, 2])
            with cols[0]:
                st.warning(
                    "TMvisDB does not collect usage statistics. However, this site uses Streamlit, which may process some data. "
                    "Please review [Streamlit's Privacy Policy](https://streamlit.io/privacy-policy) for details."
                )
            with cols[1]:
                if st.button("Dismiss"):
                    st.session_state.show_banner = False
                    st.rerun()

    # Add some vertical space after the banner
    st.markdown("<br>", unsafe_allow_html=True)


@st.cache_resource
def initialize_database_connection():
    """
    Initialize the Database connection.
    Returns the database connection if successful or None if unsuccessful.
    """
    try:
        return database.initialize_database_connection()
    except OperationalError:
        logging.exception("Failed to connect to SQLite")
        st.error(
            "Error establishing a connection to TMvisDB! Please try again later, and/or contact us here: service+tmvisdb@rostlab.org",  # noqa: E501
            icon="🚨",
        )
        return None
    except Exception:
        logging.exception("An error occurred during database initialization")
        st.error(
            "An unexpected error occurred. Please try again later.",
            icon="🚨",
        )
        return None


def initialize_session_state():
    """
    Initialize the Streamlit session state variables.
    """
    default_state = {
        "data": pd.DataFrame(),
        "user_display": "",
        "database_filter": DBFilter(),
        "visualization_filter": VizFilter(style=Style(), selected_id="Q9NVH1"),
    }
    for key, value in default_state.items():
        if key not in st.session_state:
            st.session_state[key] = value


def handle_list_tab(db_conn):
    """
    Handle the logic for the 'Database' tab.
    """
    if db_conn is None:
        db_error()
        return

    database_filter = st.session_state.database_filter

    if database_filter.random_selection:
        protein_list.display_random_data(database_filter)
    else:
        protein_list.display_filtered_data(database_filter)

    if not st.session_state.user_display == "":
        st.markdown(st.session_state.user_display)
        st.markdown("---")

    if not st.session_state.data.empty:
        protein_list.show_table(st.session_state.data, paginate=True)
        st.download_button(
            "Download selection",
            protein_list.convert_df(st.session_state.data),
            "file.csv",
            "text/csv",
            key="download-csv",
        )


def show_3d_visualization(visualization_filter: VizFilter):
    """
    Display 3D visualization of the protein.
    """
    try:
        protein_info = ProteinInfo.collect_for_id(visualization_filter.selected_id)
        protein_detail.create_visualization_for_id(
            protein_info, visualization_filter.style
        )
    except Exception:
        logging.exception(
            f"An error occurred while visualizing protein {visualization_filter.selected_id}:"
        )
        st.error(
            "Sorry, we could not visualize your selected protein. Please contact us, so we can help you with your search.",  # noqa: E501
            icon="🚨",
        )


def maintenance_mode():
    st.title("We'll be back soon!")
    st.image(
        "https://via.placeholder.com/600x200.png?text=Maintenance+Mode",
        use_column_width=True,
    )
    st.warning(
        "The TMvisDB is currently undergoing maintenance to improve your experience. Please check back later.",
        icon="🚧",
    )
    st.info(
        """
        **Why am I seeing this?**

        - We are performing essential maintenance or upgrades.
        - We apologize for the inconvenience and appreciate your patience.

        **When will the service be back?**

        - We expect the service to be restored within a couple of days.
        - Please check back tomorrow.

        **Need assistance?**

        - If you have any questions, please contact us at [service+tmvisdb@rostlab.org](mailto:service+tmvisdb@rostlab.org).
        """
    )


def _setup_logging():
    log_level = os.getenv("LOG_LEVEL", "ERROR").upper()
    numeric_log_level = getattr(logging, log_level, logging.ERROR)
    logging.basicConfig(level=numeric_log_level)


def main():
    st.set_page_config(page_title="TMvisDB", page_icon="⚛️", layout="wide")

    _setup_logging()

    maintenance_mode_enabled = os.getenv("MAINTENANCE_MODE", "false").lower()
    logging.debug(f"MAINTENANCE_MODE: {maintenance_mode_enabled}")

    # Header
    display_gdpr_banner()
    header.title()

    if maintenance_mode_enabled == "true":
        logging.debug("Entering maintenance mode")
        maintenance_mode()
        return

    try:
        db_conn = initialize_database_connection()
        initialize_session_state()

        # Sidebar
        sidebar.display_sidebar()

        tab_overview, tab_database, tab_visualization, tab_faq, tab_about = st.tabs(
            ["Overview", "Database", "Visualization", "FAQ", "About"]
        )

        # Tabs handling
        with tab_overview:
            overview.intro()

        with tab_database:
            handle_list_tab(db_conn)

            st.markdown("---")
            if not st.session_state.data.empty:
                options = st.session_state.data[
                    ["UniProt Accession", "UniProt ID"]
                ].drop_duplicates()
                options_dict = dict(
                    zip(options["UniProt ID"], options["UniProt Accession"])
                )

                local_id = st.selectbox(
                    "Choose an ID to visualize predicted transmembrane topology below",
                    options=options["UniProt ID"],
                    format_func=lambda x: f"{x} ({options_dict[x]})",
                    index=0,
                )

                filter = VizFilter(
                    style=Style(color_scheme=ColorScheme.TRANSMEMBRANE_PREDICTION),
                    selected_id=options_dict[local_id],
                )

                with st.spinner("Loading Protein Data"):
                    show_3d_visualization(filter)
            st.markdown("---")

        with tab_visualization:
            # TODO Move Sidebar filter here
            show_3d_visualization(st.session_state.visualization_filter)
            st.markdown("---")

        with tab_faq:
            faq.quest()

        with tab_about:
            about.handle_about()

    finally:
        if db_conn is not None and not db_conn.is_closed():
            db_conn.close()


if __name__ == "__main__":
    main()
