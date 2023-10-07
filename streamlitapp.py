import os
import streamlit as st
import pymongo
import logging
from pymongo.errors import ConnectionFailure
import pandas as pd

from app import faq, table, overview, visualization, about, sidebar, header
from utils import db, api
from utils.protein_info import ProteinInfo


def collect_and_display_protein_info(db_conn, selected_id):
    uniprot_acc_type = api.check_input_format(selected_id)

    if uniprot_acc_type == "unknown":
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

    if protein_info.has_annotations and protein_info.uniprot_accession is not None:
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
            st.write("Displaying protein with ID: ", protein_info.uniprot_accession)

    return protein_info


if __name__ == "__main__":
    st.set_page_config(page_title="TMvisDB", page_icon="⚛️", layout="wide")

    logging.basicConfig(level=logging.DEBUG)


# Check if the user has acknowledged the usage statistics warning
if "user_acknowledged_stats" not in st.session_state:
    st.session_state.user_acknowledged_stats = False

# Display the usage statistics warning if not acknowledged
if not st.session_state.user_acknowledged_stats:
    stats_warning_message = st.warning(
        "Welcome to TMvisDB. The authors of TMvisDB opted out of gathering any usage summary statistics.  \n"
        "However, this web application is implemented with Streamlit. "
        "Please familiarize yourself with [Streamlit's privacy policy](https://streamlit.io/privacy-policy) before proceeding. "
        "The authors of TMvisDB have no insight into or control over Streamlit's data collection process and, thus, cannot accept any liability for said process.",
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


# If the user acknowledged the warning, initialize the database connection
if st.session_state.user_acknowledged_stats:
    # Function to initialize the MongoDB connection

    client = db.init_connection()

    # Attempt to connect to the database and test the connection
    db_conn = None
    try:
        client.admin.command("ismaster")
        db_conn = client.microscope
    except pymongo.errors.ServerSelectionTimeoutError:
        logging.error("Could not connect to MongoDB. Server selection timed out.")
    except pymongo.errors.ConnectionFailure:
        logging.error("Failed to connect to MongoDB.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")

    if db_conn is None:
        st.error(
            "Error establishing a connection to TMvisDB! Please try again later, and/or contact us here: tmvisdb@rostlab.org",
            icon="🚨",
        )

    ## Initialize session ##
    if "rndm" not in st.session_state:
        st.session_state.rndm = False
    if "filt" not in st.session_state:
        st.session_state.filt = False
    if "tbl" not in st.session_state:
        st.session_state.tbl = pd.DataFrame()
    if "txt" not in st.session_state:
        st.session_state.txt = ""

    ## Sidebar ##
    [
        selected_organismid,
        selected_domain,
        selected_kingdom,
        selected_type,
        selected_sp,
        selected_limit,
        select_random,
        selected_length,
    ] = sidebar.filters()
    [selected_id, style, color_prot, spin] = sidebar.vis()
    sidebar.end()

    ## Header ##
    header.title()
    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        ["Overview", "Database", "Visualization", "FAQ", "About"]
    )

    ## Overview ##
    with tab1:
        overview.intro()

    ## Database ##
    with tab2:  #
        if db_conn is None:
            st.error(
                "Error establishing a connection to TMvisDB! Please try again later, and/or contact us here: tmvisdb@rostlab.org",
                icon="🚨",
            )
        else:
            df = pd.DataFrame()

            if st.session_state.rndm:
                st.session_state.filt = False
                st.session_state.txt = "The table below shows a random selection. To personalize your selection, use the sidebar filters. (Note: Your current random selection will not be saved after reloading.)"
                with st.spinner("Loading your data"):
                    df = db.get_random_data(db_conn, selected_limit)
                    st.session_state.tbl = df
                st.session_state.rndm = False
                st.experimental_rerun()

            if st.session_state.filt:
                st.session_state.rndm = False
                query_form = db.construct_query(
                    selected_organismid,
                    selected_domain,
                    selected_kingdom,
                    selected_type,
                    selected_sp,
                    selected_length,
                )

                with st.spinner("Loading your data"):
                    df = db.get_filtered_data(db_conn, query_form, selected_limit)

                    st.session_state.txt = (
                        "The table below shows your personalized selection: "
                        "topology (" + selected_type + "), "
                        "taxonomy (Organism ID: "
                        + selected_organismid
                        + ", Domain: "
                        + selected_domain
                        + ", Kingdom: "
                        + selected_kingdom
                        + "), length "
                        + str(selected_length)
                        + "."
                        " For a random selection use the sidebar button."
                    )
                    st.session_state.tbl = df
                    # st.session_state.txt = ''
                st.session_state.filt = False
                st.experimental_rerun()

            if len(st.session_state.txt) == 0:
                st.info("Use the sidebar to access TMvisDB.")
            else:
                # Print results
                # st.write(table.query(selected_organismid, selected_domain, selected_kingdom, selected_type, selected_sp, selected_length))
                if st.session_state.tbl is None:
                    st.error(
                        "There are no entries in TMvisDB for your selection: topology ("
                        + selected_type
                        + ") and taxonomy (Organism ID: "
                        + selected_organismid
                        + ", Domain: "
                        + selected_domain
                        + ", Kingdom: "
                        + selected_kingdom
                        + "). Please alter your search and try again, or check FAQs if you believe there is something missing.",
                        icon="🚨",
                    )
                else:
                    if (
                        not selected_organismid.isnumeric()
                        and len(selected_organismid) > 0
                    ):
                        print(selected_organismid)
                        st.warning(
                            "The organism ID you entered is not numeric: "
                            + selected_organismid
                            + ". It is not included in your database query. Please check UniProt for the organism ID you are looking for."
                        )

                    st.caption(st.session_state.txt)
                    table.show_tbl(st.session_state.tbl)
                    # Download Button
                    st.download_button(
                        "Download selection",
                        table.convert_df(st.session_state.tbl),
                        "file.csv",
                        "text/csv",
                        key="download-csv",
                    )

                    # Visualize from table
                    st.markdown("---")
                    selected_dfid = st.selectbox(
                        "Choose an ID to visualize predicted transmembrane topology below",
                        st.session_state.tbl["UniProt ID"],
                        0,
                    )
                    protein_info = collect_and_display_protein_info(
                        db_conn, selected_id
                    )
                    visualization.create_visualization_for_id(
                        protein_info, style, color_prot, spin
                    )

    ## 3D vis ##
    with tab3:
        if db_conn is None:
            st.error(
                "Error establishing a connection to TMvisDB! Please try again later, and/or contact us here: tmvisdb@rostlab.org",
                icon="🚨",
            )
        else:
            try:
                protein_info = collect_and_display_protein_info(db_conn, selected_id)
                visualization.create_visualization_for_id(
                    protein_info, style, color_prot, spin
                )
            except Exception as e:
                logging.error(f"An error occurred: {e}")
                st.error(
                    "Sorry, we could not visualize your selected protein. Please contact us, so we can help you with your search.",
                    icon="🚨",
                )
            st.markdown("---")

    ## FAQ ##
    with tab4:
        faq.quest()

    ## About ##
    with tab5:
        about.references()
        about.software()
        about.author()
        about.impr()
