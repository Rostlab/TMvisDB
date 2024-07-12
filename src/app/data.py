import logging

import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder

from utils import database
from utils.database import DBFilter


@st.cache_data
def convert_df(df):
    return df.to_csv().encode("utf-8")


def left_align(s, props="text-align: left;"):
    return props


def filter_to_markdown(db_filter: DBFilter):
    parts = []

    # Topology
    if db_filter.topology is not None:
        parts.append(
            f"**topology**: [{db_filter.topology.value} ({' not' if (not db_filter.signal_peptide and db_filter.topology is not Topology.ALL) else ''} including signal peptides )]"
        )

    # Taxonomy
    taxonomy_parts = []
    if db_filter.organism_id is not None:
        taxonomy_parts.append(f"Organism ID: {db_filter.organism_id}")
    if db_filter.domain is not None:
        taxonomy_parts.append(f"Domain: {db_filter.domain.value}")
    if db_filter.kingdom is not None:
        taxonomy_parts.append(f"Kingdom: {db_filter.kingdom.value}")

    if taxonomy_parts:
        parts.append(f"**taxonomy**: [{', '.join(taxonomy_parts)}]")

    # Sequence Lengths
    if db_filter.sequence_lengths != (16, 5500):
        parts.append(
            f"**lengths**: [{db_filter.sequence_lengths[0]}-{db_filter.sequence_lengths[1]}]"
        )
    else:
        parts.append("**lengths**: [all]")

    # Final String
    return ", ".join(parts)


def show_table(df: pd.DataFrame):

    builder = GridOptionsBuilder.from_dataframe(df, columnwidth=3)
    builder.configure_pagination(
        enabled=True, paginationAutoPageSize=False, paginationPageSize=25
    )
    builder.configure_grid_options(enableCellTextSelection=True)
    #    builder.configure_default_column(cellStyle: { textAlign: 'left' } )
    go = builder.build()
    AgGrid(
        df, gridOptions=go, fit_columns_on_grid_load=True
    )  # , update_mode='manual')#theme='alpine',


@st.cache_data(ttl=60, show_spinner=False)
def display_random_data(_db_conn, db_filter: DBFilter):
    with st.spinner("Loading random data..."):
        st.session_state.data = database.get_filtered_data(_db_conn, db_filter)
    st.session_state.user_display = "The table below shows a random selection. You can retrieve new random data every minute. Use the sidebar filters for a personalized selection."  # noqa: E501


@st.cache_data(ttl=600, show_spinner=False)
def display_filtered_data(_db_conn, db_filter: DBFilter):
    with st.spinner("Loading filtered data..."):
        st.session_state.data = database.get_filtered_data(_db_conn, db_filter)
        st.session_state.user_display = f"The table below shows your personalized selection -  {filter_to_markdown(db_filter)}. For a random selection use the sidebar button."  # noqa: E501
