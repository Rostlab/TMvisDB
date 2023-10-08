import logging

import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder

from utils import db
from utils.db import DBFilter


@st.cache_data
def convert_df(df):
    return df.to_csv().encode("utf-8")


def left_align(s, props="text-align: left;"):
    return props


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


def display_random_data(db_conn, num_elements):
    logging.debug(num_elements)
    with st.spinner("Loading random data..."):
        st.session_state.tbl = db.get_random_data(db_conn, num_elements)
    st.session_state.txt = "The table below shows a random selection. Use the sidebar filters for a personalized selection."  # noqa: E501
    st.session_state.rndm = False
    st.experimental_rerun()


def display_filtered_data(db_conn, db_filter: DBFilter):
    query_form = db_filter.construct_query()
    with st.spinner("Loading filtered data..."):
        st.session_state.tbl = db.get_filtered_data(
            db_conn, query_form, db_filter.num_sequences
        )
    st.session_state.txt = (
        f"The table below shows your personalized selection: topology ({db_filter.selected_type}),"  # noqa: E501
        f"taxonomy (Organism ID: {db_filter.selected_organismid}, Domain: {db_filter.selected_domain}, Kingdom: {db_filter.selected_kingdom}), "  # noqa: E501
        f"length {str(db_filter.selected_length)}. For a random selection use the sidebar button."  # noqa: E501
    )
    st.session_state.filt = False
    st.experimental_rerun()
