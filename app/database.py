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


@st.cache_data(ttl=60, show_spinner=False)
def display_random_data(_db_conn, db_filter: DBFilter):
    with st.spinner("Loading random data..."):
        st.session_state.data = db.get_random_data(_db_conn, db_filter.num_sequences)
    st.session_state.user_display = "The table below shows a random selection. You can retrieve new random data every minute. Use the sidebar filters for a personalized selection."  # noqa: E501


@st.cache_data(ttl=600, show_spinner=False)
def display_filtered_data(_db_conn, db_filter: DBFilter):
    query_form = db_filter.construct_query()
    with st.spinner("Loading filtered data..."):
        st.session_state.data = db.get_filtered_data(
            _db_conn, query_form, db_filter.num_sequences
        )
        st.session_state.user_display = f"The table below shows your personalized selection: {str(db_filter)}. For a random selection use the sidebar button."  # noqa: E501
