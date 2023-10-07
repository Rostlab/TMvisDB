import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder


@st.cache
def convert_df(df):
    return df.to_csv().encode("utf-8")


def left_align(s, props="text-align: left;"):
    return props


def show_tbl(df: pd.DataFrame):
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
