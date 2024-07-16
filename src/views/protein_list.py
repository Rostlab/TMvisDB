import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder, JsCode

from utils import database, api
from utils.database import DBFilter
from utils.lineage_definitions import Topology
from utils import protein_info


@st.cache_data
def convert_df(df):
    return df.to_csv().encode("utf-8")


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


def show_table(df: pd.DataFrame, paginate=True):
    if "Organism ID" in df.columns:
        df["Organism ID URL"] = df["Organism ID"].apply(api.uniprot_taxonomy_url)

    js_code = JsCode("""
        class UrlCellRenderer {
          init(params) {
            this.eGui = document.createElement('a');
            this.eGui.innerText = params.value;
            this.eGui.setAttribute('href', params['data']['Organism ID URL']);
            this.eGui.setAttribute('style', "text-decoration:none");
            this.eGui.setAttribute('target', "_blank");
          }
          getGui() {
            return this.eGui;
          }
        }
    """)

    builder = GridOptionsBuilder.from_dataframe(df)
    if paginate:
        builder.configure_pagination(
            enabled=True, paginationAutoPageSize=False, paginationPageSize=25
        )
    builder.configure_grid_options(enableCellTextSelection=True)

    if "Organism ID" in df.columns:
        builder.configure_column(
            "Organism ID",
            cellRenderer=js_code,
            dangerouslyAllowHTML=True,
        )
        builder.configure_column("Organism ID URL", hide=True)

    go = builder.build()
    AgGrid(df, gridOptions=go, fit_columns_on_grid_load=True, allow_unsafe_jscode=True)


@st.cache_data(ttl=60, show_spinner=False)
def display_random_data(db_filter: DBFilter):
    with st.spinner("Loading random data..."):
        query = database.get_sequence_data(db_filter)
        st.session_state.data = protein_info.db_to_df(query)
    st.session_state.user_display = "The table below shows a random selection. You can retrieve new random data every minute. Use the sidebar filters for a personalized selection."  # noqa: E501


@st.cache_data(ttl=600, show_spinner=False)
def display_filtered_data(db_filter: DBFilter):
    with st.spinner("Loading filtered data..."):
        query = database.get_sequence_data(db_filter)
        st.session_state.data = protein_info.db_to_df(query)
    st.session_state.user_display = f"The table below shows your personalized selection -  {filter_to_markdown(db_filter)}. For a random selection use the sidebar button."  # noqa: E501
