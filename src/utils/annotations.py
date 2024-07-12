# Copyright 2024 Tobias Olenyi.
# SPDX-License-Identifier: Apache-2.0


def construct_topdb_annotation(item: None | dict[str, str | dict[str, str]]):
    topdb_value = item.get("topdb", {}).get("TopDB_Entry", None)
    return None if not topdb_value else list(topdb_value)


def construct_membranome_annotation(item: None | dict[str, str | dict[str, str]]):
    """
    Extracts Membranome data from the given item.
    """
    if "membranomedb" not in item:
        return None

    seq_length = int(item["seq_length"]) - 1
    pos_start = int(item["membranomedb"]["tm_seq_start"]) - 1
    pos_end = int(item["membranomedb"]["tm_seq_end"]) - 1

    membdb_annotation = ["*"] * seq_length
    membdb_annotation[pos_start:pos_end] = ["AH"] * (pos_end - pos_start + 1)

    return membdb_annotation


def construct_tmembd_annotation(item: None | dict[str, str | dict[str, str]]):
    pass


ANNOTATION_HANDLER_MAP = {
    "membranome": construct_membranome_annotation,
    "topdb": construct_topdb_annotation,
    "tembed": construct_tmembd_annotation,
}
