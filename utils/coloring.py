import pandas as pd
from enum import Enum, auto
from dataclasses import dataclass


# Enum for Protein Styles
class ProteinStyle(Enum):
    CARTOON = "cartoon"
    LINE = "line"
    CROSS = "cross"
    STICK = "stick"
    SPHERE = "sphere"


# Enum for Color Schemes
class ColorScheme(Enum):
    TRANSMEMBRANE_PREDICTION = "Transmembrane Prediction"
    ALPHAFOLD_PLDDT_SCORE = "Alphafold pLDDT score"


# Enum for Color Codes
class ColorCode(Enum):
    HELIX_LIGHT = "light green"
    HELIX_DARK = "dark green"
    BETA_LIGHT = "light blue"
    BETA_DARK = "dark blue"
    INSIDE = "light grey"
    OUTSIDE = "dark grey"
    SIGNAL_PEPTIDE = "pink"
    VERY_LOW = "#FF0000"
    LOW = "#FFA500"
    CONFIDENT = "#00C900"
    VERY_HIGH = "#0000FF"


# Dataclass for Style
@dataclass
class Style:
    visualization_style: ProteinStyle = ProteinStyle.CARTOON
    color_scheme: ColorScheme = ColorScheme.TRANSMEMBRANE_PREDICTION
    spin: bool = False


ANNOTATION_LEGEND_DF = pd.DataFrame.from_records(
    [
        {
            "Topology": "Helix",
            "Abbreviation": "H",
            "Orientation": "IN-->OUT",
            "Color": ColorCode.HELIX_LIGHT.value,
        },
        {
            "Topology": "Helix",
            "Abbreviation": "h",
            "Orientation": "OUT-->IN",
            "Color": ColorCode.HELIX_DARK.value,
        },
        {
            "Topology": "Beta-Strand",
            "Abbreviation": "B",
            "Orientation": "IN-->OUT",
            "Color": ColorCode.BETA_LIGHT.value,
        },
        {
            "Topology": "Beta-Strand",
            "Abbreviation": "b",
            "Orientation": "OUT-->IN",
            "Color": ColorCode.BETA_DARK.value,
        },
        {
            "Topology": "inside",
            "Abbreviation": "i",
            "Orientation": "inside",
            "Color": ColorCode.INSIDE.value,
        },
        {
            "Topology": "outside",
            "Abbreviation": "o",
            "Orientation": "outside",
            "Color": ColorCode.OUTSIDE.value,
        },
        {
            "Topology": "Signal Peptide",
            "Abbreviation": "S",
            "Orientation": "NA",
            "Color": ColorCode.SIGNAL_PEPTIDE.value,
        },
    ]
)

ALPHAFOLD_LEGEND_DF = pd.DataFrame(
    [
        "Very low (pLDDT < 50)",
        "Low (70 > pLDDT > 50)",
        "Confident (90 > pLDDT > 70)",
        "Very high (pLDDT > 90)",
    ],
    columns=["pLDDT score"],
)


def color_prediction(annotation_table):
    tmbed_prediction = annotation_table.loc["TMbed Prediction"]
    color_map = {
        "S": ColorCode.SIGNAL_PEPTIDE.value,
        "H": ColorCode.HELIX_LIGHT.value,
        "h": ColorCode.HELIX_DARK.value,
        "B": ColorCode.BETA_LIGHT.value,
        "b": ColorCode.BETA_DARK.value,
        "i": ColorCode.INSIDE.value,
    }
    color = color_map.get(tmbed_prediction, ColorCode.OUTSIDE.value)
    return [f"background-color: {''.join(color.split())}"] * len(annotation_table.index)


def annotation_legend_coloring(value_name):
    return f"background-color: {''.join(value_name.split())}"


def alphafold_legend_coloring(value_name):
    color_map = {
        "Very low (pLDDT < 50)": ColorCode.VERY_LOW.value,
        "Low (70 > pLDDT > 50)": ColorCode.LOW.value,
        "Confident (90 > pLDDT > 70)": ColorCode.CONFIDENT.value,
        "Very high (pLDDT > 90)": ColorCode.VERY_HIGH.value,
    }
    color = color_map.get(value_name)
    return f"background-color: {color}"


def map_annotation_to_color(pred):
    atom_color = dict()
    color_map = {
        "S": ColorCode.SIGNAL_PEPTIDE.value,
        "H": ColorCode.HELIX_LIGHT.value,
        "h": ColorCode.HELIX_DARK.value,
        "B": ColorCode.BETA_LIGHT.value,
        "b": ColorCode.BETA_DARK.value,
        "i": ColorCode.INSIDE.value,
    }
    for nr, res_type in enumerate(pred):
        atom_color[nr] = color_map.get(res_type, ColorCode.OUTSIDE.value)
    return atom_color
