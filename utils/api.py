from urllib.request import urlopen
import requests
import re
from enum import Enum
import logging

## check if ID input format is correct
ACCESSION_NUMBER_RE = re.compile(
    "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}$"
)  # noqa: E501
UNIPROT_ID_RE = re.compile("^[A-Z0-9]{3,20}_[A-Z0-9]{3,20}$")


class UniprotACCType(Enum):
    UNIPROT_ID = 0
    UNIPROT_NAME = 1
    UNKNOWN = -1


def check_input_format(selected_id):
    """
    Checks if the input ID matches expected formats: uniprot_id or uniprot_acc_num.
    Returns the type of ID format or 'unknown' if the format is not recognized.
    """
    test_str = selected_id.upper()
    if UNIPROT_ID_RE.match(test_str):
        return UniprotACCType.UNIPROT_ID
    elif ACCESSION_NUMBER_RE.match(test_str):
        return UniprotACCType.UNIPROT_NAME
    else:
        return UniprotACCType.UNKNOWN


def get_af_structure(selected_id):
    """
    Fetches the AlphaFold structure for a given ID from the AlphaFold DB API.
    Returns the sequence and the associated PDB file content.
    """
    afdb_api_path = f"https://www.alphafold.ebi.ac.uk/api/prediction/{selected_id}"
    try:
        afdb_json = requests.get(afdb_api_path).json()
        seq = afdb_json[0]["uniprotSequence"]
        afdb_pdb_path = afdb_json[0]["pdbUrl"]
        afdb_file = urlopen(afdb_pdb_path).read().decode("utf-8")
        return seq, afdb_file
    except Exception as e:
        logging.error(f"Error fetching AlphaFold structure: {e}")
        return None, None


def get_uniprot_tmvec(selected_id, input_type):
    """
    Fetches the transmembrane vector information for a given ID from the UniProt database.
    Returns relevant information or default values in case of errors.
    """
    # Construct URL based on the input type
    query_prefix = {
        "uniprot_acc_num": f"accession:{selected_id}",
        "uniprot_id": f"id:{selected_id}",
        "unknown": selected_id,
    }

    url = f"https://rest.uniprot.org/uniprotkb/search?query={query_prefix.get(input_type, selected_id)} AND active:true &fields=id,accession,length,ft_transmem&format=json&size=1"  # noqa: E501

    try:
        body = requests.get(url).json()
    except Exception as e:
        logging.error(f"Error contacting {url}: \n\t{e}")
        return None, None, None, 0

    if body.get("results") and len(body["results"]) > 0:
        body = body["results"][0]
        up_acc = body["primaryAccession"]
        up_name = body["uniProtkbId"]
        seq_length = body["sequence"]["length"]
        UP_TM_vec = ["*"] * seq_length

        for entry in body.get("features", []):
            if entry["type"] == "Transmembrane":
                if "Beta" in entry["description"]:
                    annotation = "BS"
                elif "Helical" in entry["description"]:
                    annotation = "AH"
                pos_start = int(entry["location"]["start"]["value"]) - 1
                pos_end = int(entry["location"]["end"]["value"])
                UP_TM_vec[pos_start:pos_end] = [annotation] * (pos_end - pos_start)
        return up_acc, up_name, UP_TM_vec, seq_length
    else:
        return None, None, None, 0


def get_tmalphafold_annotation(up_name, seq_length):
    """
    Fetches transmembrane annotation for a given protein name from TmAlphaFold.
    Returns the annotation or a default value in case of errors.
    """
    url = f"https://tmalphafold.ttk.hu/api/tmdet/{up_name}.json"
    try:
        body = requests.get(url).json()
    except Exception as e:
        logging.error(f"Error contacting {url}: \n\t{e}")
        return None

    tmaf_tm_vec = ["*"] * seq_length

    if "CHAIN" in body:
        for entry in body["CHAIN"][0]["REGION"]:
            if entry["_attributes"]["type"] == "M":
                pos_start = int(entry["_attributes"]["seq_beg"]) - 1
                pos_end = int(entry["_attributes"]["seq_end"])
                tmaf_tm_vec[pos_start:pos_end] = ["AH"] * (pos_end - pos_start)
        return tmaf_tm_vec
    else:
        return None
