import httpx
import re
from enum import Enum
import logging
from dataclasses import dataclass

from utils.membrane_annotation import ResidueAnnotation


def _fetch_api_data(url):
    """
    Fetches data from a remote API and determines the response type based on the Content-Type header.

    Args:
        url (str): The URL to fetch data from.

    Returns:
        The response data in the appropriate format (JSON or text), or None if an error occurs.
    """
    try:
        response = httpx.get(url)
        response.raise_for_status()

        content_type = response.headers.get("Content-Type", "")
        if "application/json" in content_type:
            return response.json()
        elif "text/plain" in content_type or "application/octet-stream" in content_type:
            return response.text
        else:
            logging.error(f"Unsupported Content-Type: {content_type}")
            return None
    except Exception as e:
        logging.error(f"Error contacting {url}: \n\t{e}")
        return None


## check if ID input format is correct
ACCESSION_NUMBER_RE = re.compile(
    "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}$"
)  # noqa: E501
UNIPROT_ID_RE = re.compile("^[A-Z0-9]{3,20}_[A-Z0-9]{3,20}$")


class UniprotACCType(Enum):
    UNIPROT_ID = 0
    UNIPROT_NAME = 1
    UNKNOWN = -1


@dataclass
class UniprotResponse:
    accession: str
    name: str
    sequence_length: int
    membrane_annotations: list[ResidueAnnotation]


def uniprot_get_input_type(selected_id):
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


def uniprot_query_url(selected_id, input_type):
    """
    Constructs the URL for querying the UniProt database.
    """
    query_prefix = {
        "uniprot_acc_num": f"accession:{selected_id}",
        "uniprot_id": f"id:{selected_id}",
        "unknown": selected_id,
    }

    return f"https://rest.uniprot.org/uniprotkb/search?query={query_prefix.get(input_type, selected_id)} AND active:true&fields=id,accession,length,ft_transmem&format=json&size=1"


def uniprot_taxonomy_url(taxon_id):
    return f"https://www.uniprot.org/taxonomy/{taxon_id}"


def uniprot_parse_response(body):
    """
    Parses the UniProt API response and extracts relevant information.
    """

    annotations = []
    if body and body.get("results") and len(body["results"]) > 0:
        result = body["results"][0]
        up_acc = result["primaryAccession"]
        up_name = result["uniProtkbId"]
        seq_length = result["sequence"]["length"]

        for entry in result.get("features", []):
            if entry["type"] == "Transmembrane":
                label = "BS" if "Beta" in entry["description"] else "AH"
                pos_start = int(entry["location"]["start"]["value"])
                pos_end = int(entry["location"]["end"]["value"])
                annotations.append(ResidueAnnotation(pos_start, pos_end, label))
        return UniprotResponse(
            accession=up_acc,
            name=up_name,
            sequence_length=seq_length,
            membrane_annotations=annotations,
        )
    else:
        return None


def uniprot_fetch_annotation(selected_id):
    """
    Fetches the transmembrane vector information for a given ID from the UniProt database.
    Returns relevant information or default values in case of errors.
    """
    input_type = uniprot_get_input_type(selected_id)
    url = uniprot_query_url(selected_id, input_type)
    body = _fetch_api_data(url)
    return uniprot_parse_response(body)


def tmalphafold_query_url(up_name):
    """
    Constructs the URL for querying the TmAlphaFold database.
    """
    return f"https://tmalphafold.ttk.hu/api/tmdet/{up_name}.json"


def tmalphafold_parse_response(body):
    """
    Parses the TmAlphaFold API response and extracts relevant information.
    """
    annotations: list[ResidueAnnotation] = []

    if body and "CHAIN" in body:
        for entry in body["CHAIN"][0]["REGION"]:
            if entry["_attributes"]["type"] == "M":
                annotations.append(
                    ResidueAnnotation(
                        int(entry["_attributes"]["seq_beg"]),
                        int(entry["_attributes"]["seq_end"]),
                        "AH",
                    )
                )
        return annotations
    return annotations if len(annotations) > 0 else None


def tmalphafold_fetch_annotation(up_name):
    """
    Fetches transmembrane annotation for a given protein name from TmAlphaFold.
    Returns the annotation or a default value in case of errors.
    """
    url = tmalphafold_query_url(up_name)
    body = _fetch_api_data(url)
    tmaf_annotations = tmalphafold_parse_response(body)
    return tmaf_annotations


def alphafolddb_fetch_structure(selected_id):
    """
    Fetches the AlphaFold structure for a given ID from the AlphaFold DB API.
    Returns the sequence and the associated PDB file content.
    """
    afdb_api_path = f"https://www.alphafold.ebi.ac.uk/api/prediction/{selected_id}"
    afdb_json = _fetch_api_data(afdb_api_path)

    if not afdb_json:
        return None, None

    try:
        seq = afdb_json[0]["uniprotSequence"]
        afdb_pdb_path = afdb_json[0]["pdbUrl"]
        afdb_file = _fetch_api_data(afdb_pdb_path)
        return seq, afdb_file

    except (KeyError, IndexError) as e:
        logging.error(f"Error processing AlphaFold structure data: {e}")
        return None, None
