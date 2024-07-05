"""
This module defines the connection to the Glycosmos database.
Glycosylator can query the database using the Glycosmos API and retrieve the condensed IUPAC string of a glycan to make a model from it. 
"""

import requests

# The SPARQL parts are encoded as URL components:
# url_encoded_characters = {
#     ' ': '%20',
#     '!': '%21',
#     '"': '%22',
#     '#': '%23',
#     '$': '%24',
#     '%': '%25',
#     '&': '%26',
#     "'": '%27',
#     '(': '%28',
#     ')': '%29',
#     '*': '%2A',
#     '+': '%2B',
#     ',': '%2C',
#     '-': '%2D',
#     '.': '%2E',
#     '/': '%2F',
#     ':': '%3A',
#     ';': '%3B',
#     '<': '%3C',
#     '=': '%3D',
#     '>': '%3E',
#     '?': '%3F',
#     '@': '%40',
#     '[': '%5B',
#     '\\': '%5C',
#     ']': '%5D',
#     '^': '%5E',
#     '_': '%5F',
#     '`': '%60',
#     '{': '%7B',
#     '|': '%7C',
#     '}': '%7D',
#     '~': '%7E',
#     '\n': '%0A',  # Line feed
#     '\r': '%0D',  # Carriage return
#     '\t': '%09'   # Tab
# }
#

BASE_URL = "https://ts.glycosmos.org/sparql?default-graph-uri=&query="
"""
The base URL of the Glycosmos SPARQL endpoint.
"""

PREFIXES = "PREFIX+glycan%3A+%3Chttp%3A%2F%2Fpurl.jp%2Fbio%2F12%2Fglyco%2Fglycan%23%3E%0D%0APREFIX+dcterms%3A+%3Chttp%3A%2F%2Fpurl.org%2Fdc%2Fterms%2F%3E%0D%0A%0D%0A"
"""
The namespace prefixes
"""

# The original SPARQL query is as follows (for an example glycan G00573XU):
# ```sparql
# PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
# PREFIX dcterms: <http://purl.org/dc/terms/>

# SELECT *
# WHERE {
# ?entry dcterms:identifier "G00573XU".
# ?entry glycan:has_glycosequence ?seqval.
# ?seqval glycan:has_sequence ?seq.
# FILTER (CONTAINS(STR(?seqval), "iupac_condensed"))
# }
# ```

QUERY_FOR_IUPAC = "SELECT+%3Fseq%0D%0AWHERE+%7B%0D%0A%3Fentry+dcterms%3Aidentifier+%22{}%22.%0D%0A%3Fentry+glycan%3Ahas_glycosequence+%3Fseqval.%0D%0A%3Fseqval+glycan%3Ahas_sequence+%3Fseq.%0D%0AFILTER+%28CONTAINS%28STR%28%3Fseqval%29%2C+%22iupac_condensed%22%29%29%0D%0A%7D+&format=application%2Fsparql-results%2Bjson&timeout=0"
"""
The SPARQL query to retrieve the IUPAC string of a glycan.
This string has one placeholder for the GlyCosmos/GlyTouCan ID of the glycan.
"""

# The original SPARQL query is as follows (for the same example glycan G00573XU -> Gal(a1-4)Gal(b1-):
# ```sparql
# PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
# PREFIX dcterms: <http://purl.org/dc/terms/>

# SELECT ?entryID
# WHERE {
# ?entry dcterms:identifier ?entryID.
# ?entry glycan:has_glycosequence ?seqval.
# ?seqval glycan:has_sequence ?seq.
# FILTER (CONTAINS(STR(?seqval), "iupac_condensed") AND STR(?seq) = "Gal(a1-4)Gal(b1-")
# }
# ```

QUERY_FOR_FULL_MATCH_ID = "https://ts.glycosmos.org/sparql?default-graph-uri=&query=PREFIX+glycan%3A+%3Chttp%3A%2F%2Fpurl.jp%2Fbio%2F12%2Fglyco%2Fglycan%23%3E%0D%0APREFIX+dcterms%3A+%3Chttp%3A%2F%2Fpurl.org%2Fdc%2Fterms%2F%3E%0D%0A%0D%0ASELECT+%3FentryID%0D%0AWHERE+%7B%0D%0A%3Fentry+dcterms%3Aidentifier+%3FentryID.%0D%0A%3Fentry+glycan%3Ahas_glycosequence+%3Fseqval.%0D%0A%3Fseqval+glycan%3Ahas_sequence+%3Fseq.%0D%0AFILTER+%28CONTAINS%28STR%28%3Fseqval%29%2C+%22iupac_condensed%22%29+AND+STR%28%3Fseq%29+%3D+%22{}%22%29%0D%0A%7D+&format=application%2Fsparql-results%2Bjson&timeout=0"
"""
The SPARQL query to retrieve the GlyCosmos/GlyTouCan ID of a glycan.
This string has one placeholder for the IUPAC string of the glycan.
"""

# The original SPARQL query is as follows (for the same example):
# ```sparql
# PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
# PREFIX dcterms: <http://purl.org/dc/terms/>

# SELECT ?entryID
# WHERE {
# ?entry dcterms:identifier ?entryID.
# ?entry glycan:has_glycosequence ?seqval.
# ?seqval glycan:has_sequence ?seq.
# FILTER (CONTAINS(STR(?seqval), "iupac_condensed") AND CONTAINS(STR(?seq), "Gal(a1-4)Gal(b1-"))
# }
# ```

QUERY_FOR_PARTIAL_MATCH_ID = "SELECT+%3FentryID+%3Fseq%0D%0AWHERE+%7B%0D%0A%3Fentry+dcterms%3Aidentifier+%3FentryID.%0D%0A%3Fentry+glycan%3Ahas_glycosequence+%3Fseqval.%0D%0A%3Fseqval+glycan%3Ahas_sequence+%3Fseq.%0D%0AFILTER+%28CONTAINS%28STR%28%3Fseqval%29%2C+%22iupac_condensed%22%29+AND+CONTAINS%28STR%28%3Fseq%29%2C+%22{}%22%29%29%0D%0A%7D+&format=application%2Fsparql-results%2Bjson&timeout=0"
"""
The SPARQL query to retrieve GlyCosmos/GlyTouCan IDs of glycans that are partial matches of the query glycan.
This string has one placeholder for the IUPAC string of the glycan. 
"""


__all__ = [
    "get_iupac_from_glycosmos",
    "get_glytoucan_id_from_iupac",
    "find_glytoucan_ids_from_iupac",
]


def get_iupac_from_glycosmos(id: str) -> str:
    """
    Get the IUPAC string of a glycan from the Glycosmos database.

    Parameters
    ----------
    id : str
        The GlyCosmos/GlyTouCan ID of the glycan.

    Returns
    -------
    str
        The IUPAC string of the glycan.
    """
    url = make_query_url(QUERY_FOR_IUPAC, id)
    response = query_glycosmos(url)
    bindings = response["results"]["bindings"]
    if len(bindings) == 0:
        return None
    return bindings[0]["seq"]["value"]


def get_glytoucan_id_from_iupac(iupac: str) -> str:
    """
    Get the GlyTouCan ID from the IUPAC condensed format.

    Parameters
    ----------
    iupac : str
        The IUPAC condensed format.

    Returns
    -------
    str
        The GlyTouCan ID.
    """
    url = make_query_url(QUERY_FOR_FULL_MATCH_ID, iupac)
    response = query_glycosmos(url)
    bindings = response["results"]["bindings"]
    if len(bindings) == 0:
        return None
    return bindings[0]["entryID"]["value"]


def find_glytoucan_ids_from_iupac(iupac: str) -> list:
    """
    Find GlyTouCan IDs for glycans that are partial matches of the query glycan.

    Parameters
    ----------
    iupac : str
        The IUPAC condensed format.

    Returns
    -------
    list
        The GlyTouCan IDs and IUPAC strings of the glycans.
    """
    url = make_query_url(QUERY_FOR_PARTIAL_MATCH_ID, iupac)
    response = query_glycosmos(url)
    bindings = response["results"]["bindings"]
    if len(bindings) == 0:
        return None
    return [(b["entryID"]["value"], b["seq"]["value"]) for b in bindings]


def make_query_url(template: str, q: str) -> str:
    """
    Make a query URL from a template and a string.

    Parameters
    ----------
    template : str
        The template of the query URL.
    q : str
        The query string.

    Returns
    -------
    str
        The query URL.
    """
    return BASE_URL + PREFIXES + template.format(q)


def query_glycosmos(url: str) -> dict:
    """
    Query the Glycosmos database.

    Parameters
    ----------
    url : str
        The query URL.

    Returns
    -------
    dict
        The JSON response.
    """
    response = requests.get(url)
    return response.json()
