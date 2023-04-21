"""
Pre-defined sequons for protein glycosylation sites
"""

SEQUONS = {
    "N-linked": "(N)[^P](?=S|T)",
    "O-linked": "(S|T)[^P](?=S|T)",
}
"""
Default sequons for protein glycosylation sites
"""
