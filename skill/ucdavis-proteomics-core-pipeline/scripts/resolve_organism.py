#!/usr/bin/env python3
"""
resolve_organism.py -- turn what the user SAYS into the FASTA you actually search.

THE GAP THIS FILLS
------------------
Organism cannot be detected from a raw file, so the skill asks for it. But asking is
only half the job: the search needs a UniProt *proteome accession* (UP000000589), not
the word "mouse". Until now that accession only ever came hardcoded out of a workflow
bundle, so an organism with no bundle -- or an organism-agnostic bundle -- had nowhere
to get it, and the run either stopped or silently searched the wrong species. Searching
mouse data against a human FASTA does not error; it just quietly loses most of the
proteome, which is the worst possible failure mode.

Accepts a common name, a Latin name, a taxid, or a UP accession, and returns the
proteome to search. Built-in table for the organisms a core facility actually sees;
UniProt REST lookup for everything else.

Usage
  python3 resolve_organism.py --organism "mouse"
  python3 resolve_organism.py --organism 10090
  python3 resolve_organism.py --organism "Mus musculus" --json
  python3 resolve_organism.py --list
"""
import argparse, json, re, sys, urllib.parse, urllib.request

# taxid -> (canonical name, reference proteome, common aliases)
TABLE = {
    9606:  ("Homo sapiens",             "UP000005640", ["human", "hsapiens", "h. sapiens", "hs"]),
    10090: ("Mus musculus",             "UP000000589", ["mouse", "mice", "murine", "mmusculus", "mm"]),
    10116: ("Rattus norvegicus",        "UP000002494", ["rat", "rnorvegicus"]),
    559292:("Saccharomyces cerevisiae", "UP000002311", ["yeast", "s. cerevisiae", "scerevisiae", "budding yeast"]),
    83333: ("Escherichia coli K-12",    "UP000000625", ["e. coli", "ecoli", "e coli"]),
    7227:  ("Drosophila melanogaster",  "UP000000803", ["fly", "fruit fly", "drosophila"]),
    7955:  ("Danio rerio",              "UP000000437", ["zebrafish", "danio"]),
    3702:  ("Arabidopsis thaliana",     "UP000006548", ["arabidopsis", "thale cress"]),
    6239:  ("Caenorhabditis elegans",   "UP000001940", ["c. elegans", "celegans", "worm", "nematode"]),
    9913:  ("Bos taurus",               "UP000009136", ["cow", "bovine", "cattle"]),
    9823:  ("Sus scrofa",               "UP000008227", ["pig", "porcine", "swine"]),
    9031:  ("Gallus gallus",            "UP000000539", ["chicken", "chick"]),
    9615:  ("Canis lupus familiaris",   "UP000002254", ["dog", "canine"]),
    9541:  ("Macaca fascicularis",      "UP000233100", ["cynomolgus", "cyno", "macaque"]),
    284812:("Schizosaccharomyces pombe","UP000002485", ["fission yeast", "s. pombe", "pombe"]),
    5691:  ("Trypanosoma brucei",       "UP000008524", ["trypanosoma", "t. brucei"]),
    1773:  ("Mycobacterium tuberculosis","UP000001584", ["mtb", "m. tuberculosis", "tb"]),
    5833:  ("Plasmodium falciparum",    "UP000001450", ["plasmodium", "malaria", "p. falciparum"]),
}

# NOTE: no `fields=` here. The proteomes endpoint rejects fields=proteome_type
# ("Invalid fields parameter value"), and the full record already carries what we need.
UNIPROT_PROTEOME_SEARCH = (
    "https://rest.uniprot.org/proteomes/search"
    "?query=organism_id:{taxid}&format=json&size=25"
)
UNIPROT_TAXON_SEARCH = (
    "https://rest.uniprot.org/taxonomy/search?query={q}&format=json&size=25"
)


def _get(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "ucdavis-proteomics-pipeline"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def from_table(text):
    t = text.strip().lower()
    if t.isdigit() and int(t) in TABLE:
        tx = int(t)
        name, up, _ = TABLE[tx]
        return {"taxid": tx, "organism": name, "uniprot_proteome": up, "source": "builtin"}
    for tx, (name, up, aliases) in TABLE.items():
        if t == name.lower() or t in aliases:
            return {"taxid": tx, "organism": name, "uniprot_proteome": up, "source": "builtin"}
    return None


def _pick_taxon(hits, query):
    """Choose the right taxon. NEVER just take hits[0].

    A search for "Mus musculus" returns "Mus musculus musculus x Mus musculus
    molossinus" (taxid 3004188) ahead of plain Mus musculus (10090). Taking the first
    hit therefore searches a hybrid subspecies proteome and silently loses most of the
    identifications. Prefer an exact scientific-name match, then a plain species-rank
    hit, and only then fall back to first."""
    q = query.strip().lower()
    exact = [h for h in hits if (h.get("scientificName") or "").strip().lower() == q]
    if exact:
        return exact[0]
    common = [h for h in hits
              if any((c or "").strip().lower() == q for c in (h.get("commonName"), h.get("mnemonic")))]
    if common:
        return common[0]
    species = [h for h in hits if (h.get("rank") or "").lower() == "species"
               and " x " not in (h.get("scientificName") or "").lower()]
    if species:
        # shortest name = least sub-specific
        return sorted(species, key=lambda h: len(h.get("scientificName") or ""))[0]
    return hits[0]


def from_uniprot(text):
    """Resolve anything not in the table by asking UniProt directly."""
    t = text.strip()
    taxid = int(t) if t.isdigit() else None
    organism = None
    ambiguous = None
    if taxid is None:
        try:
            res = _get(UNIPROT_TAXON_SEARCH.format(q=urllib.parse.quote(t)))
        except Exception as e:
            return {"error": f"UniProt taxonomy lookup failed for {t!r}: {e}"}
        hits = res.get("results") or []
        if not hits:
            return {"error": f"no UniProt taxonomy match for {t!r}"}
        chosen = _pick_taxon(hits, t)
        taxid = int(chosen["taxonId"])
        organism = chosen.get("scientificName")
        if len(hits) > 1 and int(hits[0]["taxonId"]) != taxid:
            ambiguous = (f"UniProt's top taxonomy hit was {hits[0].get('scientificName')!r} "
                         f"(taxid {hits[0]['taxonId']}); selected {organism!r} "
                         f"(taxid {taxid}) as the better match for {t!r}. Confirm.")
    try:
        res = _get(UNIPROT_PROTEOME_SEARCH.format(taxid=taxid))
    except Exception as e:
        return {"error": f"UniProt proteome lookup failed for taxid {taxid}: {e}"}
    hits = res.get("results") or []
    if not hits:
        return {"error": f"no proteome in UniProt for taxid {taxid}",
                "taxid": taxid, "organism": organism}

    def rank(h):
        pt = (h.get("proteomeType") or "").lower()
        if "reference" in pt and "representative" not in pt:
            return 0
        if "representative" in pt:
            return 1
        if "other" in pt:
            return 3
        return 2
    hits.sort(key=rank)
    h = hits[0]
    out = {"taxid": taxid,
           "organism": organism or (h.get("taxonomy") or {}).get("scientificName"),
           "uniprot_proteome": h["id"],
           "proteome_type": h.get("proteomeType"),
           "source": "uniprot",
           "n_candidates": len(hits),
           "alternatives": [x["id"] for x in hits[1:4]]}
    if ambiguous:
        out["ambiguous_name"] = ambiguous
    return out


def resolve(text, offline=False):
    """UniProt is the source of truth. The built-in table is only a fallback for when
    the network is unavailable — reference proteomes get superseded (an organism's
    current UP accession is not a constant), so a hardcoded value can silently go
    stale. Anything served from the table is flagged so the caller can say so."""
    t = text.strip()
    # already an accession? take it at face value
    if re.fullmatch(r"UP\d{9}", t.upper()):
        return {"taxid": None, "organism": None, "uniprot_proteome": t.upper(),
                "source": "given-verbatim"}

    if not offline:
        live = from_uniprot(t)
        if "error" not in live:
            cached = from_table(t)
            if cached and cached["uniprot_proteome"] != live["uniprot_proteome"]:
                live["note_offline_cache_differs"] = (
                    f"offline fallback table says {cached['uniprot_proteome']}; "
                    f"UniProt currently says {live['uniprot_proteome']}. Using UniProt.")
            return live

    cached = from_table(t)
    if cached:
        cached["source"] = "offline-fallback" if not offline else "offline-table"
        cached["warning"] = ("Resolved from the built-in fallback table, not UniProt. "
                             "Reference proteome accessions change over time — verify "
                             "before relying on this for a published method.")
        return cached
    return {"error": f"could not resolve {t!r}: UniProt unreachable and no offline entry"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", help="common name, Latin name, taxid, or UPxxxxxxxxx")
    ap.add_argument("--list", action="store_true", help="show the built-in table")
    ap.add_argument("--offline", action="store_true",
                    help="skip UniProt and use the built-in fallback table only")
    a = ap.parse_args()

    if a.list:
        print(json.dumps([{"taxid": tx, "organism": n, "uniprot_proteome": up,
                           "aliases": al} for tx, (n, up, al) in sorted(TABLE.items())], indent=2))
        return
    if not a.organism:
        sys.exit("Need --organism (or --list). The organism cannot be detected from raw "
                 "files — ASK the user, then resolve it here.")

    out = resolve(a.organism, offline=a.offline)
    out["query"] = a.organism
    if "error" not in out:
        out["confirm_with_user"] = (
            f"Searching {out.get('organism') or 'this organism'} "
            f"(taxid {out.get('taxid')}) against UniProt {out['uniprot_proteome']}. "
            "Confirm before the search starts — a wrong proteome does not error, it "
            "silently loses most of the identifications.")
    print(json.dumps(out, indent=2))
    if "error" in out:
        sys.exit(2)


if __name__ == "__main__":
    main()
