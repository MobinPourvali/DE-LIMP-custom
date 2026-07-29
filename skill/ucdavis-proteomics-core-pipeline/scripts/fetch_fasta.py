#!/usr/bin/env python3
"""
fetch_fasta.py  --  Resolve the organism, then build the search FASTA.

Two modes:

  resolve  Organism name (or taxid) -> UniProt proteome candidates. The organism
           is NEVER inferred from a workflow bundle; the bundle's proteome is only
           a default to confirm. Always ask the user and resolve here.

             python3 fetch_fasta.py resolve --organism "mouse"
             python3 fetch_fasta.py resolve --taxid 10090

  fetch    Build the FASTA for a confirmed proteome ID.

             python3 fetch_fasta.py fetch --proteome UP000005640 \
                 --content one_per_gene --contaminants universal \
                 --out ./search.fasta [--hive]

Proteome resolution priority (cheapest / most-trusted first):
  1. --path override            -> used verbatim (e.g. a pre-staged proteome).
  2. UC Davis HIVE (--hive)     -> reuse /quobyte/proteomics-grp/MRS/ instead of
                                   downloading.
  3. UniProt.

CONTENT TYPE (why the default is one_per_gene)
  UniProt's REST `&onePerGene=true` is SILENTLY IGNORED -- verified 2026-07-29:
  the yeast stream returns byte-identical output (3,855,392 bytes, 6,067 entries)
  with and without it. So a REST `(proteome:X)` stream is always the FULL set:
  human UP000005640 = 147,506 sequences vs 20,652 canonical one-per-gene, a 7.1x
  larger search space. The canonical set only comes from the reference-proteome
  FTP tree, which is what DE-LIMP (R/helpers_search.R) does and what the Core has
  staged on HIVE (UP000005640_9606.fasta = 20,663). Hence: one_per_gene default,
  served from FTP, with a loud REST fallback if FTP has no file for the proteome.

CONTAMINANTS
  Sets are Cont_-tagged, which is what DIA-NN's --cont-quant-exclude Cont_ keys on
  (DIA-NN README "Contaminants"; the Linux binary ships no contaminant FASTA of its
  own, so it must be appended here). The Core's own staged database is
  proteome + Universal contaminants, so `universal` is the default.

Emits JSON on stdout and writes a `<out>.meta.json` sidecar with the same content
plus checksums, for the reproducibility bundle.
"""
import sys, os, re, json, argparse, glob, gzip, shutil, hashlib
import urllib.request, urllib.error, urllib.parse

HIVE_MRS = "/quobyte/proteomics-grp/MRS"
UNIPROT_REST = "https://rest.uniprot.org"
FTP_REF = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release"
           "/knowledgebase/reference_proteomes")

# Hao lab contaminant libraries (public, JPR 2022, doi:10.1021/acs.jproteome.2c00145).
# These are the same files DE-LIMP bundles in contaminants/ and the Core stages on
# HIVE. Headers are Cont_-tagged -> compatible with DIA-NN --cont-quant-exclude.
CONTAM_REPO = ("https://raw.githubusercontent.com/HaoGroup-ProtContLib"
               "/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/main")
CONTAM_SETS = {
    "universal": {
        "url_path": "Universal protein contaminant FASTA/0602_Universal Contaminants.fasta",
        "delimp": "Universal_Contaminants.fasta",
        "hive_tokens": ("universal", "contamin"),
        "desc": "Universal contaminants -- the default; what the Core stages on HIVE",
    },
    "cell_culture": {
        "url_path": "Sample-type specific contaminant FASTA/0602_Cell Culture Contaminants.fasta",
        "delimp": "Cell_Culture_Contaminants.fasta",
        "hive_tokens": ("cell culture", "contamin"),
        "desc": "Cell line / cell culture samples",
    },
    "mouse_tissue": {
        "url_path": "Sample-type specific contaminant FASTA/Aug2022_Mouse Tissue Contaminants.fasta",
        "delimp": "Mouse_Tissue_Contaminants.fasta",
        "hive_tokens": ("mouse tissue", "contamin"),
        "desc": "Mouse tissue samples",
    },
    "rat_tissue": {
        "url_path": "Sample-type specific contaminant FASTA/Aug2022_Rat Tissue Contaminants.fasta",
        "delimp": "Rat_Tissue_Contaminants.fasta",
        "hive_tokens": ("rat tissue", "contamin"),
        "desc": "Rat tissue samples",
    },
    "neuron_culture": {
        "url_path": "Sample-type specific contaminant FASTA/Dec2022_Neuron Culture Contaminants.fasta",
        "delimp": "Neuron_Culture_Contaminants.fasta",
        "hive_tokens": ("neuron", "contamin"),
        "desc": "Neuronal culture samples",
    },
    "stem_cell_culture": {
        "url_path": "Sample-type specific contaminant FASTA/2026_Stem Cell Culture Contaminants.fasta",
        "delimp": "Stem_Cell_Culture_Contaminants.fasta",
        "hive_tokens": ("stem cell", "contamin"),
        "desc": "Stem cell culture samples",
    },
}
CONTAM_CITATION = ("Frankenfield AM, Ni J, Ahmed M, Hao L. (2022) J Proteome Res "
                   "21(9):2104-2113. doi:10.1021/acs.jproteome.2c00145")
CONT_TAG = "Cont_"

KINGDOM_DIR = {"eukaryota": "Eukaryota", "bacteria": "Bacteria",
               "archaea": "Archaea", "viruses": "Viruses"}

# Curated organisms a core facility actually sees: taxid -> (name, offline UP, aliases).
# We resolve the proteome LIVE from the taxid rather than trusting the accession --
# reference proteomes get superseded, so a pinned UP silently goes stale, whereas an
# NCBI taxid is stable. The accession is kept only as an offline fallback and as a
# staleness cross-check. This table also fixes free-text searches that the proteomes
# endpoint answers badly: "Escherichia coli K-12" returns five Non-Reference MG1655
# assemblies and never surfaces the real reference UP000000625 at all.
ORGANISM_TAXIDS = {
    9606:  ("Homo sapiens",              "UP000005640", ["human", "hsapiens", "h. sapiens", "hs"]),
    10090: ("Mus musculus",              "UP000000589", ["mouse", "mice", "murine", "mmusculus", "mm"]),
    10116: ("Rattus norvegicus",         "UP000002494", ["rat", "rnorvegicus"]),
    559292:("Saccharomyces cerevisiae",  "UP000002311", ["yeast", "s. cerevisiae", "scerevisiae",
                                                         "budding yeast", "baker's yeast"]),
    83333: ("Escherichia coli K-12",     "UP000000625", ["e. coli", "ecoli", "e coli",
                                                         "escherichia coli", "escherichia coli k-12"]),
    7227:  ("Drosophila melanogaster",   "UP000000803", ["fly", "fruit fly", "drosophila"]),
    7955:  ("Danio rerio",               "UP000000437", ["zebrafish", "danio"]),
    3702:  ("Arabidopsis thaliana",      "UP000006548", ["arabidopsis", "thale cress"]),
    6239:  ("Caenorhabditis elegans",    "UP000001940", ["c. elegans", "celegans", "worm", "nematode"]),
    9913:  ("Bos taurus",                "UP000009136", ["cow", "bovine", "cattle"]),
    9823:  ("Sus scrofa",                "UP000008227", ["pig", "porcine", "swine"]),
    9031:  ("Gallus gallus",             "UP000000539", ["chicken", "chick"]),
    9615:  ("Canis lupus familiaris",    "UP000002254", ["dog", "canine"]),
    9541:  ("Macaca fascicularis",       "UP000233100", ["cynomolgus", "cyno", "macaque"]),
    284812:("Schizosaccharomyces pombe", "UP000002485", ["fission yeast", "s. pombe", "pombe"]),
    5691:  ("Trypanosoma brucei",        "UP000008524", ["trypanosoma", "t. brucei"]),
    1773:  ("Mycobacterium tuberculosis","UP000001584", ["mtb", "m. tuberculosis", "tb"]),
    5833:  ("Plasmodium falciparum",     "UP000001450", ["plasmodium", "malaria", "p. falciparum"]),
}

CONTENT_TYPES = ("one_per_gene", "reviewed", "reviewed_isoforms",
                 "full", "full_isoforms")


# --------------------------------------------------------------------------
# HTTP
# --------------------------------------------------------------------------
def _open(url, timeout=300):
    req = urllib.request.Request(url, headers={"User-Agent": "proteomics-pipeline-skill"})
    return urllib.request.urlopen(req, timeout=timeout)


def _get_json(url, timeout=60):
    with _open(url, timeout) as r:
        # Return the header object, not dict(r.headers): HTTP header lookup must be
        # case-insensitive ('X-UniProt-Release' vs 'x-uniprot-release').
        return json.loads(r.read().decode("utf-8")), r.headers


def _download(url, dest, timeout=600):
    with _open(url, timeout) as r, open(dest, "wb") as fh:
        shutil.copyfileobj(r, fh)


def _read_fasta_text(path):
    opn = gzip.open if path.endswith(".gz") else open
    with opn(path, "rt", errors="replace") as fh:
        return fh.read()


def _count(text):
    return sum(1 for ln in text.splitlines() if ln.startswith(">"))


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _warn(msg):
    sys.stderr.write(f"[fetch_fasta] WARNING: {msg}\n")


# --------------------------------------------------------------------------
# resolve: organism -> proteome candidates
# --------------------------------------------------------------------------
def resolve_proteomes(organism=None, taxid=None, size=25):
    """Search UniProt proteomes. Mirrors DE-LIMP search_uniprot_proteomes().

    NOTE: DE-LIMP's server-side `AND (proteome_type:1)` filter returns ZERO results
    against the current API (verified 2026-07-29) -- that filter would silently make
    every organism look unavailable. We filter on proteomeType client-side instead.
    """
    if taxid:
        query = f"taxonomy_id:{int(taxid)}"
    elif organism:
        query = organism
    else:
        raise ValueError("resolve needs --organism or --taxid")

    url = (f"{UNIPROT_REST}/proteomes/search?"
           + urllib.parse.urlencode({"query": query, "format": "json", "size": size}))
    data, _ = _get_json(url)

    rows = []
    for r in data.get("results", []):
        tax = r.get("taxonomy") or {}
        ptype = r.get("proteomeType") or ""
        rows.append({
            "proteome_id": r.get("id") or "",
            "organism": tax.get("scientificName") or "",
            "common_name": tax.get("commonName") or "",
            "taxid": int(tax.get("taxonId") or 0),
            "protein_count": int(r.get("proteinCount") or 0),
            "proteome_type": ptype,
            # Exact match, NOT a substring test: UniProt's "Non Reference proteome"
            # contains "reference", so `"reference" in ptype.lower()` ranks strain
            # assemblies as reference proteomes. (DE-LIMP's R picker has this bug --
            # helpers_search.R uses grepl("Reference", ...) -- so "baker's yeast"
            # selects UP000077179 over the real reference UP000002311.)
            "is_reference": ptype.strip().lower() == "reference proteome",
        })

    # Drop UniProt's "Excluded" (redundant/low-quality) proteomes unless they are
    # all we have -- selecting one silently would give a wrong-sized database.
    kept = [r for r in rows if r["proteome_type"].lower() != "excluded"]
    if not kept:
        kept = rows

    # A free-text organism search is a keyword match, so "mouse" also returns
    # Myotis myotis and mouse-ear cress, and "human" returns 20-odd human viruses --
    # all genuine reference proteomes. An exact name hit is what disambiguates.
    # Match the bare species name too: UniProt's scientificName carries a strain
    # qualifier ("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"), so a user
    # typing the plain scientific name never matches exactly and the sort falls through
    # to protein_count -- which put S. pastorianus above the real S. cerevisiae.
    q = _norm(organism) if organism else ""
    for r in kept:
        r["exact_name_match"] = bool(q) and q in (
            _norm(r["common_name"]), _norm(r["organism"]), _strip_strain(r["organism"]))

    kept.sort(key=lambda r: (not r["exact_name_match"], not r["is_reference"],
                             -r["protein_count"]))
    return kept


def _norm(s):
    return " ".join((s or "").strip().lower().split())


def _strip_strain(s):
    """'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)' -> 'saccharomyces cerevisiae'"""
    return _norm(re.sub(r"\s*\([^)]*\)", " ", s or ""))


def alias_taxid(text):
    """Curated name/alias/taxid -> taxid. None if we don't recognise it."""
    t = _norm(text)
    if not t:
        return None
    if t.isdigit():
        return int(t)
    for tx, (name, _up, aliases) in ORGANISM_TAXIDS.items():
        if t == _norm(name) or t in [_norm(x) for x in aliases]:
            return tx
    return None


def cmd_resolve(a):
    organism, taxid = a.organism, a.taxid
    notes = []

    # A UniProt accession typed straight in: take it at face value, but look it up so
    # the user still gets shown the organism they are about to search.
    if organism and re.fullmatch(r"UP\d{9}", organism.strip(), re.I):
        upid = organism.strip().upper()
        try:
            m = proteome_meta(upid)
        except Exception as e:
            sys.exit(f"'{upid}' is not a resolvable UniProt proteome: {e}")
        sel = {"proteome_id": upid, "organism": m["organism"], "common_name": "",
               "taxid": m["taxid"], "protein_count": m["protein_count"],
               "proteome_type": m["proteome_type"],
               "is_reference": m["proteome_type"].strip().lower() == "reference proteome",
               "exact_name_match": True}
        print(json.dumps({"query": {"organism": organism}, "n_candidates": 1,
                          "candidates": [sel], "selected": sel,
                          "needs_menu": False,
                          "notes": ["accession supplied verbatim"]}, indent=2))
        return 0

    # Curated alias -> taxid. Resolving BY TAXID is far more reliable than free text:
    # the proteomes endpoint ranks strain assemblies above the species reference.
    if organism and taxid is None:
        tx = alias_taxid(organism)
        if tx is not None:
            taxid, notes = tx, [f"'{organism}' matched the curated organism table "
                                f"(taxid {tx}); resolved by taxid, not free-text search"]

    try:
        cands = resolve_proteomes(None if taxid else organism, taxid, a.size)
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        # Offline fallback: the curated accession is better than nothing, but flag it —
        # reference proteomes are superseded over time.
        tx = taxid if taxid is not None else alias_taxid(organism or "")
        if tx in ORGANISM_TAXIDS:
            name, up, _ = ORGANISM_TAXIDS[tx]
            sel = {"proteome_id": up, "organism": name, "common_name": "", "taxid": tx,
                   "protein_count": 0, "proteome_type": "", "is_reference": True,
                   "exact_name_match": True}
            print(json.dumps({"query": {"organism": organism, "taxid": taxid},
                              "n_candidates": 1, "candidates": [sel], "selected": sel,
                              "needs_menu": True,
                              "notes": [f"UniProt unreachable ({e}); resolved from the "
                                        f"built-in fallback table. Reference proteome "
                                        f"accessions change over time — VERIFY before "
                                        f"relying on this."]}, indent=2))
            return 0
        sys.exit(f"UniProt proteome search failed: {e}")

    # Staleness cross-check against the curated table (same idea as the offline note).
    tx = taxid if taxid is not None else None
    if cands and tx in ORGANISM_TAXIDS and cands[0]["proteome_id"] != ORGANISM_TAXIDS[tx][1]:
        notes.append(f"curated table lists {ORGANISM_TAXIDS[tx][1]} for taxid {tx}; "
                     f"UniProt currently returns {cands[0]['proteome_id']}. Using UniProt.")

    out = {
        "query": {"organism": a.organism, "taxid": a.taxid},
        "notes": notes,
        "n_candidates": len(cands),
        "candidates": cands[: a.size],
        "selected": cands[0] if cands else None,
        # Auto-pick only when the answer is unambiguous: exactly one reference
        # proteome whose name the user actually named (or, for a taxid lookup,
        # exactly one reference proteome at all). Anything else -> show the menu.
        # "baker's yeast" stays a menu: two S. cerevisiae strains match exactly.
        # Branch on how we actually resolved, not on whether --organism was typed:
        # an alias that mapped to a taxid took the taxid path, where the exact-name
        # flags are all False and counting them would force a needless menu.
        "needs_menu": (not cands or not cands[0]["is_reference"] or (
            sum(1 for c in cands if c["is_reference"]) != 1
            if taxid is not None else
            sum(1 for c in cands if c["exact_name_match"] and c["is_reference"]) != 1)),
    }
    if not cands:
        out["hint"] = ("No proteome matched. Try the scientific name "
                       "(e.g. 'Danio rerio') or an NCBI taxid.")
    print(json.dumps(out, indent=2))
    return 0 if cands else 1


# --------------------------------------------------------------------------
# fetch: proteome -> FASTA text
# --------------------------------------------------------------------------
def proteome_meta(proteome):
    url = f"{UNIPROT_REST}/proteomes/{proteome}?format=json"
    data, headers = _get_json(url)
    tax = data.get("taxonomy") or {}
    return {
        "superkingdom": (data.get("superkingdom") or "").lower(),
        "taxid": int(tax.get("taxonId") or 0),
        "organism": tax.get("scientificName") or "",
        "protein_count": int(data.get("proteinCount") or 0),
        "proteome_type": data.get("proteomeType") or "",
        "uniprot_release": headers.get("x-uniprot-release", ""),
        "uniprot_release_date": headers.get("x-uniprot-release-date", ""),
    }


def download_ftp_one_per_gene(proteome, meta, workdir):
    """Canonical one-per-gene set from the reference-proteome FTP tree.

    URL: {FTP_REF}/{Kingdom}/{UPID}/{UPID}_{TAXID}.fasta.gz
    """
    kingdom = KINGDOM_DIR.get(meta["superkingdom"])
    if not kingdom:
        raise RuntimeError(
            f"cannot map superkingdom '{meta['superkingdom']}' to an FTP directory")
    if not meta["taxid"]:
        raise RuntimeError("UniProt returned no taxonomy ID for this proteome")

    url = f"{FTP_REF}/{kingdom}/{proteome}/{proteome}_{meta['taxid']}.fasta.gz"
    tmp = os.path.join(workdir, f"{proteome}.fasta.gz")
    _download(url, tmp)
    if not os.path.exists(tmp) or os.path.getsize(tmp) < 100:
        raise RuntimeError("FTP download returned an empty file")
    text = _read_fasta_text(tmp)
    os.remove(tmp)
    if _count(text) == 0:
        raise RuntimeError("FTP file decompressed to zero sequences")
    return text, url


def download_rest(proteome, content, workdir):
    """REST stream. Always the FULL set for the given filter -- see module docstring."""
    q = f"(proteome:{proteome})"
    if content in ("reviewed", "reviewed_isoforms"):
        q += " AND (reviewed:true)"
    params = {"query": q, "format": "fasta", "compressed": "false"}
    if content in ("full_isoforms", "reviewed_isoforms"):
        params["includeIsoform"] = "true"
    url = f"{UNIPROT_REST}/uniprotkb/stream?" + urllib.parse.urlencode(params)

    tmp = os.path.join(workdir, f"{proteome}.rest.fasta")
    try:
        _download(url, tmp)
        if not os.path.exists(tmp) or os.path.getsize(tmp) < 100:
            raise RuntimeError("REST download failed or returned an empty file")
        return _read_fasta_text(tmp), url
    finally:
        if os.path.exists(tmp):        # never leave a partial download next to the output
            os.remove(tmp)


def hive_proteome(proteome, content):
    """Pre-staged proteome on HIVE.

    Only matches files whose basename STARTS with the proteome ID, and skips names
    carrying a 'contam'/'plus'/'decoy' suffix -- a loose glob happily returns a
    proteome+contaminants database, which would then get contaminants appended twice.
    """
    if not os.path.isdir(HIVE_MRS):
        return None
    hits = []
    for ext in ("fasta", "fa", "fasta.gz"):
        hits += glob.glob(os.path.join(HIVE_MRS, "**", f"{proteome}*.{ext}"), recursive=True)
    clean = []
    for h in hits:
        low = os.path.basename(h).lower()
        if not low.startswith(proteome.lower()):
            continue
        if any(t in low for t in ("contam", "_plus", "decoy", "predicted")):
            continue
        clean.append(h)
    if not clean:
        return None
    # Shortest basename = the plain proteome file rather than an annotated variant.
    return sorted(clean, key=lambda p: (len(os.path.basename(p)), p))[0]


# --------------------------------------------------------------------------
# contaminants
# --------------------------------------------------------------------------
def _delimp_contaminants_dir():
    """contaminants/ from a DE-LIMP checkout, if we're running inside one."""
    here = os.path.dirname(os.path.abspath(__file__))
    cands = [os.environ.get("DELIMP_APP_DIR", ""),
             os.path.abspath(os.path.join(here, "..", "..", "..")),
             os.path.abspath(os.path.join(here, "..", "..")),
             os.getcwd()]
    for d in cands:
        if d and os.path.isdir(os.path.join(d, "contaminants")):
            return os.path.join(d, "contaminants")
    return None


def hive_contaminants(set_name):
    """Match the REQUESTED set on HIVE, not just any file with 'contaminant' in it."""
    if not os.path.isdir(HIVE_MRS):
        return None
    tokens = CONTAM_SETS[set_name]["hive_tokens"]
    hits = glob.glob(os.path.join(HIVE_MRS, "**", "*.fasta"), recursive=True)
    for h in sorted(hits):
        low = os.path.basename(h).lower()
        if all(t in low for t in tokens):
            return h
    return None


def resolve_contaminants(set_name, explicit_path, use_hive, workdir):
    """-> (text, source). Raises on total failure; the caller decides fatality."""
    if explicit_path:
        if not os.path.exists(explicit_path):
            raise RuntimeError(f"--contaminants-path not found: {explicit_path}")
        return _read_fasta_text(explicit_path), f"path:{explicit_path}"

    if use_hive:
        p = hive_contaminants(set_name)
        if p:
            return _read_fasta_text(p), f"hive:{p}"

    d = _delimp_contaminants_dir()
    if d:
        p = os.path.join(d, CONTAM_SETS[set_name]["delimp"])
        if os.path.exists(p):
            return _read_fasta_text(p), f"delimp:{p}"

    url = f"{CONTAM_REPO}/{urllib.parse.quote(CONTAM_SETS[set_name]['url_path'])}"
    tmp = os.path.join(workdir, "contaminants.fasta")
    _download(url, tmp, timeout=120)
    text = _read_fasta_text(tmp)
    os.remove(tmp)
    if _count(text) == 0:
        raise RuntimeError(f"contaminant download from {url} had zero sequences")
    return text, f"url:{url}"


# --------------------------------------------------------------------------
def cmd_fetch(a):
    if a.content not in CONTENT_TYPES:
        sys.exit(f"--content must be one of {', '.join(CONTENT_TYPES)}")
    if a.contaminants not in ("none",) and a.contaminants not in CONTAM_SETS:
        sys.exit(f"--contaminants must be 'none' or one of {', '.join(CONTAM_SETS)}")

    outdir = os.path.dirname(os.path.abspath(a.out)) or "."
    os.makedirs(outdir, exist_ok=True)

    meta = {}
    warnings = []
    base_text = source = base_url = None
    content_used = a.content

    # 1. explicit override
    if a.path:
        if not os.path.exists(a.path):
            sys.exit(f"--path given but not found: {a.path}")
        base_text = _read_fasta_text(a.path)
        source, content_used = f"override:{a.path}", "unknown"

    # 2. HIVE pre-staged
    if base_text is None and a.hive and a.proteome:
        staged = hive_proteome(a.proteome, a.content)
        if staged:
            base_text = _read_fasta_text(staged)
            source, content_used = f"hive:{staged}", "as_staged"

    # 3. UniProt
    if base_text is None:
        if not a.proteome:
            sys.exit("Need --proteome (or --path). Run `fetch_fasta.py resolve` first.")
        try:
            meta = proteome_meta(a.proteome)
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            sys.exit(f"UniProt proteome lookup failed for {a.proteome}: {e}")

        if a.content == "one_per_gene":
            try:
                base_text, base_url = download_ftp_one_per_gene(a.proteome, meta, outdir)
                source = f"uniprot_ftp:{a.proteome}"
            except Exception as e:
                # Loud, recorded fallback: the user asked for canonical and is not
                # getting it. This CHANGES the search space, so it must never be silent.
                msg = (f"no one-per-gene FTP file for {a.proteome} ({e}). "
                       f"Falling back to the REST full proteome -- this is a LARGER "
                       f"database (all isoforms + unreviewed), not the canonical set.")
                _warn(msg)
                warnings.append(msg)
                try:
                    base_text, base_url = download_rest(a.proteome, "full", outdir)
                except Exception as e2:
                    sys.exit(
                        f"Could not download {a.proteome} from UniProt.\n"
                        f"  one-per-gene (FTP): {e}\n"
                        f"  full proteome (REST): {e2}\n"
                        f"  '{meta.get('proteome_type') or 'unknown'}' proteomes are often "
                        f"not downloadable. Re-run `fetch_fasta.py resolve` and pick a "
                        f"Reference proteome, or pass --path with your own FASTA.")
                source, content_used = f"uniprot_rest:{a.proteome}", "full"
        else:
            try:
                base_text, base_url = download_rest(a.proteome, a.content, outdir)
            except (urllib.error.URLError, urllib.error.HTTPError, RuntimeError) as e:
                sys.exit(f"UniProt download failed for {a.proteome}: {e}")
            source = f"uniprot_rest:{a.proteome}"

    n_base = _count(base_text)
    if n_base == 0:
        sys.exit(f"Resolved FASTA has 0 sequences (source={source}). Refusing to proceed.")

    # Sanity: a truncated stream is the classic silent failure. UniProt's declared
    # proteinCount is the full-set count, so only compare for full-set content.
    declared = meta.get("protein_count") or 0
    if declared and content_used in ("full", "full_isoforms") and n_base < declared * 0.95:
        msg = (f"downloaded {n_base} sequences but UniProt declares {declared} for "
               f"{a.proteome} -- the stream may have been truncated.")
        _warn(msg)
        warnings.append(msg)

    # A supplied/staged database may ALREADY carry contaminants (on HIVE,
    # UP000005640_9606_plus_universal_contam.fasta does). Appending again duplicates
    # every accession, which corrupts protein inference and quant. The --hive filename
    # filter can't see this for --path, so detect it in the sequences themselves.
    n_cont_in_base = sum(1 for ln in base_text.splitlines()
                         if ln.startswith(">") and CONT_TAG in ln)
    if a.contaminants != "none" and n_cont_in_base >= 20:
        msg = (f"the supplied database already contains {n_cont_in_base} "
               f"'{CONT_TAG}'-tagged contaminant sequences (source={source}); "
               f"NOT appending the '{a.contaminants}' set again — appending would "
               f"duplicate every contaminant accession.")
        _warn(msg)
        warnings.append(msg)
        a.contaminants = "none"

    # contaminants
    n_contam, contam_text, contam_source = 0, "", None
    if a.contaminants != "none":
        try:
            contam_text, contam_source = resolve_contaminants(
                a.contaminants, a.contaminants_path, a.hive, outdir)
            n_contam = _count(contam_text)
        except Exception as e:
            # Fail loud by default. A search that silently ran WITHOUT contaminants
            # makes the pipeline's own contaminant-dominance QC check meaningless,
            # and contaminant peptides get misassigned to real proteins instead.
            if not a.allow_missing_contaminants:
                sys.exit(f"Could not obtain the '{a.contaminants}' contaminant set: {e}\n"
                         f"  Fix the source, pass --contaminants-path <file>, or "
                         f"re-run with --contaminants none (recorded as a deliberate\n"
                         f"  choice) / --allow-missing-contaminants to proceed anyway.")
            msg = f"could not fetch contaminants ({e}); proceeding WITHOUT them"
            _warn(msg)
            warnings.append(msg)

    tagged = CONT_TAG in contam_text if contam_text else False
    if contam_text and not tagged:
        msg = (f"contaminant set '{a.contaminants}' has no '{CONT_TAG}' header tag -- "
               f"DIA-NN --cont-quant-exclude {CONT_TAG} will not exclude them from quant.")
        _warn(msg)
        warnings.append(msg)

    with open(a.out, "w") as fh:
        fh.write(base_text)
        if not base_text.endswith("\n"):
            fh.write("\n")
        if contam_text:
            fh.write(contam_text)

    result = {
        "fasta": os.path.abspath(a.out),
        "sha256": _sha256(a.out),
        "source": source,
        "url": base_url,
        "proteome": a.proteome,
        "organism": meta.get("organism", ""),
        "taxid": meta.get("taxid", 0),
        # Methods text must not call a strain assembly or a user-supplied file a
        # "reference proteome" -- record what it actually is.
        "proteome_type": meta.get("proteome_type", ""),
        "content_requested": a.content,
        "content_used": content_used,
        "uniprot_release": meta.get("uniprot_release", ""),
        "uniprot_release_date": meta.get("uniprot_release_date", ""),
        "n_sequences": n_base + n_contam,
        # When the base already carried contaminants we appended none, but the search
        # database still HAS them -- report them so the counts and the DIA-NN flag stay
        # truthful instead of claiming a contaminant-free database.
        "n_proteome": n_base - n_cont_in_base,
        "n_contaminants_appended": n_contam,
        "n_contaminants_already_present": n_cont_in_base,
        "contaminant_set": (a.contaminants if n_contam
                            else ("already_in_supplied_database" if n_cont_in_base else "none")),
        "contaminant_source": contam_source or (source if n_cont_in_base else None),
        "contaminant_citation": CONTAM_CITATION if n_contam else None,
        # Pass this to DIA-NN so contaminants are identified but kept out of
        # quantification/normalisation (DIA-NN README, --cont-quant-exclude).
        "diann_cont_quant_exclude": (CONT_TAG if ((n_contam and tagged) or n_cont_in_base)
                                     else None),
        "warnings": warnings,
    }
    with open(a.out + ".meta.json", "w") as fh:
        json.dump(result, fh, indent=2)
    print(json.dumps(result, indent=2))
    return 0


# --------------------------------------------------------------------------
def main():
    # Back-compat: the old flat form (no subcommand) means `fetch`.
    argv = sys.argv[1:]
    legacy = not argv or argv[0].startswith("-")
    if legacy:
        argv = ["fetch"] + argv
    explicit_contam = any(x == "--contaminants" or x.startswith("--contaminants=")
                          for x in argv)

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("resolve", help="organism name or taxid -> proteome candidates")
    r.add_argument("--organism", help="common name, scientific name, NCBI taxid, or a "
                                      "UP accession — e.g. 'mouse', 'Mus musculus', "
                                      "10090, UP000000589")
    r.add_argument("--taxid", type=int, help="NCBI taxid, e.g. 10090")
    r.add_argument("--list", action="store_true",
                   help="print the curated organism table and exit")
    r.add_argument("--size", type=int, default=25)

    f = sub.add_parser("fetch", help="build the search FASTA")
    f.add_argument("--proteome", help="UniProt proteome ID, e.g. UP000005640")
    f.add_argument("--path", help="explicit FASTA override; used verbatim if set")
    f.add_argument("--content", default="one_per_gene", choices=CONTENT_TYPES,
                   help="default one_per_gene (canonical, via reference-proteome FTP)")
    f.add_argument("--contaminants", default="universal",
                   help="one of: " + ", ".join(list(CONTAM_SETS) + ["none"]))
    f.add_argument("--contaminants-path", help="use this contaminant FASTA verbatim")
    f.add_argument("--add-contaminants", action="store_true",
                   help="(legacy) in the old flat CLI this opted IN to contaminants; "
                        "omitting it there still means none")
    f.add_argument("--allow-missing-contaminants", action="store_true",
                   help="downgrade a contaminant-fetch failure from fatal to a warning")
    f.add_argument("--hive", action="store_true",
                   help="prefer pre-staged HIVE FASTAs (set when env is uc_davis_hive)")
    f.add_argument("--out", required=True)

    a = ap.parse_args(argv)

    # In the old flat CLI, contaminants were opt-IN via --add-contaminants. Honour that:
    # a legacy call that omitted it meant "no contaminants", and must not silently gain
    # 381 sequences now that --contaminants defaults to 'universal'.
    if legacy and not explicit_contam:
        a.contaminants = "universal" if a.add_contaminants else "none"

    if a.cmd == "resolve":
        if a.list:
            print(json.dumps([{"taxid": tx, "organism": n, "offline_fallback_proteome": up,
                               "aliases": al}
                              for tx, (n, up, al) in sorted(ORGANISM_TAXIDS.items())], indent=2))
            return 0
        if not (a.organism or a.taxid):
            ap.error("resolve needs --organism or --taxid (or --list)")
        return cmd_resolve(a)
    return cmd_fetch(a)


if __name__ == "__main__":
    sys.exit(main())
