"""Protein-level annotation for guide design: domains, cross-species conservation,
and alternative start codons.  Standard library only.

Three questions this answers for every guide:

  1. Which annotated domain does the cut site sit in?           -> UniProt features
  2. Is that part of the protein conserved to the human gene?   -> global alignment
  3. Could an alternative start codon rescue the protein?       -> in-frame ATG scan

Domains come from UniProt for the zebrafish protein itself, plus the human
ortholog's features projected onto zebrafish residue numbers through a global
BLOSUM62 alignment (``protein_align``).  A handful of domains that matter for
design but that UniProt does not annotate are curated in ``EXTRA`` below, in
human coordinates, with their source; they are projected the same way and are
always labelled as projected in the output so they are never confused with a
direct annotation.
"""
from __future__ import annotations

import json
import os
import time
import urllib.parse
import urllib.request

import protein_align as pa

UA = "chopchop-pipeline"
UNIPROT = "https://rest.uniprot.org/uniprotkb"

# Feature types worth showing on a guide row. UniProt has many more; these are the
# ones that can change a design decision.
KEEP = ("Domain", "Transmembrane", "Topological domain", "Region", "Motif",
        "Short sequence motif", "Coiled coil", "Zinc finger", "Repeat",
        "DNA binding")
# Single-residue features, shown as sites rather than spans.
SITES = ("Modified residue", "Active site", "Binding site")

# Literature-curated domains missing from UniProt, in HUMAN residue coordinates.
# (name, start, end, source, why it matters for guide design)
EXTRA: dict[str, list[tuple]] = {
    "YAP1": [
        ("TEAD-binding domain", 50, 171,
         "Li 2010 Genes Dev; Tian 2010 PNAS (YAP-TEAD co-crystal)",
         "Required for TEAD-driven co-activation, so a frameshift at or before "
         "this point abolishes the canonical transcriptional output."),
        ("PDZ-binding motif", 501, 504,
         "C-terminal -TWL motif; Oka 2010 JBC",
         "Dispensable for TEAD-driven transcription; tunes localisation only."),
        ("14-3-3 / LATS site S127", 127, 127,
         "Zhao 2007 Genes Dev",
         "The phospho-switch that holds YAP in the cytoplasm. A guide here would "
         "make a point mutant, not a knockout, which is not what an F0 frameshift "
         "screen is for."),
    ],
    "LBR": [
        ("Sterol C14 reductase domain", 212, 615,
         "Tsai 2016 eLife; ERG24/TM7SF2 homology",
         "A second, metabolic function. Any frameshift upstream of residue 212 "
         "removes it as well as the chromatin tether, so a nuclear-shape "
         "phenotype cannot be attributed to tethering without a control."),
        ("HP1 / histone tether (Tudor + basic region)", 1, 100,
         "Ye and Worman 1996 JBC; Polioudaki 2001 EMBO Rep; Hirano 2012 JBC; "
         "Solovei 2013 Cell",
         "The heterochromatin-tethering module. A frameshift here removes the "
         "tether while leaving the reductase domain intact."),
    ],
}


def _get(url: str, tries: int = 5, timeout: int = 120):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": UA, "Accept": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=timeout))
        except Exception as e:                                     # noqa: BLE001
            last = e
            if i < tries - 1:
                time.sleep(2 + 3 * i)
    raise RuntimeError(f"UniProt request failed after {tries} tries: {url}: {last}")


def _reviewed(entry: dict) -> bool:
    """True for Swiss-Prot. Written as a prefix test on purpose: the TrEMBL label
    is "UniProtKB unreviewed (TrEMBL)", which contains "reviewed"."""
    return entry.get("entryType", "").startswith("UniProtKB reviewed")


def lookup(gene: str, taxon: int, cache_dir: str = "") -> dict | None:
    """Best UniProt entry for a gene symbol in one organism. Reviewed entries win,
    then the longest sequence. Returns None when nothing matches."""
    key = (os.path.join(cache_dir, f"uniprot_{gene.lower()}_{taxon}.json")
           if cache_dir else None)
    if key and os.path.exists(key):
        try:
            return json.load(open(key, encoding="utf-8"))
        except Exception:                                          # noqa: BLE001
            pass
    q = urllib.parse.quote(f"gene_exact:{gene} AND organism_id:{taxon}")
    url = (f"{UNIPROT}/search?query={q}&format=json&size=25"
           "&fields=accession,reviewed,protein_name,length,sequence")
    res = _get(url).get("results", [])
    if not res:
        return None
    # "UniProtKB unreviewed (TrEMBL)" contains the substring "reviewed", so test
    # the prefix, never `in`. Reviewed (Swiss-Prot) entries always win; among
    # equals, the longest sequence. Swiss-Prot's canonical isoform is what the
    # curated domain boundaries in EXTRA are numbered against.
    res.sort(key=lambda d: (_reviewed(d), d.get("sequence", {}).get("length", 0)),
             reverse=True)
    acc = res[0]["primaryAccession"]
    full = _get(f"{UNIPROT}/{acc}.json")
    out = {
        "acc": acc,
        "reviewed": _reviewed(full),
        "name": (full.get("proteinDescription", {}).get("recommendedName", {})
                 .get("fullName", {}).get("value", "")),
        "seq": full["sequence"]["value"],
        "feat": [],
    }
    for f in full.get("features", []):
        t = f["type"]
        if t not in KEEP and t not in SITES:
            continue
        out["feat"].append({
            "type": t,
            "start": f["location"]["start"]["value"],
            "end": f["location"]["end"]["value"],
            "desc": f.get("description", ""),
        })
    if key:
        try:
            json.dump(out, open(key, "w", encoding="utf-8"))
        except Exception:                                          # noqa: BLE001
            pass
    return out


def _label(f: dict) -> str:
    d = (f.get("desc") or "").strip()
    t = f["type"]
    if t == "Domain":
        return d or "domain"
    if t == "Transmembrane":
        return "TM helix"
    if t == "Topological domain":
        return d or "topological"
    if t == "Coiled coil":
        return "coiled coil"
    if t == "Region":
        return d or "region"
    if d:
        return f"{t}: {d}"
    return t


def _same_feature(d: dict, lab: str, ftype: str, span) -> bool:
    """Is an existing domain the same feature as a projected span?

    Same label, or the same generic type (TM helices and WW domains carry
    numbered labels that differ between species), plus reciprocal overlap of at
    least 60% of each span."""
    if d["proj"]:
        return False
    same_name = d["label"] == lab or (d["type"] == ftype and ftype in
                                      ("Transmembrane", "Coiled coil"))
    if not same_name and not (d["label"].split()[0] == lab.split()[0]
                              and d["type"] == ftype):
        return False
    lo = max(d["start"], span[0])
    hi = min(d["end"], span[1])
    if hi < lo:
        return False
    ov = hi - lo + 1
    return (ov / (d["end"] - d["start"] + 1) >= 0.6
            and ov / (span[1] - span[0] + 1) >= 0.6)


def build(gene: str, cache_dir: str = "", zf_taxon: int = 7955,
          hs_taxon: int = 9606, window: int = 11) -> dict | None:
    """Everything the pipeline needs for one gene. Returns None if UniProt has no
    zebrafish entry, in which case the pipeline skips protein annotation."""
    zf = lookup(gene, zf_taxon, cache_dir)
    if not zf:
        return None
    hs = lookup(gene, hs_taxon, cache_dir)
    out = {
        "zf_acc": zf["acc"], "zf_len": len(zf["seq"]), "zf_name": zf["name"],
        "zf_reviewed": zf["reviewed"], "zf_seq": zf["seq"],
        "hs_acc": None, "hs_len": None, "hs_name": None,
        "ident_pct": None, "simil_pct": None, "aligned": None,
        "domains": [], "sites": [], "ident": [], "simil": [],
        "window": window, "notes": [],
    }

    # --- domains annotated directly on the zebrafish protein -----------------
    for f in zf["feat"]:
        if f["type"] in SITES and f["start"] == f["end"]:
            out["sites"].append({"aa": f["start"], "label": _label(f),
                                 "src": f"UniProt {zf['acc']}", "proj": False,
                                 "note": ""})
            continue
        out["domains"].append({
            "start": f["start"], "end": f["end"], "label": _label(f),
            "type": f["type"], "src": f"UniProt {zf['acc']}", "proj": False,
            "note": "",
        })

    if not hs:
        out["notes"].append(
            f"No human UniProt entry found for {gene}, so conservation and "
            "projected domains are unavailable.")
        out["domains"].sort(key=lambda d: (d["start"], d["end"]))
        return out

    # --- align and project ---------------------------------------------------
    a1, a2 = pa.align(hs["seq"], zf["seq"])
    h2z, _z2h = pa.residue_maps(a1, a2)
    npair = sum(1 for x, y in zip(a1, a2) if x != "-" and y != "-")
    nid = sum(1 for x, y in zip(a1, a2) if x == y and x != "-")
    nsim = sum(1 for x, y in zip(a1, a2)
               if x != "-" and y != "-" and pa.sc(x, y) > 0)
    out.update(hs_acc=hs["acc"], hs_len=len(hs["seq"]), hs_name=hs["name"],
               aligned=npair,
               ident_pct=round(100.0 * nid / npair, 1) if npair else None,
               simil_pct=round(100.0 * nsim / npair, 1) if npair else None)
    ident, simil = pa.conservation(a1, a2, len(zf["seq"]), window)
    out["ident"], out["simil"] = ident, simil

    def project(s: int, e: int):
        zz = [h2z[p] for p in range(s, e + 1) if p in h2z]
        return (min(zz), max(zz)) if zz else None

    for f in hs["feat"]:
        span = project(f["start"], f["end"])
        if not span:
            continue
        lab = _label(f)
        if f["type"] in SITES and f["start"] == f["end"]:
            out["sites"].append({
                "aa": span[0], "label": lab, "note": "",
                "src": f"UniProt {hs['acc']} (human), projected", "proj": True})
            continue
        # Drop a projected span that merely restates a zebrafish annotation: same
        # kind of feature, heavily overlapping. Without this, every TM helix and
        # WW domain appears twice (once annotated, once projected) a few residues
        # apart, which clutters the guide rows without adding information.
        if any(_same_feature(d, lab, f["type"], span) for d in out["domains"]):
            continue
        out["domains"].append({
            "start": span[0], "end": span[1], "label": lab, "type": f["type"],
            "src": f"UniProt {hs['acc']} (human), projected", "proj": True,
            "note": "",
        })
    for nm, s, e, src, why in EXTRA.get(gene.upper(), []):
        span = project(s, e)
        if not span:
            continue
        if s == e:
            out["sites"].append({"aa": span[0], "label": nm, "note": why,
                                 "src": src + ", projected", "proj": True})
            continue
        out["domains"].append({
            "start": span[0], "end": span[1], "label": nm, "type": "Curated",
            "src": src + ", projected", "proj": True, "note": why,
        })
    out["domains"].sort(key=lambda d: (d["start"], d["end"]))
    out["sites"].sort(key=lambda d: d["aa"])
    return out


# ------------------------------------------------------------------ CDS geometry
def cds_index(model):
    """Per-exon coding geometry for the reference transcript.

    Returns a list of dicts in transcript order, each giving the exon's coding
    interval in genomic coordinates, its offset into the CDS, the reading frame
    at its 5' end, and the residue range it encodes.

    Coordinates follow the rest of the pipeline and UCSC: 0-based, half-open, so
    `cstart <= pos < cend`. Exons arrive already in transcription order."""
    _t, strand, exons, cds0, cds1, _cdslen = model
    out, cum = [], 0
    for n, (s, e) in enumerate(exons, 1):
        cs, ce = max(s, cds0), min(e, cds1)
        clen = max(0, ce - cs)
        d = {"n": n, "start": s, "end": e, "len": e - s,
             "cstart": cs if clen else None, "cend": ce if clen else None,
             "coding": clen, "cds_offset": cum if clen else None,
             "phase": (cum % 3) if clen else None,
             "aa0": (cum // 3 + 1) if clen else None,
             "aa1": ((cum + clen - 1) // 3 + 1) if clen else None}
        if clen:
            cum += clen
        out.append(d)
    return out


def cds_offset_of(pos: int, idx, strand: str):
    """0-based offset of a 0-based genomic position into the CDS, or None if the
    position is not inside coding sequence."""
    for d in idx:
        if d["coding"] and d["cstart"] <= pos < d["cend"]:
            off = (d["cend"] - 1 - pos) if strand == "-" else (pos - d["cstart"])
            return d["cds_offset"] + off
    return None


def alt_starts(cds: str):
    """In-frame ATG codons after the initiator, as 1-based residue numbers.

    An in-frame downstream ATG is the documented escape route for a 5' frameshift:
    the ribosome re-initiates and makes an N-terminally truncated protein
    (Anderson 2017; Lalonde 2017; Smits 2019, cited in Kroll 2022). This is the
    reason the protocol puts "not exon 1" on its exon-priority list."""
    return [p // 3 + 1 for p in range(3, len(cds) - 2, 3) if cds[p:p + 3] == "ATG"]


def first_coding_exon(idx) -> int | None:
    for d in idx:
        if d["coding"]:
            return d["n"]
    return None
