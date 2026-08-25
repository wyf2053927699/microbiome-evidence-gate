import csv, hashlib, json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw = root / "01_raw_data" / "external_provenance"
out = root / "06_results" / "provenance_and_portability"
audit = root / "09_audit"
out.mkdir(parents=True, exist_ok=True)

edge_file = root / "06_results" / "axis_prioritisation" / "edge_level_provenance.csv"
with edge_file.open(encoding="utf-8-sig", newline="") as f:
    edges = list(csv.DictReader(f))

prov_fields = [
    "candidate_family", "source_node", "target_node", "edge_type",
    "project_measurement_state", "external_prior_state", "gutmgene_v2_state",
    "evidence_species", "identifier_state", "claim_use", "DOI", "PMID", "notes"
]
prov_rows = []
for e in edges:
    measured = "project-measured" if e.get("measured_in_project", "").upper() == "YES" else "not-project-measured"
    prior = "curated-literature-prior" if e.get("DOI") or e.get("PMID") else "missing"
    model = e.get("evidence_model", "").lower()
    species = "human" if "human" in model else ("mouse/non-human" if "mouse" in model or "pig" in model else "unspecified")
    prov_rows.append({
        "candidate_family": e.get("candidate_family", ""),
        "source_node": e.get("source_node", ""),
        "target_node": e.get("target_node", ""),
        "edge_type": e.get("edge_type", ""),
        "project_measurement_state": measured,
        "external_prior_state": prior,
        "gutmgene_v2_state": "not-verified-bulk-endpoint-unavailable-2026-08-19",
        "evidence_species": species,
        "identifier_state": "text-label-only; database identifier mapping pending",
        "claim_use": "hypothesis-generation-only",
        "DOI": e.get("DOI", ""),
        "PMID": e.get("PMID", ""),
        "notes": e.get("notes", "")
    })

with (out / "gutmgene_provenance_layer.csv").open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=prov_fields)
    w.writeheader(); w.writerows(prov_rows)

metadata = raw / "hmp2_metadata_2018-08-20.csv"
type_participants = {}
type_rows = {}
if metadata.exists():
    with metadata.open(encoding="utf-8-sig", newline="", errors="replace") as f:
        for row in csv.DictReader(f):
            dtype = (row.get("data_type") or "MISSING").strip()
            participant = (row.get("Participant ID") or "").strip()
            type_rows[dtype] = type_rows.get(dtype, 0) + 1
            if participant:
                type_participants.setdefault(dtype, set()).add(participant)

summary_rows = []
for dtype in sorted(type_rows):
    summary_rows.append({"data_type": dtype, "metadata_rows": type_rows[dtype],
                         "unique_participants": len(type_participants.get(dtype, set()))})
with (out / "hmp2_metadata_coverage.csv").open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["data_type","metadata_rows","unique_participants"])
    w.writeheader(); w.writerows(summary_rows)

manifest = {
    "generated": "2026-08-19",
    "gutmgene_v2": {
        "article_doi": "10.1093/nar/gkae1002",
        "official_site": "https://bio-computing.hrbmu.edu.cn/gutmgene/",
        "bulk_data_status": "NO-GO: official site timed out and no versioned bulk export was acquired",
        "permitted_use": "literature-prior rows only; do not label them as gutMGene-verified"
    },
    "hmp2": {
        "official_portal": "https://www.ibdmdb.org/results",
        "bioproject": "PRJNA395569 / PRJNA398089",
        "metadata_file": str(metadata.relative_to(root)) if metadata.exists() else None,
        "metadata_sha256": hashlib.sha256(metadata.read_bytes()).hexdigest() if metadata.exists() else None,
        "pilot_status": "GO for frozen acquisition/mapping pilot; NOT YET external portability validation"
    }
}
(out / "provenance_portability_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

