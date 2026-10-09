#!/usr/bin/env python3
"""Write one BED row per mature tRNA position with Sprinzl-style labels."""

import argparse
import csv
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECTS_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
DEFAULT_TRNA_DIR = os.path.join(PROJECTS_DIR, "ncRNA_biallelic_project_files", "tRNA", "hg38-tRNAs")

FIELDS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "trna_id",
    "trnascan_id",
    "amino_acid",
    "anticodon",
    "mature_pos",
    "pretrna_pos",
    "sprinzl_pos",
    "base",
    "trnascan_score",
    "origin",
    "note",
]


def read_name_map(path):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["tRNAscan-SE_id"]: row["GtRNAdb_id"] for row in reader}


def read_bed(path):
    genes = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, start0, end, name, score, strand = parts[:6]
            genes[name] = {
                "chrom": chrom,
                "start0": int(start0),
                "end": int(end),
                "score": score,
                "strand": strand,
                "mature_len": sum(int(x) for x in parts[10].rstrip(",").split(",") if x),
            }
    return genes


def read_trnascan_out(path):
    records = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("-", "Sequence", "Name")):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 15:
                continue
            begin, intron_begin, intron_end = int(parts[2]), int(parts[6]), int(parts[7])
            strand = "+" if begin <= int(parts[3]) else "-"
            intron = None
            if intron_begin and intron_end:
                intron = (intron_begin - begin + 1, intron_end - begin + 1) if strand == "+" else (begin - intron_begin + 1, begin - intron_end + 1)
                intron = tuple(sorted(intron))
            records[f"{parts[0]}.trna{parts[1]}"] = {
                "amino_acid": parts[4],
                "anticodon": parts[5],
                "score": parts[8],
                "origin": parts[14],
                "note": " ".join(parts[15:]) if len(parts) > 15 else "",
                "intron": intron,
            }
    return records


def read_ss_sequences(path):
    seqs = {}
    name = None
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("chr") and ".trna" in line:
                name = line.split()[0]
            elif line.startswith("Seq:"):
                seqs[name] = line.split(":", 1)[1].strip().upper()
    return seqs


def sprinzl_label(pos, mature_len):
    variable_end = mature_len - 25
    if pos <= 43:
        return str(pos)
    if pos <= variable_end:
        offset = pos - 44
        return str(44 + offset) if offset < 4 else f"47{chr(ord('a') + offset - 4)}"
    if pos <= mature_len - 8:
        return str(49 + pos - variable_end - 1)
    return str(66 + pos - (mature_len - 7))


def mature_to_pre_pos(pos, intron):
    if not intron:
        return pos
    intron_start, intron_end = intron
    return pos + intron_end - intron_start + 1 if pos >= intron_start else pos


def pre_to_bed(chrom_start0, chrom_end, strand, pos):
    if strand == "+":
        return chrom_start0 + pos - 1, chrom_start0 + pos
    return chrom_end - pos, chrom_end - pos + 1


def wanted(record, exclude_pseudogenes):
    return not exclude_pseudogenes or "pseudo" not in record["note"].lower()


def write_positions(trna_dir, out_path, header, exclude_pseudogenes):
    names = read_name_map(os.path.join(trna_dir, "hg38-tRNAs_name_map.txt"))
    genes = read_bed(os.path.join(trna_dir, "hg38-tRNAs.bed"))
    scan = read_trnascan_out(os.path.join(trna_dir, "hg38-tRNAs-detailed.out"))
    seqs = read_ss_sequences(os.path.join(trna_dir, "hg38-tRNAs-detailed.ss"))
    rows = []
    for scan_id, trna_id in names.items():
        if trna_id not in genes or scan_id not in scan or scan_id not in seqs or not wanted(scan[scan_id], exclude_pseudogenes):
            continue
        gene = genes[trna_id]
        rec = scan[scan_id]
        for mature_pos in range(1, gene["mature_len"] + 1):
            pre_pos = mature_to_pre_pos(mature_pos, rec["intron"])
            start0, end = pre_to_bed(gene["start0"], gene["end"], gene["strand"], pre_pos)
            sprinzl = sprinzl_label(mature_pos, gene["mature_len"])
            rows.append((
                gene["chrom"], start0, end, f"{trna_id}|sprinzl_{sprinzl}", gene["score"], gene["strand"],
                trna_id, scan_id, rec["amino_acid"], rec["anticodon"], mature_pos, pre_pos, sprinzl,
                seqs[scan_id][pre_pos - 1], rec["score"], rec["origin"], rec["note"],
            ))
    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        if header:
            writer.writerow(FIELDS)
        writer.writerows(rows)
    return len(rows), len({row[6] for row in rows})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trna-dir", default=DEFAULT_TRNA_DIR)
    parser.add_argument("--out", default="hg38-tRNA-sprinzl-positions.bed")
    parser.add_argument("--header", action="store_true")
    parser.add_argument("--exclude-pseudogenes", action="store_true")
    args = parser.parse_args()
    n_rows, n_genes = write_positions(args.trna_dir, args.out, args.header, args.exclude_pseudogenes)
    label = "non-pseudo " if args.exclude_pseudogenes else ""
    print(f"Wrote {n_rows} Sprinzl-position intervals for {n_genes} {label}tRNAs to {args.out}")


if __name__ == "__main__":
    main()
