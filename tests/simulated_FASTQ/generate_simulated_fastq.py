#!/usr/bin/env python3
"""Generate 10 single-end 35-bp positive-control FASTQ files for ViruScanNIPT.

Reads are sampled from public RefSeq-derived viral gene/region fragments.
No sequencing errors are introduced, so every read is an exact substring of
one of the source fragments (or its reverse complement).
"""
from pathlib import Path
import argparse
import csv
import random

REFS = {
    "HSV1": {
        "accession": "NC_001806.2",
        "groupname": "Herpesvirus 1/2",
        "organism": "Human herpesvirus 1 strain 17",
        "segments": [
            {
                "name": "HSV1_segment_A",
                "coordinates": "141142-143696 (sequence excerpt)",
                "sequence": """
GCATTTAAGGCGTTGTTGTGTTGACTTTGCCTCTTCTGGCGGGTTGGTGCGGTGCTGTTTGTTGGGCTCC
CATTTTACCCGAAGATCGGCTGCTATCCCCGGGACATGGATCGCGGGGCGGTGGTGGGGTTTCTTCTCGG
TGTTTGTGTTGTATCGTGCTTGGCGGGAACGCCCAAAACGTCCTGGAGACGGGTGAGTGTCGGCGAGGAC
GTTTCGTTGCTTCCAGCTCCGGGGCCTACGGGGCGCGGCCCGACCCAGAAACTACTATGGGCCGTGGAAC
CCCTGGATGGGTGCGGCCCCTTACACCCGTCGTGGGTCTCGCTGATGCCCCCCAAGCAGGTGCCCGAGAC
GGTCGTGGATGCGGCGTGCATGCGCGCTCCGGTCCCGCTGGCGATGGCGTACGCCCCCCCGGCCCCATCT
GCGACCGGGGGTCTACGAACGGACTTCGTGTGGCAGGAGCGCGCGGCCGTGGTTAACCGGAGTCTGGTTA
TTCACGGGGTCCGAGAGACGGACAGCGGCCTGTATACCCTGTCCGTGGGCGACATAAAGGACCCGGCTCG
CCAAGTGGCCTCGGTGGTCCTGGTGGTGCAACCGGCCCCAGTTCCGACCCCACCCCCGACCCCAGCCGAT
TACGACGAGGATGACAATGACGAGGGCGAGGACGAAAGTCTCGCCGGCACTCCCGCCAGCGGGACCCCCC
GGCTCCCGCCTCCCCCCGCCCCCCCGAGGTCTTGGCCCAGCGCCCCCGAAGTCTCACATGTGCGTGGGGT
GACCGTGCGTATGGAGACTCCGGAAGCTATCCTGTTTTCCCCCGGGGAGACGTTCAGCACGAACGTCTCC
""",
            },
            {
                "name": "HSV1_segment_B",
                "coordinates": "96311-98668 (sequence excerpt)",
                "sequence": """
CATGGCCCCGGGGCGGGTGGGCCTTGCCGTGGTCCTGTGGAGCCTGTTGTGGCTCGGGGCGGGGGTGTCC
GGGGGCTCGGAAACTGCCTCCACCGGGCCCACGATCACCGCGGGAGCGGTGACGAACGCGAGCGAGGCCC
CCACATCGGGGTCCCCCGGGTCAGCCGCCAGCCCGGAGGTCACCCCCACATCGACCCCAAACCCCAACAA
TGTCACACAAAACAAAACCACCCCCACCGAGCCGGCCAGCCCCCCAACAACCCCCAAGCCCACCTCCACG
""",
            },
            {
                "name": "HSV1_segment_C",
                "coordinates": "137598-141051 (sequence excerpt)",
                "sequence": """
TTAAATACATGCCCGGGACCCATGAAGCATGCGCGACTTCCGGGCCTCAGAACCCACCCGAAACGGCCAA
CGGACGTCTGAGCCAGGCCTGGCTATCCGGAGAAACAGCACACGACTTGGCGTTCTGTGTGTCGCGATGT
CTCTGCGCGCAGTCTGGCATCTGGGGCTTTTGGGAAGCCTCGTGGGGGCTGTTCTTGCCGCCACCCATCG
GGGACCTGCGGCCAACACAACGGACCCCTTAACGCACGCCCCAGTGTCCCCTCACCCCAGCCCCCTGGGG
""",
            },
        ],
        "source": "VIOLIN/Protegen gene sequence entries; RefSeq accession NC_001806.2",
    },
    "CMV": {
        "accession": "NC_006273.2",
        "groupname": "Herpesvirus 5",
        "organism": "Human herpesvirus 5 strain Merlin",
        "segments": [
            {
                "name": "CMV_gB",
                "coordinates": "82065-84788",
                "sequence": """
TTCAGACGTTCTCTTCTTCGTCAGAGTCTTTCAAGTGTCGGTAGCCGTTTTTGCGATGTCGCAGTCGGTC
TAGTAGGTTGGGCTTCTGTCCCTTGTCCTGCGTGCCAGTCCGTCCGTCCAAAGAATCTGTACCGTTCTGC
TGCGCTCGCTGCTCTGCGTCCAGACGGGCCAGGGCCAGAAGCATCTGGTAAGCCTGCTCGTTGGTGTAAG
GCGGAGCCGCCGTGGATGCATCAGACGACGGTGGTCCCGGTCCTTTGCGACCAGAATTATAAACACTTTC
CTCGTAGGAAGGCGGAGCCTGTAACGACGTGTCTTTGGTGCTGCCCGACGTCACGGTGGTCCCGTCGGCG
GACACCAGATAGGGAAAGAGGTTCTGCAGCGGCTGCGTGCACAAACGCCGCTGTCGAGTATAGATCAAAT
AAGTGATAATGACTACAGCTATGGCCACGAGGATGATGGTGAACGCTCCGAAGGGGTTTTTGAGGAAGGT
GGCAACGCCTTCGACCACGGAGGCCACCGCGCCACCCACGGCCCCAATGGCTACGCCAACGGCCTTTCCC
GCGGCGCCCAGGCCGCTCATGAGGTCGTCCAGACCCTTGAGGTAGGGCGGTAGCGGGTCGACTACCTTGT
CCTCCACGTACTTTACCCGCTGCTTGTACGAGTTGAATTCGCGCATGATCTCTTCGAGGTCAAAAACGTT
GCTGGAACGCAGCTCTTTCTGCGAGTAAAGTTCCAGTACCCTGAAGTCGGTATTTTCCAGCGGGTCGATA
TCCAGGGCGATCATGCTGTCGACGGTGGAGATACTGCTGAGGTCAATCATGCGTTTGAAGAGGTAGTCCA
CGTACTCGTAGGCCGAGTTCCCGGCGATGAAGATCTTGAGGCTGGGAAGCTGACATTCCTCAGTGCGGTG
GTTGCCCAACAGGATTTCGTTGTCCTCGCCCAGTTGACCGTACTGCACGTACGAGCTGTTGGCGAAATTA
AAGATGACCACGGGTCGTGAGTAGCAGCGTCCTGGCGACTCCTTCACGTTCATATCACGCAGCACCTTGA
CGCTGGTTTGGTTGATGGTCACGCAGCTGGCCAGGCCCAAGACATCACCCATGAAACGCGCGGCAATCGG
TTTGTTGTAAATGGCCGAGAGAATGGCTGACGGGTTGATCTTGCTGAGTTCCTTGAAGACCTCTAGGGTG
CGCCGTTGATCCACACACCAGGCTTCTGCGATTTGCGCCAGCGCCCGGTTGATGTAACCGCGCAACGTGT
CATAGGTGAACTGCAGCTGGGCGTAGACCAGATTGTGCACCGATTCCATGTTGGATAAATGAGTTGCATT
GTTGCCATCTGTACTTCTTTTGGTTCTATTATGAGTAAGATTCAGACTGGAGCGGTTGGCCAAACGTTCG
AGTTCCACCAGAGATTTTTGCTTGATACCTTGCCAGAACACTACCAAACCACCAGTGGTTTCAAAGACGG
ACACGTTTCCATATTTTTCATATGTTTGATTGTATGAAGTATTGAAAATCTGCTGTAACTTATTTATAGC
CTCATCACGTACGCAGTCCAGCGCAGAGTCGGACATGTTCACCTCTTGCTTCTTAGATAAGAAAGTGGCG
GTCATTTTGGCAGAAGAAAAGTGATACGAGTCCTCGGCTTCGGAACGAATGGTGCGTTCCGAGGCTTCCC
AGAAAGTGAGTTGACAAGTGACATTCTTTTCGTCCTGTATATCCCAGGAGATCACCGAGTCCGCACGTTC
AAGAAAAGCCACCAACCTGTGGGTCTCTAACGCAGAATTCGGTCTTCCAAAGTCGGAGACGATAGTGTAG
TTCGGAAAAATGAAAAACTTGTCGGCGTTTTCTCCAAAGTAGCTGGCATTGCGATTGGTTCCGTTGTAGA
AAGGAGAAATGTCAACCACGTCACCCGTGGAAGTGGCGAAAAAATGATAAGGATATTTGGAGCGCGCAGT
AGTGATGGTCACCATACAATTCAGATTACAGGTCTCACGATAGAGCCAGGTGCTGCCGCGGCTGTGCCAT
TGATCCTTGACCGTCACGTAACGGGTACTGTGGGTGTTGGAATAATCGTCGGGCATTAATTGCATGGTTT
TGTTTTCATAGCTGTCCCTATGATAAGCCACGAAAACCGTGCCTGCTATAACGCGGCTGTAGGAACTGTA
GCACTGACTGTGGCTGTTGATATGATGAATCTCCCACATAGGAGGCGCCACGTATTCCGTGTTGCTGCCC
AGCAGATAAGTGGTGTGGATGTAAGCGTAGCTACGACGAAACGTCAAAACCTTCTGGTAGACTCGTACCT
TAAAGGTGTGCGCGACGATGTTGCGTTTGTAGACCACCATGATGCCCTCGTCCAGGTCTTCATTGATGGG
CTTCATCGAGGTGCAGACGATATTACGTTCAAAGCGAATAAGATCCGTACCCTGGGCCATAGAACACACG
CGATAGGGGTACTTGGTGGTATTGACCCCCACCACATCTCCGTACTTGAGGGTAGTGTTGTAGATGGTCT
CGTTAACACCATGGCTGACCGTTTGGGAAGAAGTTACGCGTTGAGAGACTGAACCGGATCGAGAGTGAGC
AGCAGACGTCGTATGAGAGGAATGGTGACTGTGAGTAGCAGAAGTTCCACGAGTAGAAGATGAGGAAACC
GCAGCACCCAGACAGACGATACACAAGTTAACGCAGACTACCAGGCACCAGATCCTGGATTCCA
""",
            },
        ],
        "source": "VIOLIN/Vaxgen gene sequence entry for glycoprotein B; RefSeq accession NC_006273.2",
    },
    "HPV16": {
        "accession": "NC_001526.4",
        "groupname": "Papillomavirus",
        "organism": "Human papillomavirus type 16",
        "segments": [
            {
                "name": "HPV16_L1",
                "coordinates": "4774-6291",
                "sequence": """
GATGTCTCTTTGGCTGCCTAGTGAGGCCACTGTCTACTTGCCTCCTGTCCCAGTATCTAAGGTTGTAAGC
ACGGATGAATATGTTGCACGCACAAACATATATTATCATGCAGGAACATCCAGACTACTTGCAGTTGGAC
ATCCCTATTTTCCTATTAAAAAACCTAACAATAACAAAATATTAGTTCCTAAAGTATCAGGATTACAATA
CAGGGTATTTAGAATACATTTACCTGACCCCAATAAGTTTGGTTTTCCTGACACCTCATTTTATAATCCA
GATACACAGCGGCTGGTTTGGGCCTGTGTAGGTGTTGAGGTAGGTCGTGGTCAGCCATTAGGTGTGGGCA
TTAGTGGCCATCCTTTATTAAATAAATTGGATGACACAGAAAATGCTAGTGCTTATGCAGCAAATGCAGG
TGTGGATAATAGAGAATGTATATCTATGGATTACAAACAAACACAATTGTGTTTAATTGGTTGCAAACCA
CCTATAGGGGAACACTGGGGCAAAGGATCCCCATGTACCAATGTTGCAGTAAATCCAGGTGATTGTCCAC
CATTAGAGTTAATAAACACAGTTATTCAGGATGGTGATATGGTTGATACTGGCTTTGGTGCTATGGACTT
TACTACATTACAGGCTAACAAAAGTGAAGTTCCACTGGATATTTGTACATCTATTTGCAAATATCCAGAT
TATATTAAAATGGTGTCAGAACCATATGGCGACAGCTTATTTTTTTATTTACGAAGGGAACAAATGTTTG
TTAGACATTTATTTAATAGGGCTGGTACTGTTGGTGAAAATGTACCAGACGATTTATACATTAAAGGCTC
TGGGTCTACTGCAAATTTAGCCAGTTCAAATTATTTTCCTACACCTAGTGGTTCTATGGTTACCTCTGAT
GCCCAAATATTCAATAAACCTTATTGGTTACAACGAGCACAGGGCCACAATAATGGCATTTGTTGGGGTA
ACCAACTATTTGTTACTGTTGTTGATACTACACGCAGTACAAATATGTCATTATGTGCTGCCATATCTAC
TTCAGAAACTACATATAAAAATACTAACTTTAAGGAGTACCTACGACATGGGGAGGAATATGATTTACAG
TTTATTTTTCAACTGTGCAAAATAACCTTAACTGCAGACGTTATGACATACATACATTCTATGAATTCCA
CTATTTTGGAGGACTGGAATTTTGGTCTACAACCTCCCCCAGGAGGCACACTAGAAGATACTTATAGGTT
TGTAACATCCCAGGCAATTGCTTGTCAAAAACATACACCTCCAGCACCTAAAGAAGATCCCCTTAAAAAA
TACACTTTTTGGGAAGTAAATTTAAAGGAAAAGTTTTCTGCAGACCTAGATCAGTTTCCTTTAGGACGCA
AATTTTTACTACAAGCAGGATTGAAGGCCAAACCAAAATTTACATTAGGAAAACGAAAAGCTACACCCAC
CACCTCATCTACCTCTACAACTGCTAAACGCAAAAAACGTAAGCTGTA
""",
            },
        ],
        "source": "VIOLIN/Protegen gene sequence entry for HPV16 L1; RefSeq accession NC_001526.4",
    },
    "VZV": {
        "accession": "NC_001348.1",
        "groupname": "Herpesvirus 3",
        "organism": "Human herpesvirus 3 (Varicella-zoster virus)",
        "segments": [
            {
                "name": "VZV_gE",
                "coordinates": "115807-117720",
                "sequence": """
TATGGGGACAGTTAATAAACCTGTGGTGGGGGTATTGATGGGGTTCGGAATTATCACGGGAACGTTGCGT
ATAACGAATCCGGTCAGAGCATCCGTCTTGCGATACGATGATTTTCACACCGATGAAGACAAACTGGATA
CAAACTCCGTATATGAGCCTTACTACCATTCAGATCATGCGGAGTCTTCATGGGTAAATCGGGGAGAGTC
TTCGCGAAAAGCGTACGATCATAACTCACCTTATATATGGCCACGTAATGATTATGATGGATTTTTAGAG
AACGCACACGAACACCATGGGGTGTATAATCAGGGCCGTGGTATCGATAGCGGGGAACGGTTAATGCAAC
CCACACAAATGTCTGCACAGGAGGATCTTGGGGACGATACGGGCATCCACGTTATCCCTACGTTAAACGG
CGATGACAGACATAAAATTGTAAATGTGGACCAACGTCAATACGGTGACGTGTTTAAAGGAGATCTTAAT
CCAAAACCCCAAGGCCAAAGACTCATTGAGGTGTCAGTGGAAGAAAATCACCCGTTTACTTTACGCGCAC
CGATTCAGCGGATTTATGGAGTCCGGTACACCGAGACTTGGAGCTTTTTGCCGTCATTAACCTGTACGGG
AGACGCAGCGCCCGCCATCCAGCATATATGTTTAAAACATACAACATGCTTTCAAGACGTGGTGGTGGAT
GTGGATTGCGCGGAAAATACTAAAGAGGATCAGTTGGCCGAAATCAGTTACCGTTTTCAAGGTAAGAAGG
AAGCGGACCAACCGTGGATTGTTGTAAACACGAGCACACTGTTTGATGAACTCGAATTAGACCCCCCCGA
GATTGAACCGGGTGTCTTGAAAGTACTTCGGACAGAAAAACAATACTTGGGTGTGTACATTTGGAACATG
CGCGGCTCCGATGGTACGTCTACCTACGCCACGTTTTTGGTCACCTGGAAAGGGGATGAAAAAACAAGAA
ACCCTACGCCCGCAGTAACTCCTCAACCAAGAGGGGCTGAGTTTCATATGTGGAATTACCACTCGCATGT
ATTTTCAGTTGGTGATACGTTTAGCTTGGCAATGCATCTTCAGTATAAGATACATGAAGCGCCATTTGAT
TTGCTGTTAGAGTGGTTGTATGTCCCCATCGATCCTACATGTCAACCAATGCGGTTATATTCTACGTGTT
TGTATCATCCCAACGCACCCCAATGCCTCTCTCATATGAATTCCGGTTGTACATTTACCTCGCCACATTT
AGCCCAGCGTGTTGCAAGCACAGTGTATCAAAATTGTGAACATGCAGATAACTACACCGCATATTGTCTG
GGAATATCTCATATGGAGCCTAGCTTTGGTCTAATCTTACACGACGGGGGCACCACGTTAAAGTTTGTAG
ATACACCCGAGAGTTTGTCGGGATTATACGTTTTTGTGGTGTATTTTAACGGGCATGTTGAAGCCGTAGC
ATACACTGTTGTATCCACAGTAGATCATTTTGTAAACGCAATTGAAGAGCGTGGATTTCCGCCAACGGCC
GGTCAGCCACCGGCGACTACTAAACCCAAGGAAATTACCCCCGTAAACCCCGGAACGTCACCACTTCTAC
GATATGCCGCATGGACCGGAGGGCTTGCAGCAGTAGTACTTTTATGTCTCGTAATATTTTTAATCTGTAC
GGCTAAACGAATGAGGGTTAAAGCCTATAGGGTAGACAAGTCCCCGTATAACCAAAGCATGTATTACGCT
GGCCTTCCAGTGGACGATTTCGAGGACTCGGAATCTACGGATACGGAAGAAGAGTTTGGTAACGCGATTG
GAGGGAGTCACGGGGGTTCGAGTTACACGGTGTATATAGATAAGACCCGGTGATCACCGAACCGGGGCAA
CGCCGAGCGTGTAAATTTAAATAA
""",
            },
        ],
        "source": "VIOLIN/Vaxgen gene sequence entry for VZV glycoprotein E; RefSeq accession NC_001348.1",
    },
}

# Per-sample read counts by virus. All totals are exactly 1,000 reads.
SAMPLES = {
    "SimSample01": {"HSV1": 1000},
    "SimSample02": {"CMV": 1000},
    "SimSample03": {"HPV16": 1000},
    "SimSample04": {"VZV": 1000},
    "SimSample05": {"HSV1": 400, "CMV": 300, "HPV16": 200, "VZV": 100},
    "SimSample06": {"CMV": 500, "VZV": 300, "HSV1": 200},
    "SimSample07": {"HPV16": 500, "CMV": 300, "HSV1": 200},
    "SimSample08": {"VZV": 400, "HSV1": 300, "HPV16": 200, "CMV": 100},
    "SimSample09": {"HSV1": 250, "CMV": 250, "HPV16": 250, "VZV": 250},
    "SimSample10": {"HSV1": 300, "CMV": 300, "HPV16": 200, "VZV": 200},
}

COMP = str.maketrans("ACGT", "TGCA")
def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]

def clean_segments():
    for virus, meta in REFS.items():
        for seg in meta["segments"]:
            seq = "".join(seg["sequence"].split()).upper()
            bad = sorted(set(seq) - set("ACGT"))
            if bad:
                raise ValueError(f"Invalid bases in {virus}/{seg['name']}: {bad}")
            if len(seq) < 35:
                raise ValueError(f"Segment too short: {virus}/{seg['name']}")
            seg["sequence"] = seq

def candidates_for_virus(virus: str, read_len: int):
    candidates = []
    for seg in REFS[virus]["segments"]:
        seq = seg["sequence"]
        for start in range(len(seq) - read_len + 1):
            fwd = seq[start:start + read_len]
            rev = revcomp(fwd)
            candidates.append((fwd, seg, start + 1, "+"))
            if rev != fwd:
                candidates.append((rev, seg, start + 1, "-"))
    # De-duplicate exact sequence strings while preserving first source coordinate.
    seen = set()
    unique = []
    for item in candidates:
        if item[0] not in seen:
            seen.add(item[0])
            unique.append(item)
    return unique

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="simulated_fastq")
    ap.add_argument("--read-length", type=int, default=35)
    ap.add_argument("--reads-per-sample", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261106)
    args = ap.parse_args()
    if args.reads_per_sample != 1000:
        raise SystemExit("This bundled design is fixed at 1,000 reads per sample.")
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    clean_segments()
    rng = random.Random(args.seed)
    candidate_cache = {v: candidates_for_virus(v, args.read_length) for v in REFS}
    for virus, cand in candidate_cache.items():
        if len(cand) < max(count for comp in SAMPLES.values() for v, count in comp.items() if v == virus):
            raise SystemExit(f"Not enough unique candidate reads for {virus}: {len(cand)}")

    # Save the exact sequence fragments used to create the reads.
    fasta = outdir / "reference_segments.fasta"
    with fasta.open("w") as fh:
        for virus, meta in REFS.items():
            for seg in meta["segments"]:
                fh.write(f">{meta['accession']}|{virus}|{seg['name']}|{seg['coordinates']}\n")
                seq = seg["sequence"]
                for i in range(0, len(seq), 70):
                    fh.write(seq[i:i+70] + "\n")

    manifest_rows = []
    fq_paths = []
    for sample, composition in SAMPLES.items():
        if sum(composition.values()) != args.reads_per_sample:
            raise SystemExit(f"Bad read total for {sample}")
        sample_records = []
        for virus, n_reads in composition.items():
            pool = candidate_cache[virus]
            chosen = rng.sample(pool, n_reads)
            meta = REFS[virus]
            for seq, seg, start, strand in chosen:
                read_no = len(sample_records) + 1
                qname = (f"@{sample}:read{read_no:06d}|{meta['accession']}|"
                         f"{seg['name']}|start={start}|strand={strand}")
                sample_records.append((qname, seq, "+", "I" * args.read_length))
                manifest_rows.append({
                    "sample": sample,
                    "virus": virus,
                    "accession": meta["accession"],
                    "organism": meta["organism"],
                    "groupname": meta["groupname"],
                    "segment": seg["name"],
                    "coordinates": seg["coordinates"],
                    "reads": 1,
                })
        rng.shuffle(sample_records)
        # Re-number after shuffling to keep read IDs sequential and stable.
        fq = outdir / f"{sample}.fq"
        with fq.open("w") as fh:
            for i, (qname, seq, plus, qual) in enumerate(sample_records, 1):
                fields = qname.split("|")
                fields[0] = f"@{sample}:read{i:06d}"
                fh.write("|".join(fields) + "\n")
                fh.write(seq + "\n")
                fh.write(plus + "\n")
                fh.write(qual + "\n")
        fq_paths.append(fq.name)

    # Aggregate the manifest at sample-virus-segment level.
    agg = {}
    for row in manifest_rows:
        key = tuple(row[k] for k in ("sample", "virus", "accession", "organism", "groupname", "segment", "coordinates"))
        agg[key] = agg.get(key, 0) + 1
    with (outdir / "simulation_manifest.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "virus", "accession", "organism", "groupname", "segment", "coordinates", "read_count"])
        for key, count in sorted(agg.items()):
            writer.writerow([*key, count])

    with (outdir / "fq.list").open("w") as fh:
        for name in fq_paths:
            fh.write(name + "\n")

    with (outdir / "sample_composition.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "virus", "accession", "groupname", "read_count"])
        for sample, composition in SAMPLES.items():
            for virus, count in composition.items():
                meta = REFS[virus]
                writer.writerow([sample, virus, meta["accession"], meta["groupname"], count])

    source_rows = [
        ("HSV1", "NC_001806.2", "Herpesvirus 1/2", "https://violinet.org/protegen/gene_detail.php?c_gene_id=507; https://violinet.org/vaxgen/gene_detail.php?c_gene_id=506; https://violinet.org/protegen/gene_detail.php?c_gene_id=501"),
        ("CMV", "NC_006273.2", "Herpesvirus 5", "https://violinet.org/vaxgen/gene_detail.php?c_gene_id=1239"),
        ("HPV16", "NC_001526.4", "Papillomavirus", "https://violinet.org/vaxgen/gene_detail.php?c_gene_id=512"),
        ("VZV", "NC_001348.1", "Herpesvirus 3", "https://violinet.org/vaxgen/gene_detail.php?c_gene_id=662"),
    ]
    with (outdir / "reference_sources.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["virus", "accession", "groupname", "source_page"])
        writer.writerows(source_rows)

    readme = f'''# ViruScanNIPT simulated FASTQ positive controls\n\n'''
    readme += f'''This directory contains 10 single-end FASTQ files. Each file has exactly 1,000 reads, each read is {args.read_length} bp long. Total: 10,000 reads.\n\n'''
    readme += '''## Simulation design\n\n'''
    readme += '''- All reads are derived from DNA-virus reference gene/region fragments listed in `reference_segments.fasta`.\n- The reads are exact 35-bp substrings sampled from the source fragments, including reverse-complement orientations.\n- No sequencing errors or human-background reads are added. This is an intentionally enriched positive-control dataset to make viral BLAST hits easier to observe.\n- FASTQ quality is set to `I` at every base (Phred+33 Q40).\n- Each sample contains 1,000 reads; the virus composition varies by sample. See `simulation_manifest.tsv`.\n- Read identifiers include the reference accession, source segment, start offset within the included segment and strand.\n- Reference accessions are present in this project's `viralGroupname.tsv`: NC_001806, NC_006273, NC_001526 and NC_001348.\n\n'''
    readme += '''## Files\n\n- `SimSample01.fq` ... `SimSample10.fq`: simulated single-end FASTQ files.\n- `fq.list`: ready-to-use FASTQ list for `ViruScanNIPT.pl`; paths are relative to this directory.\n- `simulation_manifest.tsv`: read counts by sample and reference virus/segment.\n- `sample_composition.tsv`: expected per-sample read counts by virus.\n- `reference_sources.tsv`: accession, Groupname and source-page traceability.\n- `reference_segments.fasta`: exact source fragments used to generate reads.\n- `generate_simulated_fastq.py`: deterministic generator (seed 20261106).\n\n'''
    readme += '''## Reference sources and Groupname coverage\n\nThe four selected RefSeq accessions are included in the project's `database/annotation/viralGroupname.tsv`:\n\n| Virus | Accession | Groupname |\n|---|---|---|\n| HSV-1 | NC_001806.2 | Herpesvirus 1/2 |\n| Human cytomegalovirus | NC_006273.2 | Herpesvirus 5 |\n| HPV16 | NC_001526.4 | Papillomavirus |\n| Varicella-zoster virus | NC_001348.1 | Herpesvirus 3 |\n\nSource pages are recorded in `reference_sources.tsv`. The HSV-1 sequences are partial excerpts from three gene entries; the CMV gB, HPV16 L1 and VZV gE sequences correspond to the gene intervals listed in `reference_segments.fasta`.\n\n## Run with ViruScanNIPT\n\nFrom the ViruScanNIPT repository root, either copy this folder under `examples/` or pass the path directly:\n\n```bash\nperl ViruScanNIPT.pl \\\n  --config conf/ViruScanNIPT.conf \\\n  --fq-list /path/to/simulated_fastq/fq.list \\\n  --outdir ViruScanNIPT_simulation_result\n```\n\nThe pipeline performs BWA alignment to hg38 first and BLASTN on unmapped reads. Since this dataset contains only viral reads, most reads are expected to remain unmapped to hg38 and proceed to viral BLASTN; a small number may still align to human sequence by chance or due to local similarity.\n\n## Important limitations\n\nThese are synthetic positive-control reads, not real NIPT/cfDNA data. They are sampled from selected reference gene/region fragments rather than uniformly across complete viral genomes. They are not suitable for estimating clinical sensitivity, specificity, prevalence, viral load or real-world background rates. For full validation, use complete RefSeq FASTA files and add human cfDNA background, realistic fragment-length distributions, sequencing errors and contamination controls.\n'''
    (outdir / "README.md").write_text(readme)

    print(f"Output directory: {outdir}")
    print(f"FASTQ files: {len(fq_paths)}")
    print(f"Reads per file: {args.reads_per_sample}")
    print(f"Read length: {args.read_length}")
    for virus, cand in candidate_cache.items():
        print(f"Unique candidate reads {virus}: {len(cand)}")

if __name__ == "__main__":
    main()
