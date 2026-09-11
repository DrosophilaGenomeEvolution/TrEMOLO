"""Sniffles 2 records normalized to TrEMOLO's existing event-header contract.

VCF POS is the one-based preceding base: numerically it is also the zero-based
insertion breakpoint used by the downstream BAM/BED code. Do not subtract one.
"""

import re
from pathlib import Path


def is_sniffles2(path):
    with open(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if line.lower().startswith("##source=sniffles2"):
                return True
    return False


def records(path):
    seen = set()
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                raise ValueError("{}:{}: malformed VCF record".format(path, number))
            chrom, pos, identifier, ref, alt = fields[:5]
            info = dict((part.split("=", 1) if "=" in part else (part, True))
                        for part in fields[7].split(";") if part)
            kind = info.get("SVTYPE", alt.strip("<>"))
            if kind not in ("INS", "DEL"):
                continue
            if "," in alt:
                raise ValueError("Multiallelic VCF records are not supported: " + identifier)
            if identifier == "." or identifier in seen or re.search(r"[:|\s/]", identifier):
                raise ValueError("Missing, duplicate or unsafe VCF ID: " + identifier)
            seen.add(identifier)
            start = int(pos)
            length = abs(int(info["SVLEN"]))
            end = start if kind == "INS" else int(info.get("END", start + length))
            if start < 1 or end < start:
                raise ValueError("Invalid coordinates: " + identifier)
            sequence = None
            # With --reference, Sniffles 2 prefixes the preceding base. Without
            # a reference it may emit an unanchored ALT with REF=N. SVLEN
            # disambiguates these representations; matching bases alone cannot.
            allele = alt if kind == "INS" else ref
            anchor = ref if kind == "INS" else alt
            if not alt.startswith("<") and re.fullmatch(r"[ACGTURYKMSWBDHVNacgturykmswbdhvn]+", allele):
                if len(allele) == length:
                    sequence = allele
                elif len(allele) == length + len(anchor) and allele.startswith(anchor):
                    sequence = allele[len(anchor):]
                else:
                    raise ValueError("Sequence/SVLEN mismatch: " + identifier)
            yield {"chrom": chrom, "start": start, "end": end, "kind": kind,
                   "id": "sniffles.{}.{}".format(kind, identifier),
                   "support": int(info["SUPPORT"]),
                   "precision": "PRECISE" if "PRECISE" in info else "IMPRECISE",
                   "sequence": sequence, "rnames": str(info.get("RNAMES", "")).split(",")}


def selected(record, args):
    excluded = {value.strip("<>") for value in args.type.split(",")}
    return (record["kind"] not in excluded and
            any(re.search(pattern, record["chrom"]) for pattern in args.chrom.split(",")))


def extract_sequences(args):
    count = 0
    skipped = 0
    with open(args.fasta_out, "w") as output:
        for record in records(args.vcf_file):
            if not selected(record, args):
                continue
            sequence = record["sequence"]
            if sequence is None:
                skipped += 1
                continue
            sequence = sequence.upper()
            if not args.keep and set(sequence) <= {"N"}:
                continue
            if len(sequence) <= args.min_len_seq or (args.max_len_seq != -1 and len(sequence) > args.max_len_seq):
                continue
            # Preserve unknown nucleotides in the new backend.
            header = [record["chrom"], "<" + record["kind"] + ">", str(record["start"]),
                      str(record["end"]), record["id"], str(record["support"]), record["precision"]]
            output.write(">" + ":".join(header) + "\n" + sequence + "\n")
            count += 1
    print("Sniffles 2: {} sequences; {} symbolic alleles without sequence skipped".format(count, skipped))


def extract_read_names(args):
    ids = set(Path(args.id_sv).read_text().splitlines()) if args.id_sv else None
    directory = Path(args.directory_name)
    directory.mkdir(parents=True, exist_ok=True)
    for previous in directory.glob("reads_*.txt"):
        if previous.is_file():
            previous.unlink()
    for record in records(args.vcf_file):
        if not selected(record, args) or (ids is not None and record["id"] not in ids):
            continue
        if re.search(r"[/\\]", record["chrom"]):
            raise ValueError("Unsafe chromosome for read-name filename")
        path = directory / "reads_{}:{}:{}-{}.txt".format(
            record["chrom"], record["id"], record["start"], record["end"])
        path.write_text("".join(name + "\n" for name in record["rnames"] if name and name != "."))
