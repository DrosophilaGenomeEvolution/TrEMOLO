"""SAM status of the exact genomic alignment used for sequence extraction."""

ALIGNMENT_COLUMNS = ("genomic_alignment_flag", "genomic_alignment_status")


def alignment_metadata(flag=None, sequence_origin="read_sequence"):
    if sequence_origin == "vcf_reported":
        if flag not in (None, "", "."):
            raise ValueError("A VCF allele has no single genomic alignment")
        return dict(zip(ALIGNMENT_COLUMNS, (".", "not_applicable")))
    if flag in (None, "", "."):
        return dict(zip(ALIGNMENT_COLUMNS, (".", "unknown")))
    number = int(flag)
    if not 0 <= number <= 65535:
        raise ValueError("Invalid SAM flag")
    if number & 4:
        status = "unmapped"
    elif number & 256 and number & 2048:
        status = "secondary_and_supplementary"
    elif number & 256:
        status = "secondary"
    elif number & 2048:
        status = "supplementary"
    else:
        status = "primary"
    return dict(zip(ALIGNMENT_COLUMNS, (str(number), status)))
