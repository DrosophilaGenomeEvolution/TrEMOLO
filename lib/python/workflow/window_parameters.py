"""Validated genomic distances shared by the migrated workflow."""

DEFAULT_WINDOWS = {
    "OUTSIDER_INSERTION_CLUSTER_DISTANCE": 50,
    "OUTSIDER_SOFT_CLIP_WINDOW": 50,
    "OUTSIDER_HARD_CLIP_CLUSTER_DISTANCE": 100,
    "OUTSIDER_EVIDENCE_CLUSTER_DISTANCE": 100,
    "OUTSIDER_EVIDENCE_MATCH_WINDOW": 100,
    "OUTSIDER_SNIFFLES_READ_MATCH_WINDOW": 30,
    "OUTSIDER_SNIFFLES_READ_FETCH_WINDOW": 1,
    "OUTSIDER_SNIFFLES_READ_FLANK_SIZE": 30,
    "OUTSIDER_CLIPPED_SEQUENCE_PADDING": 100,
    "OUTSIDER_FREQUENCY_WINDOW": 50,
    "OUTSIDER_FREQUENCY_REGION_MERGE_DISTANCE": 1,
    "INSIDER_FREQUENCY_OUTER_MARGIN": 100,
    "INSIDER_FREQUENCY_INNER_MARGIN": 5,
    "INSIDER_FREQUENCY_DEPTH_MARGIN": 30,
    "INSIDER_FREQUENCY_SIZE_WINDOW": 30,
    "INSIDER_FREQUENCY_BREAKPOINT_DISTANCE": 30,
    "TE_INFOS_DELETION_EXCLUSION_WINDOW": 30,
}


def resolve_windows(overrides=None):
    if overrides is None:
        overrides = {}
    if not isinstance(overrides, dict):
        raise ValueError("PARAMS.WINDOWS must be a mapping")
    unknown = set(overrides) - set(DEFAULT_WINDOWS)
    if unknown:
        raise ValueError("Unknown PARAMS.WINDOWS option(s): " + ", ".join(sorted(map(str, unknown))))
    values = dict(DEFAULT_WINDOWS, **overrides)
    for name, value in values.items():
        if type(value) is not int or value < 0:
            raise ValueError("PARAMS.WINDOWS.{} must be a non-negative integer (bp)".format(name))
    return values
