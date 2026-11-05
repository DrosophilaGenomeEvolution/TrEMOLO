"""Attach the depth actually used by reported frequency estimates; never infer it."""
import math

from build_te_infos import index_frequency, read_rows, valid_percentage
from restore_te_names import read_names


def attach_frequency_depth(calls, outsider=None, precise=None, insider=None, te_headers=None):
    """Join internal sidecars to public calls, mirroring TE_INFOS fallback rules.

    OUTSIDER depth is an adjusted frequency denominator, not raw BAM coverage.
    INSIDER depth is the integer mean of two sampled flanking base depths.
    """
    def frequency_index(path):
        return index_frequency(read_rows(path), str(path)) if path else {}

    raw = frequency_index(outsider)
    corrected = frequency_index(precise)
    names = read_names(te_headers) if te_headers else {}
    resident = {}
    if insider:
        for row in read_rows(insider, skip_header=True):
            if len(row) != 8:
                raise ValueError("{}: expected 8 INSIDER depth columns".format(insider))
            family, separator, event = row[7].partition("|")
            if not separator:
                raise ValueError("{}: invalid INSIDER TE|event".format(insider))
            # Mirror the first-row choice used by build_te_infos.
            resident.setdefault((names.get(family, family), event), row)

    def number(value):
        try:
            value = float(value)
            return value if math.isfinite(value) else None
        except (TypeError, ValueError):
            return None

    for call in calls:
        call.update(frequency_depth=None, frequency_support_reads=None,
                    frequency_depth_method=None, frequency_depth_table=None,
                    frequency_depth_status="unavailable")
        frequency = call['display_frequency']
        if frequency is None:
            call['frequency_depth_status'] = 'frequency_unavailable'
            continue
        row = None
        if call['source'] == 'OUTSIDER':
            row = corrected.get(call['event_id'])
            table = 'precise'
            if row is None or not (valid_percentage(row[8]) and valid_percentage(row[9])):
                row = raw.get(call['event_id'])
                table = 'raw'
            if row is not None:
                with_clipped = call['frequency_with_clipped'] is not None
                depth, support, estimate = (row[7], row[6], row[9]) if with_clipped else (row[5], row[3], row[8])
                # DEL frequency always measures TE presence.
                if len(row) > 10 and row[10] == 'DEL':
                    support = row[6]
                chrom = row[0].rpartition(':')[0]
                method = 'adjusted_read_denominator'
        elif call['source'] == 'INSIDER':
            row = resident.get((call['family'], call['event_id']))
            if row is not None:
                chrom, depth, support, estimate = row[0], row[2], row[4], row[5]
                table, method = 'insider_flanks', 'mean_flanking_depth'
        if row is None:
            continue
        if chrom != call['chrom']:
            call['frequency_depth_status'] = 'chromosome_mismatch'
            continue
        depth, support, estimate = number(depth), number(support), number(estimate)
        if (estimate is None or not math.isclose(estimate, frequency, rel_tol=1e-5, abs_tol=1e-4)):
            call['frequency_depth_status'] = 'frequency_mismatch'
            continue
        call.update(frequency_depth=depth, frequency_support_reads=support,
                    frequency_depth_method=method, frequency_depth_table=table)
        if (depth is None or support is None or depth < 0 or support < 0
                or not depth.is_integer() or not support.is_integer() or support > depth):
            status = 'invalid_counts'
        elif depth == 0:
            status = 'zero_depth'
        elif not 0 <= frequency <= 100:
            status = 'invalid_frequency'
        else:
            status = 'matched'
        call['frequency_depth_status'] = status
