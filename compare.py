import collections
import json

Pair = collections.namedtuple('Pair', ['start', 'stop'])

def closest_indices(contig_expected):

    def _closest(start, stop):

        if contig_expected:
            start_ = min(contig_expected, key=lambda p: abs(p.start - start))
            stop_ = min(contig_expected, key=lambda p: abs(p.stop - stop))

            out = Pair(start_.start, stop_.stop)
        else:
            out = Pair(None, None)

        return out

    return _closest

def call_type(locus_info):

    cov = locus_info['indices']['coverage']

    if locus_info['contig_proportion'] >= 0.95:
        out = 'contamination'

    elif cov is None:
        out = 'false_positive'
    elif any((cov > 1.1, cov < 0.9,
              None in locus_info['indices']['closest_expected'])):
        out = 'false_positive'

    else:
        out = 'foreign'

    return out

def coverage(found, closest):

    if not found or None in found or None in closest:
        out = None
    else:
        out = (found.stop - found.start) / (closest.stop - closest.start)
    return out

def compare_each_genome(expecteds: dict, results: dict):

    out = {}
    for key in results:

        out[key] = [compare(expect, json.loads(res))
                    for expect, res in zip(expecteds[key], results[key])]

    return out

def compare(expected: dict, report: dict) -> dict:

    comparison = {'calls': {'n_loci': 0,
                            'foreign': 0,
                            'contamination': 0,
                            'false_positive': 0,
                            'false_negative': 0},
                  'loci': collections.defaultdict(list)}

    for contig in filter(None, set(report) | set(expected['contigs'])):

        if contig in report:

            if contig in expected['contigs']:
                comparison = compare_hit(contig, expected, report, comparison)

            else:
                comparison['calls']['n_loci'] += 1
                comparison['calls']['false_positive'] += 1

        else:
            comparison['calls']['n_loci'] += 1
            comparison['calls']['false_negative'] += 1

    return comparison

def compare_hit(contig: str, expected: dict, report: dict, comparison: dict):

    find_closest = closest_indices(expected['contigs'][contig])

    for locus in report[contig]:

        found_index = Pair(locus['index']['start'], locus['index']['stop'])

        closest = find_closest(found_index.start, found_index.stop)

        indices = {'observed': locus['index'],
                   'closest_expected': {'start': closest.start,
                                        'stop': closest.stop},
                   'coverage': coverage(found_index, closest)}

        organism = {'expected': expected['organism'],
                    'blast': locus['blast_hits'],
                    'kraken': locus['read_classification']}

        locus_info = {'indices': indices, 'organism': organism,
                      'contig_proportion': locus['length']['ratio']}

        locus_info['call'] = call_type(locus_info)

        comparison['loci'][contig].append(locus_info)

        comparison['calls']['n_loci'] += 1
        comparison['calls'][locus_info['call']] += 1

    return comparison

def create_expected(results: list, organism: str) -> list:
    def parse_metadata(metadata: list) -> dict:

        expected = {'organism': organism, 'contigs': collections.defaultdict(list)}

        for m in metadata:

            try:
                p = Pair(m.start, m.start + m.length)
            except TypeError:
                p = Pair(None, None)

            expected['contigs'][m.contig].append(p)

        return expected
    print('results:', type(results))
    print('inside_results:', [type(x) for x in results])
    return [parse_metadata(result.metadata) for result in results]

def metacompare(comparisons: dict):

    metacomparison = collections.defaultdict(int)

    for comparison in comparisons:
        for genome in comparisons[comparison]:
            metacomparison['genomes'] += 1
            for key in genome['calls']:
                metacomparison[key] += genome['calls'][key]

    return metacomparison
