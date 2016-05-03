import collections

Pair = collections.namedtuple('Pair', ['start', 'stop'])

def closest_indices(contig_expected):

    indices = lambda locus: (contig_expected[locus]['index']['start'],
                             contig_expected[locus]['index']['stop'])

    pairs = [Pair(*indices(locus)) for locus in contig_expected]

    def _closest(start, stop):

        best_start = min(pairs, key = lambda p: abs(p.start - start))
        best_stop = min(pairs, key = lambda p: abs(p.stop - stop))

        return Pair(best_start, best_stop)

    return _closest

def call_type(locus):

    pass

def compare(expected: dict, report: dict) -> dict:

    comparison = {'calls': {'n_loci': 0,
                            'foreign': 0,
                            'contamination': 0,
                            'false_positive': 0},
                  'loci': collections.defaultdict(list)}

    for contig in report:

        find_closest = closest_indices(expected[contig])

        for locus in report[contig]:

            comparison['calls']['n_loci'] += 1

            found_index = Pair(locus['index']['start'], locus['index']['stop'])

            closest = find_closest(found_index.start, found_index.stop)

            coverage = len(range(*found_index)) / len(range(*closest))

            indices = {'observed': locus['index'],
                       'closest_expected': closest,
                       'coverage': coverage}

            organism = {'expected': expected['organism'],
                        'blast': locus['blast_hits'],
                        'kraken': locus['read_classification']}

            locus_info = {'indices': indices, 'organism': organism}
            locus_info['call'] = call_type(locus_info)
