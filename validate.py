#!/usr/bin/env python3

import argparse
import compare
import os
import utilities
import validators

def arguments():

    docs = {'description': 'A script for validating the output of fngr.py',
            'epilog': 'Output is returned in JSON format to stdout'}

    parser = argparse.ArgumentParser(**docs)

    data = parser.add_argument_group(title='Data',
                                     description='Input FASTA files')

    data.add_argument('--ingroup', nargs='+', type=utilities.fasta,
                      metavar='FASTAs', required=True,
                      help='FASTA file(s) belonging to the target species')

    data.add_argument('--outgroup', nargs='+', type=utilities.fasta,
                      metavar='FASTAs', required=True,
                      help='FASTA file(s) serving as the source \
                              of foreign sequence to the ingroup')

    parms = parser.add_argument_group(title='Analysis Parameters',
                                      description='')

    parms.add_argument('--organism', type=str.lower,
                       metavar='TAXONOMY', required=True,
                       help='Taxonomic name of target species')

    parms.add_argument('--foreign', type=str.lower, metavar='TAXONOMY',
                       required=True, help='Taxonomic name of outgroup')

    parms.add_argument('--contig-mean', type=float, metavar='NUM',
                       required=True, help='Mean size of synthetic contigs')

    parms.add_argument('--contig-stdev', type=float,
                       metavar='NUM', required=True,
                       help='Standard deviation in synthetic contig lengths')

    parms.add_argument('--insertion-mean', type=float,
                       metavar='NUM', required=True,
                       help='Mean length of foreign sequence insertions')

    parms.add_argument('--insertion-stdev', type=float,
                       metavar='NUM', required=True,
                       help='Standard deviation in foreign \
                            sequence insertions')

    parms.add_argument('--iterations', type=int, metavar='INT', default=10,
                       help='Iterations of each test to perform [10]')

    parms.add_argument('--fragment', type=int, metavar='INT', default=250,
                       help='Pseudoread size (bp) into which Fngr divides \
                            assemblies for kraken [250]')

    parms.add_argument('--threshold', type=int, metavar='INT', default=100,
                       help='Number of consecutive reads required to \
                             call a sequence foreign [100]')

    resources = parser.add_argument_group(title='Resources',
                                          description='External resources \
                                                      used by Fngr')

    resources.add_argument('--fngr', required=True, metavar='PATH',
                           help='Path to fngr.py')

    resources.add_argument('--kraken-database', required=True, metavar='PATH',
                           help='Path to Kraken database')

    resources.add_argument('--nt-database', default=None, metavar='PATH',
                           help='Path to NCBI `nt` database')

    resources.add_argument('--cores', type=int, default=None, metavar='INT',
                           help='Number of CPU cores to use [all]')

    return parser.parse_args()

def main():

    args = arguments()

    fngr = utilities.Fngr(prog=args.fngr, organism=args.organism,
                          fragment=args.fragment, threshold=args.threshold,
                          kraken_db=args.kraken_database,
                          nt_db=args.nt_database,
                          cores=args.cores)

    results = validators.validate(sources=args.outgroup,
                                  recipients=args.ingroup,
                                  gene_mean=args.insertion_mean,
                                  gene_stdev=args.insertion_stdev,
                                  contig_mean=args.contig_mean,
                                  contig_stdev=args.contig_stdev,
                                  iterations=args.iterations,
                                  fngr=fngr)


    utilities.user_msg('')
    expected = {result: compare.create_expected(results[result], args.foreign)
                for result in results}

    fngr_output = {res: [r.result for r  in results[res]] for res in results}

    comparisons = compare.compare_each_genome(expected, fngr_output)
    import pprint
    pp = pprint.PrettyPrinter(indent =2)
    for i in comparisons:
        pp.pprint(i)
        pp.pprint(comparisons[i])
        print('')

    metacomparison = compare.metacompare(comparisons)
    pp.pprint(metacomparison)
if __name__ == '__main__':
    main()
