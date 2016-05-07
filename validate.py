#!/usr/bin/env python3

from Bio import SeqIO
from itertools import dropwhile
from utilities import fasta, Fngr
import argparse
import collections
import compare
import os
import random
import sys

def arguments():

    docs = {'description': 'A script for validating the output of fngr.py',
            'epilog': 'Output is returned in JSON format to stdout'}

    parser = argparse.ArgumentParser(**docs)

    data = parser.add_argument_group(title='Data',
                                     description='Input FASTA files')

    data.add_argument('--ingroup', nargs='+', type=fasta,
                      metavar='FASTAs', required=True,
                      help='FASTA file(s) belonging to the target species')

    data.add_argument('--outgroup', nargs='+', type=fasta,
                      metavar='FASTAs', required=True,
                      help='FASTA file(s) serving as the source \
                              of foreign sequence to the ingroup')

    parms = parser.add_argument_group(title='Analysis Parameters',
                                      description='')

    parms.add_argument('--organism', type=str.lower,
                       metavar='TAXONOMY', required=True,
                       help='Taxonomic name of target species')

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

    parms.add_argument('--contaminants', type=int,
                       metavar='INT', required=True,
                       help='Number of contaminating contigs \
                            to add to each genome')

    parms.add_argument('--integrations', type=int,
                       metavar='INT', required=True,
                       help='Number of transposon-like integrations of \
                            foreign nucleotide sequence per recipient genome')

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

def load_genomes(paths: list) -> dict:

    def load_genome(filepath: str) -> dict:
        """Parse FASTA formatted string"""

        with open(filepath) as f:
            g = {nucl.id: str(nucl.seq) for nucl in SeqIO.parse(f, 'fasta')}
        return g

    return (load_genome(p) for p in paths if '.f' in p)


def validate_group(group: list, contig_mean: float,
                   contig_stdev: float) -> list:
    """Use a set of high quality genomes to look for any spurious
    identification of 'foreign' sequence
    """

    def fngr_contigs(genome, seed):

        random.seed(seed)
        return fngr.fngr(contigify(genome, contig_mean, contig_stdev))

    return [fngr_contigs(g, s) for s, g in enumerate(load_genomes(group))]

def validate_insertions(sources: list, recipients: list,
                        mean: float, stdev: float, iterations: int) -> list:

    def fngr_insertion(source: dict, recipient: dict) -> (dict, str, int, int):

        source_contig = random.choice(source.values())
        recipient_contig_name = random.choice(recipient.keys())

        recipient_contig = recipient[recipient_contig_name]

        transposon, source_entry = select_subsequence(sequence=source_contig,
                                                      mean=mean, stdev=stdev)

        pivot = random.randint(0, len(recipient_contig))

        recipient[recipient_contig_name] = integrate(transposon,
                                                     recipient_contig, pivot)

        return recipient, recipient_contig_name, pivot, len(transposon)

    def prepare_insert_func(sources: list, mean: float, stdev: float,
                            iterations: int) -> 'function':

        def _func(recipient: dict, seed: int) -> (dict, (str, int, int)):

            random.seed(seed)

            source = random.choice(sources)

            insertion_metadata = []

            Metadata = collections.namedtuple('Metadata',
                                              ['contig', 'pivot', 'length'])

            for _ in range(iterations):
                recipient, contig, pivot, length = fngr_insertion(source,
                                                                  recipient)

                insertion_metadata.append(Metadata(contig, pivot, length))

            return recipient, insertion_metadata

        return _func

    # eliminate side effects of dict alterations being passed back
    # up to other areas of the script
    genomes = [dict(recipient.items()) for recipient in recipients]

    insert = prepare_insert_func(sources, mean, stdev, iterations)

    return [insert(recipient, seed) for seed, recipient in enumerate(genomes)]

def contigify(genome: dict, mean: float, stdev: float) -> dict:
    """Cut genomes into artificial contigs based on
    empirical distribution contig sizes in draft assemblies
    """

    def chunk_genome():

        counter = 0

        for sequence in genome.values():

            processed = 0

            while processed < len(sequence):

                counter += 1

                n = int(random.gauss(mean, stdev))

                chunk = sequence[processed:processed + n]

                processed += n

                yield 'contig_{:04d}'.format(counter), chunk

    return dict(chunk_genome())

def select_subsequence(sequence: str, mean: float, stdev: float) -> (str, int):
    """Return a gene-like subsequence from a source genome"""

    subseq_length = int(random.gauss(mean, stdev))
    entry = random.randint(0, len(sequence) - subseq_length - 1)

    return sequence[entry:entry + subseq_length], entry

def integrate(transposon: str, contig: str, pivot: int) -> str:
    """Take a gene-like subsequence and integrate it into a target contig"""

    first, last = contig[:pivot], contig[pivot:]
    return ''.join((first, transposon, last))

def contaminate_genomes(sources: list, recipients: list, contig_mean: float,
                        contig_stdev: float, interations: int):

    def contaminate_func(sources: list, mean: float,
                         stdev: float, iterations: int) -> 'function':

        def _func(recipient: dict, seed: int) -> (dict, (str, int)):

            random.seed(seed)

            source = random.choice(sources)

            contamination_metadata = []

            Metadata = collections.namedtuple('Metadata', ['contig', 'length'])

            for _ in range(iterations):

                contaminant = random.choice(source.values())

                recipient, contig, length = contaminate(contaminant, recipient)

                contamination_metadata.append(Metadata(contig, length))

            return recipient, contamination_metadata
        return _func

    genomes = [dict(recipient.items()) for recipient in recipients]

    add_contamination = contaminate_func(sources, mean, stdev, interations)

    return [add_contamination(r, s) for s, r in enumerate(recipients)]

def contaminate(contaminant: str, genome: dict) -> (dict, str, int):
    """Add a contamination contig to genome"""

    def suffix_max(dictionary: dict) -> int:

        k = dictionary.keys()

        not_int = lambda x: x not in map(str, range(10))

        suffix = lambda z: int(''.join(dropwhile(not_int, z)) or 0)

        return max([suffix(key) for key in k] or [0])

    # not strictly necessary, but makes the program easier to follow
    contaminated_genome = dict(genome.items())

    name = 'contamination_{}'.format(suffix_max(contaminated_genome) + 1)
    contaminated_genome[name] = contaminant

    return contaminated_genome, name, len(contaminant)

def compare_to_expected():

    pass

def iterate():
    """Repeat a treatment with a random variable"""
    pass

def main():

    args = arguments()

    global fngr
    fngr = Fngr(prog=args.fngr, organism=args.organism, cores=args.cores,
                kraken_db=args.kraken_database, nt_db=args.nt_database)

    print(validate_ingroup(args.ingroup, args.contig_mean, args.contig_stdev))

if __name__ == '__main__':
    main()
