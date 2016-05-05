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

    parser = argparse.ArgumentParser()

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

    return (load_genome for p in paths if '.f' in paths)

def validate_ingroup(ingroup: list):
    """Use a set of high quality genomes of the target species
    to look for any spurious identification of 'foreign' sequence
    """

    for genome in load_genomes(ingroup):

        random.seed(1)  # reproducibility

        # run contigified genomes; expected results are no foreign loci
        genome_contigs = contigify(genome, mean, stdev)


def validate_outgroup(outgroup: list):
    """Use a set of high quality genomes known not to be of the target
    species to validate that fngr correctly identifies foreign sequence
    """
    for genome in load_genome(outgroup):
        # run as-is; expected results are nearly all foreign
        pass
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

def select_subsequence(sequence: str, mean: float, stdev: float) -> str:
    """Return a gene-like subsequence from a source genome"""

    subseq_length = int(random.gauss(mean, stdev))
    entry = random.randint(0, len(sequence) - subseq_length - 1)

    return sequence[entry:entry + subseq_length]

def integrate(transposon: str, contig: str, breakpoint: int) -> str:
    """Take a gene-like subsequence and integrate it into a target contig"""

    first, last = contig[:breakpoint], contig[breakpoint:]
    return ''.join((first, transposon, last))

def contaminate(contaminant: str, genome: dict) -> dict:
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

    return contaminated_genome

def iterate():
    """Repeat a treatment with a random variable"""
    pass

def main():

    args = arguments()
if __name__ == '__main__':
    main()
