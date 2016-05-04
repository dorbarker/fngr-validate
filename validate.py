#!/usr/bin/env python3

from Bio import SeqIO
from itertools import dropwhile
import argparse
import collections
import compare
import os
import random
import sys

def arguments():

    def fasta(f):
        if '.f' not in f:
            msg = 'Requires FASTA files (*.fasta, *.f, *.fna, etc)'
            raise argparse.ArgumentTypeError(msg)
        return f

    parser = argparse.ArgumentParser()

    data = parser.add_argument_group('Data', 'Input FASTA files')

    data.add_argument('--ingroup', nargs='+', type=fasta, metavar='FASTAs',
                      help='FASTA file(s) belonging to the target species')

    data.add_argument('--outgroup', nargs='+', type=fasta, metavar='FASTAs',
                      help='FASTA file(s) belonging to a different species')

    return parser.parse_args()

def load_genomes(paths: list) -> dict:

    def load_genome(filepath: str) -> dict:
        """Parse FASTA formatted string"""

        with open(filepath) as f:
            g = {contig.id: str(contig.seq) for contig in SeqIO.parse(f, 'fasta')}
        return g

    return (load_genome for p in paths if '.f' in paths)

def validate_ingroup(ingroup: list):
    """Use a set of high quality genomes of the target species
    to look for any spurious identification of 'foreign' sequence
    """

    for genome in load_genomes(ingroup):

        random.seed(1)
        # run contigified genomes; expected results are no foreign
        genome_contigs = contigify(genome, mean, stdev)

def validate_outgroup(outgroup: list):
    """Use a set of high quality genomes known not to be of the target
    species to validate that fngr correctly identifies foreign sequence
    """
    for genome in load_genome(outgroup):
        # run as-is; expected results are nearly all foreign

def contigify(sequence: str, mean: float, stdev: float) -> dict:
    """Cut single-sequence genomes into artificial contigs based on
    empirical distribution contig sizes in draft assemblies
    """

    def chunk_genome():

        processed = 0
        counter = 0

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
