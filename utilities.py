#!/usr/bin/env python3

from Bio import SeqIO
from functools import partial, wraps
from multiprocessing import cpu_count
import argparse
import collections
import glob
import random
import subprocess
import sys

Metadata = collections.namedtuple('Metadata', ['contig', 'start', 'length'])

Result = collections.namedtuple('Results', ['result', 'metadata'])

user_msg = partial(print, file=sys.stderr)

def fasta(f):
    if '.f' not in f:
        msg = 'Requires FASTA files (*.fasta, *.f, *.fna, etc)'
        raise argparse.ArgumentTypeError(msg)
    return f

def load_genomes(paths: list) -> 'generator':

    def load_genome(filepath: str) -> dict:
        """Parse FASTA formatted string"""

        with open(filepath) as f:
            g = {nucl.id: str(nucl.seq) for nucl in SeqIO.parse(f, 'fasta')}
        return g

    return (load_genome(p) for p in paths if '.f' in p)

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

def prepare_genomes(func, contig_mean: float, contig_stdev: float):
    @wraps(func)
    def wrapper(sources: list, recipients: list, *args):

        contigulate = partial(contigify, mean=contig_mean, stdev=contig_stdev)

        sources_ = [contigulate(genome)
                    for genome in load_genomes(sources)]

        recipients_ = [contigulate(genome)
                       for genome in load_genomes(recipients)]

        return func(sources_, recipients_, *args)
    return wrapper

class Fngr(object):

    def __init__(self, prog, organism, kraken_db, nt_db,
                 threshold = 100, fragment = 250, cores = None):

        self.prog = prog
        self.organism = organism
        self.kraken_db = kraken_db
        self.nt_db = nt_db
        self.fragment = fragment
        self.threshold = threshold
        self.cores = cores or cpu_count()

    def fngr(self, assembly):

        assem = self._format_assembly(assembly)

        cmd = ('python3', self.prog,
               '--organism', self.organism,
               '--kraken-database', self.kraken_db,
               '--nt-database' if self.nt_db else '', self.nt_db or '',
               '--threshold', str(self.threshold),
               '--fragment', str(self.fragment),
               '--cores', str(self.cores),
               '-')

        out = subprocess.check_output([x for x in cmd if x], input = assem,
                                      universal_newlines = True)

        return out  # JSON

    def _format_assembly(self, assembly):

        def dict_to_fasta(assem_dict):

            out = []
            for name, seq in assem_dict.items():

                out.extend([">{}".format(name), seq])

            return '\n'.join(out)

        if type(assembly) is dict:
            assem = dict_to_fasta(assembly)
        elif type(assembly) is str:
            assem = assembly
        else:
            raise TypeError('assembly is not of type str or dict')

        return assem
