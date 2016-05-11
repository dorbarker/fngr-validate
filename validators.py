from functools import partial
from itertools import dropwhile
from utilities import prepare_genomes, Metadata, Result, user_msg
import random

def validate(sources: list, recipients: list,
             gene_mean: float, gene_stdev: float,
             contig_mean: float, contig_stdev: float,
             iterations: int, fngr: object):

    duplicate_and_contigify = partial(prepare_genomes,
                                      contig_mean=contig_mean,
                                      contig_stdev=contig_stdev)

    @duplicate_and_contigify
    def validate_groups(sources: list, recipients: list) -> list:
        """Use a set of high quality genomes to look for any spurious
        identification of 'foreign' sequence
        """

        metadata = [Metadata(None, None, None)]

        src_validated = [Result(fngr.fngr(src), metadata)
                         for src in sources]

        recip_validated = [Result(fngr.fngr(rec), metadata)
                           for rec in recipients]

        return src_validated, recip_validated

    @duplicate_and_contigify
    def validate_insertions(sources: list, recipients: list, mean: float,
                            stdev: float, iterations: int) -> list:

        def fngr_insertion(source: dict, recipient: dict) -> (dict, str,
                                                              int, int):

            source_contig = random.choice(list(source.values()))
            recipient_contig_name = random.choice(list(recipient.keys()))

            recipient_contig = recipient[recipient_contig_name]

            transposon, source_entry = select_subsequence(source_contig,
                                                          mean, stdev)

            pivot = random.randint(0, len(recipient_contig))

            recipient[recipient_contig_name] = integrate(transposon,
                                                         recipient_contig,
                                                         pivot)

            return recipient, recipient_contig_name, pivot, len(transposon)

        def prepare_insert_func(sources: list, mean: float, stdev: float,
                                iterations: int) -> 'function':

            def _func(recipient: dict, seed: int) -> (dict, (str, int, int)):

                random.seed(seed)

                source = random.choice(sources)

                insertion_metadata = []

                for _ in range(iterations):
                    recipient, *metadata = fngr_insertion(source, recipient)

                    insertion_metadata.append(Metadata(*metadata))

                return Result(fngr.fngr(recipient), insertion_metadata)

            return _func

        insert = prepare_insert_func(sources, mean, stdev, iterations)

        return [insert(recipient, seed)
                for seed, recipient in enumerate(recipients)]

    def select_subsequence(sequence: str,
                           mean: float, stdev: float) -> (str, int):
        """Return a gene-like subsequence from a source genome"""

        subseq_length = int(random.gauss(mean, stdev))
        entry = random.randint(0, len(sequence) - subseq_length - 1)

        return sequence[entry:entry + subseq_length], entry

    def integrate(transposon: str, contig: str, pivot: int) -> str:
        """Take a gene-like subsequence and
        integrate it into a target contig
        """

        first, last = contig[:pivot], contig[pivot:]
        return ''.join((first, transposon, last))

    @duplicate_and_contigify
    def validate_contamination(sources: list, recipients: list,
                               iterations: int) -> list:

        def contaminate_func(sources: list, iterations: int) -> 'function':

            def _func(recipient: dict, seed: int) -> (dict, (str, int)):

                random.seed(seed)

                source = random.choice(sources)

                contamination_metadata = []

                for _ in range(iterations):

                    contaminant = random.choice(list(source.values()))

                    recipient, contig, length = contaminate(contaminant,
                                                            recipient)

                    contamination_metadata.append(Metadata(contig, 0, length))

                return Result(fngr.fngr(recipient), contamination_metadata)
            return _func

        add_contamination = contaminate_func(sources, iterations)

        return [add_contamination(r, s) for s, r in enumerate(recipients)]

    def contaminate(contaminant: str, genome: dict) -> (dict, str, int):
        """Add a contamination contig to genome"""

        def suffix_max(dictionary: dict) -> int:

            k = dictionary.keys()

            not_int = lambda x: x not in map(str, range(10))

            suffix = lambda z: int(''.join(dropwhile(not_int, z)) or 0)

            return max([suffix(key) for key in k] or [0])

        name = 'contamination_{}'.format(suffix_max(genome) + 1)
        genome[name] = contaminant

        return genome, name, len(contaminant)

    user_msg('Validating unmodified genomes')
    source_validated, recipient_validated = validate_groups(sources,
                                                            recipients)
    user_msg('Validating insertions')
    insertion_validated = validate_insertions(sources, recipients,
                                               gene_mean, gene_stdev,
                                               iterations)

    user_msg('Validating contamination')
    contamination_validated = validate_contamination(sources, recipients,
                                                      iterations)

    return {'sources_clean': source_validated,
            'recipients_clean': recipient_validated,
            'insertions': insertion_validated,
            'contamination': contamination_validated}
