from random import choices
from collections.abc import Generator, Iterator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ReadSequencer:
    """ReadSequencer class

    Args:
        fasta: path fasta file
        output: path output fasta file(s)
        read_length: read length, defaults to 150.
        chunk_size: batch size used for memory efficient processing,
         only used when number of sequences greater
         than number of passed sequences. Defaults to 10000.

    Returns:
        None
    """

    def __init__(
        self,
        fasta: str = None,
        output: str = None,
        read_length: int = 150,
        chunk_size: int = 10000,
    ) -> None:

        self.fasta = fasta
        self.output = output
        self.read_length = read_length
        self.chunk_size = chunk_size
        self.random = False
        self.bases = ("A", "T", "C", "G")
        self.n_sequences = None

    def get_n_sequences(self) -> None:
        """
        Helper function to detect number of sequences present in set fasta file.

        Returns:
             None
        """
        self.n_sequences = len(list(SeqIO.parse(self.fasta, "fasta")))

    def define_random_sequences(self, n_seq: int) -> None:
        """
        Defines random sequences.

        Args:
             n_seq: number of random sequences to be generated

        Returns:
            None
        """
        self.random = True
        self.n_sequences = n_seq

    def generate_random_sequence(self, length: int) -> Seq:
        """
        Generates random sequence.

        Args:
            length: length of sequence

        Returns:
            random sequence of length n
        """
        seq = choices(self.bases, k=length)
        seq = Seq("".join(seq))
        return seq

    def resize_sequence(self, record: SeqRecord) -> SeqRecord:
        """Resizes sequence

        Resizes sequence according to set read length. If sequence is
         shorter than read length, fills up with random nucleotides.

        Args:
            record: SeqRecord

        Returns:
            resized SeqRecord
        """
        if (len(record)) >= self.read_length:
            record.seq = record.seq[0:self.read_length - 1]
        else:
            n_add = self.read_length - len(record)
            add_seq = self.generate_random_sequence(n_add)
            record.seq = record.seq + add_seq
        return record.seq

    def batch_iterator(self, iterator: Iterator, batch_size: int) -> Generator:
        """Generates batch iterator.

        This is a generator function, and it returns lists of the
        entries from the supplied iterator.  Each list will have
        batch_size entries, although the final list may be shorter.

        Args:
            iterator: iterator object generated with Bio.SeqIO.parse()
            batch_size: batch size to use for the generator

        Returns:
            list of entries from supplied iterator according to batch_size
        """
        batch = []
        for entry in iterator:
            batch.append(entry)
            if len(batch) == batch_size:
                yield batch
                batch = []

    def run_sequencing(self) -> None:
        """Runs sequencing.

        Runs read sequencing of specified sequences from input fasta file or
         generates random sequences for a given read length. If number of
         sequences exceeds chunk-size, it will switch to batch processing mode.

        Returns:
             Writes processed sequences to output fasta file(s).
        """
        if self.random:
            if self.n_sequences <= self.chunk_size:
                with open(self.output, "w") as output_handle:
                    for i in range(self.n_sequences):
                        record = SeqRecord(
                            self.generate_random_sequence(self.read_length),
                            id="random_seq: " + str(i + 1),
                        )
                        SeqIO.write(record, output_handle, "fasta")
            else:
                batch_generator = self.batch_iterator(
                    range(self.n_sequences), self.chunk_size
                )
                for i, batch in enumerate(batch_generator):
                    filename = self.output.replace(".fasta", "") + "_chunk_%i.fasta" % (
                        i + 1
                    )
                    with open(filename, "w") as output_handle:
                        for j, k in enumerate(batch):
                            record = SeqRecord(
                                self.generate_random_sequence(self.read_length),
                                id="random_seq: " + str(j + 1),
                            )
                            SeqIO.write(record, output_handle, "fasta")
        else:
            if self.n_sequences <= self.chunk_size:
                with open(self.fasta) as input_handle, open(
                    self.output, "w"
                ) as output_handle:
                    for record in SeqIO.parse(input_handle, "fasta"):
                        record.seq = self.resize_sequence(record)
                        SeqIO.write(record, output_handle, "fasta")

            else:
                record_iter = SeqIO.parse(open(self.fasta), "fasta")
                for i, batch in enumerate(
                    self.batch_iterator(record_iter, self.chunk_size)
                ):
                    filename = self.output.replace(".fasta", "") + "_chunk_%i.fasta" % (i + 1)
                    for j, record in enumerate(batch):
                        batch[j].seq = self.resize_sequence(record)
                    with open(filename, "w") as handle:
                        SeqIO.write(batch, handle, "fasta")
