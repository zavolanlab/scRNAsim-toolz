"""Main module for read sequencer."""
from random import choices
from collections.abc import Generator, Iterator
from Bio import SeqIO  # type: ignore
from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore


class ReadSequencer:
    """ReadSequencer class.

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
        fasta=None,
        output=None,
        read_length: int = 150,
        chunk_size: int = 10000,
    ) -> None:
        """Initialise class."""
        self.fasta = fasta
        self.output = output
        self.read_length = read_length
        self.chunk_size = chunk_size
        self.random = False
        self.bases = ("A", "T", "C", "G")
        self.n_sequences: int

    def get_n_sequences(self) -> None:
        """Detect number of sequences present in set fasta file.

        Returns:
             None
        """
        self.n_sequences = len(list(SeqIO.parse(self.fasta, "fasta")))

    def define_random_sequences(self, n_seq: int) -> None:
        """Define random sequences.

        Args:
             n_seq: number of random sequences to be generated

        Returns:
            None
        """
        self.random = True
        self.n_sequences = n_seq

    def generate_random_sequence(self, length: int) -> Seq:
        """Generate random sequence.

        Args:
            length: length of sequence

        Returns:
            random sequence of length n
        """
        seq = choices(self.bases, k=length)
        seq = Seq("".join(seq))
        return seq

    def resize_sequence(self, record: SeqRecord) -> SeqRecord:
        """Resize sequence.

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
        """Generate batch iterator.

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
        """Run sequencing.

        Runs read sequencing of specified sequences from input fasta file or
         generates random sequences for a given read length. If number of
         sequences exceeds chunk-size, it will switch to batch processing mode.

        Returns:
             Writes processed sequences to output fasta file(s).
        """
        if self.random:
            if self.n_sequences <= self.chunk_size:
                with open(self.output, "w", encoding="utf-8") as output_handle:
                    for i in range(self.n_sequences):
                        record = SeqRecord(
                            self.generate_random_sequence(self.read_length),
                            id="random_seq: " + str(i + 1),
                        )
                        SeqIO.write(record, output_handle, "fasta")
            else:
                batch_generator = self.batch_iterator(
                    iter(range(self.n_sequences)), self.chunk_size
                )
                for i, batch in enumerate(batch_generator):
                    filename = (
                        self.output.replace(".fasta", "") +
                        f"_chunk_{i + 1}.fasta"
                    )
                    with open(
                        filename, "w", encoding="utf-8"
                    ) as output_handle:
                        for j, _ in enumerate(batch):
                            record = SeqRecord(
                                self.generate_random_sequence(
                                    self.read_length
                                ),
                                id="random_seq: " + str(j + 1),
                            )
                            SeqIO.write(record, output_handle, "fasta")
        else:
            if self.n_sequences <= self.chunk_size:
                with open(self.fasta, encoding="utf-8") as input_handle, open(
                    self.output, "w", encoding="utf-8"
                ) as output_handle:
                    for record in SeqIO.parse(input_handle, "fasta"):
                        record.seq = self.resize_sequence(record)
                        SeqIO.write(record, output_handle, "fasta")

            else:
                with open(self.fasta, encoding="utf-8") as file:
                    record_iter = SeqIO.parse(file, "fasta")
                    for i, batch in enumerate(
                        self.batch_iterator(record_iter, self.chunk_size)
                    ):
                        filename = (
                            self.output.replace(".fasta", "") +
                            f"_chunk_{i + 1}.fasta"
                        )
                        for j, record in enumerate(batch):
                            record.seq = self.resize_sequence(record)
                        with open(filename, "w", encoding="utf-8") as handle:
                            SeqIO.write(batch, handle, "fasta")
