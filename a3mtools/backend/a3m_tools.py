from pathlib import Path
import re
import string

class ProteinSequence:
    def __init__(self, header: str, sequence: str):
        self.header = header
        self.seq_str = sequence

    def __str__(self):
        return f">{self.header}\n{self.seq_str}"

    def __repr__(self):
        return f">{self.header}\n{self.seq_str}"

    def __len__(self):
        return len(self.seq_str)

    def __getitem__(self, index):
        """This should allow for slicing of the sequence"""
        return ProteinSequence(self.header, self.seq_str[index])

    def __add__(self, other):
        """This should allow for concatenation of the sequence"""
        if isinstance(other, ProteinSequence):
            return ProteinSequence(self.header, self.seq_str + other.seq_str)
        elif isinstance(other, str):
            return ProteinSequence(self.header, self.seq_str + other)
        else:
            raise TypeError(f"Cannot concatenate ProteinSequence with {type(other)}")

    def __radd__(self, other):
        """This should allow for concatenation of the sequence"""
        if isinstance(other, ProteinSequence):
            return ProteinSequence(self.header, other.seq_str + self.seq_str)
        elif isinstance(other, str):
            return ProteinSequence(self.header, other + self.seq_str)
        else:
            raise TypeError(f"Cannot concatenate ProteinSequence with {type(other)}")


def parse_header(header_line: str):
    header = header_line.rstrip()
    assert header.startswith(
        ">"
    ), f"Expected header to start with >, instead saw: {header_line}"
    return header[1:]


def import_a3m(filepath):
    with open(filepath) as f:
        info_line = f.readline().rstrip()
        if not info_line.startswith("#"):
            raise ValueError(
                f"{filepath} should begin with #, instead saw: {info_line}"
            )
        lines = [i.strip() for i in f.readlines()]
        query_header = parse_header(lines[0])
        query_sequence = lines[1].rstrip()
        query = ProteinSequence(query_header, query_sequence)
        sequences = []
        for c in range(2, len(lines), 2):
            name = parse_header(lines[c])
            seq = lines[c + 1].rstrip()
            sequences.append(ProteinSequence(name, seq))
    return info_line, query, sequences


def parse_info_line(info_line: str):
    """ """
    info = info_line.split("\t")
    assert (
        len(info) == 2
    ), f"Expected 2 fields separated by a tab in the info line, instead saw: {info_line}"
    query_seq_lengths = info[0][1:].split(",")
    query_seq_cardinality = info[1].split(",")
    return query_seq_lengths, query_seq_cardinality


def diagonal_concat(msa_1, msa_2):
    """
    Concatenates two MSAs in the unpaired format
    """
    new_msa_1 = MSAa3m(info_line=msa_1.info_line,
                       query=msa_1.query,
                       sequences=[]
                      )
    new_msa_2 = MSAa3m(info_line=msa_2.info_line,
                       query=msa_2.query,
                       sequences=[]
                      )

    msa_1_no_gaps = []
    for seq in msa_1.sequences:
        if seq.seq_str == len(seq.seq_str) * '-':
            continue
        elif re.search("^10[1-9]$", seq.header): # query, again
            continue
        else:
            msa_1_no_gaps.append(seq)
    new_msa_1.sequences = msa_1_no_gaps

    msa_2_no_gaps = []
    for seq in msa_2.sequences:
        if seq.seq_str == len(seq.seq_str) * '-':
            continue
        elif re.search("^10[1-9]$", seq.header): # query, again
            continue
        else:
            msa_2_no_gaps.append(seq)
    new_msa_2.sequences = msa_2_no_gaps

    # Make msa_1 gap sequences
    for seq in new_msa_2.sequences:
        gap_seq = ProteinSequence(seq.header, '-'*len(new_msa_1.query.seq_str))
        new_msa_1.sequences.append(gap_seq)

    # Make msa_2 gap sequences
    for seq in reversed(new_msa_1.sequences):
        gap_seq = ProteinSequence(seq.header, '-'*len(new_msa_2.query.seq_str))
        new_msa_2.sequences.insert(0, gap_seq)

    return new_msa_1 + new_msa_2


class MSAa3m:

    def __init__(
        self, info_line: str, query: ProteinSequence, sequences: list[ProteinSequence]
    ):
        self.info_line = info_line
        self.query = query
        assert "-" not in self.query.seq_str, "Query sequence should not contain gaps"
        assert all(
            [i.isupper() for i in self.query.seq_str]
        ), "Query sequence should not contain insertions"
        self.sequences = sequences
        self.query_seq_lengths, self.query_seq_cardinality = self._parse_info_line(
            self.info_line
        )

    def __str__(self):
        msa_str = f"{self.info_line}\n{self.query}\n"
        for seq in self.sequences:
            msa_str += f"{seq}\n"
        return msa_str

    def __repr__(self):
        msa_str = f"{self.info_line}\n{self.query}\n"
        for seq in self.sequences:
            msa_str += f"{seq}\n"
        return msa_str

    @classmethod
    def from_a3m_file(cls, file_path: Path | str):
        """Create an MSAa3m object from an a3m file

        Parameters
        ----------
        file_path : Path | str
            path to the a3m file

        Returns
        -------
        MSAa3m
            a a3mtools.MSAa3m object
        """        
        info_line, query, sequences = import_a3m(file_path)
        return cls(info_line, query, sequences)

    @staticmethod
    def _parse_info_line(info_line: str):
        """
        #157,235	1,1
        #398	1
        """
        info = info_line.split("\t")
        assert (
            len(info) == 2
        ), f"Expected 2 fields separated by a tab in the info line, instead saw: {info_line}"
        query_seq_lengths = info[0][1:].split(",")
        query_seq_cardinality = info[1].split(",")
        return query_seq_lengths, query_seq_cardinality

    def _slice_alignment(self, start, end):
        if len(self.query_seq_lengths) > 1:
            raise NotImplementedError(
                "Cannot slice alignment with multiple query sequences yet"
            )
        # handle cases where the start/end are none (i.e. [:-1] or [1:])
        if start is None:
            start = 0
        if end is None:
            end = len(self.query)
        # handle negative indices
        if start < 0:
            start = len(self.query) + start
        if end < 0:
            end = len(self.query) + end
        sliced_query = ProteinSequence(self.query.header, self.query.seq_str[start:end])
        new_info_line = f"#{len(sliced_query.seq_str)}\t{self.query_seq_cardinality[0]}"
        sliced_sequences = []
        for seq in self.sequences:
            core_index = 0
            sliced_seq = ""
            for char in seq.seq_str:
                if core_index >= end:
                    break
                if char.isupper():
                    if core_index >= start:
                        sliced_seq += char
                    core_index += 1
                elif char.islower():
                    if core_index > start:
                        sliced_seq += char
                elif char == "-":
                    if core_index >= start:
                        sliced_seq += char
                    core_index += 1
            sliced_sequences.append(ProteinSequence(seq.header, sliced_seq))
        return MSAa3m(new_info_line, sliced_query, sliced_sequences)

    def __getitem__(self, index):
        """This should allow for slicing of the sequence"""
        if isinstance(index, slice):
            # raise error if the slice is not a step of 1
            if index.step is not None and index.step != 1:
                raise NotImplementedError("Slicing with a step other than 1 is not supported")
            return self._slice_alignment(index.start, index.stop)
        elif isinstance(index, int):
            return self._slice_alignment(index, index + 1)
        else:
            raise TypeError("Only integers and slices are valid indices for A3M_MSA")

    def __add__(self, other):
        """This should allow for concatenation of the msas."""
        if isinstance(other, MSAa3m):
            new_info_line = f"#{','.join(self.query_seq_lengths)},{','.join(other.query_seq_lengths)}\t{','.join(self.query_seq_cardinality)},{','.join(other.query_seq_cardinality)}"

            # change other query header to a new index.
            # create a new object to avoid changing the original
            last_self_int = [int(i) for i in self.query.header.split("\t")][-1]
            other_ints = [int(i) for i in other.query.header.split("\t")]
            c = last_self_int + 1
            other_int_map = {}
            new_other_int_list = []
            for i in other_ints:
                other_int_map[str(i)] = str(c)
                new_other_int_list.append(str(c))
                c += 1
            new_other_header = "\t".join(new_other_int_list)
            other_query = ProteinSequence(new_other_header, other.query.seq_str)

            new_query = ProteinSequence(
                f"{self.query.header}\t{other_query.header}",
                f"{self.query.seq_str}{other_query.seq_str}",
            )
            new_sequences = []
            if len(self.query_seq_lengths) == 1:
                new_sequences.append(self.query + "-" * len(other_query))
            for seq in self.sequences:
                new_seq = seq + "-" * len(other_query)
                if not all([i == "-" for i in new_seq.seq_str]):
                    new_sequences.append(new_seq)
            if len(other.query_seq_lengths) == 1:
                new_sequences.append("-" * len(self.query) + other_query)
            for seq in other.sequences:
                new_seq = "-" * len(self.query) + seq
                if new_seq.header in other_int_map.keys():
                    new_seq.header = other_int_map[new_seq.header]
                if not all([i == "-" for i in new_seq.seq_str]):
                    new_sequences.append(new_seq)
            # add the original 2 queries back to the MSA but only
            # if the MSA was not previously concatenated
            return MSAa3m(new_info_line, new_query, new_sequences)
        else:
            raise TypeError(f"Cannot concatenate MSAa3m with {type(other)}")
    
    def save(self, filepath):
        with open(filepath, "w") as f:
            f.write(str(self))


class PairedMSAa3m(MSAa3m):
    def __init__(self, 
                 info_line: str, 
                 query: ProteinSequence, 
                 sequences: list[ProteinSequence],
                 paired: bool = False
                ):
        super().__init__(info_line, query, sequences)

        self.query_sequences = []
        self.query_idxs = []
        self.paired = paired

        for _, i in zip(string.ascii_uppercase, range(len(self.query_seq_lengths))):
            if i == 0:
                end = int(self.query_seq_lengths[i])
                seq = self.query.seq_str[:end]
                self.query_idxs.append((0, end))
            else:
                start = int(self.query_seq_lengths[i-1])
                end = start + int(self.query_seq_lengths[i])
                seq = self.query.seq_str[start:end]
                self.query_idxs.append((start, end))
            self.query_sequences.append(ProteinSequence(f'10{i+1}', seq))

    
    def __str__(self):
        if self.paired:
            msa_str = f"{self.info_line}\n{self.query}\n{self.query}\n"
        else:
            msa_str = f"{self.info_line}\n{self.query}\n"
        # Write paired region
        for seq in self.sequences:
            msa_str += f"{seq}\n"

        if self.paired:
            # Write unpaired region
            for i, query_seq in enumerate(self.query_sequences):
                msa_str += f'>{query_seq.header}\n'
                pre_padding = ''
                for j in range(0, i):
                    pre_padding += '-'*self.query_seq_lengths[j]
                post_padding = ''
                for j in range(i+1, len(self.query_sequences)):
                    post_padding += '-'*self.query_seq_lengths[j]
                query_seq_msa_row = f"{pre_padding}{query_seq.seq_str}{post_padding}"
                assert len(query_seq_msa_row) == len(self.query.seq_str)
                msa_str += f'{query_seq_msa_row}\n'
        return msa_str


    def _slice_alignment(self, start, end):
        sliced_query = ProteinSequence('101', self.query.seq_str[start:end])
        new_info_line = f"#{len(sliced_query.seq_str)}\t{self.query_seq_cardinality[0]}"
        sliced_sequences = []
        for seq in self.sequences:
            core_index = 0
            sliced_seq = ""
            for char in seq.seq_str:
                if core_index >= end:
                    break
                if char.isupper():
                    if core_index >= start:
                        sliced_seq += char
                    core_index += 1
                elif char.islower():
                    if core_index > start:
                        sliced_seq += char
                elif char == "-":
                    if core_index >= start:
                        sliced_seq += char
                    core_index += 1
            sliced_sequences.append(ProteinSequence(seq.header, sliced_seq))
        return MSAa3m(new_info_line, sliced_query, sliced_sequences)


    def get_msa_by_chain(self, chain_idx: int):
        """
        Returns the MSA associated with the chain index (int)
        """
        return self._slice_alignment(*self.query_idxs[chain_idx])
