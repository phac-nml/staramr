import math

import Bio.Seq


class NucleotideMutationPosition:

    def __init__(self, match_position, database_string, query_string, database_start, database_frame, query_frame):
        self._preconditions(match_position, database_string, query_string, database_start, database_frame, query_frame)

        self._database_start = database_start
        self._database_frame = database_frame
        self._query_frame = query_frame

        if database_frame == 1:
            self._nucleotide_position_database = database_start + match_position
        else:
            self._nucleotide_position_database = database_start - match_position

        self._codon_start_database = math.ceil(self._nucleotide_position_database / 3)
        frame_shift = (self._nucleotide_position_database - 1) % 3

        self._database_codon = self._find_codon(database_string, match_position, database_frame, frame_shift)
        self._query_codon = self._find_codon(query_string, match_position, database_frame, frame_shift)

    def _find_codon(self, string, match_position, frame, frame_shift):
        if frame == 1:
            codon_start_index = match_position - frame_shift
            return string[codon_start_index:(codon_start_index + 3)].upper()
        else:
            codon_end_index = match_position + frame_shift
            return Bio.Seq.reverse_complement(string[(codon_end_index - 3 + 1):(codon_end_index + 1)].upper())

    def _preconditions(self, match_position, database_string, query_string, database_start, database_frame,
                       query_frame):
        self._check_frame(database_frame)
        self._check_frame(query_frame)

    def _check_frame(self, frame):
        if frame not in [1, -1]:
            raise Exception("Error, frame=" + frame + " not in [1, -1].")

    def get_nucleotide_position(self):
        return self._nucleotide_position_database

    def get_codon_start(self):
        return self._codon_start_database

    def get_database_codon(self):
        return self._database_codon

    def get_database_amino_acid(self):
        if '-' in self.get_database_codon():
            return 'X'
        else:
            return Bio.Seq.translate(self.get_database_codon(), table='Standard')

    def get_query_amino_acid(self):
        if '-' in self.get_query_codon():
            return 'X'
        else:
            return Bio.Seq.translate(self.get_query_codon(), table='Standard')

    def get_query_codon(self):
        return self._query_codon

    def __repr__(self):
        return "[database_start=" + str(self._database_start) + ", database_frame=" + str(
            self._database_frame) + ", query_frame=" + str(self._query_frame) + ", nucleotide_position=" \
               + str(self._nucleotide_position_database) + ", codon_start=" + str(self._codon_start_database) \
               + ", codon=" + self._database_codon + "]"
