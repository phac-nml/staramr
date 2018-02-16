import abc
import os
import sys

from Bio.Blast import NCBIXML

class BlastResultsParser:

    def __init__(self, file_blast_map, pid_threshold, plength_threshold):
        __metaclass__ = abc.ABCMeta
        self._file_blast_map = file_blast_map
        self._pid_threshold = pid_threshold
        self._plength_threshold = plength_threshold

    def _parse_blast_results(self, file_handle):
        for file in self._file_blast_map:
            databases = self._file_blast_map[file]
            for database_name, blast_out in databases.items():
                if (not os.path.exists(blast_out)):
                    raise Exception("Blast output [" + blast_out + "] does not exist")
                self._handle_blast_hit(file, database_name, blast_out, file_handle)

    def _handle_blast_hit(self, in_file, database_name, blast_file, file_handle):
        blast_handle = open(blast_file)
        blast_records = NCBIXML.parse(blast_handle)
        for blast_record in blast_records:
            hits = []
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    hit = self._create_hit(in_file,alignment,hsp)
                    if (hit.get_pid() > self._pid_threshold and hit.get_plength() > self._plength_threshold):
                        hits.append(hit)
            # sort by pid and then by plength
            hits.sort(key=lambda x: (x.get_pid(), x.get_plength()), reverse=True)
            if (len(hits) >= 1):
                hit = hits[0]
                self._print_hit(hit, file_handle)
        blast_handle.close()

    @abc.abstractmethod
    def _create_hit(self, alignment, hsp):
        pass

    @abc.abstractmethod
    def _print_hit(self, hit, file_handle):
        pass

    def print_to_file(self, file=None):
        file_handle = sys.stdout

        if not file:
            file_handle = open(file, 'w')

        file_handle.write("File\tResistance gene\tPhenotype\t% Identity\t% Overlap\tHSP/Alignment\tContig\n")
        self._parse_blast_results(file_handle)

        if not file:
            file_handle.close()
