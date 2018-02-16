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
        results = []

        for file in self._file_blast_map:
            databases = self._file_blast_map[file]
            for database_name, blast_out in databases.items():
                if (not os.path.exists(blast_out)):
                    raise Exception("Blast output [" + blast_out + "] does not exist")
                self._handle_blast_hit(file, database_name, blast_out, results, file_handle)

        return self._create_data_frame(results)

    def _handle_blast_hit(self, in_file, database_name, blast_file, results, file_handle):
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
            if len(hits) >= 1:
                hit = hits[0]
                self._append_results_to(hit, results)
        blast_handle.close()

    @abc.abstractmethod
    def _create_data_frame(self, results):
        pass

    @abc.abstractmethod
    def _create_hit(self, alignment, hsp):
        pass

    @abc.abstractmethod
    def _append_results_to(self, hit, results):
        pass

    def print_to_file(self, file=None):
        file_handle = sys.stdout

        #if not file:
            #file_handle = open(file, 'w')

        data_frame = self._parse_blast_results(file_handle)
        data_frame.to_csv(file_handle, sep="\t", index=False, float_format="%0.2f")

        #if not file:
         #   file_handle.close()
