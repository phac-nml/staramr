import abc
import logging
import os

from Bio.Blast import NCBIXML

from staramr.blast.results.BlastHitPartitions import BlastHitPartitions

logger = logging.getLogger('BlastResultsParser')

"""
Class for parsing BLAST results.
"""


class BlastResultsParser:

    def __init__(self, file_blast_map, blast_database, pid_threshold, plength_threshold):
        """
        Creates a new class for parsing BLAST results.
        :param file_blast_map: A map/dictionary linking input files to BLAST results files.
        :param blast_database: The particular staramr.blast.AbstractBlastDatabase to use.
        :param pid_threshold: A percent identity threshold for BLAST results.
        :param plength_threshold: A percent length threshold for results.
        """
        __metaclass__ = abc.ABCMeta
        self._file_blast_map = file_blast_map
        self._blast_database = blast_database
        self._pid_threshold = pid_threshold
        self._plength_threshold = plength_threshold

    def parse_results(self):
        """
        Parses the BLAST files passed to this particular object.
        :return: A pandas.DataFrame containing the AMR matches from BLAST.
        """
        results = []

        for file in self._file_blast_map:
            databases = self._file_blast_map[file]
            for database_name, blast_out in databases.items():
                logger.debug(str(blast_out))
                if (not os.path.exists(blast_out)):
                    raise Exception("Blast output [" + blast_out + "] does not exist")
                self._handle_blast_hit(file, database_name, blast_out, results)

        return self._create_data_frame(results)

    def _handle_blast_hit(self, in_file, database_name, blast_file, results):
        blast_handle = open(blast_file)
        blast_records = NCBIXML.parse(blast_handle)
        for blast_record in blast_records:
            partitions = BlastHitPartitions()
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    hit = self._create_hit(in_file, database_name, blast_record, alignment, hsp)
                    if hit.get_pid() > self._pid_threshold and hit.get_plength() > self._plength_threshold:
                        partitions.append(hit)
            for hits_non_overlapping in partitions.get_hits_nonoverlapping_regions():
                # sort by pid and then by plength
                hits_non_overlapping.sort(key=lambda x: (x.get_pid(), x.get_plength()), reverse=True)
                if len(hits_non_overlapping) >= 1:
                    hit = hits_non_overlapping[0]
                    self._append_results_to(hit, database_name, results)
        blast_handle.close()

    @abc.abstractmethod
    def _create_data_frame(self, results):
        pass

    @abc.abstractmethod
    def _create_hit(self, file, database_name, blast_record, alignment, hsp):
        pass

    @abc.abstractmethod
    def _append_results_to(self, hit, database_name, results):
        pass
