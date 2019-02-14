import abc
import logging
import os
from typing import List

import Bio.SeqIO
import numpy as np
import pandas as pd

from staramr.blast.BlastHandler import BlastHandler
from staramr.blast.results.BlastHitPartitions import BlastHitPartitions

logger = logging.getLogger('BlastResultsParser')

"""
Class for parsing BLAST results.
"""


class BlastResultsParser:
    INDEX = 'Isolate ID'  # type: str
    COLUMNS = []  # type: List[str]
    SORT_COLUMNS = []  # type: List[str]
    BLAST_SORT_COLUMNS = [x.strip() for x in '''
    plength
    pident
    sstart
    '''.strip().split('\n')]  # type: List[str]

    def __init__(self, file_blast_map, blast_database, pid_threshold, plength_threshold, report_all=False,
                 output_dir=None, genes_to_exclude=[]):
        """
        Creates a new class for parsing BLAST results.
        :param file_blast_map: A map/dictionary linking input files to BLAST results files.
        :param blast_database: The particular staramr.blast.AbstractBlastDatabase to use.
        :param pid_threshold: A percent identity threshold for BLAST results.
        :param plength_threshold: A percent length threshold for results.
        :param report_all: Whether or not to report all blast hits.
        :param output_dir: The directory where output files are being written.
        :param genes_to_exclude: A list of gene IDs to exclude from the results.
        """
        __metaclass__ = abc.ABCMeta
        self._file_blast_map = file_blast_map
        self._blast_database = blast_database
        self._pid_threshold = pid_threshold
        self._plength_threshold = plength_threshold
        self._report_all = report_all
        self._output_dir = output_dir
        self._genes_to_exclude = genes_to_exclude

    def parse_results(self):
        """
        Parses the BLAST files passed to this particular object.
        :return: A pd.DataFrame containing the AMR matches from BLAST.
        """
        results = []

        for file in self._file_blast_map:
            databases = self._file_blast_map[file]
            hit_seq_records = []
            for database_name, blast_out in sorted(databases.items()):
                logger.debug(str(blast_out))
                if (not os.path.exists(blast_out)):
                    raise Exception("Blast output [" + blast_out + "] does not exist")
                self._handle_blast_hit(file, database_name, blast_out, results, hit_seq_records)

            if self._output_dir:
                out_file = self._get_out_file_name(file)
                if hit_seq_records:
                    logger.debug("Writting hits to %s", out_file)
                    Bio.SeqIO.write(hit_seq_records, out_file, 'fasta')
                else:
                    logger.debug("No hits found, skipping writing output file to %s", out_file)
            else:
                logger.debug("No output directory defined for blast hits, skipping writing file")

        return pd.DataFrame(results, columns=self.COLUMNS).sort_values(by=self.SORT_COLUMNS).set_index(self.INDEX)

    @abc.abstractmethod
    def _get_out_file_name(self, in_file):
        """
        Gets hits output file name given input file.
        :param in_file: The input file name.
        :return: The output file name.
        """
        pass

    def _handle_blast_hit(self, in_file, database_name, blast_file, results, hit_seq_records):
        blast_table = pd.read_csv(blast_file, sep='\t', header=None, names=BlastHandler.BLAST_COLUMNS,
                                  index_col=False).astype(
            dtype={'qseqid': np.unicode_, 'sseqid': np.unicode_})
        partitions = BlastHitPartitions()

        blast_table['plength'] = (blast_table.length / blast_table.qlen) * 100.0
        blast_table = blast_table[
            (blast_table.pident >= self._pid_threshold) & (blast_table.plength >= self._plength_threshold) &
            ~blast_table.qseqid.isin(self._genes_to_exclude)]
        blast_table.sort_values(by=self.BLAST_SORT_COLUMNS, inplace=True)
        for index, blast_record in blast_table.iterrows():
            partitions.append(self._create_hit(in_file, database_name, blast_record))

        for hits_non_overlapping in partitions.get_hits_nonoverlapping_regions():
            for hit in self._select_hits_to_include(hits_non_overlapping):
                blast_results = self._get_result_rows(hit, database_name)
                if blast_results is not None:
                    logger.debug("record = %s", blast_results)
                    results.extend(blast_results)
                    hit_seq_records.append(hit.get_seq_record())

    def _select_hits_to_include(self, hits):
        hits_to_include = []

        if len(hits) >= 1:
            sorted_hits_pid_first = sorted(hits, key=lambda x: (
                x.get_pid(), x.get_plength(), x.get_amr_gene_length(), x.get_amr_gene_id()), reverse=True)
            sorted_hits_length_first = sorted(hits, key=lambda x: (
                x.get_amr_gene_length(), x.get_pid(), x.get_plength(), x.get_amr_gene_id()), reverse=True)

            if self._report_all:
                hits_to_include = sorted_hits_pid_first
            else:
                first_hit_pid = sorted_hits_pid_first[0]
                first_hit_length = sorted_hits_length_first[0]

                if first_hit_pid == first_hit_length:
                    hits_to_include.append(first_hit_length)
                # if the top length hit is significantly longer, and the pid is not too much below the top pid hit (nor percent overlap too much below top pid hit), use the longer hit
                elif (first_hit_length.get_amr_gene_length() - first_hit_pid.get_amr_gene_length()) > 10 and (
                        first_hit_length.get_pid() - first_hit_pid.get_pid()) > -1 and (
                        first_hit_length.get_plength() - first_hit_pid.get_plength()) > -1:
                    hits_to_include.append(first_hit_length)
                # otherwise, prefer the top pid hit, even if it's shorter than the longest hit
                else:
                    hits_to_include.append(first_hit_pid)

        return hits_to_include

    @abc.abstractmethod
    def _create_hit(self, file, database_name, blast_record):
        pass

    @abc.abstractmethod
    def _get_result_rows(self, hit, database_name):
        pass
