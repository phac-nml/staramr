import abc
import logging
import os

import Bio.SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from staramr.blast.results.BlastHitPartitions import BlastHitPartitions

logger = logging.getLogger('BlastResultsParser')

"""
Class for parsing BLAST results.
"""


class BlastResultsParser:

    def __init__(self, file_blast_map, blast_database, pid_threshold, plength_threshold, report_all=False,
                 output_dir=None):
        """
        Creates a new class for parsing BLAST results.
        :param file_blast_map: A map/dictionary linking input files to BLAST results files.
        :param blast_database: The particular staramr.blast.AbstractBlastDatabase to use.
        :param pid_threshold: A percent identity threshold for BLAST results.
        :param plength_threshold: A percent length threshold for results.
        :param report_all: Whether or not to report all blast hits.
        :param output_dir: The directory where output files are being written.
        """
        __metaclass__ = abc.ABCMeta
        self._file_blast_map = file_blast_map
        self._blast_database = blast_database
        self._pid_threshold = pid_threshold
        self._plength_threshold = plength_threshold
        self._report_all = report_all
        self._output_dir = output_dir

    def parse_results(self):
        """
        Parses the BLAST files passed to this particular object.
        :return: A pandas.DataFrame containing the AMR matches from BLAST.
        """
        results = []

        for file in self._file_blast_map:
            databases = self._file_blast_map[file]
            hit_seq_records = []
            for database_name, blast_out in databases.items():
                logger.debug(str(blast_out))
                if (not os.path.exists(blast_out)):
                    raise Exception("Blast output [" + blast_out + "] does not exist")
                self._handle_blast_hit(file, database_name, blast_out, results, hit_seq_records)

            if self._output_dir:
                out_file = self._get_out_file_name(file)
                if hit_seq_records:
                    logger.debug("Writting hits to " + out_file)
                    Bio.SeqIO.write(hit_seq_records, out_file, 'fasta')
                else:
                    logger.debug("No hits found, skipping writing output file to " + out_file)
            else:
                logger.debug("No output directory defined for blast hits, skipping writing file")

        return self._create_data_frame(results)

    @abc.abstractmethod
    def _get_out_file_name(self, in_file):
        """
        Gets hits output file name given input file.
        :param in_file: The input file name.
        :return: The output file name.
        """
        pass

    def _handle_blast_hit(self, in_file, database_name, blast_file, results, hit_seq_records):
        blast_handle = open(blast_file)
        blast_records = NCBIXML.parse(blast_handle)
        for blast_record in blast_records:
            partitions = BlastHitPartitions()
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    hit = self._create_hit(in_file, database_name, blast_record, alignment, hsp)
                    if hit.get_pid() >= self._pid_threshold and hit.get_plength() >= self._plength_threshold:
                        partitions.append(hit)
            for hits_non_overlapping in partitions.get_hits_nonoverlapping_regions():
                for hit in self._select_hits_to_include(hits_non_overlapping):
                    self._append_results_to(hit, database_name, results, hit_seq_records)
        blast_handle.close()

    def _select_hits_to_include(self, hits):
        hits_to_include = []

        if len(hits) >= 1:
            sorted_hits_pid_first = sorted(hits, key=lambda x: (
                x.get_pid(), x.get_plength(), x.get_alignment_length(), x.get_hit_id()), reverse=True)
            sorted_hits_length_first = sorted(hits, key=lambda x: (
                x.get_alignment_length(), x.get_pid(), x.get_plength(), x.get_hit_id()), reverse=True)

            if self._report_all:
                hits_to_include = sorted_hits_pid_first
            else:
                first_hit_pid = sorted_hits_pid_first[0]
                first_hit_length = sorted_hits_length_first[0]

                if first_hit_pid == first_hit_length:
                    hits_to_include.append(first_hit_length)
                # if the top length hit is significantly longer, and the pid is not too much below the top pid hit (nor percent overlap too much below top pid hit), use the longer hit
                elif (first_hit_length.get_alignment_length() - first_hit_pid.get_alignment_length()) > 10 and (
                        first_hit_length.get_pid() - first_hit_pid.get_pid()) > -1 and (
                        first_hit_length.get_plength() - first_hit_pid.get_plength()) > -1:
                    hits_to_include.append(first_hit_length)
                # otherwise, prefer the top pid hit, even if it's shorter than the longest hit
                else:
                    hits_to_include.append(first_hit_pid)

        return hits_to_include

    @abc.abstractmethod
    def _create_data_frame(self, results):
        pass

    @abc.abstractmethod
    def _create_hit(self, file, database_name, blast_record, alignment, hsp):
        pass

    @abc.abstractmethod
    def _append_results_to(self, hit, database_name, results, hit_seq_records):
        pass

    def _append_seqrecords_to(self, hit, hit_seq_records):
        seq_record = SeqRecord(Seq(hit.get_hsp_query_proper()), id=hit.get_hit_id(),
                               description='isolate: ' + hit.get_isolate_id() +
                                           ', contig: ' + hit.get_contig() +
                                           ', contig_start: ' + str(hit.get_contig_start()) +
                                           ', contig_end: ' + str(hit.get_contig_end()) +
                                           ', resistance_gene_start: ' + str(hit.get_resistance_gene_start()) +
                                           ', resistance_gene_end: ' + str(hit.get_resistance_gene_end()) +
                                           ', hsp/length: ' + str(hit.get_hsp_alignment_length()) + '/' + str(
                                   hit.get_alignment_length()) +
                                           ', pid: ' + str("%0.2f%%" % hit.get_pid()) +
                                           ', plength: ' + str("%0.2f%%" % hit.get_plength()))
        logger.debug("seq_record=" + repr(seq_record))
        hit_seq_records.append(seq_record)
