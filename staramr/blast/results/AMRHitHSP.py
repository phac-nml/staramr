import os
import re

import Bio.Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
Class used to store/parse AMR BLAST hits/hsps.
"""


class AMRHitHSP:

    def __init__(self, file, blast_record):
        """
        Creates a new AMRHitHSP.
        :param file: The particular file this BLAST hit came from.
        :param blast_record: The Bio.Blast.Record this hit came from.
        """
        self._file = file

        if blast_record is not None:
            self._blast_record = blast_record.to_dict()

    def get_alignment_length(self):
        """
        Gets the BLAST alignment length.
        :return: The BLAST alignment length.
        """
        return self._blast_record['slen']

    def get_hsp_alignment_length(self):
        """
        Gets the BLAST HSP length.
        :return: The BLAST HSP length.
        """
        return self._blast_record['length']

    def get_pid(self):
        """
        Gets the percent identity of the HSP.
        :return: The HSP percent identity.
        """
        return self._blast_record['pident']

    def get_plength(self):
        """
        Gets the percent length of the HSP.
        :return: The percent length of the HSP.
        """
        return (self.get_hsp_alignment_length() / self.get_alignment_length()) * 100

    def get_hit_id(self):
        """
        Gets the hit id.
        :return: The hit id.
        """
        return self._blast_record['sseqid']

    def get_file(self):
        """
        Gets the corresponding input file.
        :return: The corresponding input file.
        """
        return self._file

    def get_isolate_id(self):
        """
        Gets isolate id for the file.
        :return: The isolate id for the file
        """
        return os.path.splitext(self._file)[0]

    def get_contig(self):
        """
        Gets the particular contig id this HSP came from in the input file.
        :return: The contig id.
        """
        re_search = re.search(r'^(\S+)', self._blast_record['qseqid'])
        return re_search.group(1)

    def get_contig_start(self):
        """
        Gets the start of the HSP in the contig in the input file.
        :return: The start of the HSP.
        """
        return self._blast_record['qstart']

    def get_contig_end(self):
        """
        Gets the end of the HSP in the contig from the input file.
        :return: The end of the HSP.
        """
        return self._blast_record['qend']

    def get_resistance_gene_start(self):
        """
        Gets the start of the hsp to the resistance gene.
        :return: The start of the resistance gene hsp.
        """
        return self._blast_record['sstart']

    def get_resistance_gene_end(self):
        """
        Gets the end of the hsp to the resistance gene.
        :return: The end of the resistance gene hsp.
        """
        return self._blast_record['send']

    def get_query_seq(self):
        """
        Gets the query sequence from the HSP.
        :return: The query sequence (as a string) from the HSP.
        """
        return self._blast_record['qseq']

    def get_query_seq_in_database_strand(self):
        """
        Gets the query sequence from the HSP.
        :return: The query sequence (as a string) from the HSP.
        """
        if self.get_database_strand() == 'plus':
            return self.get_query_seq()
        else:
            return Bio.Seq.reverse_complement(self.get_query_seq())

    def get_database_strand(self):
        """
        Gets the database (subject) strand for the BLAST hit.
        :return: The database (subject) strand for the BLAST hit.
        """
        return self._blast_record['sstrand']

    def get_seq_record(self):
        """
        Gets a SeqRecord for this hit.
        :return: A SeqRecord for this hit.
        """
        return SeqRecord(Seq(self.get_query_seq_in_database_strand()), id=self.get_hit_id(),
                         description='isolate: ' + self.get_isolate_id() +
                                     ', contig: ' + self.get_contig() +
                                     ', contig_start: ' + str(self.get_contig_start()) +
                                     ', contig_end: ' + str(self.get_contig_end()) +
                                     ', resistance_gene_start: ' + str(self.get_resistance_gene_start()) +
                                     ', resistance_gene_end: ' + str(self.get_resistance_gene_end()) +
                                     ', hsp/length: ' + str(self.get_hsp_alignment_length()) + '/' + str(
                             self.get_alignment_length()) +
                                     ', pid: ' + str("%0.2f%%" % self.get_pid()) +
                                     ', plength: ' + str("%0.2f%%" % self.get_plength()))
