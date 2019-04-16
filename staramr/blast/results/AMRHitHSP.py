import abc
import os
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from staramr.exceptions.InvalidPositionException import InvalidPositionException

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
        __metaclass__ = abc.ABCMeta
        self._file = file

        if blast_record is not None:
            self._blast_record = blast_record

            if self.get_genome_contig_start() > self.get_genome_contig_end() and self.get_genome_contig_strand() != 'minus':
                raise InvalidPositionException(
                    "contig start = {} > contig end = {} and strand is {}".format(self.get_genome_contig_start(),
                                                                                  self.get_genome_contig_end(),
                                                                                  self.get_genome_contig_strand()))
            elif self.get_amr_gene_start() > self.get_amr_gene_end():
                raise InvalidPositionException(
                    "amr gene start = {} > amr gene end = {}".format(self.get_amr_gene_start(),
                                                                     self.get_amr_gene_end()))

    def get_amr_gene_length(self):
        """
        Gets the amr gene length.
        :return: The amr gene length.
        """
        return self._blast_record['qlen']

    def get_hsp_length(self):
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
        Gets the percent length of the HSP to the AMR gene.
        :return: The percent length of the HSP to the AMR gene.
        """
        return self._blast_record['plength']

    def get_amr_gene_id(self):
        """
        Gets the hit id.
        :return: The hit id.
        """
        return self._blast_record['qseqid']

    @abc.abstractmethod
    def get_amr_gene_name(self):
        """
        Gets the gene name for the amr gene.
        :return: The gene name.
        """
        pass

    def get_file(self):
        """
        Gets the corresponding input file.
        :return: The corresponding input file.
        """
        return self._file

    def get_genome_id(self):
        """
        Gets genome id for the file.
        :return: The genome id for the file
        """
        return os.path.splitext(self._file)[0]

    def get_genome_contig_id(self):
        """
        Gets the particular id from the genome input file.
        :return: The contig id.
        """
        re_search = re.search(r'^(\S+)', self._blast_record['sseqid'])
        return re_search.group(1)

    def get_genome_contig_start(self) -> int:
        """
        Gets the start of the HSP in the genome input file.
        :return: The start of the HSP.
        """
        return self._blast_record['sstart']

    def get_genome_contig_end(self) -> int:
        """
        Gets the end of the HSP in the genome input file.
        :return: The end of the HSP.
        """
        return self._blast_record['send']

    def get_amr_gene_start(self):
        """
        Gets the start of the hsp to the resistance gene.
        :return: The start of the resistance gene hsp.
        """
        return self._blast_record['qstart']

    def get_amr_gene_end(self):
        """
        Gets the end of the hsp to the resistance gene.
        :return: The end of the resistance gene hsp.
        """
        return self._blast_record['qend']

    def get_amr_gene_seq(self):
        """
        Gets the amr gene from the HSP.
        :return: The amr gene (as a string) from the HSP.
        """
        return self._blast_record['qseq']

    def get_genome_contig_hsp_seq(self):
        """
        Gets the genome sequence from the HSP.
        :return: The genome sequence (as a string) from the HSP.
        """
        return self._blast_record['sseq']

    def get_genome_seq_in_amr_gene_strand(self):
        """
        Gets the query sequence from the HSP.
        :return: The query sequence (as a string) from the HSP.
        """
        return self.get_genome_contig_hsp_seq()

    def get_genome_contig_strand(self):
        """
        Gets the genome contig strand for the BLAST hit.
        :return: The genome contig strand for the BLAST hit.
        """
        return self._blast_record['sstrand']

    def get_seq_record(self):
        """
        Gets a SeqRecord for this hit.
        :return: A SeqRecord for this hit.
        """
        return SeqRecord(Seq(self.get_genome_contig_hsp_seq()), id=self.get_amr_gene_id(),
                         description=(
                             'isolate: {}, contig: {}, contig_start: {}, contig_end: {}, database_gene_start: {},'
                             ' database_gene_end: {}, hsp/length: {}/{}, pid: {:0.2f}%, plength: {:0.2f}%').format(
                             self.get_genome_id(),
                             self.get_genome_contig_id(),
                             self.get_genome_contig_start(),
                             self.get_genome_contig_end(),
                             self.get_amr_gene_start(),
                             self.get_amr_gene_end(),
                             self.get_hsp_length(),
                             self.get_amr_gene_length(),
                             self.get_pid(),
                             self.get_plength()))
