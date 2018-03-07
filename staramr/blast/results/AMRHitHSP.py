import re

"""
Class used to store/parse AMR BLAST hits/hsps.
"""


class AMRHitHSP:

    def __init__(self, file, blast_record, hit, hsp):
        """
        Creates a new AMRHitHSP.
        :param file: The particular file this BLAST hit came from.
        :param blast_record: The Bio.Blast.Record this hit came from.
        :param hit: The particular Bio.Blast.Record.Alignment.
        :param hsp: The particular Bio.Blast.Record.HSP.
        """
        self._file = file
        self._blast_record = blast_record
        self.hit = hit
        self.hsp = hsp

    def get_alignment_length(self):
        """
        Gets the BLAST alignment length.
        :return: The BLAST alignment length.
        """
        return self.hit.length

    def get_hsp_alignment_length(self):
        """
        Gets the BLAST HSP length.
        :return: The BLAST HSP length.
        """
        return self.hsp.align_length

    def get_pid(self):
        """
        Gets the percent identity of the HSP.
        :return: The HSP percent identity.
        """
        return (self.hsp.identities / self.hsp.align_length) * 100

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
        return self.hit.hit_id

    def get_file(self):
        """
        Gets the corresponding input file.
        :return: The corresponding input file.
        """
        return self._file

    def get_contig(self):
        """
        Gets the particular contig id this HSP came from in the input file.
        :return: The contig id.
        """
        re_search = re.search("^(\S+)", self._blast_record.query)
        return re_search.group(1)

    def get_contig_start(self):
        """
        Gets the start of the HSP in the contig in the input file.
        :return: The start of the HSP.
        """
        return self.hsp.query_start

    def get_contig_end(self):
        """
        Gets the end of the HSP in the contig from the input file.
        :return: The end of the HSP.
        """
        return self.hsp.query_end
