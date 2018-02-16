class AMRHitHSP:

    def __init__(self, file, blast_record, hit, hsp):
        self._file = file
        self._blast_record = blast_record
        self.hit = hit
        self.hsp = hsp

    def get_alignment_length(self):
        return self.hit.length

    def get_hsp_alignment_length(self):
        return self.hsp.align_length

    def get_pid(self):
        return (self.hsp.identities / self.hsp.align_length) * 100

    def get_plength(self):
        return (self.get_hsp_alignment_length() / self.get_alignment_length()) * 100

    def get_hit_id(self):
        return self.hit.hit_id

    def get_file(self):
        return self._file

    def get_contig(self):
        return self._blast_record.query

    def get_contig_start(self):
        return self.hsp.query_start

    def get_contig_end(self):
        return self.hsp.query_end
