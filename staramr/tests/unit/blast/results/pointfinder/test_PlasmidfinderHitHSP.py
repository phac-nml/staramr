import unittest

from staramr.blast.results.plasmidfinder.PlasmidfinderHitHSP import PlasmidfinderHitHSP


class PlasmidfinderHitHSPTest(unittest.TestCase):

    def testParseSequenceId1(self):
        test_blast_record = {"sstart": 20, "send": 30, "sstrand": "ABC", "qstart": 1, "qend": 10,
                             'qseqid': 'RepA_1_pKPC-CAV1321_CP011611'}

        plasmid_hit_hsp = PlasmidfinderHitHSP('test_file', test_blast_record)

        self.assertEqual('RepA', plasmid_hit_hsp.get_amr_gene_name(), 'Did not parse correct gene name')
        self.assertEqual('1', plasmid_hit_hsp.get_amr_gene_variant(), 'Did not parse correct gene variant')
        self.assertEqual('RepA_1', plasmid_hit_hsp.get_amr_gene_name_with_variant(),
                         'Did not parse correct gene name variant')
        self.assertEqual('CP011611', plasmid_hit_hsp.get_amr_gene_accession(),
                         'Did not parse correct gene name variant')
        self.assertEqual('RepA_1_CP011611', plasmid_hit_hsp.get_amr_gene_variant_accession(),
                         'Did not parse correct gene name variant accession')

    def testParseSequenceId2(self):
        test_blast_record = {"sstart": 20, "send": 30, "sstrand": "ABC", "qstart": 1, "qend": 10,
                             'qseqid': 'IncHI2_1__BX664015'}

        plasmid_hit_hsp = PlasmidfinderHitHSP('test_file', test_blast_record)

        self.assertEqual('IncHI2', plasmid_hit_hsp.get_amr_gene_name(), 'Did not parse correct gene name')
        self.assertEqual('1', plasmid_hit_hsp.get_amr_gene_variant(), 'Did not parse correct gene variant')
        self.assertEqual('IncHI2_1', plasmid_hit_hsp.get_amr_gene_name_with_variant(),
                         'Did not parse correct gene name variant')
        self.assertEqual('BX664015', plasmid_hit_hsp.get_amr_gene_accession(),
                         'Did not parse correct gene name variant')
        self.assertEqual('IncHI2_1_BX664015', plasmid_hit_hsp.get_amr_gene_variant_accession(),
                         'Did not parse correct gene name variant accession')

    def testParseSequenceId3(self):
        test_blast_record = {"sstart": 20, "send": 30, "sstrand": "ABC", "qstart": 1, "qend": 10,
                             'qseqid': 'IncB/O/K/Z_1__CU928147'}

        plasmid_hit_hsp = PlasmidfinderHitHSP('test_file', test_blast_record)

        self.assertEqual('IncB/O/K/Z', plasmid_hit_hsp.get_amr_gene_name(), 'Did not parse correct gene name')
        self.assertEqual('1', plasmid_hit_hsp.get_amr_gene_variant(), 'Did not parse correct gene variant')
        self.assertEqual('IncB/O/K/Z_1', plasmid_hit_hsp.get_amr_gene_name_with_variant(),
                         'Did not parse correct gene name variant')
        self.assertEqual('CU928147', plasmid_hit_hsp.get_amr_gene_accession(),
                         'Did not parse correct gene name variant')
        self.assertEqual('IncB/O/K/Z_1_CU928147', plasmid_hit_hsp.get_amr_gene_variant_accession(),
                         'Did not parse correct gene name variant accession')

    def testParseSequenceId4(self):
        test_blast_record = {"sstart": 20, "send": 30, "sstrand": "ABC", "qstart": 1, "qend": 10,
                             'qseqid': 'IncFII(Serratia)_1_Serratia_NC_009829'}

        plasmid_hit_hsp = PlasmidfinderHitHSP('test_file', test_blast_record)

        self.assertEqual('IncFII(Serratia)', plasmid_hit_hsp.get_amr_gene_name(), 'Did not parse correct gene name')
        self.assertEqual('1', plasmid_hit_hsp.get_amr_gene_variant(), 'Did not parse correct gene variant')
        self.assertEqual('IncFII(Serratia)_1', plasmid_hit_hsp.get_amr_gene_name_with_variant(),
                         'Did not parse correct gene name variant')
        self.assertEqual('NC_009829', plasmid_hit_hsp.get_amr_gene_accession(),
                         'Did not parse correct gene name variant')
        self.assertEqual('IncFII(Serratia)_1_NC_009829', plasmid_hit_hsp.get_amr_gene_variant_accession(),
                         'Did not parse correct gene name variant accession')

    def testParseUnderscoresBrackets(self):
        # Tests to ensure that PlasmidfinderHitHSP's init function can properly parse FASTA record IDs that have
        # underscores within brackets. For example: "rep21_24_rep(CN1_plasmid2)_NC_022227". The "rep(CN1_plasmid2)"
        # needs to be parsed as a single element.

        test_blast_record = {"sstart": 20, "send": 30, "sstrand": "ABC", "qstart": 1, "qend": 10,
                             'qseqid': 'rep21_24_rep(CN1_plasmid2)_NC_022227'}

        plasmid_hit_hsp = PlasmidfinderHitHSP('test_file', test_blast_record)

        self.assertEqual('rep21', plasmid_hit_hsp.get_amr_gene_name(), 'Did not parse correct gene name')
        self.assertEqual('24', plasmid_hit_hsp.get_amr_gene_variant(), 'Did not parse correct gene variant')
        self.assertEqual('rep21_24', plasmid_hit_hsp.get_amr_gene_name_with_variant(),
                         'Did not parse correct gene name variant')
        self.assertEqual('NC_022227', plasmid_hit_hsp.get_amr_gene_accession(),
                         'Did not parse correct gene name variant')
        self.assertEqual('rep21_24_NC_022227', plasmid_hit_hsp.get_amr_gene_variant_accession(),
                         'Did not parse correct gene name variant accession')
