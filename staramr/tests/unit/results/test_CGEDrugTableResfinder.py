import logging
import unittest

from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from staramr.databases.AMRDatabasesManager import AMRDatabasesManager
from staramr.databases.resistance.cge.CGEDrugTableResfinder import CGEDrugTableResfinder

logger = logging.getLogger('CGEDrugTableResfinderTest')


class CGEDrugTableResfinderTest(unittest.TestCase):

    def setUp(self):
        blast_databases_repositories = AMRDatabasesManager.create_default_manager().get_database_repos()
        self.resfinder_dir = blast_databases_repositories.get_repo_dir('resfinder')
        self.resfinder_database = ResfinderBlastDatabase(self.resfinder_dir)
        self.cge_drug_table = CGEDrugTableResfinder(self.resfinder_database.get_phenotypes_file())

    def testGetDrug_aac_1(self):

        # As of 2023-04-21, the CGE Resfinder phenotypes.txt contains duplicate entries for:
        # aac(6')-Ib-cr_1_DQ303918
        # This test ensures that the duplicate entries are handled correctly.

        drug_class = "aminoglycoside, quinolone"
        gene_plus_variant = "aac(6')-Ib-cr_1"
        accession = "DQ303918"

        drug = self.cge_drug_table.get_drug(drug_class, gene_plus_variant, accession)

        self.assertEqual(drug,
                         "Tobramycin, Dibekacin, Amikacin, Sisomicin, Netilmicin, Fluoroquinolone, Ciprofloxacin",
                         "The drug does not match.")

    def testGetDrug_aac_2(self):

        # As of 2023-04-21, the CGE Resfinder phenotypes.txt contains duplicate entries for:
        # aac(6')-Ib-cr_2_DQ303918
        # This test ensures that the duplicate entries are handled correctly.

        drug_class = "aminoglycoside, quinolone"
        gene_plus_variant = "aac(6')-Ib-cr_2"
        accession = "EF636461"

        drug = self.cge_drug_table.get_drug(drug_class, gene_plus_variant, accession)

        self.assertEqual(drug,
                         "Tobramycin, Dibekacin, Amikacin, Sisomicin, Netilmicin, Fluoroquinolone, Ciprofloxacin",
                         "The drug does not match.")

    def test_get_cge_notes(self):

        # ant(2'')-Ia_1_X04555	[...]	Alternative name aadB
        gene_plus_variant = "ant(2'')-Ia_1"
        accession = "X04555"

        notes = self.cge_drug_table.get_cge_notes(gene_plus_variant, accession)
        self.assertEqual(notes, "Alternative name aadB", "The notes do not match.")

        # ant(2'')-Ia_10_HM367617	[...]	[NO NOTES]
        gene_plus_variant = "ant(2'')-Ia_10"
        accession = "HM367617"

        notes = self.cge_drug_table.get_cge_notes(gene_plus_variant, accession)
        self.assertEqual(notes, "", "The notes do not match.")

        # aac(6')-Ib-cr_1_DQ303918	[...]	[MIC of ciprofloxacin does not always ...]
        # (duplicate entry)
        gene_plus_variant = "aac(6')-Ib-cr_1"
        accession = "DQ303918"

        notes = self.cge_drug_table.get_cge_notes(gene_plus_variant, accession)
        print(notes)
        self.assertEqual(notes,
                         "MIC of ciprofloxacin does not always increase above ECOFF PMID 16369542;MIC of ciprofloxacin does not always increase above ECOFF PMID 16369542",
                         "The notes do not match.")
