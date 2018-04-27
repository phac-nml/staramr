from os import path

from staramr.detection.AMRDetectionFactory import AMRDetectionFactory
from staramr.detection.AMRDetectionPlus import AMRDetectionPlus
from staramr.dr.phenotype.pointfinder.ARGDrugTablePointfinder import ARGDrugTablePointfinder
from staramr.dr.phenotype.resfinder.ARGDrugTableResfinder import ARGDrugTableResfinder


class AMRDetectionFactoryPlus(AMRDetectionFactory):

    def __init__(self):
        super().__init__()

    def build(self, resfinder_database, blast_handler, pointfinder_database, include_negatives, output_dir=None):
        data_dir = path.join(path.dirname(__file__), 'data')

        arg_drug_table_resfinder_file = path.join(data_dir, 'ARG_drug_key_resfinder.tsv')
        arg_drug_table_resfinder = ARGDrugTableResfinder(arg_drug_table_resfinder_file)

        arg_drug_table_pointfinder_file = path.join(data_dir, 'ARG_drug_key_pointfinder.tsv')
        arg_drug_table_pointfinder = ARGDrugTablePointfinder(arg_drug_table_pointfinder_file)

        return AMRDetectionPlus(resfinder_database, arg_drug_table_resfinder, blast_handler, arg_drug_table_pointfinder,
                                pointfinder_database, include_negatives, output_dir=output_dir)
