from staramr.detection.AMRDetectionFactory import AMRDetectionFactory
from staramr.detection.AMRDetectionPlus import AMRDetectionPlus
from staramr.databases.dr.phenotype.pointfinder.ARGDrugTablePointfinder import ARGDrugTablePointfinder
from staramr.databases.dr.phenotype.resfinder.ARGDrugTableResfinder import ARGDrugTableResfinder


class AMRDetectionFactoryPlus(AMRDetectionFactory):

    def __init__(self):
        super().__init__()

    def build(self, resfinder_database, blast_handler, pointfinder_database, include_negatives, output_dir=None):

        arg_drug_table_resfinder = ARGDrugTableResfinder()
        arg_drug_table_pointfinder = ARGDrugTablePointfinder()

        return AMRDetectionPlus(resfinder_database, arg_drug_table_resfinder, blast_handler, arg_drug_table_pointfinder,
                                pointfinder_database, include_negatives, output_dir=output_dir)
