from staramr.databases.resistance.pointfinder.ARGDrugTablePointfinder import ARGDrugTablePointfinder
from staramr.databases.resistance.resfinder.ARGDrugTableResfinder import ARGDrugTableResfinder
from staramr.detection.AMRDetection import AMRDetection
from staramr.detection.AMRDetectionResistance import AMRDetectionResistance

"""
A Class used to construct a particular staramr.detection.AMRDetection object.
"""


class AMRDetectionFactory:

    def __init__(self):
        pass

    def build(self, plasmidfinder_database, resfinder_database, blast_handler, pointfinder_database, include_negatives,
              include_resistances=False, output_dir=None, genes_to_exclude=[]):
        """
        Builds a new AMRDetection object.
        :param plasmidfinder_database: The staramr.blast.plasmidfinder.PlasmidfinderBlastDatabase to use for the particular PlasmidFinder database.
        :param resfinder_database: The staramr.blast.resfinder.ResfinderBlastDatabase for the particular ResFinder database.
        :param blast_handler: The staramr.blast.BlastHandler to use for scheduling BLAST jobs.
        :param pointfinder_database: The staramr.blast.pointfinder.PointfinderBlastDatabase to use for the particular PointFinder database.
        :param include_negatives:  If True, include files lacking AMR genes in the resulting summary table.
        :param include_resistances: If True, include predicted drug resistances in output.
        :param output_dir: The directory where output files are being written.
        :param genes_to_exclude: A list of gene IDs to exclude from the results.
        :return: A new AMRDetection object.
        """

        if include_resistances:
            return AMRDetectionResistance(resfinder_database, ARGDrugTableResfinder(), blast_handler,
                                          ARGDrugTablePointfinder(), pointfinder_database, include_negatives,
                                          output_dir=output_dir, genes_to_exclude=genes_to_exclude,
                                          plasmidfinder_database=plasmidfinder_database)
        else:
            return AMRDetection(resfinder_database, blast_handler, pointfinder_database, include_negatives,
                                output_dir=output_dir, genes_to_exclude=genes_to_exclude,
                                plasmidfinder_database=plasmidfinder_database)
