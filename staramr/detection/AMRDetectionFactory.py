from staramr.detection.AMRDetection import AMRDetection

"""
A Class used to construct a particular staramr.detection.AMRDetection object.
"""


class AMRDetectionFactory:

    def __init__(self):
        pass

    def build(self, resfinder_database, blast_handler, pointfinder_database, include_negatives, output_dir=None):
        """
        Builds a new AMRDetection object.
        :param resfinder_database: The staramr.blast.resfinder.ResfinderBlastDatabase for the particular ResFinder database.
        :param blast_handler: The staramr.blast.BlastHandler to use for scheduling BLAST jobs.
        :param pointfinder_database: The staramr.blast.pointfinder.PointfinderBlastDatabase to use for the particular PointFinder database.
        :param include_negatives:  If True, include files lacking AMR genes in the resulting summary table.
        :param output_dir: The directory where output files are being written.
        :return: A new AMRDetection object.
        """
        return AMRDetection(resfinder_database, blast_handler, pointfinder_database, include_negatives,
                            output_dir=output_dir)
