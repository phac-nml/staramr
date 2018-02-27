from staramr.detection.AMRDetection import AMRDetection

class AMRDetectionFactory:

    def __init__(self):
        pass

    def build(self, resfinder_database, blast_handler, pointfinder_database, include_negatives):
        return AMRDetection(resfinder_database, blast_handler, pointfinder_database, include_negatives)
