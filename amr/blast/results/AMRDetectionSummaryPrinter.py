import sys

class AMRDetectionSummaryPrinter:

    def __init__(self, resfinder_hits, pointfinder_hits = None):
        self._resfinder_hits = resfinder_hits
        self._pointfinder_hits = pointfinder_hits

        if pointfinder_hits:
            self._has_pointfinder_hits = True
        else:
            self._has_pointfinder_hits = False

    def print_summary(self, file = None):
        file_handle = sys.stdout

        if not file:
            file_handle = open(file, 'w')

        for file, hits in self._resfinder_hits.items():
            for hit in hits:
                if (hit.gene):
                    return