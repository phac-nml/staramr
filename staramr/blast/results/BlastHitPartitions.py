import logging

logger = logging.getLogger('BlastHits')

"""
Class for partitioning up blast hits into non-overlapping regions.
"""


class BlastHitPartitions:

    def __init__(self):
        """
        Creates a new object to store BLAST hit partitions.
        """
        self._partitions = []

    def append(self, hit):
        """
        Adds a new blast hit to the set of partitions.
        :param hit: The hit to add.
        :return: None
        """
        if hit.get_contig_start() > hit.get_contig_end():
            raise Exception(
                "Unsupported condition: contig start > contig end for hit (contig=" + hit.get_contig() + ", start=" +
                str(hit.get_contig_start()) + ", end=" + str(hit.get_contig_end()) + ")")
        if hit.get_query_frame() != 1:
            raise Exception("Unsupported condition: query frame is not 1 for hit (contig=" + hit.get_contig() +
                            ", start=" + str(hit.get_contig_start()) + ", end=" + str(
                hit.get_contig_end()) + ", query_frame=" + str(hit.get_query_frame()) + ")")

        partition = self._find_parition(hit)
        if (partition is None):
            self._partitions.append(self._create_new_parition(hit))
        else:
            self._add_hit_partition(hit, partition)

    def _add_hit_partition(self, hit, partition):
        if partition['name'] != hit.get_contig():
            raise Exception("Cannot add hit with different contig name to partition")

        if hit.get_contig_start() < partition['start']:
            partition['start'] = hit.get_contig_start()

        if hit.get_contig_end() > partition['end']:
            partition['end'] = hit.get_contig_end()

        partition['hits'].append(hit)

    def _find_parition(self, hit):
        for partition in self._partitions:
            if self._hit_in_parition(hit, partition):
                return partition
        return None

    def _hit_in_parition(self, hit, partition):
        if partition['name'] == hit.get_contig():
            start = hit.get_contig_start()
            end = hit.get_contig_end()
            if end > partition['start'] and end < partition['end']:
                return True
            elif start < partition['end'] and start > partition['start']:
                return True
            elif start <= partition['start'] and end >= partition['end']:
                return True
        return False

    def _create_new_parition(self, hit):
        return {
            'name': hit.get_contig(),
            'start': hit.get_contig_start(),
            'end': hit.get_contig_end(),
            'hits': [hit]
        }

    def get_hits_nonoverlapping_regions(self):
        """
        Gets BLAST hits divided up into separate lists for non-overlapping regions..
        :return: A list of BLAST hits divided up into non-overlapping regions.
        """
        partitions_list = []
        for partition in self._partitions:
            logger.debug("Partition " + repr(partition))
            partitions_list.append(partition['hits'])
        logger.debug("Done printing partitions")

        return partitions_list
