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
        self._partitions = {}

    def append(self, hit):
        """
        Adds a new blast hit to the set of partitions.
        :param hit: The hit to add.
        :return: None
        """
        if hit.get_genome_contig_start() > hit.get_genome_contig_end():
            raise Exception(
                "Unsupported condition: contig start > contig end for hit (contig=" + hit.get_genome_contig_id() + ", start=" +
                str(hit.get_genome_contig_start()) + ", end=" + str(hit.get_genome_contig_end()) + ")")

        partition = self._get_existing_partition(hit)
        if (partition is None):
            self._create_new_parition(hit)
        else:
            self._add_hit_partition(hit, partition)

    def _add_hit_partition(self, hit, partition):
        if hit.get_genome_contig_start() < partition['start']:
            partition['start'] = hit.get_genome_contig_start()

        if hit.get_genome_contig_end() > partition['end']:
            partition['end'] = hit.get_genome_contig_end()

        partition['hits'].append(hit)

    def _get_existing_partition(self, hit):
        partition_name = hit.get_genome_contig_id()

        if partition_name in self._partitions:
            contig_partitions_list = self._partitions[partition_name]
            for partition in contig_partitions_list:
                if self._hit_in_parition(hit, partition):
                    return partition

        return None

    def _hit_in_parition(self, hit, partition):
        pstart, pend = partition['start'], partition['end']
        start, end = hit.get_genome_contig_start(), hit.get_genome_contig_end()

        return (pstart < start < pend) or (pstart < end < pend) or (start <= pstart and end >= pend)

    def _create_new_parition(self, hit):
        contig_name = hit.get_genome_contig_id()
        partition =  {
            'start': hit.get_genome_contig_start(),
            'end': hit.get_genome_contig_end(),
            'hits': [hit]
        }

        if contig_name in self._partitions:
            self._partitions[contig_name].append(partition)
        else:
            self._partitions[contig_name] = [partition]

    def get_hits_nonoverlapping_regions(self):
        """
        Gets BLAST hits divided up into separate lists for non-overlapping regions..
        :return: A list of BLAST hits divided up into non-overlapping regions.
        """
        partitions_list = []
        for name in self._partitions:
            for partition in self._partitions[name]:
                partitions_list.append(partition['hits'])

        return partitions_list
