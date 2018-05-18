import logging
from typing import Dict
from typing import Union
from typing import List

logger = logging.getLogger('BlastHits')

from staramr.blast.results.AMRHitHSP import AMRHitHSP

"""
Class for partitioning up blast hits into non-overlapping regions.
"""


class BlastHitPartitions:

    def __init__(self):
        """
        Creates a new object to store BLAST hit partitions.
        """
        self._partitions = {}

    def append(self, hit: AMRHitHSP) -> None:
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

    def _add_hit_partition(self, hit: AMRHitHSP, partition: Dict) -> None:
        if hit.get_genome_contig_start() < partition['start']:
            partition['start'] = hit.get_genome_contig_start()

        if hit.get_genome_contig_end() > partition['end']:
            partition['end'] = hit.get_genome_contig_end()

        partition['hits'].append(hit)

    def _get_existing_partition(self, hit: AMRHitHSP) -> Dict:
        partition_name = hit.get_genome_contig_id()

        if partition_name in self._partitions:
            contig_partitions_list = self._partitions[partition_name]
            for partition in contig_partitions_list:
                if self._hit_in_parition(hit, partition):
                    return partition

        return None

    def _hit_in_parition(self, hit: AMRHitHSP, partition: Dict) -> bool:
        pstart, pend = partition['start'], partition['end']
        start, end = hit.get_genome_contig_start(), hit.get_genome_contig_end()

        return (pstart < start < pend) or (pstart < end < pend) or (start <= pstart and end >= pend)

    def _create_new_parition(self, hit: AMRHitHSP) -> None:
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

    def get_hits_nonoverlapping_regions(self) -> List[AMRHitHSP]:
        """
        Gets BLAST hits divided up into separate lists for non-overlapping regions..
        :return: A list of BLAST hits divided up into non-overlapping regions.
        """
        return [p['hits'] for name in self._partitions for p in self._partitions[name]]
