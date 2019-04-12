import logging
from collections import OrderedDict
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple

logger = logging.getLogger('BlastHits')

from staramr.blast.results.AMRHitHSP import AMRHitHSP
from staramr.exceptions.InvalidPositionException import InvalidPositionException

"""
Class for partitioning up blast hits into non-overlapping regions.
"""


class BlastHitPartitions:

    def __init__(self):
        """
        Creates a new object to store BLAST hit partitions.
        """
        self._partitions = OrderedDict()

    def append(self, hit: AMRHitHSP) -> None:
        """
        Adds a new blast hit to the set of partitions.
        :param hit: The hit to add.
        :return: None
        """
        if hit.get_genome_contig_start() > hit.get_genome_contig_end() and hit.get_genome_contig_strand() == 'plus':
            raise InvalidPositionException(
                "Unsupported condition: strand=plus and contig start > contig end for hit (contig=" + hit.get_genome_contig_id() + ", start=" +
                str(hit.get_genome_contig_start()) + ", end=" + str(hit.get_genome_contig_end()) + ")")

        partition = self._get_existing_partition(hit)
        if (partition is None):
            self._create_new_parition(hit)
        else:
            self._add_hit_partition(hit, partition)

    def _add_hit_partition(self, hit: AMRHitHSP, partition: Dict[str, Any]) -> None:
        start, end = self._stranded_ends(hit)

        if start < partition['start']:
            partition['start'] = start

        if end > partition['end']:
            partition['end'] = end

        partition['hits'].append(hit)

    def _get_existing_partition(self, hit: AMRHitHSP) -> Optional[Dict[str, Any]]:
        partition_name = hit.get_genome_contig_id()

        if partition_name in self._partitions:
            contig_partitions_list = self._partitions[partition_name]
            for partition in contig_partitions_list:
                if self._hit_in_parition(hit, partition):
                    return partition

        return None

    def _hit_in_parition(self, hit: AMRHitHSP, partition: Dict[str, Any]) -> bool:
        pstart, pend = partition['start'], partition['end']
        start, end = self._stranded_ends(hit)

        return (pstart < start < pend) or (pstart < end < pend) or (start <= pstart and end >= pend)

    def _create_new_parition(self, hit: AMRHitHSP) -> None:
        start, end = self._stranded_ends(hit)

        contig_name = hit.get_genome_contig_id()
        partition = {
            'start': start,
            'end': end,
            'hits': [hit]
        }

        if contig_name in self._partitions:
            self._partitions[contig_name].append(partition)
        else:
            self._partitions[contig_name] = [partition]

    def get_hits_nonoverlapping_regions(self) -> List[List[AMRHitHSP]]:
        """
        Gets BLAST hits divided up into separate lists for non-overlapping regions..
        :return: A list of BLAST hits divided up into non-overlapping regions.
        """
        return [p['hits'] for name in self._partitions for p in self._partitions[name]]

    def _stranded_ends(self, hit: AMRHitHSP) -> Tuple[int, int]:
        """
        Gets the start/end coordinates, taking into account the strand.
        :param hit: The hit.
        :return: The (start,end) as a tuple.
        """
        start = hit.get_genome_contig_start() if hit.get_genome_contig_strand() == 'plus' else hit.get_genome_contig_end()
        end = hit.get_genome_contig_end() if hit.get_genome_contig_strand() == 'plus' else hit.get_genome_contig_start()
        return start, end
