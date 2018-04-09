import unittest
from unittest.mock import MagicMock

from staramr.blast.results.AMRHitHSP import AMRHitHSP
from staramr.blast.results.BlastHitPartitions import BlastHitPartitions


class BlastHitPartitionsTest(unittest.TestCase):

    def testSinglePartition(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=1)
        hit1.get_contig_end = MagicMock(return_value=10)

        parts = BlastHitPartitions()

        parts.append(hit1)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(1, len(return_list), "Should only be one partition")
        self.assertEqual(1, len(return_list[0]), "Should only be one hit")
        self.assertEqual('contig1', return_list[0][0].get_contig(), "Should have correct contig name")
        self.assertEqual(1, return_list[0][0].get_contig_start(), "Should have correct contig start")
        self.assertEqual(10, return_list[0][0].get_contig_end(), "Should have correct contig start")

    def testSinglePartitionIdenticalHits(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=1)
        hit1.get_contig_end = MagicMock(return_value=10)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit1)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(1, len(return_list), "Should only be one partition")
        self.assertEqual(2, len(return_list[0]), "Should be two hits")
        self.assertEqual(['contig1', 'contig1'], [x.get_contig() for x in return_list[0]],
                         "Should have correct contig names")
        self.assertEqual([1, 1], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10, 10], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

    def testSinglePartitionHitEndGreater(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=1)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=1)
        hit2.get_contig_end = MagicMock(return_value=11)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(1, len(return_list), "Should only be one partition")
        self.assertEqual(2, len(return_list[0]), "Should be two hits")
        self.assertEqual(['contig1', 'contig1'], [x.get_contig() for x in return_list[0]],
                         "Should have correct contig names")
        self.assertEqual([1, 1], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10, 11], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

    def testSinglePartitionHitEndGreaterStartGreater(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=1)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=2)
        hit2.get_contig_end = MagicMock(return_value=11)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(1, len(return_list), "Should only be one partition")
        self.assertEqual(2, len(return_list[0]), "Should be two hits")
        self.assertEqual(['contig1', 'contig1'], [x.get_contig() for x in return_list[0]],
                         "Should have correct contig names")
        self.assertEqual([1, 2], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10, 11], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

    def testSinglePartitionHitEndLesser(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=1)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=1)
        hit2.get_contig_end = MagicMock(return_value=9)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(1, len(return_list), "Should only be one partition")
        self.assertEqual(2, len(return_list[0]), "Should be two hits")
        self.assertEqual(['contig1', 'contig1'], [x.get_contig() for x in return_list[0]],
                         "Should have correct contig names")
        self.assertEqual([1, 1], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10, 9], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

    def testSinglePartitionHitEndLesserStartLesser(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=2)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=1)
        hit2.get_contig_end = MagicMock(return_value=9)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(1, len(return_list), "Should only be one partition")
        self.assertEqual(2, len(return_list[0]), "Should be two hits")
        self.assertEqual(['contig1', 'contig1'], [x.get_contig() for x in return_list[0]],
                         "Should have correct contig names")
        self.assertEqual([2, 1], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10, 9], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

    def testSinglePartitionHitEndGreaterStartLesser(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=2)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=1)
        hit2.get_contig_end = MagicMock(return_value=11)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(1, len(return_list), "Should only be one partition")
        self.assertEqual(2, len(return_list[0]), "Should be two hits")
        self.assertEqual(['contig1', 'contig1'], [x.get_contig() for x in return_list[0]],
                         "Should have correct contig names")
        self.assertEqual([2, 1], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10, 11], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

    def testSinglePartitionHit2InHit1(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=1)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=2)
        hit2.get_contig_end = MagicMock(return_value=9)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(1, len(return_list), "Should only be one partition")
        self.assertEqual(2, len(return_list[0]), "Should be two hits")
        self.assertEqual(['contig1', 'contig1'], [x.get_contig() for x in return_list[0]],
                         "Should have correct contig names")
        self.assertEqual([1, 2], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10, 9], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

    def testSinglePartitionHit2EdgeWithinHit1Lesser(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=5)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=1)
        hit2.get_contig_end = MagicMock(return_value=6)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(1, len(return_list), "Should only be one partition")
        self.assertEqual(2, len(return_list[0]), "Should be two hits")
        self.assertEqual(['contig1', 'contig1'], [x.get_contig() for x in return_list[0]],
                         "Should have correct contig names")
        self.assertEqual([5, 1], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10, 6], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

    def testSinglePartitionHit2EdgeSameHit1Lesser(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=5)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=1)
        hit2.get_contig_end = MagicMock(return_value=5)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(2, len(return_list), "Should be two partitions")
        self.assertEqual(1, len(return_list[0]), "Partition 1 should have 1 hit")
        self.assertEqual(['contig1'], [x.get_contig() for x in return_list[0]], "Should have correct contig names")
        self.assertEqual([5], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

        self.assertEqual(1, len(return_list[1]), "Partition 2 should have 1 hit")
        self.assertEqual(['contig1'], [x.get_contig() for x in return_list[1]], "Should have correct contig names")
        self.assertEqual([1], [x.get_contig_start() for x in return_list[1]], "Should have correct contig starts")
        self.assertEqual([5], [x.get_contig_end() for x in return_list[1]], "Should have correct contig ends")

    def testSinglePartitionHit2EdgeWithinHit1Greater(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=5)
        hit1.get_contig_end = MagicMock(return_value=11)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=10)
        hit2.get_contig_end = MagicMock(return_value=15)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(1, len(return_list), "Should only be one partition")
        self.assertEqual(2, len(return_list[0]), "Should be two hits")
        self.assertEqual(['contig1', 'contig1'], [x.get_contig() for x in return_list[0]],
                         "Should have correct contig names")
        self.assertEqual([5, 10], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([11, 15], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

    def testSinglePartitionHit2EdgeSameHit1Greater(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=5)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=10)
        hit2.get_contig_end = MagicMock(return_value=15)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(2, len(return_list), "Should be two partitions")
        self.assertEqual(1, len(return_list[0]), "Partition 1 should have 1 hit")
        self.assertEqual(['contig1'], [x.get_contig() for x in return_list[0]], "Should have correct contig names")
        self.assertEqual([5], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

        self.assertEqual(1, len(return_list[1]), "Partition 2 should have 1 hit")
        self.assertEqual(['contig1'], [x.get_contig() for x in return_list[1]], "Should have correct contig names")
        self.assertEqual([10], [x.get_contig_start() for x in return_list[1]], "Should have correct contig starts")
        self.assertEqual([15], [x.get_contig_end() for x in return_list[1]], "Should have correct contig ends")

    def testTwoPartitionsHit2EdgeHit1Lesser(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=5)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=1)
        hit2.get_contig_end = MagicMock(return_value=4)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(2, len(return_list), "Should be two partitions")
        self.assertEqual(1, len(return_list[0]), "Partition 1 should have 1 hit")
        self.assertEqual(['contig1'], [x.get_contig() for x in return_list[0]], "Should have correct contig names")
        self.assertEqual([5], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

        self.assertEqual(1, len(return_list[1]), "Partition 2 should have 1 hit")
        self.assertEqual(['contig1'], [x.get_contig() for x in return_list[1]], "Should have correct contig names")
        self.assertEqual([1], [x.get_contig_start() for x in return_list[1]], "Should have correct contig starts")
        self.assertEqual([4], [x.get_contig_end() for x in return_list[1]], "Should have correct contig ends")

    def testTwoPartitionsHit2EdgeHit1Greater(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=5)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig1")
        hit2.get_contig_start = MagicMock(return_value=11)
        hit2.get_contig_end = MagicMock(return_value=15)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(2, len(return_list), "Should be two partitions")
        self.assertEqual(1, len(return_list[0]), "Partition 1 should have 1 hit")
        self.assertEqual(['contig1'], [x.get_contig() for x in return_list[0]], "Should have correct contig names")
        self.assertEqual([5], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

        self.assertEqual(1, len(return_list[1]), "Partition 2 should have 1 hit")
        self.assertEqual(['contig1'], [x.get_contig() for x in return_list[1]], "Should have correct contig names")
        self.assertEqual([11], [x.get_contig_start() for x in return_list[1]], "Should have correct contig starts")
        self.assertEqual([15], [x.get_contig_end() for x in return_list[1]], "Should have correct contig ends")

    def testTwoPartitionsDifferentContigNames(self):
        hit1 = AMRHitHSP(None, None, None, None)
        hit1.get_query_frame = MagicMock(return_value=1)
        hit1.get_contig = MagicMock(return_value="contig1")
        hit1.get_contig_start = MagicMock(return_value=1)
        hit1.get_contig_end = MagicMock(return_value=10)

        hit2 = AMRHitHSP(None, None, None, None)
        hit2.get_query_frame = MagicMock(return_value=1)
        hit2.get_contig = MagicMock(return_value="contig2")
        hit2.get_contig_start = MagicMock(return_value=1)
        hit2.get_contig_end = MagicMock(return_value=10)

        parts = BlastHitPartitions()

        parts.append(hit1)
        parts.append(hit2)

        return_list = parts.get_hits_nonoverlapping_regions()
        self.assertEqual(2, len(return_list), "Should be two partitions")
        self.assertEqual(1, len(return_list[0]), "Partition 1 should have 1 hit")
        self.assertEqual(['contig1'], [x.get_contig() for x in return_list[0]], "Should have correct contig names")
        self.assertEqual([1], [x.get_contig_start() for x in return_list[0]], "Should have correct contig starts")
        self.assertEqual([10], [x.get_contig_end() for x in return_list[0]], "Should have correct contig ends")

        self.assertEqual(1, len(return_list[1]), "Partition 2 should have 1 hit")
        self.assertEqual(['contig2'], [x.get_contig() for x in return_list[1]], "Should have correct contig names")
        self.assertEqual([1], [x.get_contig_start() for x in return_list[1]], "Should have correct contig starts")
        self.assertEqual([10], [x.get_contig_end() for x in return_list[1]], "Should have correct contig ends")
