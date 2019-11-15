import logging
import unittest

import pandas as pd
from os import path

from staramr.results.QualityModule import QualityModule

logger = logging.getLogger('QualityModuleTest')

class QualityModuleTest(unittest.TestCase):

    def setUp(self):
        self.genome_size_lower_bound = 400
        self.genome_size_upper_bound = 600
        self.minimum_N50_value = 100
        self.minimum_contig_length = 300
        self.unacceptable_num_contigs = 10
        self.test_data_dir = path.join(path.dirname(__file__), '..', 'data')

    def testN50ExactlyMinimumValue(self):
        file = path.join(self.test_data_dir, "test-N50-Exactly-Minimum-Value.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-N50-Exactly-Minimum-Value', quality_module.index[0], 'File name not equal')
        self.assertEqual(100, quality_module['N50 value'].iloc[0], 'N50 vlaue not equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}] ; N50 value is not greater than the specified minimum value [{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testN50OneBPLargerThanMinimumValue(self):
        file = path.join(self.test_data_dir, "test-N50-One-BP-Larger-Than-Minimum-Value.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-N50-One-BP-Larger-Than-Minimum-Value', quality_module.index[0], 'File name not equal')
        self.assertEqual(101, quality_module['N50 value'].iloc[0], 'N50 vlaue not equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testN50SmallerThanMinimumValue(self):
        file = path.join(self.test_data_dir, "test-N50-Smaller-Than-Minimum-Value.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-N50-Smaller-Than-Minimum-Value', quality_module.index[0], 'File name not equal')
        self.assertEqual(20, quality_module['N50 value'].iloc[0], 'N50 vlaue not equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}] ; N50 value is not greater than the specified minimum value [{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testN50MuchLargerThanMinimumValue(self):
        file = path.join(self.test_data_dir, "test-N50-Much-Larger-Than-Minimum-Value.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-N50-Much-Larger-Than-Minimum-Value', quality_module.index[0], 'File name not equal')
        self.assertEqual(1000, quality_module['N50 value'].iloc[0], 'N50 vlaue not equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testN50Calculation(self):
        #tests to make sure N50 contig length +all contig lengths greater than it >= half of genome length, here we are testing the = part
        file = path.join(self.test_data_dir, "test-N50-Calculation.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-N50-Calculation', quality_module.index[0], 'File name not equal')
        self.assertEqual(102, quality_module['N50 value'].iloc[0], 'N50 vlaue not equal')
        self.assertEqual('Passed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('', quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testN50UnaffectedByEmptyContigs(self):
        file = path.join(self.test_data_dir, "test-N50-Unaffected-By-Empty-Contigs.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-N50-Unaffected-By-Empty-Contigs', quality_module.index[0], 'File name not equal')
        self.assertEqual(101, quality_module['N50 value'].iloc[0], 'N50 vlaue not equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testGenomeSizeExactlyMinimum(self):
        file = path.join(self.test_data_dir, "test-Genome-Size-Exactly-Minimum.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-Genome-Size-Exactly-Minimum', quality_module.index[0], 'File name not equal')
        self.assertEqual(400, quality_module['Genome Length'].iloc[0], 'Genome length not equal')
        self.assertEqual('Passed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('', quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testGenomeSizeExactlyMaximum(self):
        file = path.join(self.test_data_dir, "test-Genome-Size-Exactly-Maximum.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-Genome-Size-Exactly-Maximum', quality_module.index[0], 'File name not equal')
        self.assertEqual(600, quality_module['Genome Length'].iloc[0], 'Genome length not equal')
        self.assertEqual('Passed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('', quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testGenomeSizeWithinAcceptedRange(self):
        file = path.join(self.test_data_dir, "test-Genome-Size-Within-Accepted-Range.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-Genome-Size-Within-Accepted-Range', quality_module.index[0], 'File name not equal')
        self.assertEqual(500, quality_module['Genome Length'].iloc[0], 'Genome length not equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('N50 value is not greater than the specified minimum value [{}]'.format(self.minimum_N50_value), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testGenomeSizeSmallerThanMinimum(self):
        file = path.join(self.test_data_dir, "test-Genome-Size-Smaller-Than-Minimum.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-Genome-Size-Smaller-Than-Minimum', quality_module.index[0], 'File name not equal')
        self.assertEqual(200, quality_module['Genome Length'].iloc[0], 'Genome length not equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}] ; N50 value is not greater than the specified minimum value [{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testGenomeSizeLargerThanMaximum(self):
        file = path.join(self.test_data_dir, "test-Genome-Size-Larger-Than-Maximum.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-Genome-Size-Larger-Than-Maximum', quality_module.index[0], 'File name not equal')
        self.assertEqual(700, quality_module['Genome Length'].iloc[0], 'Genome length not equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}] ; N50 value is not greater than the specified minimum value [{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testGenomeSizeUnaffectedByEmptyContigs(self):
        file = path.join(self.test_data_dir, "test-Genome-Size-Unaffected-By-Empty-Contigs.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-Genome-Size-Unaffected-By-Empty-Contigs', quality_module.index[0], 'File name not equal')
        self.assertEqual(600, quality_module['Genome Length'].iloc[0], 'Genome length not equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('N50 value is not greater than the specified minimum value [{}]'.format(self.minimum_N50_value), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')
    
    #NCO=Number of contigs with length above minimum contig length	NCU=Number of contigs with length under minimum contig length	NC=Number of contigs whose length is exactly minimum contig length
    def testNCExactlyUnacceptable(self):
        file = path.join(self.test_data_dir, "test-NC-Exactly-Unacceptable.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-NC-Exactly-Unacceptable', quality_module.index[0], 'File name not equal')
        self.assertEqual(10, quality_module['Number of Contigs Greater Than Or Equal To '+ str(self.minimum_contig_length) +' bp'].iloc[0], 'Number of Contigs Greater Than Or Equal To Our Minimum Contig Length Not Equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}] ; Number of Contigs with a length greater than or equal to the minimum Contig length [{}] exceeds the acceptable number [{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_contig_length,self.unacceptable_num_contigs), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testNCUExactlyUnacceptable(self):
        file = path.join(self.test_data_dir, "test-NCU-Exactly-Unacceptable.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-NCU-Exactly-Unacceptable', quality_module.index[0], 'File name not equal')
        self.assertEqual(0, quality_module['Number of Contigs Greater Than Or Equal To '+ str(self.minimum_contig_length) +' bp'].iloc[0], 'Number of Contigs Greater Than Or Equal To Our Minimum Contig Length Not Equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testNCOneUnderUnacceptable(self):
        file = path.join(self.test_data_dir, "test-NC-One-Under-Unacceptable.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-NC-One-Under-Unacceptable', quality_module.index[0], 'File name not equal')
        self.assertEqual(9, quality_module['Number of Contigs Greater Than Or Equal To '+ str(self.minimum_contig_length) +' bp'].iloc[0], 'Number of Contigs Greater Than Or Equal To Our Minimum Contig Length Not Equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testNCOOneUnderUnacceptable(self):
        file = path.join(self.test_data_dir, "test-NCO-One-Under-Unacceptable.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-NCO-One-Under-Unacceptable', quality_module.index[0], 'File name not equal')
        self.assertEqual(9, quality_module['Number of Contigs Greater Than Or Equal To '+ str(self.minimum_contig_length) +' bp'].iloc[0], 'Number of Contigs Greater Than Or Equal To Our Minimum Contig Length Not Equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testNCUnacceptableByEmptyContigs(self):
        file = path.join(self.test_data_dir, "test-NC-Unacceptable-By-Empty-Contigs.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-NC-Unacceptable-By-Empty-Contigs', quality_module.index[0], 'File name not equal')
        self.assertEqual(9, quality_module['Number of Contigs Greater Than Or Equal To '+ str(self.minimum_contig_length) +' bp'].iloc[0], 'Number of Contigs Greater Than Or Equal To Our Minimum Contig Length Not Equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testNCOMuchLowerThanUnacceptable(self):
        file = path.join(self.test_data_dir, "test-NCO-Much-Lower-Than-Unacceptable.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-NCO-Much-Lower-Than-Unacceptable', quality_module.index[0], 'File name not equal')
        self.assertEqual(1, quality_module['Number of Contigs Greater Than Or Equal To '+ str(self.minimum_contig_length) +' bp'].iloc[0], 'Number of Contigs Greater Than Or Equal To Our Minimum Contig Length Not Equal')
        self.assertEqual('Passed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('', quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testNCUMuchHigherThanUnacceptable(self):
        file = path.join(self.test_data_dir, "test-NCU-Much-Higher-Than-Unacceptable.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-NCU-Much-Higher-Than-Unacceptable', quality_module.index[0], 'File name not equal')
        self.assertEqual(0, quality_module['Number of Contigs Greater Than Or Equal To '+ str(self.minimum_contig_length) +' bp'].iloc[0], 'Number of Contigs Greater Than Or Equal To Our Minimum Contig Length Not Equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}] ; N50 value is not greater than the specified minimum value [{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')

    def testNCOMuchHigherThanUnacceptable(self):
        file = path.join(self.test_data_dir, "test-NCO-Much-Higher-Than-Unacceptable.fasta")
        files = [file]
        quality_module = QualityModule(files,self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_N50_value,self.minimum_contig_length,self.unacceptable_num_contigs)
        quality_module=quality_module.create_quality_module_dataframe()
        self.assertEqual(1, len(quality_module.index), 'Invalid number of rows in results')
        self.assertEqual('test-NCO-Much-Higher-Than-Unacceptable', quality_module.index[0], 'File name not equal')
        self.assertEqual(110, quality_module['Number of Contigs Greater Than Or Equal To '+ str(self.minimum_contig_length) +' bp'].iloc[0], 'Number of Contigs Greater Than Or Equal To Our Minimum Contig Length Not Equal')
        self.assertEqual('Failed', quality_module['Quality Module'].iloc[0], 'Quality result not equal')
        self.assertEqual('Genome length is not within the acceptable length range [{},{}] ; Number of Contigs with a length greater than or equal to the minimum Contig length [{}] exceeds the acceptable number [{}]'.format(self.genome_size_lower_bound,self.genome_size_upper_bound,self.minimum_contig_length,self.unacceptable_num_contigs), quality_module['Quality Module Feedback'].iloc[0], 'Quality feedback not equal')