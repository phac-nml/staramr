import logging
from os import path
from typing import Set

import pandas as pd
from pandas import DataFrame

logger = logging.getLogger("AMRDetectionSummary")

"""
Summarizes both ResFinder, PointFinder, and PlasmidFinder database results into a single table.
"""


class AMRDetectionSummary:
    SEPARATOR = ','
    FLOAT_DECIMALS = 2

    def __init__(self, files, resfinder_dataframe: DataFrame, pointfinder_dataframe=None,
                 plasmidfinder_dataframe=None, mlst_dataframe=None) -> None:
        """
        Constructs an object for summarizing AMR detection results.
        :param files: The list of genome files we have scanned against.
        :param resfinder_dataframe: The pd.DataFrame containing the ResFinder results.
        :param pointfinder_dataframe: The pd.DataFrame containing the PointFinder results.
        """
        self._names = [path.splitext(path.basename(x))[0] for x in files]
        self._resfinder_dataframe = resfinder_dataframe
        self._plasmidfinder_dataframe = plasmidfinder_dataframe
        self._mlst_dataframe = mlst_dataframe
        self.files=files
        self.num_files = len(self.files)

        if pointfinder_dataframe is not None:
            self._has_pointfinder = True
            self._pointfinder_dataframe = pointfinder_dataframe
        else:
            self._has_pointfinder = False
    def parse_fasta(self,filepath):
        #This solution was taken directly from https://github.com/phac-nml/sistr_cmd and was in no way desogned specifically designed for starAMR
        '''
        Parse a fasta file returning a generator yielding tuples of fasta headers to sequences.
        Note:
            This function should give equivalent results to SeqIO from BioPython
            .. code-block:: python
                from Bio import SeqIO
            # biopython to dict of header-seq
            hseqs_bio = {r.description:str(r.seq) for r in SeqIO.parse(fasta_path, 'fasta')}
            # this func to dict of header-seq
            hseqs = {header:seq for header, seq in parse_fasta(fasta_path)}
            # both methods should return the same dict
            assert hseqs == hseqs_bio
        Args:
            filepath (str): Fasta file path
        Returns:
            generator: yields tuples of (<fasta header>, <fasta sequence>)
        '''
        with open(filepath, 'r') as f:
            seqs = []
            header = ''
            for line in f:
                line = line.strip()
                if line == '':
                    continue
                if line[0] == '>':
                    if header == '':
                        header = line.replace('>','')
                    else:
                        yield header, ''.join(seqs)
                        seqs = []
                        header = line.replace('>','')
                else:
                    seqs.append(line)
            yield header, ''.join(seqs)
    
    def _get_files_contigs_lengths(self):
        #Goes through each file in self.files and for each file determines the length of each contig. 
        #Returns an array where each element represents a file and is itself an array of the contig lengths 
        files_contigs_lengths =[]
        for filepath in self.files:
            with open(filepath,'r') as g:
                contig_lengths = []
                length = 0
                for line in g:
                    line = line.strip()
                    if line == '':
                        continue
                    if line[0] == '>':
                        if length == 0:
                            continue
                        else:
                            contig_lengths.append(length)
                            length = 0        
                    else:
                        length = length + len(line)
            contig_lengths.append(length)
            files_contigs_lengths.append(contig_lengths)
        return files_contigs_lengths

    def _get_genome_lengths(self):
        #Goes through each file in self.files and for each file determines the length of the genome. 
        #Returns an array where each element is the genome length of the corresponding file
        files_genome_lengths=[]
        for myFile in self.files:
            parsedFile=self.parse_fasta(myFile)
            genome_size = sum([len(s) for h, s in parsedFile])
            files_genome_lengths.append(genome_size)
        return files_genome_lengths

    def _get_genome_length_feedback(self,files_genome_lengths):
        #Takes as input an array where each element is the genome length of a corresponding file. 
        #Returns an array where each elements is the feedback of whether or not the genome length for the corresponding file is between 4 Mbp and 6 Mbp
        lb_gsize= 4000000
        ub_gsize = 6000000
        feedback=[]
        for genome_length in files_genome_lengths:
            is_gsize_acceptable = (genome_length >= lb_gsize and genome_length <= ub_gsize)
            if is_gsize_acceptable is True:
                feedback.append("Passed as Genome size is "+str(genome_length/1000000)+" Mbp which is between 4 Mbp and 6 Mbp")
            else:
                feedback.append("Failed as Genome size is "+str(genome_length/1000000)+" Mbp which is not between 4 Mbp and 6 Mbp")
        return feedback

    def _get_N50(self,files_contigs_lengths,files_genome_lengths):
        #Takes as input an array where each element represents a file and is itself an array where each element corresponds to the length
        #of a corresponding contig. It also takes as input an array where each element is the genome length of a corresponding file. 
        #Returns an array where each element is the N50 value for a corresponding file
        #For information on what N50 is and how it is calculated see https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics
    
        file_index =0
        file_N50=[]
        while file_index < self.num_files:
            contig_lengths = files_contigs_lengths[file_index]
            half_length=(files_genome_lengths[file_index])/2
            contig_lengths.sort()
            contig_num=len(contig_lengths)
            contig_index = 1
            sum_of_largest_contigs=0
            while contig_index < contig_num:
                if sum_of_largest_contigs+contig_lengths[contig_num-contig_index] >=half_length:
                    break
                else: 
                    sum_of_largest_contigs=sum_of_largest_contigs+contig_lengths[contig_num-contig_index]
                    contig_index=contig_index+1
            file_N50.append(contig_lengths[contig_num-contig_index])
            file_index = file_index +1
        return file_N50

    def _get_N50_feedback(self,N50_values:[]):
        #Takes as input an array where each element is the N50 value of a corresponding file. 
        #Returns an array where each elements is the feedback of whether or not the N50 length for the corresponding file is over 10 000
        N50_feedback = []
        for file_N50_value in N50_values:
            if file_N50_value > 10000:
                N50_feedback.append("Passed as N50 value is "+str(file_N50_value)+" which is greater than 10000")
            else:
                N50_feedback.append("Failed as N50 value is "+str(file_N50_value)+" which is less than or equal to 10000")
        return N50_feedback
    
    def _get_num_contigs_under_minimum_bp(self,files_contigs_lengths,minimum_contig_length):
        #Takes as input an array where each element represents a file and is itself an array where each element corresponds to the length
        #of a corresponding contig. It also take in as input the predefined minimum length that a contig can be without raising concern. 
        #Returns an array where each element denotes the number of contigs below this predefined minimum contig length for the corresponding file
        file_num_contigs=[]
        file_index = 0
        while file_index < self.num_files:
            num_contigs = 0
            for contig in files_contigs_lengths[file_index]:
                if contig < minimum_contig_length:
                    num_contigs = num_contigs+1
            file_num_contigs.append(num_contigs)
            file_index=file_index+1
        return file_num_contigs
    
    def _get_num_contigs_under_minimum_bp_feedback(self,num_contigs_under_minimum_bp,minimum_contig_length,unacceptable_num_contigs_under_minimum_bp):
        #Takes as input an array where each element denotes the number of contigs below our predefined minimum contig length for the corresponding file.
        #It also takes as input the predefined minimum length that a contig can be without raising concern. Finally, it takes as input the minimum
        #number of contigs which raise concern in order for us to determine that a file fails. 
        #Returns an array where each elements is the feedback of whether or not the corresponding file has too many contigs that are smaller than our predefined 
        #minimum contig length
        contigs_under_minimum_bp_feedback=[]
        for file_num_contigs_under_minimum_bp in num_contigs_under_minimum_bp:
            if file_num_contigs_under_minimum_bp >= unacceptable_num_contigs_under_minimum_bp:
                contigs_under_minimum_bp_feedback.append("Failed as the number of contigs less than "+str(minimum_contig_length)+" bp is "+str(file_num_contigs_under_minimum_bp)+" which is greater than or equal to "+str(unacceptable_num_contigs_under_minimum_bp))
            else:
                contigs_under_minimum_bp_feedback.append("Passed as the number of contigs less than "+str(minimum_contig_length)+" bp is "+str(file_num_contigs_under_minimum_bp)+" which is less than "+str(unacceptable_num_contigs_under_minimum_bp))
        return contigs_under_minimum_bp_feedback

    def _compile_results(self, resistance_frame: DataFrame) -> DataFrame:
        df_summary = resistance_frame.sort_values(by=['Gene']).groupby(['Isolate ID']).aggregate(
            lambda x: {'Gene': (self.SEPARATOR + ' ').join(x['Gene'])})
        return df_summary[['Gene']]

    def _compile_plasmids(self, plasmid_frame: DataFrame) -> DataFrame:
        ds_summary = plasmid_frame.sort_values(by=['Gene']).groupby(['Isolate ID']).aggregate(
            lambda x: {'Gene': (self.SEPARATOR + ' ').join(x['Gene'])})

        ds_frame = ds_summary[['Gene']]

        plasmid_frame = ds_frame.rename(columns={'Gene': 'Plasmid'})

        return plasmid_frame

    def _include_negatives(self, resistance_frame: DataFrame) -> DataFrame:
        result_names_set = set(resistance_frame.index.tolist())
        names_set = set(self._names)

        negative_names_set = names_set - result_names_set
        negative_entries = pd.DataFrame([[x, 'None'] for x in negative_names_set],
                                        columns=('Isolate ID', 'Gene')).set_index('Isolate ID')

        return resistance_frame.append(negative_entries, sort=True)

    def _get_detailed_negative_columns(self):
        return ['Isolate ID', 'Gene', 'Start', 'End']

    def _include_phenotype(self):
        return False

    def _get_negative_resistance_entries(self, names_set: Set, resistance_frame: DataFrame) -> DataFrame:
        resfinder_names_set = set(resistance_frame.index.tolist())
        negative_res_names_set = names_set - resfinder_names_set
        negative_columns = self._get_detailed_negative_columns()

        negative_resistance_entries = None

        if len(negative_res_names_set) != len(names_set) or resistance_frame.empty:

            if self._include_phenotype():
                negative_resistance_entries = pd.DataFrame(
                    [[x, 'None', 'Sensitive', '', ''] for x in negative_res_names_set],
                    columns=negative_columns).set_index('Isolate ID')
            else:
                negative_resistance_entries = pd.DataFrame([[x, 'None', '', ''] for x in negative_res_names_set],
                                                           columns=negative_columns).set_index('Isolate ID')

            negative_resistance_entries['Data Type'] = 'Resistance'

        return negative_resistance_entries

    def _include_detailed_negatives(self, resistance_frame: DataFrame, plasmid_frame: DataFrame = None) -> DataFrame:
        names_set = set(self._names)
        set_used = names_set

        negative_entries = self._get_negative_resistance_entries(names_set, resistance_frame)

        if plasmid_frame is not None:
            plasmid_frame = self._compile_plasmids(plasmid_frame)
            plasmidfinder_names_set = set(plasmid_frame.index.tolist())
            negative_plasmid_names_set = names_set - plasmidfinder_names_set

            if not plasmid_frame.empty:
                set_used = negative_plasmid_names_set

        negative_plasmid_entries = pd.DataFrame([[x, 'None'] for x in set_used],
                                                columns=('Isolate ID', 'Gene')).set_index('Isolate ID')
        negative_plasmid_entries['Data Type'] = 'Plasmid'

        if negative_entries is None:
            negative_entries = negative_plasmid_entries
        else:
            negative_entries = negative_entries.append(negative_plasmid_entries, sort=True)

        return resistance_frame.append(negative_entries, sort=True)

    def _get_summary_empty_values(self):
        return {'Genotype': 'None'}

    def _get_summary_resistance_columns(self):
        return ['Genotype', 'Plasmid']

    def create_summary(self, include_negatives: bool = False) -> DataFrame:
        """
        Constructs a summary pd.DataFrame for all ResFinder/PointFinder/PlasmidFinder results.
        :param include_negatives: If True, include files with no ResFinder/PointFinder/PlasmidFinder results.
        :return: A pd.DataFrame summarizing the results.
        """
        resistance_frame = self._resfinder_dataframe
        plasmid_frame = self._plasmidfinder_dataframe
        mlst_frame = self._mlst_dataframe

        if self._has_pointfinder:
            resistance_frame = resistance_frame.append(self._pointfinder_dataframe, sort=True)

        resistance_frame = self._compile_results(resistance_frame)

        if include_negatives:
            resistance_frame = self._include_negatives(resistance_frame)

        resistance_frame.rename(columns={'Gene': 'Genotype'}, inplace=True)

        fill_values = self._get_summary_empty_values()
        resistance_columns = self._get_summary_resistance_columns()

        if plasmid_frame is not None:
            plasmid_frame = self._compile_plasmids(plasmid_frame)

            if resistance_frame.empty:
                resistance_frame = resistance_frame.append(plasmid_frame)
            else:
                resistance_frame = resistance_frame.merge(plasmid_frame, on='Isolate ID', how='left').fillna(
                    value={'Plasmid': 'None'})

            resistance_frame = resistance_frame.fillna(value=fill_values)
            resistance_frame = resistance_frame.reindex(columns=resistance_columns)

        if mlst_frame is not None:
            mlst_merging_frame = mlst_frame[['Scheme', 'Sequence Type']]
            resistance_frame = resistance_frame.merge(mlst_merging_frame, on='Isolate ID', how='left')

        name_set=[]
        for myFile in self.files:
            name_set.append(path.splitext(path.basename(myFile))[0])
        files_genome_lengths = self._get_genome_lengths()
        files_genome_length_feedback = self._get_genome_length_feedback(files_genome_lengths)
        files_contigs_lengths=self._get_files_contigs_lengths()
        files_N50_values=self._get_N50(files_contigs_lengths,files_genome_lengths)
        files_N50_feedback=self._get_N50_feedback(files_N50_values)
        minimum_contig_length=1000
        files_contigs_under_minimum_bp= self._get_num_contigs_under_minimum_bp(files_contigs_lengths,minimum_contig_length)
        unacceptable_num_contigs_under_minimum_bp= 1000
        file_num_contigs_under_minimum_bp_feedback= self._get_num_contigs_under_minimum_bp_feedback(files_contigs_under_minimum_bp,minimum_contig_length,unacceptable_num_contigs_under_minimum_bp)

        quality_module_frame=pd.DataFrame([[t,u,v,w,x,y,z] for t,u,v,w,x,y,z in zip(name_set,files_genome_lengths,files_genome_length_feedback,files_N50_values, files_N50_feedback,files_contigs_under_minimum_bp,file_num_contigs_under_minimum_bp_feedback)],
            columns=('Isolate ID', 'Genome Length','Genome Length Feedback','N50 value','N50 Feedback','Number of Contigs Under '+str(minimum_contig_length)+' bp','Number of Contigs Under '+str(minimum_contig_length)+' bp Feedback')).set_index('Isolate ID')

        resistance_frame = resistance_frame.merge(quality_module_frame, on='Isolate ID', how='left')
        return resistance_frame.sort_index()

    def _get_detailed_summary_columns(self):
        return ['Gene', 'Data Type', '%Identity', '%Overlap', 'HSP Length/Total Length', 'Contig', 'Start', 'End', 'Accession']

    def create_detailed_summary(self, include_negatives: bool = True) -> DataFrame:
        mlst_merging_frame = None

        if self._mlst_dataframe is None:
            mlst_frame = None
        else:
            mlst_frame = self._mlst_dataframe.copy()
            mlst_merging_frame = mlst_frame.loc[:, ('Scheme', 'Sequence Type')]
            mlst_merging_frame.loc[:, ('Gene')] = 'ST' + mlst_frame['Sequence Type'].map(str) + ' (' + mlst_frame['Scheme'] + ')'
            mlst_merging_frame.loc[:, ('Data Type')] = 'MLST'
            mlst_merging_frame = mlst_merging_frame.loc[:, ('Gene', 'Data Type')]

        if self._resfinder_dataframe is None:
            resistance_frame = None
        else:
            resistance_frame = self._resfinder_dataframe.copy()
            resistance_frame['Data Type'] = 'Resistance'
            resistance_frame = resistance_frame.round(
                {'%Identity': self.FLOAT_DECIMALS, '%Overlap': self.FLOAT_DECIMALS})

        if self._plasmidfinder_dataframe is None:
            plasmid_frame = None
        else:
            plasmid_frame = self._plasmidfinder_dataframe.copy()

        column_names = self._get_detailed_summary_columns()

        if self._has_pointfinder:
            if self._pointfinder_dataframe is None:
                point_frame = None
            else:
                point_frame = self._pointfinder_dataframe.copy()
                point_frame['Data Type'] = 'Resistance'
                point_frame = point_frame.round({'%Identity': self.FLOAT_DECIMALS, '%Overlap': self.FLOAT_DECIMALS})
                point_frame = point_frame.reindex(columns=column_names)

            if resistance_frame is not None:
                resistance_frame = resistance_frame.append(point_frame, sort=True)

        if include_negatives:
            if plasmid_frame is not None:
                plasmid_frame = plasmid_frame.reindex(columns=column_names)
                resistance_frame = self._include_detailed_negatives(resistance_frame, plasmid_frame)
            else:
                resistance_frame = self._include_detailed_negatives(resistance_frame)
            resistance_frame = resistance_frame.reindex(columns=column_names)

        if plasmid_frame is not None:
            plasmid_frame['Data Type'] = 'Plasmid'

            if self._include_phenotype():
                plasmid_frame['Predicted Phenotype'] = ''

            plasmid_frame = plasmid_frame.round({'%Identity': self.FLOAT_DECIMALS, '%Overlap': self.FLOAT_DECIMALS})

            if resistance_frame is not None:
                resistance_frame = resistance_frame.append(plasmid_frame, sort=True)
                resistance_frame = resistance_frame.reindex(columns=column_names)
                resistance_frame = resistance_frame.sort_values(['Isolate ID', 'Data Type', 'Gene'])

        if mlst_merging_frame is not None and resistance_frame is not None:
            resistance_frame = resistance_frame.append(mlst_merging_frame, sort=True)
            resistance_frame = resistance_frame.reindex(columns=column_names)
            resistance_frame = resistance_frame.sort_values(['Isolate ID', 'Data Type', 'Gene'])

        if resistance_frame is not None:
            resistance_frame = resistance_frame.fillna("")

        return resistance_frame
