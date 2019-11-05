import logging
from os import path
from typing import Set

import pandas as pd
from pandas import DataFrame

logger = logging.getLogger("QualityModule")

"""
Summarizes the quality metric feedback into a single table.
"""


class QualityModule:

    def __init__(self, files, genome_size_lower_bound,genome_size_upper_bound,minimum_N50_value,minimum_contig_length,unacceptable_num_contigs) -> None:
        """
        Constructs an object for summarizing our quality module.
        :param files: The list of genome files we have scanned against.
        :param genome_size_lower_bound: The lower bound for the genome size as defined by the user for quality metrics
        :param genome_size_upper_bound: The upper bound for the genome size as defined by the user for quality metrics
        :param minimum_N50_value: The minimum N50 value as defined by the user for quality metrics
        :param minimum_contig_length: The minimum contig length as defined by the user for quality metrics
        :param unacceptable_num_contigs: The number of contigs in a file, equal to or above our minimum contig length, for which to raise a flag as defined by the user for quality metrics
        """
        self._files = files
        self._genome_size_lower_bound = genome_size_lower_bound
        self._genome_size_upper_bound = genome_size_upper_bound
        self._minimum_N50_value = minimum_N50_value
        self._minimum_contig_length = minimum_contig_length
        self._unacceptable_num_contigs = unacceptable_num_contigs
        
    
    def _create_quality_module_dataframe(self):
        """
        Goes through the files and creates a dataframe consisting of the file's genome length, N50 value and the number of contigs greater or equal to the minimum contig length as
        specified by the quality metrics. It also consists of the feedback for whether or not the file passed the quality metrics and if it didn't feedback on why it failed
        :return: A pd.dataframe containing the genome size, N50 value, number of contigs equal to or above our user defined minimum contig length
        as well as the results of our quality metrics (pass or fail) and the corresponding feedback
        """
        name_set=[]
        for myFile in self._files:
            name_set.append(path.splitext(path.basename(myFile))[0])

        files_contigs_and_genomes_lengths=self._get_files_contigs_and_genomes_lengths(self._files)

        files_genome_length_feedback = self._get_genome_length_feedback(files_contigs_and_genomes_lengths[1],self._genome_size_lower_bound,self._genome_size_upper_bound)
        
        files_N50_value_feedback=self._get_N50_feedback(files_contigs_and_genomes_lengths[0],files_contigs_and_genomes_lengths[1],self._minimum_N50_value)
        
        file_num_contigs_over_minimum_bp_feedback= self._get_num_contigs_over_minimum_bp_feedback(files_contigs_and_genomes_lengths[0],self._minimum_contig_length,self._unacceptable_num_contigs)
        
        quality_module = self._get_quality_module(files_genome_length_feedback,files_N50_value_feedback[1],file_num_contigs_over_minimum_bp_feedback[1])      
        quality_module_feedback = quality_module[0]
        quality_module_result = quality_module[1]
        
        quality_module_frame=pd.DataFrame([[t,u,v,w,x,y] for t,u,v,w,x,y in zip(name_set,files_contigs_and_genomes_lengths[1],files_N50_value_feedback[0],file_num_contigs_over_minimum_bp_feedback[0],quality_module_result,quality_module_feedback)],
            columns=('Isolate ID', 'Genome Length','N50 value','Number of Contigs Greater Than Or Equal To '+str(self._minimum_contig_length)+' bp','Quality Module','Quality Module Feedback')).set_index('Isolate ID')
        
        return quality_module_frame

    def _get_files_contigs_and_genomes_lengths(self,files):
        #This solution was taken almost directly from https://github.com/phac-nml/sistr_cmd and was in no way specifically designed for starAMR
        """
        Goes through the files and determines their genome length as well as the length of each contig
        :param files: The files for which we wish to determine the genome length as well as the length of each contig
        :return: An array where the first element is itself an array where each element represents the corresponding 
        file and is itself an array where each element is the length for the corresponding contig inside of this file.
        The second element is itself an array where each element is the genome length for the corresponding file
        """
        #This array will contain the contig lengths for all the files to be used for our quality metrics, in particular to determine if the files pass our N50
        #metric, and to determine if the files pass the number of contigs equal to or above the minimum contig length metric
        files_contigs_lengths =[]
        #This array will contain the genome lengths for all the files to be used for our quality metrics, in particualr to determine if the files pass our genome 
        #length metric, and to calculate our N50 value 
        files_genomes_lengths =[]
        #This array is what we will return and will contain our files genome lengths and their contig lengths
        feedback = []
        for filepath in files:
            genome_length = 0

            with open(filepath,'r') as g:
                #This array will hold the contig lengths for all contigs in this file to be used for our quality metrics
                contig_lengths = []
                #This variable will hold the genome length for this file to be used for our quality metrics
                length = 0  

                for line in g:
                    line = line.strip()

                    if line == '':
                        continue

                    if line[0] == '>':
                        if length == 0:
                            continue

                        else:
                            genome_length = genome_length + length
                            contig_lengths.append(length)
                            length = 0    

                    else:
                        length = length + len(line)

            contig_lengths.append(length)
            files_contigs_lengths.append(contig_lengths)
            files_genomes_lengths.append(genome_length+length)


        feedback.append(files_contigs_lengths)
        feedback.append(files_genomes_lengths)
        return feedback

    def _get_genome_length_feedback(self,files_genome_lengths,lb_gsize,ub_gsize):
        """
        Goes through the files and determines whether or not they pass the quality metrics for genome length
        :param files_genome_lengths: An array where each element is the genome length of the corresponding file
        :param lb_gsize: The lower bound for the genome size as defined by the user for quality metrics
        :param ub_gsize: The upper bound for the genome size as defined by the user for quality metrics
        :return: An array where each element corresponds to the feedback (true or false) for the corresponding file in regards to the
        genome size quality metric
        """
        #The array will contain the feedback of either True or false for whether or not the corresponding files pass the genome length quality metric, and
        #this feedback will be used to construc our quality module dataframe
        files_genome_feedback=[]
        for genome_length in files_genome_lengths:
            if genome_length >= lb_gsize and genome_length <= ub_gsize:
                files_genome_feedback.append(True)

            else:
                files_genome_feedback.append(False)


        return files_genome_feedback

    def _get_N50_feedback(self,files_contigs_lengths,files_genome_lengths,minimum_N50):
        #For information on what N50 is and how it is calculated see https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics
        """
        Goes through the files and determines whether or not they pass the quality metrics for N50 value
        :param files_contigs_lengths: The lengths of the contigs for the files
        :param files_genome_lengths: An array where each element is the genome length of the corresponding file
        :param minimum_N50_value: The minimum N50 value as defined by the user for quality metrics
        :return: An array where the first element is itself an array where each element is the N50 value for 
        the corresponding file. The second element is itself an array where each element is the feedback (true or false) 
        for whether the corresponding file passes the N50 quality metrics
        """
        #This array is what we will return and will contain the files N50 values and their feedback for whether or not they pass the N50 quality metric
        #and will be used to construct our quality module dataframe
        feedback = []
        #This array will contain the N50 values for the files, it will be returned as the first element of feedback
        files_N50=[]
        #This array will contain the feedback of either True or false for whether or not the corresponding files pass the N50 quality metric, it will be
        #returned as the second element of feedback
        N50_feedback = []
        for file_genome_length,file_contigs_lengths in zip(files_genome_lengths,files_contigs_lengths):
            half_length=(file_genome_length)/2
            file_contigs_lengths.sort()
            contig_num=len(file_contigs_lengths)
            contig_index = 1
            sum_of_largest_contigs=0

            while contig_index < contig_num:
                if sum_of_largest_contigs+file_contigs_lengths[contig_num-contig_index] >=half_length:
                    break

                else: 
                    sum_of_largest_contigs=sum_of_largest_contigs+file_contigs_lengths[contig_num-contig_index]
                    contig_index=contig_index+1

            files_N50.append(file_contigs_lengths[contig_num-contig_index])


        for file_N50_value in files_N50:
            if file_N50_value > minimum_N50:
                N50_feedback.append(True)

            else:
                N50_feedback.append(False)


        feedback.append(files_N50)
        feedback.append(N50_feedback)

        return feedback
    
    def _get_num_contigs_over_minimum_bp_feedback(self,files_contigs_lengths,minimum_contig_length,unacceptable_num_contigs_over_minimum_bp):
        """
        Goes through the files and determines whether or not they pass the quality metrics for the acceptable number of contigs equal to or above the minimum contig length
        :param files_contigs_lengths: The lengths of the contigs for the files
        :param minimum_contig_length: The minimum contig length as defined by the user for quality metrics
        :param unacceptable_num_contigs: The number of contigs in a file, equal to or above our minimum contig length, for which to raise a flag as defined by the user for quality metrics
        :return: An array where the first element is itself an array where each element is the number of contigs equal to or above the minimum contig length for
        the corresponding file. The second element is itself an array where each element is the feedback (True or False) 
        for whether the corresponding file passes the acceptable number of contigs equal to or above the minimum contig length quality metric
        """
        #This array is what we will return and will contain for each of the files the number of contigs of length greater than or equal to the minimum contig length 
        #as well as the feedback for whether or not this number of contigs is greater than or equal to the unacceptable number of contigs and thus wether it passes or fails the quality metric
        #this array will be used to construct our quality module dataframe
        feedback=[]
        #This array will contain the number of contigs of length greater than or equal to the minimum contig length for the files, it will be returned as the first element of feedback
        file_num_contigs=[]
        #This array will contain the feedback of either True or False for whether or not the corresponding files pass the quality metrics for the acceptable number of contigs equal to 
        #or above the minimum contig length, it will be returned as the second element of feedback
        contigs_over_minimum_bp_feedback=[]
        for file_contigs_lengths in files_contigs_lengths:
            num_contigs = 0

            for contig in file_contigs_lengths:
                if contig >= minimum_contig_length:
                    num_contigs = num_contigs+1

            file_num_contigs.append(num_contigs)


        for file_num_contigs_over_minimum_bp in file_num_contigs:
            if file_num_contigs_over_minimum_bp >= unacceptable_num_contigs_over_minimum_bp:
                contigs_over_minimum_bp_feedback.append(False)

            else:
                contigs_over_minimum_bp_feedback.append(True)


        feedback.append(file_num_contigs)
        feedback.append(contigs_over_minimum_bp_feedback)

        return feedback

    def _get_quality_module(self,genome_length_feedback,N50_feedback,contigs_over_minimum_bp_feedback):
        """
        Goes through the files and for each provides detailed feedback for why they failed the quality metrics
        :param genome_length_feedback: An array where each element is the feedback (true or false) for the corresponding file in regards to the
        genome length quality metric
        :param N50_feedback: An array where each element is the feedback (true or false) for the corresponding file in regards to the
        the N50 quality metric
        :param contigs_over_minimum_bp_feedback: An array where each element is the feedback (true or false) for the corresponding file in regards to the
        the acceptable number of contigs equal to or above the minimum contig length quality metric
        :return: An array where the first element is itself an array where each element is the detailed quality metric feedback for
        the corresponding file. The second element is itself an array where each element is the feedback (true or false) 
        for whether the corresponding file passes all of the quality metrics
        """
        #This array is what we will return and will contain the files feedback for whether or not they pass the quality metrics as well as the specific metrics
        #the files fail on, if they fail any of the metrics and will be used to construct our quality module dataframe
        feedback = []
        #This array contains the feedback of either True or False for whether the corresponding file passes all the quality metrics 
        #it will be returned as the first element of feedback
        quality_parameter = []
        #This array contains the detailed feeback or either nothing if the corresponding file passes all the quality metrics or the metrics which the corresponding file
        #fails on and it will be returned as the second element of feedback
        quality_parameter_feedback = []
        for file_genome_length_feedback,file_N50_feedback,file_contigs_over_minimum_bp_feedback in zip(genome_length_feedback,N50_feedback,contigs_over_minimum_bp_feedback):
            if file_genome_length_feedback == True & file_N50_feedback == True & file_contigs_over_minimum_bp_feedback == True:
                quality_parameter_feedback_for_file=""
                quality_parameter.append("Passed")

            else:
                quality_parameter_feedback_for_file=""
                quality_parameter.append("Failed")
                if file_genome_length_feedback == False:
                    quality_parameter_feedback_for_file = quality_parameter_feedback_for_file + "Genome length is not within the acceptable length range. "

                if file_N50_feedback == False:
                    quality_parameter_feedback_for_file = quality_parameter_feedback_for_file + "N50 value is not greater than the specified minimum value. "

                if file_contigs_over_minimum_bp_feedback == False:
                    quality_parameter_feedback_for_file = quality_parameter_feedback_for_file + "Number of Contigs with a length greater than or equal to the minimum Contig length exceeds the acceptable number. "

            quality_parameter_feedback.append(quality_parameter_feedback_for_file)


        feedback.append(quality_parameter_feedback)
        feedback.append(quality_parameter)

        return feedback

