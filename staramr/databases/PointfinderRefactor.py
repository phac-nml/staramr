import os
import timeit
from pathlib import Path
import numpy as np
import pandas as pd
import glob
import logging
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

logger = logging.getLogger('PointfinderMutationsRefactor')

"""
Class to compare fasta file list  and resistance gene list into resistens_overview.txt from pointfinder database
"""

class PointfinderOverviewCheck:

    def __init__(self,poinfinder_dir):
        """
        create the pointfinder path
        :param pointfinder_database_path: the path of pointfinder database
        :param pointfinder_database_name: the name of pointfinder database
        """
        # path to pointfinder DB
        self._poinfinder_dir = poinfinder_dir

    # lambda to extract the name of directories
    def _extract_dir_path(self,dir_name):
        return os.path.dirname(dir_name)
    # lambda to extract the filename without extension
    def _extract_filename(self,file_name):
        return os.path.splitext(os.path.basename(file_name))[0]
        # lambda to check available fasta in repertories
    def _extract_fasta_path(self,fasta_path):
        return glob.glob(fasta_path + "/*.f[a-z]a")

    def _resistance_overview_list(self):
        """
        Extract the list of resistens-overview files in each species subdirectories
        :return: the list of overview files in pandas dataframe format
        """
        resistens_overview_list = list(Path(self._poinfinder_dir).rglob("resistens-overview.txt"))
        resistens_overview_df = pd.DataFrame(resistens_overview_list, columns=["overview_file_path"], dtype="string")
        resistens_overview_df["overview_file_path"] = resistens_overview_df["overview_file_path"].apply(self._extract_dir_path)
        resistens_overview_df["overview_dirname"] = resistens_overview_df["overview_file_path"].apply(self._extract_filename)
        resistens_overview_df["overview_fasta_path_list"] = resistens_overview_df["overview_file_path"].apply(self._extract_fasta_path)
        return resistens_overview_df

    def get_gene_list_df(self):
        """
        :return: the dataframe of resistens overview with file path, species name, a list of fasta file and all available genes
        """
        return self._extract_fasta_list()

    def _extract_fasta_list(self):
        """
        From the panda dataframe, loop on each row, extract gene list and make a dict before to merge by species
        :return: the final dataframe with gene list
        """
        overview_gene_from_fasta_dic = {}
        for overview_key in self._resistance_overview_list()["overview_dirname"].index:
            extract_fasta_list = self._extract_gene_list(overview_key)
            species_name = self._resistance_overview_list()["overview_dirname"].iloc[overview_key]
            overview_gene_from_fasta_dic[species_name] = set(extract_fasta_list)

        overview_gene_from_fasta_df = pd.DataFrame(list(overview_gene_from_fasta_dic.items()),columns=["overview_dirname", "overview_fasta_gene_list"])
        resistens_overview_df = pd.merge(self._resistance_overview_list(), overview_gene_from_fasta_df,on="overview_dirname")
        return resistens_overview_df

    def _extract_gene_list(self,overview_key):
        """
        A loop to extract only filename without file extension
        :param overview_key: index localisation in panda dataframe
        :return: a list of gene determined from the fasta file list
        """
        tempory_gene_list = []
        [tempory_gene_list.append(self._extract_filename(fasta_name)) for fasta_name in self._resistance_overview_list()["overview_fasta_path_list"].iloc[overview_key]]
        return tempory_gene_list

    def _compare_overview_and_fasta_gene(self):
        """
        Extract from the fasta file names, a set of genes for each species directories in PointFinder database
        :return: a new dataframe with a set of each genes in repositories 
        """
        overview_gene_from_overview_dic = {}
        for overview_dir in self.get_gene_list_df()["overview_file_path"]:
            overview_path = overview_dir + "/resistens-overview.txt"
            read_overview_txt = pd.read_csv(overview_path, sep="\t")
            gene_in_overview = set(read_overview_txt["#Gene_ID"][read_overview_txt["#Gene_ID"].str.contains("#") == False].unique())
            overview_gene_from_overview_dic[overview_dir] = gene_in_overview

        overview_gene_from_overview_df = pd.DataFrame(list(overview_gene_from_overview_dic.items()),columns=["overview_file_path", "overview_gene_list"])
        overview_gene_all = pd.merge(self.get_gene_list_df(), overview_gene_from_overview_df,on="overview_file_path")
        return overview_gene_all

    def pandas_set_diff(self):
        """
        Compare the set of fasta gene files to the set of gene in resistens_overview.txt files in each species directories
        :return: final dataframe of all available informations from each species directories 
        """
        gene_list_difference = []
        compare_overview_df = self._compare_overview_and_fasta_gene()
        for row_index in compare_overview_df.index:
            gene_list_from_fasta = compare_overview_df.iloc[row_index]["overview_fasta_gene_list"]
            gene_list_from_file = compare_overview_df.iloc[row_index]["overview_gene_list"]
            gene_list_difference.append(gene_list_from_fasta.symmetric_difference(gene_list_from_file))

        compare_overview_df["fasta_file_difference"] = gene_list_difference
        gene_differences = compare_overview_df['fasta_file_difference']
        compare_overview_df['fasta_file_difference'] = np.where(gene_differences.values == set(), "",gene_differences.values)
        return compare_overview_df

    def fasta_overview_differences(self):
        file_ambiguity = self._pandas_set_diff()[~(self._pandas_set_diff()["fasta_file_difference"] == '')]
        missing_species = str(file_ambiguity['overview_dirname'].values)
        missing_genes = list(file_ambiguity['fasta_file_difference'].values)
        if len(missing_genes) > 0:
            logger.warning("Error for " + missing_species + " species. No correspondance between fasta and overview list for " + str(missing_genes))

"""
Class to control each described mutation in resistens_overview files comparing to fasta sequences
"""

class PointFinderMutationRefactor1(PointfinderOverviewCheck):
    super()__init__(self, poinfinder_dir):

    @classmethod
    def fasta_reader(cls,file):
        fasta_df = pd.read_csv(file, sep='>', lineterminator='>', header=None)
        fasta_df[['Accession', 'Sequence']] = fasta_df[0].str.split(r'\r|\n', 1, expand=True)
        fasta_df.drop(0, axis=1, inplace=True)
        fasta_df['Sequence'] = fasta_df['Sequence'].replace('\n', '', regex=True)
        return(fasta_df)

    @classmethod
    def fasta_to_pandas_df(cls,dataframe_of_fasta_path):
        whole_fasta_df = pd.DataFrame()
        whole_fasta_df = whole_fasta_df.append([PointFinderMutationRefactor.partial_fasta_to_pandas_df(fasta_paths) for fasta_paths in dataframe_of_fasta_path["overview_fasta_path_list"]])
        return whole_fasta_df

    @classmethod
    def partial_fasta_to_pandas_df(cls,partial_dataframe_of_fasta_path):
        tmp_df_each_species = pd.DataFrame()
        tmp_df_each_species = tmp_df_each_species.append([PointFinderMutationRefactor.fasta_reader(one_fasta_path) for one_fasta_path in partial_dataframe_of_fasta_path])
        return tmp_df_each_species







class encours(PointfinderOverviewCheck):

    super()__init__(self, poinfinder_dir):

    def overview_mutation_df(self):
        all_pointfinder_mutation_df = pd.DataFrame()
        overview_file_list = self._resistance_overview_list()
        for overview_index in overview_file_list.index:
            overview_files = overview_file_list.iloc[overview_index]["overview_file_path"] + "/resistens-overview.txt"
            overview_file_df = pd.read_csv(overview_files,sep="\t")
            overview_file_df = overview_file_df[overview_file_df["#Gene_ID"].str.contains("#") == False]
            overview_file_df["Species"] = overview_file_list.iloc[overview_index]["overview_dirname"]
            all_pointfinder_mutation_df = all_pointfinder_mutation_df.append(overview_file_df)
        all_pointfinder_mutation_df["group"] = all_pointfinder_mutation_df["Species"] + "_" + all_pointfinder_mutation_df["#Gene_ID"].str.lower()
        return all_pointfinder_mutation_df

    def get_all_pointfinder_mutations(self):
        return self.all_fasta_df()

    def get_all_pointfinder_mutations(self):
        return self.overview_mutation_df()


    def all_fasta_df(self):
        all_pointfinder_fasta_df = pd.DataFrame()
        fasta_file_list = self._resistance_overview_list()
        for overview_index in fasta_file_list.index:
            overview_files = fasta_file_list.iloc[overview_index]["overview_fasta_path_list"]
            allfasta_df = pd.DataFrame()
            for fasta_file in overview_files:
                final_liste = []
                for fasta_seq in SeqIO.parse(fasta_file, "fasta"):
                    listeo = [fasta_seq.id, str(fasta_seq.translate(table="11", id=True).seq), str(fasta_seq.seq)]
                    final_liste.append(listeo)

                fasta_pandas = pd.DataFrame(final_liste, columns=['Name', 'amino_sequence', 'nucleotide_sequence'])
                allfasta_df = allfasta_df.append(fasta_pandas)
            allfasta_df["Species"] = fasta_file_list.iloc[overview_index]["overview_dirname"]
            all_pointfinder_fasta_df = all_pointfinder_fasta_df.append(allfasta_df)

        all_pointfinder_fasta_df["group"] = all_pointfinder_fasta_df["Species"] + "_" + all_pointfinder_fasta_df["Name"].str.lower()
        return all_pointfinder_fasta_df

    def encore_unpeu(self):
        all_pointfinder_fasta_df
        all_pointfinder_mutation_df





    def read_fasta(fasta_file):
        final_liste = []
        for fasta_seq in SeqIO.parse(fasta_file,"fasta"):
            listeo = [fasta_seq.id,str(fasta_seq.translate(table="11",id=True).seq),str(fasta_seq.seq)]
            final_liste.append(listeo)

        fasta_pandas = pd.DataFrame(final_liste,columns=['Name','amino_sequence','nucleotide_sequence'])
        return fasta_pandas




data = PointfinderOverviewCheck(overview)

file="/home/pierre/PycharmProjects/pythonProject/staramr_update/staramr/staramr/databases/data/update/pointfinder/mycobacterium_tuberculosis/ahpC_promoter_size_180bp.fsa"
file="/home/pierre/PycharmProjects/pythonProject/staramr_update/staramr/staramr/databases/data/update/pointfinder/enterococcus_faecium/resistens-overview.txt"
overview = "/home/pierre/PycharmProjects/pythonProject/staramr_update/staramr/staramr/databases/data/update/pointfinder/"
file="/home/pierre/PycharmProjects/pythonProject/staramr_update/staramr/staramr/databases/data/update/pointfinder/mycobacterium_tuberculosis/resistens-overview.txt"
file = "/home/pierre/PycharmProjects/pythonProject/staramr_update/staramr/staramr/databases/data/update/pointfinder/mycobacterium_tuberculosis/gyrA.fsa"



    def refactor(self):
        print(self.pointfinder_database_path,"___", self.pointfinder_database_name)
        chemin = os.path.join(self.pointfinder_database_path, self.pointfinder_database_name)
        for file in os.listdir(chemin):
            print("\nespèce :", file)
            path = os.path.join(chemin, file)
            if os.path.isdir(path):
                print("fzef")
                for file2 in os.listdir(path):
                    print("espèce")
                    print(file2)

                    for line in reader:
                        gene = line[gene]


                # traitement
folp = Seq("ATGAAACTCTTTGCCCAGGGTACTTCACTGGACCTTAGCCATCCTCACGTAATGGGGATCCTCAACGTCACGCCTGATTCCTTTTCGGATGGTGGCACGCATAACTCGCTGATAGATGCGGTGAAACATGCGAATCTGATGATCAACGCTGGCGCGACGATCATTGACGTTGGTGGCGAGTCCACGCGCCCAGGGGCGGCGGAAGTTAGCGTTGAAGAAGAGTTGCAACGTGTTATTCCTGTGGTTGAGGCAATTGCTCAACGCTTCGAAGTCTGGATCTCAGTCGATACATCCAAACCAGAAGTCATCCGTGAGTCAGCGAAAGTTGGCGCTCACATTATTAATGATATCCGCTCCCTTTCCGAACCTGGCGCTCTGGAGGCGGCTGCAGAAACCGGTTTACCGGTTTGTCTGATGCATATGCAGGGAAATCCAAAAACCATGCAGGAAGCTCCGAAGTATGACGATGTCTTTGCAGAAGTGAATCGCTACTTTATTGAGCAAATAGCACGTTGCGAGCAGGCGGGTATCGCAAAAGAGAAATTGTTGCTCGACCCCGGATTCGGTTTCGGTAAAAATCTCTCCCATAACTATTCATTACTGGCGCGCCTGGCTGAATTTCACCATTTCAACCTGCCGCTGTTGGTGGGTATGTCACGAAAATCGATGATTGGGCAGCTGCTGAACGTGGGGCCGTCCGAGCGCCTGAGCGGTAGTCTGGCCTGTGCGGTCATTGCCGCAATGCAAGGCGCGCACATCATTCGTGTTCATGACGTCAAAGAAACCGTAGAAGCGATGCGGGTGGTGGAAGCCACTCTGTCTGCAAAGGAAAACAAACGCTATGAGTAA")
