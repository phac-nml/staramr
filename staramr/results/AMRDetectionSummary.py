from os import path

import pandas as pd
from pandas import DataFrame

from typing import List

import logging
logger = logging.getLogger("AMRDetectionSummary")

"""
Summarizes both ResFinder, PointFinder, and PlasmidFinder database results into a single table.
"""


class AMRDetectionSummary:
    SEPARATOR = ','

    def __init__(self, files, resfinder_dataframe: DataFrame, pointfinder_dataframe: DataFrame=None, plasmidfinder_dataframe: DataFrame=None) -> None:
        """
        Constructs an object for summarizing AMR detection results.
        :param files: The list of genome files we have scanned against.
        :param resfinder_dataframe: The pd.DataFrame containing the ResFinder results.
        :param pointfinder_dataframe: The pd.DataFrame containing the PointFinder results.
        """
        self._names = [path.splitext(path.basename(x))[0] for x in files]
        self._resfinder_dataframe = resfinder_dataframe
        self._plasmidfinder_dataframe = plasmidfinder_dataframe

        if pointfinder_dataframe is not None:
            self._has_pointfinder = True
            self._pointfinder_dataframe = pointfinder_dataframe
        else:
            self._has_pointfinder = False

    def _compile_results(self, resistance_frame: DataFrame) -> DataFrame:
        df_summary = resistance_frame.sort_values(by=['Gene']).groupby(['Isolate ID']).aggregate(
            lambda x: {'Gene': (self.SEPARATOR + ' ').join(x['Gene'])})
        return df_summary[['Gene']]

    def _compile_plasmids(self, plasmid_frame: DataFrame) -> DataFrame:
        ds_summary = plasmid_frame.sort_values(by=['Gene']).groupby(['Isolate ID']).aggregate(
            lambda x: {'Gene': (self.SEPARATOR + ' ').join(x['Gene'])})

        ds_frame = ds_summary[['Gene']]

        plasmid_frame = ds_frame.rename(columns={'Gene': 'Plasmid Genes'})

        return plasmid_frame

    def _include_negatives(self, resistance_frame: DataFrame) -> DataFrame:
        result_names_set = set(resistance_frame.index.tolist())
        names_set = set(self._names)

        negative_names_set = names_set - result_names_set
        negative_entries = pd.DataFrame([[x, 'None'] for x in negative_names_set],
                                        columns=('Isolate ID', 'Gene')).set_index('Isolate ID')
        
        return resistance_frame.append(negative_entries, sort=True)

    def _include_detailed_negatives(self, resistance_frame: DataFrame, plasmid_frame: DataFrame) -> DataFrame:
        resfinder_names_set = set(resistance_frame.index.tolist())
        plasmid_frame = self._compile_plasmids(plasmid_frame)
        plasmidfinder_names_set = set(plasmid_frame.index.tolist())

        names_set = set(self._names)

        negative_res_names_set = names_set - resfinder_names_set
        negative_plasmid_names_set = names_set - plasmidfinder_names_set
        set_used = None

        negative_resistance_entries = pd.DataFrame([[x, 'None', 'Sensitive'] for x in negative_res_names_set],
                                        columns=('Isolate ID', 'Gene', 'Predicted Phenotype')).set_index('Isolate ID')
        negative_resistance_entries['Data Type']='Resistance'
        negative_entries = negative_resistance_entries

        if plasmid_frame.empty:
            set_used = names_set
        else:
            set_used = negative_plasmid_names_set

        negative_plasmid_entries = pd.DataFrame([[x, 'None'] for x in set_used],
                                        columns=('Isolate ID', 'Gene')).set_index('Isolate ID')
        negative_plasmid_entries['Data Type']='Plasmid'
        negative_entries = negative_entries.append(negative_plasmid_entries, sort=True)

        return resistance_frame.append(negative_entries, sort=True)

    def create_summary(self, include_negatives: bool=False) -> DataFrame:
        """
        Constructs a summary pd.DataFrame for all ResFinder/PointFinder/PlasmidFinder results.
        :param include_negatives: If True, include files with no ResFinder/PointFinder/PlasmidFinder results.
        :return: A pd.DataFrame summarizing the results.
        """
        resistance_frame = self._resfinder_dataframe
        plasmid_frame = self._plasmidfinder_dataframe

        if self._has_pointfinder:
            resistance_frame = resistance_frame.append(self._pointfinder_dataframe, sort=True)

        resistance_frame = self._compile_results(resistance_frame)

        if include_negatives:
            resistance_frame = self._include_negatives(resistance_frame)

        resistance_frame.rename(columns={'Gene': 'Genotype'}, inplace=True)

        if plasmid_frame is not None:
            plasmid_frame = self._compile_plasmids(plasmid_frame)

            resistance_frame = resistance_frame.merge(plasmid_frame, on='Isolate ID', how='left').fillna(value={'Plasmid Genes': 'None'})

        return resistance_frame.sort_index()

    def create_detailed_summary(self, include_negatives: bool=True) -> DataFrame:
        resistance_frame = self._resfinder_dataframe
        resistance_frame['Data Type']='Resistance'
        plasmid_frame = self._plasmidfinder_dataframe

        column_names = ['Gene', 'Predicted Phenotype','%Identity', '%Overlap', 'HSP Length/Total Length','Contig', 'Start', 'End', 'Accession', 'Data Type']

        plasmid_frame = plasmid_frame.reindex(columns=column_names)

        if self._has_pointfinder:
            point_frame = self._pointfinder_dataframe
            point_frame['Data Type']='Resistance'
            point_frame = point_frame.reindex(columns=column_names)
            resistance_frame = resistance_frame.append(point_frame, sort=True)

        if include_negatives:
            resistance_frame = self._include_detailed_negatives(resistance_frame, plasmid_frame)
            resistance_frame = resistance_frame.fillna(" ")

        if plasmid_frame is not None:
            plasmid_frame['Data Type']='Plasmid'
            plasmid_frame['Predicted Phenotype']=''
            resistance_frame = resistance_frame.append(plasmid_frame, sort=True)
            resistance_frame = resistance_frame.reindex(columns=column_names)
            resistance_frame = resistance_frame.sort_values(['Isolate ID', 'Data Type', 'Gene'])

        return resistance_frame
