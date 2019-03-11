from os import path

import pandas as pd

import logging

logger = logging.getLogger("AMRDetectionSummary")

"""
Summarizes both ResFinder and PointFinder database results into a single table.
"""


class AMRDetectionSummary:
    SEPARATOR = ','

    def __init__(self, files, resfinder_dataframe, pointfinder_dataframe=None, plasmidfinder_dataframe=None):
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

    def _compile_results(self, df):
        df_summary = df.sort_values(by=['Gene']).groupby(['Isolate ID']).aggregate(
            lambda x: {'Gene': (self.SEPARATOR + ' ').join(x['Gene'])})
        return df_summary[['Gene']]

    def _compile_plasmids(self, ds):
        ds_summary = ds.sort_values(by=['Gene']).groupby(['Isolate ID']).aggregate(
            lambda x: {'Gene': (self.SEPARATOR + ' ').join(x['Gene'])})

        ds_frame = ds_summary[['Gene']]

        plasmid_frame = ds_frame.rename(columns={'Gene': 'Plasmid Genes'})

        return plasmid_frame

    def _include_negatives(self, df):
        result_names_set = set(df.index.tolist())
        names_set = set(self._names)

        negative_names_set = names_set - result_names_set
        negative_entries = pd.DataFrame([[x, 'None'] for x in negative_names_set],
                                        columns=('Isolate ID', 'Gene')).set_index('Isolate ID')
        
        return df.append(negative_entries, sort=True)

    def _include_detailed_negatives(self, df, ds):
        resfinder_names_set = set(df.index.tolist())
        ds = self._compile_plasmids(ds)
        plasmidfinder_names_set = set(ds.index.tolist())

        names_set = set(self._names)

        negative_res_names_set = names_set - resfinder_names_set
        negative_plasmid_names_set = names_set - plasmidfinder_names_set
        set_used = None

        negative_resistance_entries = pd.DataFrame([[x, 'None', 'Sensitive'] for x in negative_res_names_set],
                                        columns=('Isolate ID', 'Gene', 'Predicted Phenotype')).set_index('Isolate ID')
        negative_resistance_entries['Data Type']='Resistance'
        negative_entries = negative_resistance_entries

        if ds.empty:
            set_used = names_set
        else:
            set_used = negative_plasmid_names_set

        negative_plasmid_entries = pd.DataFrame([[x, 'None'] for x in set_used],
                                        columns=('Isolate ID', 'Gene')).set_index('Isolate ID')
        negative_plasmid_entries['Data Type']='Plasmid'
        negative_entries = negative_entries.append(negative_plasmid_entries, sort=True)

        return df.append(negative_entries, sort=True)

    def create_summary(self, include_negatives=False):
        """
        Constructs a summary pd.DataFrame for all ResFinder/PointFinder/PlasmidFinder results.
        :param include_negatives: If True, include files with no ResFinder/PointFinder/PlasmidFinder results.
        :return: A pd.DataFrame summarizing the results.
        """
        df = self._resfinder_dataframe
        ds = self._plasmidfinder_dataframe

        if self._has_pointfinder:
            df = df.append(self._pointfinder_dataframe, sort=True)

        df = self._compile_results(df)

        if include_negatives:
            df = self._include_negatives(df)

        df.rename(columns={'Gene': 'Genotype'}, inplace=True)

        if ds is not None:
            ds = self._compile_plasmids(ds)

            df = df.merge(ds, on='Isolate ID', how='left').fillna(value={'Plasmid Genes': 'None'})

        return df.sort_index()

    def create_detailed_summary(self, include_negatives=True):
        df = self._resfinder_dataframe
        df['Data Type']='Resistance'
        ds = self._plasmidfinder_dataframe

        column_names = ['Gene', 'Predicted Phenotype','%Identity', '%Overlap', 'HSP Length/Total Length','Contig', 'Start', 'End', 'Accession', 'Data Type']

        ds = ds.reindex(columns=column_names)

        if self._has_pointfinder:
            dt = self._pointfinder_dataframe
            dt['Data Type']='Resistance'
            dt = dt.reindex(columns=column_names)
            df = df.append(dt, sort=True)

        if include_negatives:
            df = self._include_detailed_negatives(df, ds)
            df = df.fillna(" ")

        if ds is not None:
            ds['Data Type']='Plasmid'
            ds['Predicted Phenotype']=''
            df = df.append(ds, sort=True)
            df = df.reindex(columns=column_names)
            df = df.sort_values(['Isolate ID', 'Data Type', 'Gene'])

        return df
