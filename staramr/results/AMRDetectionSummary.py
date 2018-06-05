from os import path

import pandas as pd

"""
Summarizes both ResFinder and PointFinder database results into a single table.
"""


class AMRDetectionSummary:
    SEPARATOR = ','

    def __init__(self, files, resfinder_dataframe, pointfinder_dataframe=None):
        """
        Constructs an object for summarizing AMR detection results.
        :param files: The list of genome files we have scanned against.
        :param resfinder_dataframe: The pd.DataFrame containing the ResFinder results.
        :param pointfinder_dataframe: The pd.DataFrame containing the PointFinder results.
        """
        self._names = [path.splitext(path.basename(x))[0] for x in files]
        self._resfinder_dataframe = resfinder_dataframe

        if pointfinder_dataframe is not None:
            self._has_pointfinder = True
            self._pointfinder_dataframe = pointfinder_dataframe
        else:
            self._has_pointfinder = False

    def _compile_results(self, df):
        df_summary = df.sort_values(by=['Gene']).groupby(['Isolate ID']).aggregate(
            lambda x: {'Gene': (self.SEPARATOR + ' ').join(x['Gene'])})
        return df_summary[['Gene']]

    def _include_negatives(self, df):
        result_names_set = set(df.index.tolist())
        names_set = set(self._names)

        negative_names_set = names_set - result_names_set
        negative_entries = pd.DataFrame([[x, 'None'] for x in negative_names_set],
                                        columns=('Isolate ID', 'Gene')).set_index('Isolate ID')
        return df.append(negative_entries, sort=True)

    def create_summary(self, include_negatives=False):
        """
        Constructs a summary pd.DataFrame for all ResFinder/PointFinder results.
        :param include_negatives: If True, include files with no ResFinder/PointFinder results.
        :return: A pd.DataFrame summarizing the results.
        """
        df = self._resfinder_dataframe

        if self._has_pointfinder:
            df = df.append(self._pointfinder_dataframe, sort=True)

        df = self._compile_results(df)

        if include_negatives:
            df = self._include_negatives(df)

        df.rename(columns={'Gene': 'Genotype'}, inplace=True)

        return df.sort_index()
