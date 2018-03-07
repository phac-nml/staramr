from os import path

import pandas

"""
Summarizes both ResFinder and PointFinder database results into a single table.
"""


class AMRDetectionSummary:

    def __init__(self, files, phenotype_column_name, resfinder_dataframe, pointfinder_dataframe=None):
        """
        Constructs an object for summarizing AMR detection results.
        :param files: The list of genome files we have scanned against.
        :param phenotype_column_name: The name of the column in the generated tables containing the phenotype information.
        :param resfinder_dataframe: The pandas.DataFrame containing the ResFinder results.
        :param pointfinder_dataframe: The pandas.DataFrame containing the PointFinder results.
        """
        self._names = [path.basename(x) for x in files]
        self._phenotype_column_name = phenotype_column_name
        self._resfinder_dataframe = resfinder_dataframe

        if pointfinder_dataframe is not None:
            self._has_pointfinder = True
            self._pointfinder_dataframe = pointfinder_dataframe
        else:
            self._has_pointfinder = False

    def _compile_results(self, df):
        df_summary = df.groupby(['FILE']).aggregate(lambda x: {'GENE': "%s" % ', '.join(x['GENE']),
                                                               self._phenotype_column_name: "%s" % ', '.join(
                                                                   x[self._phenotype_column_name])})
        return df_summary[['GENE', self._phenotype_column_name]]

    def _include_negatives(self, df):
        result_names_set = set(df.index.tolist())
        names_set = set(self._names)

        negative_names_set = names_set - result_names_set
        negative_entries = pandas.DataFrame([[x, 'None', 'Sensitive'] for x in negative_names_set],
                                            columns=('FILE', 'GENE', self._phenotype_column_name)).set_index('FILE')
        return df.append(negative_entries)

    def create_summary(self, include_negatives=False):
        """
        Constructs a summary pandas.DataFrame for all ResFinder/PointFinder results.
        :param include_negatives: If True, include files with no ResFinder/PointFinder results.
        :return: A pandas.DataFrame summarizing the results.
        """
        df = None

        if self._has_pointfinder:
            if include_negatives:
                df = self._include_negatives(
                    self._compile_results(self._resfinder_dataframe.append(self._pointfinder_dataframe)))
            else:
                df = self._compile_results(self._resfinder_dataframe.append(self._pointfinder_dataframe))
        else:
            if include_negatives:
                df = self._include_negatives(self._compile_results(self._resfinder_dataframe))
            else:
                df = self._compile_results(self._resfinder_dataframe)

        return df.sort_index()
