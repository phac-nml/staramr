import numpy
import pandas
from collections import OrderedDict

from staramr.results.AMRDetectionSummary import AMRDetectionSummary

"""
Summarizes both ResFinder and PointFinder database results into a single table.
"""


class AMRDetectionSummaryResistance(AMRDetectionSummary):
    blank = '-'

    def __init__(self, files, resfinder_dataframe, pointfinder_dataframe=None):
        super().__init__(files, resfinder_dataframe, pointfinder_dataframe)

    def _aggregate_gene_phenotype(self, dataframe):
        uniq_phenotype = OrderedDict.fromkeys(dataframe['Predicted Phenotype'])

        return {'Gene': "%s" % ', '.join(dataframe['Gene']),
                'Predicted Phenotype': "%s" % ', '.join(list(uniq_phenotype))
                }

    def _compile_results(self, df):
        df_summary = df.replace(numpy.nan, self.blank)

        # Used to sort by gene names, ignoring case
        df_summary['Gene.Lower'] = df['Gene'].str.lower()

        df_summary = df_summary.sort_values(by=['Gene.Lower']).groupby(['Isolate ID']).aggregate(
            self._aggregate_gene_phenotype).replace({'Predicted Phenotype': {self.blank: 'Sensitive'}}).replace(
            {'Predicted Phenotype': {', ' + self.blank: '', self.blank + ', ': ''}}, regex=True)
        return df_summary[['Gene', 'Predicted Phenotype']]

    def _include_negatives(self, df):
        result_names_set = set(df.index.tolist())
        names_set = set(self._names)

        negative_names_set = names_set - result_names_set
        negative_entries = pandas.DataFrame([[x, 'None', 'Sensitive'] for x in negative_names_set],
                                            columns=('Isolate ID', 'Gene', 'Predicted Phenotype')).set_index(
            'Isolate ID')
        return df.append(negative_entries)
