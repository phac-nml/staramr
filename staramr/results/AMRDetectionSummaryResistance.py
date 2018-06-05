from collections import OrderedDict

import pandas as pd

from staramr.results.AMRDetectionSummary import AMRDetectionSummary

"""
Summarizes both ResFinder and PointFinder database results into a single table.
"""


class AMRDetectionSummaryResistance(AMRDetectionSummary):

    def __init__(self, files, resfinder_dataframe, pointfinder_dataframe=None):
        """
        Creates a new AMRDetectionSummaryResistance.
        :param files: The list of genome files we have scanned against.
        :param resfinder_dataframe: The pd.DataFrame containing the ResFinder results.
        :param pointfinder_dataframe: The pd.DataFrame containing the PointFinder results.
        """
        super().__init__(files, resfinder_dataframe, pointfinder_dataframe)

    def _aggregate_gene_phenotype(self, dataframe):
        flattened_phenotype_list = [y.strip() for x in dataframe['Predicted Phenotype'].tolist() for y in
                                    x.split(self.SEPARATOR)]
        uniq_phenotype = OrderedDict.fromkeys(flattened_phenotype_list)

        return {'Gene': (self.SEPARATOR + ' ').join(dataframe['Gene']),
                'Predicted Phenotype': (self.SEPARATOR + ' ').join(list(uniq_phenotype))
                }

    def _compile_results(self, df):
        df_summary = df.copy()

        # Used to sort by gene names, ignoring case
        df_summary['Gene.Lower'] = df['Gene'].str.lower()

        # Compiles the gene/phenotype results into a single entry per isolate (groupby)
        df_summary = df_summary \
            .sort_values(by=['Gene.Lower']) \
            .groupby(['Isolate ID'], sort=True) \
            .aggregate(self._aggregate_gene_phenotype)
        return df_summary[['Gene', 'Predicted Phenotype']]

    def _include_negatives(self, df):
        result_names_set = set(df.index.tolist())
        names_set = set(self._names)

        negative_names_set = names_set - result_names_set
        negative_entries = pd.DataFrame([[x, 'None', 'Sensitive'] for x in negative_names_set],
                                        columns=('Isolate ID', 'Gene', 'Predicted Phenotype')).set_index(
            'Isolate ID')
        return df.append(negative_entries, sort=True)
