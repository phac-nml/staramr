from collections import OrderedDict

import numpy
import pandas as pd

from staramr.results.AMRDetectionSummary import AMRDetectionSummary

"""
Summarizes both ResFinder and PointFinder database results into a single table.
"""


class AMRDetectionSummaryResistance(AMRDetectionSummary):
    BLANK = '-'

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

        return {'Gene': "%s" % (self.SEPARATOR + ' ').join(dataframe['Gene']),
                'Predicted Phenotype': "%s" % (self.SEPARATOR + ' ').join(list(uniq_phenotype))
                }

    def _compile_results(self, df):
        df_summary = df.replace(numpy.nan, self.BLANK)

        # Used to sort by gene names, ignoring case
        df_summary['Gene.Lower'] = df['Gene'].str.lower()

        # Compiles the gene/phenotype results into a single entry per isolate (groupby)
        # Removes any any blank character ('-') results in the predicted phenotype unless there is no other phenotypes
        # in which case, replace with 'Sensitive'
        df_summary = df_summary \
            .sort_values(by=['Gene.Lower']) \
            .groupby(['Isolate ID']) \
            .aggregate(self._aggregate_gene_phenotype) \
            .replace({'Predicted Phenotype': {self.BLANK: 'Sensitive'}}) \
            .replace({'Predicted Phenotype': {(self.SEPARATOR + ' ') + self.BLANK: '',
                                              self.BLANK + (self.SEPARATOR + ' '): ''}},
                     regex=True)
        return df_summary[['Gene', 'Predicted Phenotype']]

    def _include_negatives(self, df):
        result_names_set = set(df.index.tolist())
        names_set = set(self._names)

        negative_names_set = names_set - result_names_set
        negative_entries = pd.DataFrame([[x, 'None', 'Sensitive'] for x in negative_names_set],
                                        columns=('Isolate ID', 'Gene', 'Predicted Phenotype')).set_index(
            'Isolate ID')
        return df.append(negative_entries)
