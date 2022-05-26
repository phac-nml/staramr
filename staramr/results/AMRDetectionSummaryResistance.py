from collections import OrderedDict

from staramr.results.AMRDetectionSummary import AMRDetectionSummary

"""
Summarizes both ResFinder and PointFinder database results into a single table.
"""


class AMRDetectionSummaryResistance(AMRDetectionSummary):

    def __init__(self, files, resfinder_dataframe, quality_module_dataframe,pointfinder_dataframe=None, plasmidfinder_dataframe=None, mlst_dataframe=None):
        """
        Creates a new AMRDetectionSummaryResistance.
        :param files: The list of genome files we have scanned against.
        :param resfinder_dataframe: The pd.DataFrame containing the ResFinder results.
        :param quality_module_dataframe: The pd.DataFrame containing the genome size, N50 value, number of contigs under our user defined minimum length
        as well as the results of our quality metrics (pass or fail) and the corresponding feedback
        :param pointfinder_dataframe: The pd.DataFrame containing the PointFinder results.
        :param plasmidfinder_dataframe: The pd.DataFrame containing the PlasmidFinder results.
        """
        super().__init__(files, resfinder_dataframe,quality_module_dataframe, pointfinder_dataframe, plasmidfinder_dataframe, mlst_dataframe)

    def _aggregate_phenotype(self, phenotype_series):
        flattened_phenotype_list = [y.strip() for x in list(phenotype_series) for y in
                                    x.split(self.SEPARATOR)]

        # Only remove None if there is more than one entry in this list
        if len(flattened_phenotype_list) > 1 and 'None' in flattened_phenotype_list:
            flattened_phenotype_list.remove('None')

        uniq_phenotype = OrderedDict.fromkeys(flattened_phenotype_list)

        return (self.SEPARATOR + ' ').join(list(uniq_phenotype))

    def _aggregate_gene(self, gene_series):
        return (self.SEPARATOR + ' ').join(list(gene_series))

    def _compile_results(self, df):
        df_summary = df.copy()

        # Used to sort by gene names, ignoring case
        df_summary['Gene.Lower'] = df['Gene'].str.lower()

        # Compiles the gene/phenotype results into a single entry per isolate (groupby)
        df_summary = df_summary \
            .sort_values(by=['Gene.Lower']) \
            .groupby(['Isolate ID'], sort=True) \
            .aggregate({'Gene': self._aggregate_gene,
                'Predicted Phenotype': self._aggregate_phenotype})
        return df_summary[['Gene', 'Predicted Phenotype']]

    def _get_detailed_negative_columns(self):
        return ['Isolate ID', 'Gene', 'Predicted Phenotype', 'Start', 'End']

    def _get_summary_empty_values(self):
        return {'Genotype': 'None', 'Predicted Phenotype': 'Sensitive'}

    def _get_summary_resistance_columns(self):
        return ['Genotype', 'Predicted Phenotype', 'Plasmid']

    def _get_detailed_summary_columns(self):
        return ['Gene', 'Data Type', 'Predicted Phenotype', '%Identity', '%Overlap', 'HSP Length/Total Length', 'Contig', 'Start',
                'End', 'Accession']

    def _include_phenotype(self):
        return True
