class AMRDetectionSummary:

    def __init__(self, resfinder_dataframe, pointfinder_dataframe=None):
        self._resfinder_dataframe = resfinder_dataframe
        self._pointfinder_dataframe = pointfinder_dataframe

        if pointfinder_dataframe:
            self._has_pointfinder = True
        else:
            self._has_pointfinder = False

    def create_summary(self):
        res_df = self._resfinder_dataframe
        res_df_summary = res_df.groupby(['FILE']).aggregate(lambda x: {'GENE': "%s" % ', '.join(x['GENE']),
                                                                       'RESFINDER_PHENOTYPE': "%s" % ', '.join(
                                                                           x['RESFINDER_PHENOTYPE'])})
        res_df_summary = res_df_summary.iloc[:, 0:2]
        res_df_summary['FILE'] = res_df_summary.index

        return res_df_summary[['FILE', 'GENE', 'RESFINDER_PHENOTYPE']]
