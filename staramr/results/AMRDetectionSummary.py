class AMRDetectionSummary:

    def __init__(self, resfinder_dataframe, pointfinder_dataframe=None):
        self._resfinder_dataframe = resfinder_dataframe

        if pointfinder_dataframe is not None:
            self._has_pointfinder = True
            self._pointfinder_dataframe = pointfinder_dataframe
        else:
            self._has_pointfinder = False

    def _compile_results(self, df):
        df_summary = df.groupby(['FILE']).aggregate(lambda x: {'GENE': "%s" % ', '.join(x['GENE']),
                                                               'RESFINDER_PHENOTYPE': "%s" % ', '.join(
                                                                   x['RESFINDER_PHENOTYPE'])})
        return df_summary[['GENE', 'RESFINDER_PHENOTYPE']]

    def create_summary(self):
        if self._has_pointfinder:
            res_point_df = self._resfinder_dataframe.append(self._pointfinder_dataframe)
            return self._compile_results(res_point_df)
        else:
            return self._compile_results(self._resfinder_dataframe)
