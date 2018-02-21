import pandas

class AMRDetectionSummary:

    def __init__(self, resfinder_dataframe, pointfinder_dataframe=None):
        self._resfinder_dataframe = resfinder_dataframe
        self._pointfinder_dataframe = pointfinder_dataframe

        if pointfinder_dataframe is not None:
            self._has_pointfinder = True
        else:
            self._has_pointfinder = False

    def _compile_resfinder_results(self, df):
        df_summary = df.groupby(['FILE']).aggregate(lambda x: {'GENE': "%s" % ', '.join(x['GENE']),
                                                                       'RESFINDER_PHENOTYPE': "%s" % ', '.join(
                                                                           x['RESFINDER_PHENOTYPE'])})
        df_summary = df_summary.iloc[:, 0:2]
        df_summary['FILE'] = df_summary.index

        return df_summary[['FILE', 'GENE', 'RESFINDER_PHENOTYPE']]


    def _compile_pointfinder_results(self, df):
        df_summary = df.groupby(['FILE']).aggregate(lambda x: {'GENE': "%s" % ', '.join(x['GENE']),
                                                                       'POINTFINDER_PHENOTYPE': "%s" % ', '.join(
                                                                           x['POINTFINDER_PHENOTYPE'])})
        df_summary = df_summary.iloc[:, 0:2]
        df_summary['FILE'] = df_summary.index

        return df_summary[['FILE', 'GENE', 'POINTFINDER_PHENOTYPE']]


    def create_summary(self):
        res_df = self._compile_resfinder_results(self._resfinder_dataframe)

        if self._has_pointfinder:
            point_df = self._compile_pointfinder_results(self._pointfinder_dataframe)
            return pandas.concat([res_df, point_df])
        else:
            return res_df


