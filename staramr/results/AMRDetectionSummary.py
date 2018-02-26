import pandas
from os import path


class AMRDetectionSummary:

    def __init__(self, files, resfinder_dataframe, pointfinder_dataframe=None):
        self._names = [path.basename(x) for x in files]
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

    def _include_negatives(self, df):
        result_names_set = set(df.index.tolist())
        names_set = set(self._names)

        negative_names_set = names_set - result_names_set
        negative_entries = pandas.DataFrame([[x, 'None', 'Sensitive'] for x in negative_names_set],
                                            columns=('FILE', 'GENE', 'RESFINDER_PHENOTYPE')).set_index('FILE')
        return df.append(negative_entries)

    def create_summary(self, include_negatives=False):
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
