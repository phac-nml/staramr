import re

from amr.blast.results.AMRHitHSP import AMRHitHSP


class ResfinderHitHSP(AMRHitHSP):

    def __init__(self, file, hit, hsp):
        super().__init__(file, hit, hsp)

        re_search = re.search("(.+)_([^_]+)$", hit.hit_id)
        if not re_search:
            raise Exception("Could not split up seq name for [" + hit.hit_id + "]")
        self.gene = re_search.group(1)
        self.accession = re_search.group(2)

    def get_gene(self):
        return self.gene

    def get_accession(self):
        return self.accession

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            return self.__dict__ == other.__dict__
