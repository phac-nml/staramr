import re
import logging

from amr.blast.results.AMRHitHSP import AMRHitHSP

logger = logging.getLogger('ResfinderHitHSP')

class ResfinderHitHSP(AMRHitHSP):

    def __init__(self, file, blast_record, hit, hsp):
        super().__init__(file, blast_record, hit, hsp)

        re_search = re.search("([^_]+)_([^_]+)_([^_\s]+)$", hit.hit_id)
        if not re_search:
            raise Exception("Could not split up seq name for [" + hit.hit_id + "]")
        self._gene = re_search.group(1)
        self._gene_variant = re_search.group(2)
        self._accession = re_search.group(3)

        logger.debug("hit_id="+str(hit.hit_id))

    def get_gene(self):
        return self._gene

    def get_accession(self):
        return self._accession

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            return self.__dict__ == other.__dict__
