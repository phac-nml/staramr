# Version 0.2.0 (in development)

* Major
  * Inclusion of predicted resistances to antimicrobial drugs thanks to gene/drug mappings from the NARMS/CIPARS Molecular Working Group. Resistance predictions are microbiological resistances and not clinical resistances. 
* Minor
  * Switched to including negative AMR results (samples with no AMR genes) by default.  Must now use parameter `--exclude-negatives` to exclude them.
  * Only print 2 decimals in excel output.

# Version 0.1.0

* Initial release.  Supports batch scanning against the ResFinder and PointFinder databases.
