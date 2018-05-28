# Version 0.2.0 (in development)

* Major
    * Inclusion of predicted resistances to antimicrobial drugs thanks to gene/drug mappings from the NARMS/CIPARS Molecular Working Group. Resistance predictions are microbiological resistances and not clinical resistances (issue #4).
    * Adding a `staramr db restore-default` command to restore the default `staramr` database (issue #3).
    * Switched to using BLAST Tabular data + pandas to read BLAST results (issue #10).
    * Inverted direction of BLAST (we now BLAST the AMR gene files against the input genomes).
* Minor
    * Less verbose messages when encountering errors parsing the command-line options.
    * Able to support adding options after a list of files (e.g., `staramr search *.fasta -h` will print help docs).
    * Many other coding improvements (issue #11, #13 and others).

# Version 0.1.0

* Initial release.  Supports batch scanning against the ResFinder and PointFinder databases.
