# Version 0.2.0 (in development)

* Major
    * Inclusion of predicted resistances to antimicrobial drugs thanks to gene/drug mappings from the NARMS/CIPARS Molecular Working Group. Resistance predictions are microbiological resistances and not clinical resistances (issue #4, #6).
    * Adding a `staramr db restore` command to restore the default `staramr` database in case an error was encountered (issue #3).
        * Also sets a database to an **error** state in case there was an error downloading/formatting files.
* Minor
    * Less verbose messages when encountering errors parsing the command-line options.
    * Able to support adding options after a list of files (e.g., `staramr search *.fasta -h` will print help docs).
    * Switched to including negative AMR results (samples with no AMR genes) by default.  Must now use parameter `--exclude-negatives` to exclude them (issue #2).
    * Only print 2 decimals in Excel output (issue #5).
    * Automatically adjust Excel cells to better fit text (issue #7).
    * Many other coding improvements (issue #11 and others).

# Version 0.1.0

* Initial release.  Supports batch scanning against the ResFinder and PointFinder databases.
