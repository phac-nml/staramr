# Version 0.9.0

* Updates to PointFinder database handling
    * Adds the ability to handle promoters (regions with both promoter nucleotide information and non-promoter codon information)
    * Adds handling of insertions and deletions in nucleotide and codon sequence
    * Updates list of supported PointFinder species to `salmonella`, `campylobacter`, `enterococcus_faecalis`, `enterococcus_faecium`, `escherichia_coli`, `helicobacter_pylori`.
* Switch name `e.coli` to `escherichia_coli` in PointFinder gene-drug key to match organism name in PointFinder database.

# Version 0.8.0

* Fixed issue when using older version of pandas (#136) (0.8.0.dev0).
* Fixed up some Python warnings related to pandas (0.8.0.dev0).
* Adjusted `mlst` tests to account for differences in results for newer versions (0.8.0.dev0).
* Drop support for Python 3.5 as it leads to issues with managing dependency versions (0.8.0.dev0).
* Switched from disallowing to generating a warning when the PointFinder organism is not one of the validated organisms (0.8.0.dev1).
* Updated ResFinder and PointFinder gene drug key to `072621` (0.8.0.dev2).
* Updated default database commits to those corresponding to dates used by ResFinder (2022-05-24), PointFinder (2021-02-01), and PlasmidFinder (2021-11-29) (0.8.0.dev2).
* Fixed issue when multiple matches for same amino acid change to try and select the most specific amino acid change (0.8.0.dev2).
* Fixed up PlasmidFinder database handling to account for changes to PlasmidFinder database structure (0.8.0.dev2).

# Version 0.7.2.zenodo0

* Identical code to version `0.7.2` but made mainly to upload this version into Zenodo.
* Migrated integration tests from TravisCI to GitHub Actions
* A few small fixes to `README.md`.

# Version 0.7.2

* Fixed `KeyError` issue with later versions of pandas (#115, thanks @javiertognarelli).

# Version 0.7.1

* Fix a bug so that the Sequence column in resfinder.tsv uses the isolate sequence instead of the reference sequence

# Version 0.7.0

* Added quality module that adds PASS/Fail column and detail information in Summary.tsv
* Added following new optional arguments for Search.py
  - --genome-size-lower-bound
  - --genome-size-upper-bound
  - --minimum-N50-value
  - --minimum-contig-length
  - --unacceptable-number-contigs
* Add DNA column in Resfinder report

# Version 0.6.1

* Added --output-mlst in Search.py

# Version 0.6.0

* Added [coloredlogs](https://pypi.org/project/coloredlogs/) library to format the output
* Added support for [MLST](https://github.com/tseemann/mlst)

# Version 0.5.1

* Renamed the following columns for clarification:
    - `Plasmid Genes` to `Plasmid` in Summary table.
    - `Gene` to `Plasmid` in PlasmidFinder table.
    - `Gene` to `Gene/Plasmid` in Detailed Summary table.

# Version 0.5.0

* Add support for scanning against the PlasmidFinder database.
* Upgraded the testing package to use [Green test runner](https://github.com/CleanCut/green).
* Added Detailed_Summary table which combines results from Resfinder, Pointfinder (optional), and Plasmidfinder.
* Added `--ignore-invalid-files` command and check for duplicate sequence ids.

# Version 0.4.0

* Add support for campylobacter from PointFinder database.
* Fix `read_table` deprecation warnings by replacing `read_table` with `read_csv`.
* Handling issue with name of `16S` gene in PointFinder database for salmonella.
* Refactoring and simplifying some of the git ResFinder/PointFinder database code.
* Added automated type checking with [mypy](https://mypy.readthedocs.io).

# Version 0.3.0

* Exclusion of `aac(6')-Iaa` from results by default. Added ability to override this with `--no-exclude-genes` or pass a custom list of genes to exclude from results with `--exclude-genes-file`.

# Version 0.2.2

* Fix issue where `staramr` crashes if an input contig id is a number.

# Version 0.2.1

* Minor
    * Updating default ResFinder/PointFinder databases to version from July 2018.
    * Fix regex extracting gene/variant/accession values from ResFinder/PointFinder databases.
    * Fixing a few entries in table mapping genes to phenotypes.
    * Print stderr for errors with `makeblastdb`

# Version 0.2.0

* Major
    * Inclusion of predicted resistances to antimicrobial drugs thanks to gene/drug mappings from the NARMS/CIPARS Molecular Working Group. Resistance predictions are microbiological resistances and not clinical resistances (issue #4, #6).
    * Adding a `staramr db restore-default` command to restore the default `staramr` database (issue #3).
    * Switched to using BLAST Tabular data + pandas to read BLAST results (issue #10).
    * Inverted direction of BLAST (we now BLAST the AMR gene files against the input genomes).
* Minor
    * Less verbose messages when encountering errors parsing the command-line options.
    * Able to support adding options after a list of files (e.g., `staramr search *.fasta -h` will print help docs).
    * Switched to including negative AMR results (samples with no AMR genes) by default.  Must now use parameter `--exclude-negatives` to exclude them (issue #2).
    * Only print 2 decimals in Excel output (issue #5).
    * Automatically adjust Excel cells to better fit text (issue #7).
    * Many other coding improvements (issue #11, #13 and others).

# Version 0.1.0

* Initial release.  Supports batch scanning against the ResFinder and PointFinder databases.
