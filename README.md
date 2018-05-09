[![Build Status](https://travis-ci.org/phac-nml/staramr.svg?branch=development)](https://travis-ci.org/phac-nml/staramr)
[![pypi](https://badge.fury.io/py/staramr.svg)](https://pypi.python.org/pypi/staramr/)
[![conda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://anaconda.org/bioconda/staramr)

# `staramr`

`staramr` (*AMR) scans bacterial genome contigs against both the [ResFinder][resfinder-db] and [PointFinder][pointfinder-db] databases (used by the [ResFinder webservice][resfinder-web]) and compiles a summary report of detected antimicrobial resistance genes.

**Note: The predicted phenotypes/drug resistances are for microbiological resistance and *not* clinical resistance. This is an experimental feature provided with support from the NARMS/CIPARS Molecular Working Group and is continually being improved. We welcome any feedback or suggestions.**

For example:

```
staramr search -o out --pointfinder-organism salmonella *.fasta
```

**out/summary.tsv**:

| Isolate ID | Genotype                                                  | Predicted Phenotype                                                                                       |
|------------|-----------------------------------------------------------|-----------------------------------------------------------------------------------------------------------|
| SRR1952908 | aadA1, aadA2, blaTEM-57, cmlA1, gyrA (S83Y), sul3, tet(A) | streptomycin, ampicillin, chloramphenicol, ciprofloxacin I/R, nalidixic acid, sulfisoxazole, tetracycline |
| SRR1952926 | blaTEM-57, gyrA (S83Y), tet(A)                            | ampicillin, ciprofloxacin I/R, nalidixic acid, tetracycline                                               |

**out/resfinder.tsv**:

| Isolate ID | Gene      | Predicted Phenotype | %Identity | %Overlap | HSP Length/Total Length | Contig      | Start | End  | Accession |
|------------|-----------|---------------------|-----------|----------|-------------------------|-------------|-------|------|-----------|
| SRR1952908 | sul3      | sulfisoxazole       | 100.00    | 100.00   | 792/792                 | contig00030 | 2091  | 2882 | AJ459418  |
| SRR1952908 | tet(A)    | tetracycline        | 99.92     | 97.80    | 1247/1275               | contig00032 | 1476  | 2722 | AF534183  |
| SRR1952908 | cmlA1     | chloramphenicol     | 99.92     | 100.00   | 1260/1260               | contig00030 | 5448  | 6707 | M64556    |
| SRR1952908 | aadA1     | streptomycin        | 100.00    | 100.00   | 792/792                 | contig00030 | 4564  | 5355 | JQ414041  |
| SRR1952908 | aadA2     | streptomycin        | 99.75     | 100.00   | 792/792                 | contig00030 | 6969  | 7760 | JQ364967  |
| SRR1952908 | blaTEM-57 | ampicillin          | 99.88     | 100.00   | 861/861                 | contig00032 | 5387  | 6247 | FJ405211  |
| SRR1952926 | tet(A)    | tetracycline        | 99.92     | 97.80    | 1247/1275               | contig00027 | 1405  | 2651 | AF534183  |
| SRR1952926 | blaTEM-57 | ampicillin          | 99.88     | 100.00   | 861/861                 | contig00027 | 5316  | 6176 | FJ405211  |

**out/pointfinder.tsv**:

| Isolate ID | Gene        | Predicted Phenotype               | Type  | Position | Mutation            | %Identity | %Overlap | HSP Length/Total Length | Contig      | Start  | End    |
|------------|-------------|-----------------------------------|-------|----------|---------------------|-----------|----------|-------------------------|-------------|--------|--------|
| SRR1952908 | gyrA (S83Y) | ciprofloxacin I/R, nalidixic acid | codon | 83       | TCC -> TAC (S -> Y) | 99.96     | 100.00   | 2637/2637               | contig00008 | 20165  | 22801  |
| SRR1952926 | gyrA (S83Y) | ciprofloxacin I/R, nalidixic acid | codon | 83       | TCC -> TAC (S -> Y) | 99.96     | 100.00   | 2637/2637               | contig00011 | 157768 | 160404 |

# Table of Contents

   * [Quick Usage](#quick-usage)
      * [Search contigs](#search-contigs)
      * [Database Info](#database-info)
      * [Update Database](#update-database)
      * [Restore Database](#restore-database)
   * [Installation](#installation)
      * [Bioconda](#bioconda)
      * [PyPI/Pip](#pypipip)
      * [Latest Code](#latest-code)
      * [Dependencies](#dependencies)
   * [Output](#output)
   * [Usage](#usage)
      * [Main Command](#main-command)
      * [Search](#search)
      * [Database Build](#database-build)
      * [Database Update](#database-update)
      * [Database Info](#database-info-1)
      * [Databae Restore](#databae-restore)
   * [Caveats](#caveats)
   * [Acknowledgements](#acknowledgements)
   * [Citations](#citations)
   * [Legal](#legal)

# Quick Usage

## Search contigs

To search a list of contigs (in **fasta** format) for AMR genes using ResFinder please run:

```bash
staramr search -o out *.fasta
```

Output files will be located in the directory `out/`.

To include acquired point-mutation resistances using PointFinder, please run:

 ```bash
staramr search --pointfinder-organism salmonella -o out *.fasta
```

Where `--pointfinder-organism` is the specific organism you are interested in (currently only *salmonella* is supported).


## Database Info

To print information about the installed databases, please run:

```
staramr db info
```

## Update Database

If you wish to update to the latest ResFinder and PointFinder databases, you may run:

```bash
staramr db update --update-default
```

If you wish to switch to specific git commits of the ResFinder and PointFinder databases you may also pass `--resfinder-commit [COMMIT]` and `--pointfinder-commit [COMMIT]`.

## Restore Database

If you find that, somehow, the database gets messed up, say a message like:

```
subprocess.CalledProcessError: Command '['makeblastdb', '-in', 'resfinder/macrolide.fsa', '-dbtype', 'nucl', '-parse_seqids']' returned non-zero exit status 1
```

Then don't worry, you're in luck. You can restore the default database with:

```
staramr db restore
``` 

# Installation

## Bioconda

The easiest way to install `staramr` is through [Bioconda][bioconda].

```bash
conda install -c bioconda staramr
```

This will install the `staramr` Python package as well as all necessary dependencies and databases.  You can now run:

```bash
staramr --help
```

If you wish to use `staramr` in an isolated environment (in case dependencies conflict) you may alternatively install with:

```bash
conda create -c bioconda --name staramr staramr
```

To run `staramr` in this case, you must first activate the environment.  That is:

```bash
source activate staramr
staramr --help
```

## PyPI/Pip

You can also install `staramr` from [PyPI][pypi-staramr] using `pip`:

```
pip install staramr
```

However, you will have to install the external dependencies (listed below) separately.

## Latest Code

If you wish to make use of the latest in-development version of `staramr`, you may update directly from GitHub using `pip`:

```bash
pip install git+https://github.com/phac-nml/staramr
```

This will only install the Python code, you will still have to install the dependencies listed below (or run the `pip` command from the previously installed Bioconda environment).

Alternatively, if you wish to do development with `staramr` you can use a Python virtual environment (you must still install the non-Python dependencies separately).

```bash
# Clone code
git clone https://github.com/phac-nml/staramr.git
cd staramr

# Setup virtual environment
virtualenv -p /path/to/python-bin .venv
source .venv/bin/activate

# Install staramr. Use '-e' to update the install on code changes.
pip install -e .

# Now run `starmr`
starmr 
```

Due to the way I package the ResFinder/PointFinder databases, the development code will not come with a default database.  You must first build the database before usage. E.g.

```
staramr db build --resfinder-commit dc33e2f9ec2c420f99f77c5c33ae3faa79c999f2 --pointfinder-commit ba65c4d175decdc841a0bef9f9be1c1589c0070a
```

## Dependencies

* Python 3
* BLAST+
* Git

# Output

There are 5 different output files produced by `staramr`:

1. `summary.tsv`:  A summary of all detected AMR genes/mutations in each genome, one genome per line.
2. `resfinder.tsv`: A tabular file of each AMR gene and additional BLAST information from the **ResFinder** database, one gene per line.
3. `pointfinder.tsv`: A tabular file of each AMR point mutation and additional BLAST information from the **PointFinder** database, one gene per line.
4. `settings.txt`: The command-line, database versions, and other settings used to run `staramr`.
5. `results.xlsx`: An Excel spreadsheet containing the previous 4 files as separate worksheets.

In addition, the directory `hits/` stores fasta files of the specific blast hits. 

# Usage

## Main Command

Main `staramr` command. Can be used to set global options (primarily `--verbose`).

```
usage: staramr [-h] [--verbose] [-V] {search,db} ...

Do AMR detection for genes and point mutations

positional arguments:
  {search,db}    Subcommand for AMR detection.
    search       Search for AMR genes
    db           Download ResFinder/PointFinder databases

optional arguments:
  -h, --help     show this help message and exit
  --verbose      Turn on verbose logging [False].
  -V, --version  show program's version number and exit
```

## Search

Searches input FASTA files for AMR genes.

```
usage: staramr search [-h] [-n NPROCS] [--pid-threshold PID_THRESHOLD]
                      [--percent-length-overlap-resfinder PLENGTH_THRESHOLD_RESFINDER]
                      [--percent-length-overlap-pointfinder PLENGTH_THRESHOLD_POINTFINDER]
                      [--pointfinder-organism POINTFINDER_ORGANISM]
                      [--exclude-negatives] [--report-all-blast]
                      [--exclude-resistance-phenotypes] [-d DATABASE]
                      [-o OUTPUT_DIR]
                      files [files ...]

positional arguments:
  files

optional arguments:
  -h, --help            show this help message and exit
  -n NPROCS, --nprocs NPROCS
                        The number of processing cores to use [16].
  --pid-threshold PID_THRESHOLD
                        The percent identity threshold [98.0].
  --percent-length-overlap-resfinder PLENGTH_THRESHOLD_RESFINDER
                        The percent length overlap for resfinder results [60.0].
  --percent-length-overlap-pointfinder PLENGTH_THRESHOLD_POINTFINDER
                        The percent length overlap for pointfinder results [95.0].
  --pointfinder-organism POINTFINDER_ORGANISM
                        The organism to use for pointfinder {salmonella} [None].
  --exclude-negatives   Exclude negative results (those sensitive to antimicrobials) [False].
  --report-all-blast    Report all blast hits (vs. only top blast hits) [False].
  --exclude-resistance-phenotypes
                        Exclude predicted antimicrobial resistances [False].
  -d DATABASE, --database DATABASE
                        The directory containing the resfinder/pointfinder databases [staramr/databases/data].
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        The output directory for results.  If unset prints all results to stdout.

Example:
        staramr search --output-dir out *.fasta
                Searches the files *.fasta for AMR genes using only the ResFinder database,
                storing results in the out/ directory.

        staramr search --pointfinder-organism salmonella --output-dir out *.fasta
                Searches *.fasta for AMR genes using ResFinder and PointFinder database with the passed organism,
                storing results in out/.
```

## Database Build

Downloads and builds the ResFinder and PointFinder databases.

```
usage: staramr db build [-h] [--dir DESTINATION]
                        [--resfinder-commit RESFINDER_COMMIT]
                        [--pointfinder-commit POINTFINDER_COMMIT]

optional arguments:
  -h, --help            show this help message and exit
  --dir DESTINATION     The directory to download the databases into [staramr/databases/data].
  --resfinder-commit RESFINDER_COMMIT
                        The specific git commit for the resfinder database [latest].
  --pointfinder-commit POINTFINDER_COMMIT
                        The specific git commit for the pointfinder database [latest].

Example:
        staramr db build
                Builds a new ResFinder/PointFinder database under staramr/databases if it does not exist

        staramr db build --dir databases
                Builds a new ResFinder/PointFinder database under databases/
```

## Database Update

Updates an existing download of the ResFinder and PointFinder databases.

```
usage: staramr db update [-h] [-d] [--resfinder-commit RESFINDER_COMMIT]
                         [--pointfinder-commit POINTFINDER_COMMIT]
                         ...

positional arguments:
  directories

optional arguments:
  -h, --help            show this help message and exit
  -d, --update-default  Updates default database directory (staramr/databases/data).
  --resfinder-commit RESFINDER_COMMIT
                        The specific git commit for the resfinder database [latest].
  --pointfinder-commit POINTFINDER_COMMIT
                        The specific git commit for the pointfinder database [latest].

Example:
        staramr db update databases/
                Updates the ResFinder/PointFinder database under databases/

        staramr db update -d
                Updates the default ResFinder/PointFinder database under staramr/databases
```

## Database Info

Prints information about an existing build of the ResFinder/PointFinder databases.

```
usage: staramr db info [-h] ...

positional arguments:
  directories

optional arguments:
  -h, --help   show this help message and exit

Example:
        staramr db info
                Prints information about the default database in staramr/databases

        staramr db info databases
                Prints information on the database stored in databases/
```

## Databae Restore

Restores the default database for `staramr`.

```
usage: staramr db restore [-h] [-f]

optional arguments:
  -h, --help   show this help message and exit
  -f, --force  Force restore without asking for confirmation.

Example:
        staramr db restore/
                Restores the default ResFinder/PointFinder database
```

# Caveats

This software is still a work-in-progress.  In particular, not all organisms stored in the PointFinder database are supported (only *salmonella* is currently supported). Additionally, the predicted phenotypes are for microbiological resistance and *not* clinical resistance. Phenotype/drug resistance predictions are an experimental feature which is continually being improved.

`staramr` only works on assembled genomes and not directly on reads. A quick genome assembler you could use is [Shovill][shovill]. Or, you may also wish to try out the [ResFinder webservice][resfinder-web],  or the command-line tools [rgi][] or [ariba][] which will work on sequence reads as well as genome assemblies.  You may also wish to check out the [CARD webservice][card-web]. 

# Acknowledgements

Some ideas for the software were derived from the [ResFinder][resfinder-git] and [PointFinder][pointfinder-git] command-line software, as well as from [ABRicate][abricate].

Phenotype/drug resistance predictions are provided with support from the NARMS/CIPARS Molecular Working Group. 

# Citations

If you find `staramr` useful, please consider citing this GitHub repository (https://github.com/phac-nml/staramr) as well as the original ResFinder and PointFinder publications.

> **Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup FM, Larsen MV**. 2012. Identification of acquired antimicrobial resistance genes. J. Antimicrob. Chemother. 67:2640–2644. doi: [10.1093/jac/dks261][resfinder-cite]

> **Zankari E, Allesøe R, Joensen KG, Cavaco LM, Lund O, Aarestrup F**. PointFinder: a novel web tool for WGS-based detection of antimicrobial resistance associated with chromosomal point mutations in bacterial pathogens. J Antimicrob Chemother. 2017; 72(10): 2764–8. doi: [10.1093/jac/dkx217][pointfinder-cite]

# Legal

Copyright 2018 Government of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

[resfinder-db]: https://bitbucket.org/genomicepidemiology/resfinder_db
[pointfinder-db]: https://bitbucket.org/genomicepidemiology/pointfinder_db
[resfinder-web]: https://cge.cbs.dtu.dk/services/ResFinder/
[resfinder-cite]: https://dx.doi.org/10.1093/jac/dks261
[pointfinder-cite]: https://doi.org/10.1093/jac/dkx217
[Bioconda]: https://bioconda.github.io/
[requirements.txt]: requirements.txt
[resfinder-git]: https://bitbucket.org/genomicepidemiology/resfinder
[pointfinder-git]: https://bitbucket.org/genomicepidemiology/pointfinder-3.0
[abricate]: https://github.com/tseemann/abricate
[shovill]: https://github.com/tseemann/shovill
[ariba]: https://github.com/sanger-pathogens/ariba
[rgi]: https://github.com/arpcard/rgi
[pypi-staramr]: https://pypi.org/project/staramr/
[bioconda]: https://bioconda.github.io/
[card-web]: https://card.mcmaster.ca/
