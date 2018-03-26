[![Build Status](https://travis-ci.org/phac-nml/staramr.svg?branch=development)](https://travis-ci.org/phac-nml/staramr)

# `staramr`

`staramr` (*AMR) scans bacterial genome contigs against both the [ResFinder][resfinder-db] and [PointFinder][pointfinder-db] databases and complies a summary report of detected antimicrobial resistance genes.

For example:

```
staramr.py --output out --pointfinder-organism salmonella *.fasta
```

**out/summary.tsv**:

| Isolate ID | Genotype                           |
|------------|------------------------------------|
| file1      | gyrA (D87N)                        |
| file2      | fosA7                              |
| file3      | sul2, tet(D), aph(6)-Id, blaTEM-1B |
| file4      | fosA7                              |

# Quick Usage

## Build Database

To construct the ResFinder and PointFinder database please run:

```bash
staramr.py db build
```

## Search contigs

To search a list of contigs (in **fasta** format) for AMR genes using ResFinder please run:

```bash
staramr.py search -o out *.fasta
```

Output files will be located in the directory `out/`.

To include acquired point-mutation resistances using PointFinder, please run:

 ```bash
staramr.py search --pointfinder-organism salmonella -o out *.fasta
```

Where `--pointfinder-organism` is the specific organism you are interested in (currently only *salmonella* is supported).

# Installation

To get the latest code, please run:

```
git clone https://github.com/phac-nml/staramr
```

`staramr` requires the dependencies listed below.  The easiest way to install these is through a [Bioconda][] environment (a full conda package will be made when code is stable).  Assuming you have `conda` installed, you may run:

```bash
conda create --name staramr --file conda-packages.txt

# Activate environment
source activate staramr
```

Now, you may run `staramr`:

```
./staramr.py
```

## Dependencies

* Python 3
* BLAST+
* Git
* Python packages in [requirements.txt][]

# Tests

To run the test suite, please run:

```
./run-tests.sh
```

# Output

There are 5 different output files produced by `staramr`:

1. `summary.tsv`:  A summary of all detected AMR genes/mutations in each genome, one genome per line.
2. `resfinder.tsv`: A tabular file of each AMR gene and additional BLAST information from the **ResFinder** database, one gene per line.
3. `pointfinder.tsv`: A tabular file of each AMR point mutation and additional BLAST information from the **PointFinder** database, one gene per line.
4. `settings.txt`: The command-line, database versions, and other settings used to run `staramr`.
5. `results.xlsx`: An Excel spreadsheet containing the previous 4 files as separate worksheets.

# Usage

## Search

Searches input FASTA files for AMR genes.

```
usage: staramr.py search [-h] [-n NPROCS] [--pid-threshold PID_THRESHOLD]
                         [--percent-length-overlap PLENGTH_THRESHOLD]
                         [--pointfinder-organism POINTFINDER_ORGANISM]
                         [--include-negatives] [-d DATABASE] [-o OUTPUT_DIR]
                         ...

positional arguments:
  files

optional arguments:
  -h, --help            show this help message and exit
  -n NPROCS, --nprocs NPROCS
                        The number of processing cores to use [16].
  --pid-threshold PID_THRESHOLD
                        The percent identity threshold [98.0].
  --percent-length-overlap PLENGTH_THRESHOLD
                        The percent length overlap [60.0].
  --pointfinder-organism POINTFINDER_ORGANISM
                        The organism to use for pointfinder {salmonella} [None].
  --include-negatives   Inclue negative results (those sensitive to antimicrobials) [False].
  -d DATABASE, --database DATABASE
                        The directory containing the resfinder/pointfinder databases [staramr/databases].
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        The output directory for results.  If unset prints all results to stdout.

Example:
        staramr.py search --output-dir out *.fasta
                Searches the files *.fasta for AMR genes using only the ResFinder database,
                storing results in the out/ directory.

        staramr.py search --pointfinder-organism salmonella --output-dir out *.fasta
                Searches *.fasta for AMR genes using ResFinder and PointFinder database with the passed organism,
                storing results in out/.
```

## Database Build

Downloads and builds the ResFinder and PointFinder databases.

```
usage: staramr.py db build [-h] [--dir DESTINATION]

optional arguments:
  -h, --help         show this help message and exit
  --dir DESTINATION  The directory to download the databases into [staramr/databases].

Example:
        staramr.py db build
                Builds a new ResFinder/PointFinder database under staramr/databases if it does not exist

        staramr.py db build --dir databases
                Builds a new ResFinder/PointFinder database under databases/
```

## Database Update

Updates an existing download of the ResFinder and PointFinder databases.

```
usage: staramr.py db update [-h] [-d] ...

positional arguments:
  directories

optional arguments:
  -h, --help            show this help message and exit
  -d, --update-default  Updates default database directory (staramr/databases).

Example:
        staramr.py db update databases/
                Updates the ResFinder/PointFinder database under databases/

        staramr.py db update -d
                Updates the default ResFinder/PointFinder database under staramr/databases
```

## Database Info

Prints information about an existing build of the ResFinder/PointFinder databases.

```
usage: staramr.py db info [-h] ...

positional arguments:
  directories

optional arguments:
  -h, --help   show this help message and exit

Example:
        staramr.py db info
                Prints information about the default database in staramr/databases

        staramr.py db info databases
                Prints information on the database stored in databases/
```

# Caveats

This software is still a work-in-progress.  In particular, not all point mutations stored in the PointFinder database are supported.

# Acknowledgements

Some ideas for the software were derived from the [ResFinder][resfinder-git] and [PointFinder][pointfinder-git] command-line software, as well as from [ABRicate][abricate].

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
[resfinder-cite]: https://dx.doi.org/10.1093/jac/dks261
[pointfinder-cite]: https://doi.org/10.1093/jac/dkx217
[Bioconda]: https://bioconda.github.io/
[requirements.txt]: requirements.txt
[resfinder-git]: https://bitbucket.org/genomicepidemiology/resfinder
[pointfinder-git]: https://bitbucket.org/genomicepidemiology/pointfinder-3.0
[abricate]: https://github.com/tseemann/abricate
