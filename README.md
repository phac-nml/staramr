[![Build Status](https://travis-ci.org/phac-nml/staramr.svg?branch=development)](https://travis-ci.org/phac-nml/staramr)
[![pypi](https://badge.fury.io/py/staramr.svg)](https://pypi.python.org/pypi/staramr/)
[![conda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://anaconda.org/bioconda/staramr)

# `staramr`

`staramr` (*AMR) scans bacterial genome contigs against both the [ResFinder][resfinder-db] and [PointFinder][pointfinder-db] databases (used by the [ResFinder webservice][resfinder-web]) and compiles a summary report of detected antimicrobial resistance genes.

**Note: The predicted phenotypes/drug resistances are for microbiological resistance and *not* clinical resistance. This is provided with support from the NARMS/CIPARS Molecular Working Group and is continually being improved. A small comparison between phenotype/drug resistance predictions produced by `staramr` and those available from NCBI can be found in the [tutorial][tutorial]. We welcome any feedback or suggestions.**

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

| Isolate ID | Gene       | Predicted Phenotype  | %Identity  | %Overlap  | HSP Length/Total Length  | Contig       | Start  | End   | Accession |
|------------|------------|----------------------|------------|-----------|--------------------------|--------------|--------|-------|-----------|
| SRR1952908 | sul3       | sulfisoxazole        | 100.00     | 100.00    | 792/792                  | contig00030  | 2091   | 2882  | AJ459418  |
| SRR1952908 | tet(A)     | tetracycline         | 99.92      | 100.00    | 1200/1200                | contig00032  | 1551   | 2750  | AJ517790  |
| SRR1952908 | cmlA1      | chloramphenicol      | 99.92      | 100.00    | 1260/1260                | contig00030  | 6707   | 5448  | M64556    |
| SRR1952908 | aadA1      | streptomycin         | 100.00     | 100.00    | 792/792                  | contig00030  | 5355   | 4564  | JQ414041  |
| SRR1952908 | aadA2      | streptomycin         | 99.75      | 100.00    | 792/792                  | contig00030  | 7760   | 6969  | JQ364967  |
| SRR1952908 | blaTEM-57  | ampicillin           | 99.88      | 100.00    | 861/861                  | contig00032  | 6247   | 5387  | FJ405211  |
| SRR1952926 | tet(A)     | tetracycline         | 99.92      | 100.00    | 1200/1200                | contig00027  | 1480   | 2679  | AJ517790  |
| SRR1952926 | blaTEM-57  | ampicillin           | 99.88      | 100.00    | 861/861                  | contig00027  | 6176   | 5316  | FJ405211  |

**out/pointfinder.tsv**:

| Isolate ID  | Gene         | Predicted Phenotype                | Type   | Position  | Mutation             | %Identity  | %Overlap  | HSP Length/Total Length  | Contig       | Start   | End    |
|-------------|--------------|------------------------------------|--------|-----------|----------------------|------------|-----------|--------------------------|--------------|---------|--------|
| SRR1952908  | gyrA (S83Y)  | ciprofloxacin I/R, nalidixic acid  | codon  | 83        | TCC -> TAC (S -> Y)  | 99.96      | 100.00    | 2637/2637                | contig00008  | 22801   | 20165  |
| SRR1952926  | gyrA (S83Y)  | ciprofloxacin I/R, nalidixic acid  | codon  | 83        | TCC -> TAC (S -> Y)  | 99.96      | 100.00    | 2637/2637                | contig00011  | 157768  | 160404 |

# Table of Contents

- [Quick Usage](#quick-usage)
  * [Search contigs](#search-contigs)
  * [Database Info](#database-info)
  * [Update Database](#update-database)
  * [Restore Database](#restore-database)
- [Installation](#installation)
  * [Bioconda](#bioconda)
  * [PyPI/Pip](#pypipip)
  * [Latest Code](#latest-code)
  * [Dependencies](#dependencies)
- [Output](#output)
  * [summary.tsv](#summarytsv)
  * [resfinder.tsv](#resfindertsv)
  * [pointfinder.tsv](#pointfindertsv)
  * [settings.txt](#settingstxt)
  * [hits/](#hits)
- [Tutorial](#tutorial)
- [Usage](#usage)
  * [Main Command](#main-command)
  * [Search](#search)
  * [Database Build](#database-build)
  * [Database Update](#database-update)
  * [Database Info](#database-info-1)
  * [Database Restore Default](#database-restore-default)
- [Caveats](#caveats)
- [Acknowledgements](#acknowledgements)
- [Citations](#citations)
- [Legal](#legal)

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

If you have updated the ResFinder/PointFinder databases and wish to restore to the default version, you may run:

```
staramr db restore-default
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

# Now run `staramr`
staramr 
```

Due to the way I package the ResFinder/PointFinder databases, the development code will not come with a default database.  You must first build the database before usage. E.g.

```
staramr db restore-default
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

## summary.tsv

The **summary.tsv** output file generated by `staramr` contains the following columns:

* __Isolate ID__: The id of the isolate/genome file(s) passed to `staramr`.
* __Genotype__: The AMR genotype of the isolate.
* __Predicted Phenotype__: The predicted AMR phenotype (drug resistances) for the isolate.

### Example

| Isolate ID | Genotype                                                  | Predicted Phenotype                                                                                       |
|------------|-----------------------------------------------------------|-----------------------------------------------------------------------------------------------------------|
| SRR1952908 | aadA1, aadA2, blaTEM-57, cmlA1, gyrA (S83Y), sul3, tet(A) | streptomycin, ampicillin, chloramphenicol, ciprofloxacin I/R, nalidixic acid, sulfisoxazole, tetracycline |
| SRR1952926 | blaTEM-57, gyrA (S83Y), tet(A)                            | ampicillin, ciprofloxacin I/R, nalidixic acid, tetracycline                                               |

## resfinder.tsv

The **resfinder.tsv** output file generated by `staramr` contains the following columns:

* __Isolate ID__: The id of the isolate/genome file(s) passed to `staramr`.
* __Gene__: The particular AMR gene detected.
* __Predicted Phenotype__: The predicted AMR phenotype (drug resistances) for this gene.
* __%Identity__: The % identity of the top BLAST HSP to the AMR gene.
* __%Overlap__: THe % overlap of the top BLAST HSP to the AMR gene (calculated as __hsp length/total length * 100__).
* __HSP Length/Total Length__ The top BLAST HSP length over the AMR gene total length (nucleotides).
* __Contig__: The contig id containing this AMR gene.
* __Start__: The start of the AMR gene (will be greater than __End__ if on minus strand).
* __End__: The end of the AMR gene.
* __Accession__: The accession of the AMR gene in the ResFinder database.

### Example

| Isolate ID | Gene       | Predicted Phenotype  | %Identity  | %Overlap  | HSP Length/Total Length  | Contig       | Start  | End   | Accession |
|------------|------------|----------------------|------------|-----------|--------------------------|--------------|--------|-------|-----------|
| SRR1952908 | sul3       | sulfisoxazole        | 100.00     | 100.00    | 792/792                  | contig00030  | 2091   | 2882  | AJ459418  |
| SRR1952908 | tet(A)     | tetracycline         | 99.92      | 100.00    | 1200/1200                | contig00032  | 1551   | 2750  | AJ517790  |

## pointfinder.tsv

The **pointfinder.tsv** output file generated by `staramr` contains the following columns:

* __Isolate ID__: The id of the isolate/genome file(s) passed to `staramr`.
* __Gene__: The particular AMR gene detected, with the point mutation within *()*.
* __Predicted Phenotype__: The predicted AMR phenotype (drug resistances) for this gene.
* __Type__: The type of this mutation from PointFinder (either **codon** or **nucleotide**).
* __Position__: The position of the mutation. For **codon** type, the position is the codon number in the gene, for **nucleotide** type it is the nucleotide number.
* __Mutation__: The particular mutation. For **codon** type lists the codon mutation, for **nucleotide** type lists the single nucleotide mutation.
* __%Identity__: The % identity of the top BLAST HSP to the AMR gene.
* __%Overlap__: The % overlap of the top BLAST HSP to the AMR gene (calculated as __hsp length/total length * 100__).
* __HSP Length/Total Length__ The top BLAST HSP length over the AMR gene total length (nucleotides).
* __Contig__: The contig id containing this AMR gene.
* __Start__: The start of the AMR gene (will be greater than __End__ if on minus strand).
* __End__: The end of the AMR gene.

### Example

| Isolate ID  | Gene         | Predicted Phenotype                | Type   | Position  | Mutation             | %Identity  | %Overlap  | HSP Length/Total Length  | Contig       | Start   | End    |
|-------------|--------------|------------------------------------|--------|-----------|----------------------|------------|-----------|--------------------------|--------------|---------|--------|
| SRR1952908  | gyrA (S83Y)  | ciprofloxacin I/R, nalidixic acid  | codon  | 83        | TCC -> TAC (S -> Y)  | 99.96      | 100.00    | 2637/2637                | contig00008  | 22801   | 20165  |
| SRR1952926  | gyrA (S83Y)  | ciprofloxacin I/R, nalidixic acid  | codon  | 83        | TCC -> TAC (S -> Y)  | 99.96      | 100.00    | 2637/2637                | contig00011  | 157768  | 160404 |

## settings.txt

The **settings.txt** file contains the particular settings used to run `staramr`.

* __command_line__: The command line used to run `staramr`.
* __version__: The version of `staramr`.
* __start_time__,__end_time__,__total_minutes__: The start, end, and duration for running `staramr`.
* __resfinder_db_dir__, __pointfinder_db_dir__: The directory containing the ResFinder and PointFinder databases.
* __resfinder_db_url__, __pointfinder_db_url__: The URL to the git repository for the ResFinder and PointFinder databases.
* __resfinder_db_commit__, __pointfinder_db_commit__: The git commit ids for the ResFinder and PointFinder databases.
* __resfinder_db_date__, __pointfinder_db_date__: The date of the git commits of the ResFinder and PointFinder databases.
* __pointfinder_gene_drug_version__, __resfinder_gene_drug_version__: A version identifier for the gene/drug mapping table used by `staramr`.

### Example

```
command_line                  = staramr search -o out --pointfinder-organism salmonella SRR1952908.fasta SRR1952926.fasta
version                       = 0.2.0
start_time                    = 2018-06-08 10:28:47
end_time                      = 2018-06-08 10:28:59
total_minutes                 = 0.20
resfinder_db_dir              = staramr/databases/data/dist/resfinder
resfinder_db_url              = https://bitbucket.org/genomicepidemiology/resfinder_db.git
resfinder_db_commit           = dc33e2f9ec2c420f99f77c5c33ae3faa79c999f2
resfinder_db_date             = Tue, 20 Mar 2018 16:49
pointfinder_db_dir            = staramr/databases/data/dist/pointfinder
pointfinder_db_url            = https://bitbucket.org/genomicepidemiology/pointfinder_db.git
pointfinder_db_commit         = ba65c4d175decdc841a0bef9f9be1c1589c0070a
pointfinder_db_date           = Fri, 06 Apr 2018 09:02
pointfinder_gene_drug_version = 050218
resfinder_gene_drug_version   = 050218
```

## hits/

The **hits/** directory contains the BLAST HSP nucleotides for the entries listed in the **resfinder.tsv** and **pointfinder.tsv** files. There are up to two files per input genome, one for ResFinder and one for PointFinder.

For example, with an input genome named **SRR1952908.fasta** there would be two files `hits/resfinder_SRR1952908.fasta` and `hits/pointfinder_SRR1952908.fasta`. These files contain mostly the same information as in the **resfinder.tsv** and **pointfinder.tsv** files. Additional information is the **resistance_gene_start** and **resistance_gene_end** listing the start/end of the BLAST HSP on the AMR resistance gene from the ResFinder/PointFinder databases. 

### Example

```
>aadA1_3_JQ414041 isolate: SRR1952908, contig: contig00030, contig_start: 5355, contig_end: 4564, resistance_gene_start: 1, resistance_gene_end: 792, hsp/length: 792/792, pid: 100.00%, plength: 100.00%
ATGAGGGAAGCGGTGATCGCCGAAGTATCGACTCAACTATCAGAGGTAGTTGGCGTCATC
GAGCGCCATCTCGAACCGACGTTGCTGGCCGTACATTTGTACGGCTCCGCAGTGGATGGC
...
```

# Tutorial

A tutorial guiding you though the usage of `staramr`, interpreting the results, and comparing with antimicrobial resistances available on NCBI can be found at [staramr tutorial][tutorial].

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
usage: staramr search [-h] [--pointfinder-organism POINTFINDER_ORGANISM]
                      [-d DATABASE] [-n NPROCS]
                      [--pid-threshold PID_THRESHOLD]
                      [--percent-length-overlap-resfinder PLENGTH_THRESHOLD_RESFINDER]
                      [--percent-length-overlap-pointfinder PLENGTH_THRESHOLD_POINTFINDER]
                      [--exclude-negatives] [--exclude-resistance-phenotypes]
                      [--report-all-blast] [-o OUTPUT_DIR]
                      [--output-summary OUTPUT_SUMMARY]
                      [--output-resfinder OUTPUT_RESFINDER]
                      [--output-pointfinder OUTPUT_POINTFINDER]
                      [--output-settings OUTPUT_SETTINGS]
                      [--output-excel OUTPUT_EXCEL]
                      [--output-hits-dir HITS_OUTPUT_DIR]
                      files [files ...]

positional arguments:
  files

optional arguments:
  -h, --help            show this help message and exit
  --pointfinder-organism POINTFINDER_ORGANISM
                        The organism to use for pointfinder {salmonella}. Defaults to disabling search for point mutations. [None].
  -d DATABASE, --database DATABASE
                        The directory containing the resfinder/pointfinder databases [staramr/databases/data].
  -n NPROCS, --nprocs NPROCS
                        The number of processing cores to use [MAX CPU CORES].

BLAST Thesholds:
  --pid-threshold PID_THRESHOLD
                        The percent identity threshold [98.0].
  --percent-length-overlap-resfinder PLENGTH_THRESHOLD_RESFINDER
                        The percent length overlap for resfinder results [60.0].
  --percent-length-overlap-pointfinder PLENGTH_THRESHOLD_POINTFINDER
                        The percent length overlap for pointfinder results [95.0].

Reporting options:
  --exclude-negatives   Exclude negative results (those sensitive to antimicrobials) [False].
  --exclude-resistance-phenotypes
                        Exclude predicted antimicrobial resistances [False].
  --report-all-blast    Report all blast hits (vs. only top blast hits) [False].

Output:
  Use either --output-dir or specify individual output files

  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        The output directory for results [None].
  --output-summary OUTPUT_SUMMARY
                        The name of the output file containing the summary results. Not be be used with '--output-dir'. [None]
  --output-resfinder OUTPUT_RESFINDER
                        The name of the output file containing the resfinder results. Not be be used with '--output-dir'. [None]
  --output-pointfinder OUTPUT_POINTFINDER
                        The name of the output file containing the pointfinder results. Not be be used with '--output-dir'. [None]
  --output-settings OUTPUT_SETTINGS
                        The name of the output file containing the settings. Not be be used with '--output-dir'. [None]
  --output-excel OUTPUT_EXCEL
                        The name of the output file containing the excel results. Not be be used with '--output-dir'. [None]
  --output-hits-dir HITS_OUTPUT_DIR
                        The name of the directory to contain the BLAST hit files. Not be be used with '--output-dir'. [None]

Example:
	staramr search -o out *.fasta
		Searches the files *.fasta for AMR genes using only the ResFinder database, storing results in the out/ directory.

	staramr search --pointfinder-organism salmonella --output-excel results.xlsx *.fasta
		Searches *.fasta for AMR genes using ResFinder and PointFinder database with the passed organism, storing results in results.xlsx.
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
                Builds a new ResFinder/PointFinder database under staramr/databases/data if it does not exist

        staramr db build --dir databases
                Builds a new ResFinder/PointFinder database under databases/
```

## Database Update

Updates an existing download of the ResFinder and PointFinder databases.

```
usage: staramr db update [-h] [-d] [--resfinder-commit RESFINDER_COMMIT]
                         [--pointfinder-commit POINTFINDER_COMMIT]
                         [directories [directories ...]]

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
                Updates the default ResFinder/PointFinder database under staramr/databases/data
```

## Database Info

Prints information about an existing build of the ResFinder/PointFinder databases.

```
usage: staramr db info [-h] [directories [directories ...]]

positional arguments:
  directories

optional arguments:
  -h, --help   show this help message and exit

Example:
        staramr db info
                Prints information about the default database in staramr/databases/data

        staramr db info databases
                Prints information on the database stored in databases/
```

## Database Restore Default

Restores the default database for `staramr`.

```
usage: staramr db restore-default [-h] [-f]

optional arguments:
  -h, --help   show this help message and exit
  -f, --force  Force restore without asking for confirmation.

Example:
        staramr db restore-default
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
[tutorial]: doc/tutorial/staramr-tutorial.ipynb
