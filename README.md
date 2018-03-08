# `staramr`

`staramr` (*AMR) scans bacterial genome contigs against both the [ResFinder][resfinder-db] and [PointFinder][pointfinder-db] databases and complies a summary report of detected antimicrobial resistance genes and matching phenotypes.

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
staramr.py search --pointfinder-organism -o out *.fasta
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

# Citations

If you find `staramr` useful, please consider citing this GitHub repository (https://github.com/phac-nml/staramr) as well as the original ResFinder and PointFinder publications.

> **Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup FM, Larsen MV**. 2012. Identification of acquired antimicrobial resistance genes. J. Antimicrob. Chemother. 67:2640–2644. doi: [10.1093/jac/dks261](resfinder-cite)

> **Zankari E, Allesøe R, Joensen KG, Cavaco LM, Lund O, Aarestrup F**. PointFinder: a novel web tool for WGS-based detection of antimicrobial resistance associated with chromosomal point mutations in bacterial pathogens. J Antimicrob Chemother. 2017; 72(10): 2764–8. doi: [10.1093/jac/dkx217][pointfinder-cite]

[resfinder-db]: https://bitbucket.org/genomicepidemiology/resfinder_db
[pointfinder-db]: https://bitbucket.org/genomicepidemiology/pointfinder_db
[resfinder-cite]: https://dx.doi.org/10.1093/jac/dks261
[pointfinder-cite]: https://doi.org/10.1093/jac/dkx217
[Bioconda]: https://bioconda.github.io/
[requirements.txt]: requirements.txt