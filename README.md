# PROmiRNA [![Build Status](https://travis-ci.org/marsicoLab/PROmiRNA.svg?branch=master)](https://travis-ci.org/marsicoLab/PROmiRNA)

[PROmiRNA][1] is a program for annotating miRNA promoters in human, as well as other species. It uses deepCAGE data and integrated cage tag counts and other promoter features, such as CpG content, conservation and TATA box affinity, to score the potential of a candidate region to be a promoter. This is the reimplemented, faster version of PROmiRNA (old python version available [here][2]) now capable of fast-processing of large CAGE libraries such as available from FANTOM4, FANTOM5 and ENCODE.

## Prerequisites

* GCC ≥ 7
* CMake ≥ 3.0
* R with packages RInside and Rcpp installed

## Installation

Clone project with:

```
$ git clone https://github.com/marsicoLab/PROmiRNA
$ cd PROmiRNA
$ git submodule update --init
```

Build the project with:

```
$ mkdir build && cd build
$ cmake ../src/
$ make
```

Run PROmiRNA with:

```
$ ./PROmiRNA [options]
       
  -g, --genome       Path to reference genome (.fa, .fasta)
  -c, --coding       Path to coding/gene regions (.gtf)
  -s, --starts       Path to gene start regions (.gff)
  -r, --repeats      Path to repeat regions (.bed)
  -a, --annotation   Path to pre-miRNA annotation file from miRBase(<species>.gff3)
  -m, --mirna        Path to miRBase database file (mirna.txt)
  -n, --context      Path to miRBase database file (mirna_context.txt)
  -p, --psem         Path to position weight matrix (.psem)
  -w, --wig          Path to phastcons file (.wig, .txt)
  -i, --cage         Path to directory containing CAGE libraries (files need to be BED files)
  -t, --threads      Number of threads to use
```

R packages can be installed within R environment:

```
$ R
> install.packages("Rcpp")
> install.packages("RInside")
```

### Important to note:

The folder ```external_data``` contains already the current up-to-date ```mirna.txt``` and ```mirna_context.txt``` files from the newest miRBase version (release 22). In addition, the pre-miRNA annotation for the human species ```hsa.gff3``` is also already shipped with this repository. In case other species are needed please download the annotation file from http://www.mirbase.org/ftp.shtml. 

The PSEM matrix is also included in the ```external_data``` folder (```TATA_box_jaspar.psem```).

All input files need to match the reference genome that is used.

The conservation file with phastcons values can be downloaded from UCSC. It is required that in the same location as the ```.wig``` file the corresponding ```.wib``` file with same prefix/file name is located. Files can be found here (example hg19):
* http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/phastCons46way.txt.gz (need to have ```.txt``` or ```.wig``` ending for running PROmiRNA)
* http://hgdownload.cse.ucsc.edu/gbdb/hg19/multiz46way/phastCons46way.wib

The genes file can be downloaded from Ensembl (```.gtf``` ending required), for example for hg19 at ftp://ftp.ensembl.org/pub/grch37/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh37.92.gtf.gz.


[1]:https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-8-r84
[2]:http://promirna.molgen.mpg.de/
