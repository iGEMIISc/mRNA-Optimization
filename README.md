This repository uses **RiboTree-mRNA** to optimize our mRNA sequences. However, the files related to RiboTree have been omitted from the repository since the tool is distributed with a license. This is work in progress, and results will soon be added.

## Some additional things to install

```bash
$ sudo apt-get install gfortran
$ sudo apt-get install liblapack-dev liblapacke-dev
$ pip install setuptools networkx matplotlib scipy numpy pandas
```

## Arnie File

```
#Path to EternaFold
eternafold: /path/to/mRNA-Optimization/eternafold/src

# Vienna RNAfold 2 Linux build from source:
vienna: /path/to/mRNA-Optimization/ViennaRNA-2.6.3/src/bin
vienna_2: /path/to/mRNA-Optimization/ViennaRNA-2.6.3/src/bin

# CONTRAfold build
contrafold_2: /path/to/mRNA-Optimization/contrafold-se/src

# LinearFold build
linearfold: /path/to/mRNA-Optimization/LinearFold/bin

#LinearPartition build
linearpartition: /path/to/mRNA-Optimization/LinearPartition/bin

#directory to write temp files
TMP: /tmp
```
---

### Open-source Contributions 

#### Open PRs
1. [Fixed InferenceEngine.hpp](https://github.com/csfoo/contrafold-se/pull/4) in `contrafold-se`

#### Merged PRs
1. [Modify documentation for adding Vienna to arnie file ](https://github.com/DasLab/arnie/pull/30) in `arnie`
2. [Fixed error with file deletion on first run](https://github.com/DasLab/arnie/pull/31) in `arnie`

#### Other bugs reports

1. Deprecated pandas method in `RiboTree-mRNA`
