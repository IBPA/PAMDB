# PAMDB
A Parts Database with Consensus Parameter Estimation for Synthetic Circuit Design

### What is PAMDB?

DeepPep, is a **protein identification** software which uses deep-convolutional neural network to predict the protein set from a proteomics mixture, given the sequence universe of possible proteins and a target peptide profile.

### Dependencies
* [torch7](http://torch.ch/docs/getting-started.html)
* luarocks install cephes
* luarocks install csv
* [SparseNN](https://github.com/ameenetemady/SparseNN/)
* python3.4 or above
* [biopython](http://biopython.org/wiki/Download)



### Installation
```
git clone https://github.com/ameenetemady/MyCommon.git
git clone https://github.com/DeepPep/DeepPep.git
```

### Running
* Step1: prepare a directory containing your input files (with exact names):

  * ```identification.tsv```: tab-delimeted file:  **column1**: peptide, **column2**: protein name, **column3**: identification probability
  * ```db.fasta```: reference protein database in fasta format.

* Step2: ```python run.py directoryName```

Upon completion, ```pred.csv``` will contain the predicted protein identification probabilities.

### Benchmark Datasets
There are [7 example datasets](https://github.com/DeepPep/public/tree/master/data) (used for benchmarking in the paper). Each dataset is generated from MS/MS raw files using TPP pipeline. For example, to run the [18Mix benchmark dataset](https://github.com/DeepPep/public/tree/master/data/18mix), simply run the following:

```
python run.py data/18Mix
```
### Support

If you have any questions about PAMDB, please contact Linh Huynh (huynh@ucdavis.edu) or Ilias Tagkopoulos (itagkopoulos@ucdavis.edu).

### Citation
 L. Huynh, and I. Tagkopoulos, “A Parts Database with Consensus Parameter Estimation for Synthetic Circuit Design”, ACS Synthetic Biology (2016) [\[link\]](https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00205)

### Licence
See the [LICENSE](./LICENSE) file for license rights and limitations (Apache2.0).

### Acknowledgement
This work was supported by the NSF CAREER grant #1254205 to IT.

