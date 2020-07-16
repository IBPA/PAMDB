# PAMDB
A Parts Database with Consensus Parameter Estimation for Synthetic Circuit Design

### What is PAMDB?
PAMDB aims to bridge the gap between published circuits and lack of universal quantitive parameters for their parts. PAMDB includes parameter values inferred from hundreds of papers under a single unified model. PAMDB can be used to rank parts based on their activity, acquire parameter values and quantify their uncertainty, simulate circuits and visualize the part universe for the published synthetic biology work so far

### Dependencies
* Gcc
* Matlab


### Installation
```
git clone https://github.com/IBPA/PAMDB.git
./compile.sh
```

### Running
* Step1: Generate the optimization problem for MATLAB 
```
./src/sbrome
```
* Step2: Use MATLAB to estimate parameter values by running file 
```
MATLAB/ACS_2016/real_benchmark_fitting.m
```
This will print all parameters, their value, and confidence interval into file 
```
ACS_2016/simultaneous_parameter_CI.dat
```

### Support

If you have any questions about PAMDB, please contact Linh Huynh (huynh@ucdavis.edu) or Ilias Tagkopoulos (itagkopoulos@ucdavis.edu).

### Citation
 L. Huynh, and I. Tagkopoulos, “A Parts Database with Consensus Parameter Estimation for Synthetic Circuit Design”, ACS Synthetic Biology (2016) [\[link\]](https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00205)

### Licence
See the [LICENSE](./LICENSE) file for license rights and limitations (Apache2.0).

### Acknowledgement
This work was supported by the NSF CAREER grant #1254205 to IT.

