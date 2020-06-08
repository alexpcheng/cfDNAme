# ddcfDNA

## Installing ddcfDNA

* install miniconda3 (https://conda.io/miniconda.html)
* ~/miniconda3/bin/conda config --add channels r
* ~/miniconda3/bin/conda config --add channels bioconda
* ~/miniconda3/bin/conda create --name cfDNA -c bioconda python=3.5 --file workflow/requirements3.txt       
* install autograd using pip: 
source ~/miniconda3/bin/activate cfDNA; 
pip install --user autograd;

*clone this github.

## Installing WASP

* download the scripts in the mapping directory from WASP github https://github.com/gmcvicker/WASP
* copy the scripts from WASP mapping to {cfDNA scripts path}/WASP 
* install anaconda2 from https://www.continuum.io/downloads
* create anaconda environment for wasp (this is a different environment for cfDNA since it is python 2.7)
  * ~/anaconda2/bin/conda config --add channels r
  * ~/anaconda2/bin/conda config --add channels bioconda
  * ~/anaconda2/bin/conda create -n wasp python=2.7 --file workflow/requirements2.txt

## Full workflow 
The workflow download the required files, creates the input file for the inference algorithm, runs the inference for each sample and then collect the results.

See README under the workflow directory

## Running the inference step

The inference step executble will run only on a linux machine (linux-x86_64-3.5). Please contact us, if you need to run on a different machine. 

1. Create an input file using the workflow.
2. Run ddcfDNA inference:
  Before running the executable, you will need to set LD_LIBRARY_PATH to point to your miniconda installation using:
  export LD_LIBRARY_PATH=/home/<YOUR USER NAME>/miniconda3/envs/cfDNA/lib
  export PYTHONHOME=/home/<YOUR USER NAME>/miniconda3/envs/cfDNA
  You need to run these lines before using the executable.
  After running the executable return these variable values to their original values
  (usually using 
  export LD_LIBRARY_PATH=;
  and 
  export PYTHONHOME=;)

  so a run shell may look like that:

  source ~/miniconda3/bin/activate cfDNA;
  export LD_LIBRARY_PATH=/home/<YOUR USER NAME>/miniconda3/envs/cfDNA/lib; 
  export PYTHONHOME=/home/<YOUR USER NAME>/miniconda3/envs/cfDNA; 
  ddcfDNA/cfDNA_infer_donor_fraction.exe --help;
  export LD_LIBRARY_PATH=; 
  export PYTHONHOME=; 
  
3. If multiple samples are run you can use this script to colelect the results: cfDNA_collect_inference_results.py
  Otherwise, the output is a table. Each row represent one possible donor population. You should use the results in the row that is marked as selected by the algorithm. 

## Citation

Please cite http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005629 (PMID: 28771616)
