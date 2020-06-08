

export LD_LIBRARY_PATH=/home/eilon/miniconda3/lib
export PYTHONHOME=/home/eilon/miniconda3/envs/cfDNA/
gcc -I/home/eilon/miniconda3/envs/cfDNA/include/python3.5m -I/home/eilon/miniconda3/envs/cfDNA/lib/python3.5/site-packages/pandas/ -I/home/eilon/miniconda3/envs/cfDNA/lib/python3.5/site-packages/numpy/core/include -L/home/eilon/miniconda3/envs/cfDNA/lib -lpython3.5m cfDNA_infer_donor_fraction.c -ocfDNA_infer_donor_fraction.exe

