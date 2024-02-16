# PET2PHA
Repository with the resources for the in silico analysis for biodegradation of PET and accumulation of PHA by P. putida
## Citation

[![DOI](DOI)](https://doi.org/10.6084/m9.figshare.23760606.v1)

Nogales, Juan; Manoli, Maria-Tsampika; Gargantilla Becerra, Alvaro; del Cerro Sánchez, Carlos; Rivero, Virginia; Prieto, Auxiliadora (2023). PET2PHA: Source code and resources for in silico metabolic analysis in P. putida to accumulate PHA with PET as carbon source. figshare. Software. https://doi.org/10.6084/m9.figshare.23760606.v1

```
PET2PHA
  |---code
  |     |---utils
  |
  |---data
  |     |---cameo
  |     |---CometsResults  
  |
  |---figures
  |
  |---models

```
The scripts used for the search of Growth-Coupled strategies described in the paper consist on the sequential execution of `code/kt_simplification.py` and `code/gc_strategy_search.py`. For doing this, you need to have `python`, `conda` and `matlab>=MATLAB 2017` installed in your computer. Once all of them are installed, you can create the conda environment holding all the needed python packages to run the workflow. This can be done from the `environment.yml` in the root with the following command:
```
conda env create -f environment.yml
```
Then activate the environment with:
```
conda activate strain-design
``` 
Eventually run the two main scripts by executing them from shell:
```
nohup python code/kt_simplification.py && python code/gc_strategy_search.py $
```
To analyse the results, just execute the notebooks within the `code` folder.
