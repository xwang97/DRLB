# DRLB
This is the source code of DRLB, a deep learning-based Boolean matrix bias removal method published in UAI 2024.

## Folder architecture
1. The .py files in the root directory contain all the models and functions of our method.
2. The R file evaluate.R is used for model evaluation, i.e., visualization and metric computation
3. data: six sample datasets for the simulated scenarios and two real world data
4. comparison: BMF methods used in this paper


## Dependencies
torch 1.13  
pandas 1.5  
numpy 1.24  

## R dependencies (not required if you only need the debiasing step by DRLB)
Matrix  
gplots  


## Usage
1. Run main_gpu.py to debias your data with DRLB:
   ```
   python main_gpu.py -dataname data -batch_size size -alpha alpha
   ```
   For example:
   ```
   python main_gpu.py -dataname data/simu1/data.csv -batch_size 8 -alpha 3
   ```
3. After running the code, a result named debiased.csv that contains the debiased data will be saved in the root dir
   
 
 ## A demo for running on a simulated dataset
 1. In the ternimal, use the following command
    ```
    python main_gpu.py -dataname data/simu1/data.csv -batch_size 8 -alpha 3
    ```
    After running, the debiased file of simulated data 1 will be saved.
    Then you can follow the next steps for evaluation and comparison.
 2. In R, use command 'source("evaluate.R")' to compile the evaluation functions;  
    use 'source("BIND.R")' to compile the functions for BIND;  
    use 'source("comparison/MEBF.R"), source("comparison/ASSO.R"), source("comparison/PANDA.R")' to compile the BMF methods.
 4. Open "run.R", run the code row by row in R to see the evaluation and comparison results:  
    (1) Load the dataset. The data, pattern and bias matrix for simulated data 1 will be loaded.  
    (2) Run BMF without debiasing. A BMF method will run on the data (MEBF as default, you can switch to ASSO or PANDA by uncommenting the following rows of MEBF). When finished, you can print out the reconstruction         error without debiasing.  
    (3) Debias with BIND. The data will be debiased by BIND, the resulted matrix will be visualized. Then the same BMF method used in (2) will be run here, and the reconstruction error will be printed (some
        improvement may be seen here).  
    (4) Debias with DRLB. The debiased matrix of DRLB will be loaded and visualized. The same BMF method will be run again and the reconstruction error will be printed (significant improvement will be seen).
 6. Running CG will require some installation. We recommend to follow the instructions in their github page (refer to the paper for the link) to install it. They already provided functions for running BMF and
    checking the reconstruction error. You can use the data before/after debiasing to compare the results.
