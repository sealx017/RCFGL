# C/Python implementation of Rapid Condition adaptive Fused Graphical Lasso (RCFGL)

### Jupyter notebook with walkthrough
- The Jupyter notebook entitled "RCFGL_workflow_notebook.ipynb" provides a thorough guide on how to use the package on an example gene-expression dataset with three conditions. If you have problem viewing the notebook, please click on this link instead, https://bit.ly/3IFkdon.


### Overview of the main functions
- The package provides implementation of one existing joint network estimation model entitled Fused Multiple Graphical Lasso (Yang et al. (2015)) and our proposed method named Rapid Condition adaptive Fused Graphical Lasso (RCFGL). The former is implemented as a function named RFGL and the latter as a function named RCFGL. 

- Both the models can jointly estimate co-expression networks under multiple conditions. RFGL is a non-condition adaptive method, whereas RCFGL is a condition adaptive method meaning that the latter is better suited to capture condition-specific patterns in the networks across the conditions. 


### Python modules required to be installed
- The package requires the following Python modules to be pre-installed,
  1. igl, install using: "conda install -c conda-forge igl"  (https://libigl.github.io/libigl-python-bindings/)
  2. pywt, install using: "conda install pywavelets"  (https://pywavelets.readthedocs.io/en/latest/)
  3. matplotlib, install using: "conda install -c conda-forge matplotlib"  (https://matplotlib.org/)
  4. networkx, install using: "conda install -c conda-forge networkx"  (https://networkx.org/)
  5. cppyy, install using: "python -m pip install cppyy"  (https://cppyy.readthedocs.io/en/latest/)
  6. prox_tv, install using: "pip install prox_tv" (https://pypi.org/project/prox_tv/)

* We recommend using Anaconda (https://www.anaconda.com/products/individual) and Python version > 3.9. The package can only be used on a Mac or Linux system not on Windows as the "prox_tv" module is not available for the latter.


### References

1. Friedman, Jerome, et al. "Pathwise coordinate optimization." Annals of applied statistics 1.2 (2007): 302-332.
2.  Condat, Laurent. "A direct algorithm for 1-D total variation denoising." IEEE Signal Processing Letters 20.11 (2013): 1054-1057.
3. Yang, Sen, et al. "Fused multiple graphical lasso." SIAM Journal on Optimization 25.2 (2015): 916-943.
4. Barbero, Alvaro, and Suvrit Sra. "Modular proximal optimization for multidimensional total-variation regularization." The Journal of Machine Learning Research 19.1 (2018): 2232-2313.
5. Lyu, Yafei, et al. "Condition-adaptive fused graphical lasso (CFGL): An adaptive procedure for inferring condition-specific gene co-expression network." PLoS computational biology 14.9 (2018): e1006436.



