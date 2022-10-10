# Software implementation of Rapid Condition adaptive Fused Graphical Lasso (RCFGL)

### Overview of the main functions
- The package provides implementation of one existing joint network estimation model named Fused Multiple Graphical Lasso and one proposed method named Rapid Condition adaptive Fused Graphical Lasso. The former is implemented as a function named RFGL and the latter as a function named RCFGL. 

- Both the models can jointly estimate networks of multiple conditions. RFGL is a non-condition adaptive method, whereas RCFGL is a condition adaptive method meaning that the latter is better suited to capture condition-specific patterns in the networks across the conditions. 

### Jupyter notebook with walkthrough
- Refer to the Jupyter notebook entitled "RCFGL_workflow_notebook.ipynb" to get a thorough walk-through on how to use the package on an example dataset with three conditions. If you have problem viewing the notebook, please click on this link instead, https://bit.ly/3IFkdon.

### Python modules required to be installed
- The notebook requires the following Python modules to be pre-installed,
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
3. Danaher, Patrick, Pei Wang, and Daniela M. Witten. "The joint graphical lasso for inverse covariance estimation across multiple classes." Journal of the Royal Statistical Society. Series B, Statistical methodology 76.2 (2014): 373.
4. Yang, Sen, et al. "Fused multiple graphical lasso." SIAM Journal on Optimization 25.2 (2015): 916-943.
5. Barbero, Alvaro, and Suvrit Sra. "Modular proximal optimization for multidimensional total-variation regularization." The Journal of Machine Learning Research 19.1 (2018): 2232-2313.
6. Lyu, Yafei, et al. "Condition-adaptive fused graphical lasso (CFGL): An adaptive procedure for inferring condition-specific gene co-expression network." PLoS computational biology 14.9 (2018): e1006436.



