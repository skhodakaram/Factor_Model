- This repository contains MATLAB files with which the numerical results of the paper **A Saddle Point Algorithm for Robust Data-Driven Factor Model Problems** 
  [(link)](https://arxiv.org/abs/2506.09776) have been produced, including a main script file named `FM_main.m` and the corresponding functions. 

- The factor model problem focuses on the decomposition of a covariance matrix $\Sigma$ to low-rank and diagnonal positive semidefinite matrices. In practice, $\Sigma$ is 
  often not available, and only its empirical counterpart, $\widehat{\Sigma}$, is available. To robustify to this approximation error, a common practice is to consider a 
  family of covariance matrices in the vicinity of $\widehat{\Sigma}$ defined as
  
  $$ B_{\varepsilon}^{d}(\widehat{\Sigma}) := \\{ \Sigma \succeq 0 ~ \colon ~ d(\Sigma,\widehat{\Sigma}) \leq \varepsilon \\}, $$
  
  where $d$ is a generic distance function in the space of matrices, and $\varepsilon$ is the radius (size) of the set. Considering the target decomposition and the 
  uncertainty set $B_{\varepsilon}^{d_{+}}(\widehat{\Sigma})$, our robust data-driven factor model problem can be formulated as the optimization problem 
  
  $$ \begin{array}{rcl} 
		& \min\limits_{L, D } & \textrm{Tr}(L) \\
		&\textrm{s.t.} & L \in S_{+}, D \in D_{+}, \\
            	&& L + D \in B_{\varepsilon}^{d}(\widehat{\Sigma})
     \end{array} $$

- The numerical results have been produced as follows.

  **1. Convergence:** After setting the parameters values based on the paper, in the `FM_MaxMin` function, the stopping condition has been commented out. After
     running the code, the function `plotErrorObjVal` was used to plot the corresponding diagram. 

  **2. Estimation of the ground-truth $\Sigma_{\textrm{True}}$:** After fixing the parameters values based on the paper and setting the value of the parameter *SweetSpot* 
     to 1, the code was executed.

  **3. Execution time:** The `FM_main` file was executed for different dimensions in the set {50, 100, ..., 300} with the parameters set based on the paper.


# Requirements:

  YALMIP, MOSEK, MATLAB Optimization, and MATLAB Statistics and Machine Learning toolboxes are used for comparison purposes. However, they are not used in the algorithm
  and you may have to remove the related code if you do not have them.


# Instruction:

  Run `FM_main.m` file. You may set the values of *Flags* based on their description.


# Acknowledgments:

  This research is supported by the ERC Starting Grant TRUST-949796 and was conducted in Delft University of Technology.

# Citing:  

  Please cite the following paper if you use Factor_Model:

```
@article{khodakaramzadeh2025saddle,
  title={A Saddle Point Algorithm for Robust Data-Driven Factor Model Problems},
  author={Khodakaramzadeh, Shabnam and Shafiee, Soroosh and de Albuquerque Gleizer, Gabriel and Mohajerin Esfahani, Peyman},
  journal={preprint available at arXiv:2506.09776},
  year={2025}
}

```

- This work is licensed under a BSD 3-Clause OSS licence.


Shabnam Khodakaramzadeh






