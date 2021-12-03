###############
# Description #
###############
This matlab package contains files implementing the dual subgradient method for finding the center of the Grassmannian minimum enclosing ball for a finite collection of subspaces of differing dimensions, and the order-selection rules for the dimension of that center, as described in the paper:
"On a minimum enclosing ball of a collection of linear subspaces" by TIM MARRINAN, P.-A. ABSIL, and NICOLAS GILLIS.
Linear Algebra and its Applications 625 (2021): 248-278.

############
# Contents #
############
The matlab package includes the following directories and files:

01. README.txt
02. demo.m

\src\
03. GMEB_DataGen.m
04. GMEB_DualSubgrad.m
05. GMEB_OrderSelection.m
06. main.m

\src\helper\
07. aboxplot.m
08. diffusionKernel.m
09. GMEB_Primal.m
10. GrDist.m
11. GrExp.m
12. GrLog.m
13. KarcherMean.m
14. knee_pt.m
15. SimplexProj.m
16. quartile.m

\examples\
17. GMEB_Exp001.m
18. GMEB_Exp002.m
19. GMEB_Exp003.m
20. GMEB_IllustrativeExample.m
21. GMEB_ScenarioSpecification.m

############
# Abstract #
############
This paper concerns the minimax center of a collection of linear subspaces.  When the subspaces are $k$-dimensional subspaces of $\mathbb{R}^n$, this can be cast as a problem of minimizing the maximum distance on a Grassmann manifold, Gr$(k,n)$.  For subspaces of different dimension, the setting becomes a disjoint union of Grassmannians rather than a single manifold, and the problem is no longer well-defined. However, natural geometric maps exist between these manifolds with a well-defined notion of distance for the images of the subspaces under the mappings. Solving the initial problem in this context leads to a candidate minimax center on each of the constituent manifolds, but does not inherently provide intuition about which candidate is the best representation of the data.  Additionally, the solutions of different rank are generally not nested so a deflationary approach will not suffice, and the problem must be solved independently on each manifold.  We propose and solve an optimization problem parametrized by the rank of the minimax center.  The solution is computed using a subgradient algorithm on the dual. By scaling the objective and penalizing the information lost by the rank-$k$ minimax center, we jointly recover an optimal dimension, $k^*$, and a central subspace, $\*U^* \in$ Gr$(k^*,n)$ at the center of the minimum enclosing ball, that best represents the data.

##############
# File Usage #
##############
The file demo.m can be run with no modifications to add the folders to the path and run a small synthetic example. It computes the GMEB center of the data, selects the ideal order using the proposed order-selection rules, and visualizes an embedding of the subspace data and center into 2D or 3D.

Alternatively, the file main.m can be run with no modifications to recreate the figures from the associated paper.  This will call each of the functions in the 'examples' directory, whose parameters are set by GMEB_ScenarioSpecification.m. The variables 'queued' and 'nRuns' can be modified to control which experiments are run and for how many iterations.

To re-run the experiments and generate the plots from the paper:
01. Run main.m

In the paper the experiments are each run for nRuns = 100 Monte Carlo trials, however this takes numerous hours to complete. The slowest part is Exp_003.m where the order-selection rules are applied for each value of the independent variable and require max_i{dim(X_i)}+1 minimum enclosing balls to be computed each time. Alternatively, the behavior of the methods can be assessed from a much smaller number of independent trials to avoid the long computation time. 

On a laptop with a 2.6 GHz Intel Core i7-8850H processor and 16 GB of RAM, the experiments from the paper run for nRuns = 5 take approximately 25 minutes. 

To run your own experiments:
01. In the file GMEB_ScenarioSpecification.m, modify the parameters of the 'custom' scenario as desired.
02. Depending on your desired outcome, look to the files demo.m, GMEB_IllustrativeExample.m, GMEB_Exp001.m, GMEB_Exp002.m, or GMEB_Exp003.m for examples of how to use the data generation function, minimum enclosing ball dual subgradient method, and order-selection rules.


###########
# Contact #
###########

In case of questions, suggestions, problems etc. please send an email.

Tim Marrinan:
timothy.marrinan@umons.ac.be

This matlab package is hosted at:
https://sites.google.com/site/nicolasgillis/code
