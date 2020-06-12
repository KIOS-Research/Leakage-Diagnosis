# Leakage Detection and Localization in Water Distribution Systems: a Model Invalidation Approach

Abstract:
Model-based methodologies can assist in addressing the challenging problem of leakage detection and localization in water distribution systems. However, this is not trivial due to inherent non-linearities and parametric uncertainties. Most importantly, due to the small number of available sensor measurements compared to the number of system states, the inverse problem for estimating leakages is highly under-determined. In this work we propose the utilization of a priori available information about the system to formulate a hydraulic model of the system in its non-linear form, in which uncertainties are modeled by intervals defined by a lower and upper bound. A novel optimization-based methodology then utilizes pressure and flow measurements to perform leakage detection through model-invalidation. A modification of the optimization algorithm is activated in the case of a detection to refine possible leak locations and retain only the ones that can be explained by the interval model and available measurements from multiple time-steps. The proposed methodology is demonstrated on a benchmark network and evaluated using a leakage diagnosis benchmark dataset.

# Requirements
This work uses the EPANET-Matlab-Toolkit which is a Matlab class for EPANET libraries.
Please install the toolkit before use.
For more information see https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit#EPANET-MATLAB-Toolkit .

This work also uses Gurobi Optmization, LLC (http://www.gurobi.com).
For obtaining an academic licence please visit https://www.gurobi.com/academia/academic-program-and-licenses/.

# Contributors
Stelios Vrachimis, KIOS Research Center for Intelligent Systems and Networks, University of Cyprus