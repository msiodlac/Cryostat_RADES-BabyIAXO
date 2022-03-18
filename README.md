# Cryostat_RADES-BabyIAXO
Supporting code for the master thesis "Thermodynamic Design of a Cryostat for a Radiofrequency Cavity Detector in Search of Dark Matter".
See this thesis for more information and to understand the code: https://cds.cern.ch/record/2804052?ln=en

The folder structure includes three big components:

HeatLoadEstimation
Code used to estimate the heat load in the cryostat.
2D thermal model for coaxial cables, 1D thermal model for the tuning rods, and multiple smaller calculations.

ThermodynamicModel
Code for a thermodynamic model that represents the cooling system in the present cryostat.
Additionally there is a file for the modelling of the experimental setup at the CERN Cryolab at the status of 02/2022.
The experimental setup is used to validate the present model.

ExperimentPostProcessing
Code for the post processing of the data gathered from the experimental setup at the CERN Cryolab in 02/2022.
That includes the calculation of the average value of the steady state measuring point, the uncertainty, mass flow rate, and effectiveness.

python_shared
Here, the HEPAK wrapper is located to use HEPAK in python, which is needed in all the previous described codes.
supcode_py collects all small functions that are used in the post processing of the experimental data.
