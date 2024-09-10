# Quantifying and estimating dependence via sensitivity of conditional distributions

The paper is available at https://arxiv.org/abs/2308.06168

## How to Use
* Core_Functions.R contains all functions needed to calculate an estimate of the functional Lambda_phi. If you want to use our method on your own data you will first have to apply the function ECBC from the qad package[1] to your data and then use the function Lambda_phi from this file on the output.
* Simulation.R runs the Simulation. You will need to create a Folder "Simulation_Results" first.
* Figures.R creates the Figures for the Simulation. You will need to create a Folder "Figures" first.
* Data_Example.R creates the Figures and Table from the Real Data Example. You will need to create a Folder "Figures" first. Data[2] was taken from physionet.org[3] https://physionet.org/content/maternal-visceral-adipose/1.0.0/

[1] Kasper T, Griessenberger F, Junker R, Petzel V, Trutschnig W (2022). _qad: Quantification of Asymmetric Dependence_. R package version 1.0.3, <https://CRAN.R-project.org/package=qad>.

[2] Rocha, A. d. S., von Diemen, L., Kretzer, D., Matos, S., Rombaldi Bernardi, J., & Magalhães, J. A. (2020). Visceral adipose tissue measurements during pregnancy (version 1.0.0). PhysioNet. https://doi.org/10.13026/p729-7p53.

[3] Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.
