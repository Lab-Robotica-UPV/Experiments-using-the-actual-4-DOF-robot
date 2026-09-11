# Linear Algebra-Based Control (LABC) for a 4-DOF Parallel Robot
This repository contains the experimental datasets, MATLAB analysis scripts, and simulation models associated with the manuscript:
Scaglia, G., Vallés, M., Díaz-Rodríguez, M., Pulloquinga, J.L., Zamora-Ortiz, P., Valera, A. "Linear algebra-based trajectory control of a 4-DOF parallel robot with stability analysis and experimental validation" Mechanics Based Design of Structures and Machines, 2025.

# Data Format
Each run_X/ folder contains TSV files with 10 columns:
Columns	Content
1–2	Encoder timestamp (sec + nanosec)
3–6	Init/Current time (internal, not used in analysis)
7–10	Joint data: q₁,₃ · q₂,₃ · q₃,₃ · q₄,₂ [rad]
•	Sampling frequency: 100 Hz (dt = 0.01 s)
•	Duration per run: 52.1 s (5211 samples)
•	N = 12 paired runs per controller (PID, CTC, LABC)
•	Execution order: randomized to avoid systematic bias
The x_mocap.txt files (LABC only) contain OptiTrack end-effector measurements:
Columns	Content
1–2	Encoder timestamp
3–6	Init/Current time (internal)
7	x [m]
8	z [m]
9	θ [rad]
10	ψ [rad]
Note: run_7 (marker occlusion) and run_8 (duration mismatch, 67.8 s vs 52.1 s) are excluded from the MoCap analysis. The scripts handle this automatically via SKIP_RUNS = [7, 8].

________________________________________
# Requirements
MATLAB Version
•	MATLAB R2020b or later
Required Toolboxes
•	Statistics and Machine Learning Toolbox — required for: 
o	signrank (Wilcoxon signed-rank test)
o	friedman (Friedman test)
o	bootstrp (bootstrap confidence intervals)
o	swtest (Shapiro-Wilk normality test — included as a local function)

________________________________________
# How to Run the Analysis
Step 1 — Unzip the files containing the experimental data from the PID, Computed torque and Linear algebra-based controllers executions




Step 2 — Joint-space statistical analysis
analysis_complete.m

This script performs the complete statistical comparison (PID vs CTC vs LABC) following the methodology in Nehmzow (2006):
1.	Loads all 12 runs per controller from data/
2.	Computes scalar metrics per run: RMSE, MaxError, IAE, TVE, TVC
3.	Runs Friedman test (global, k=3 controllers)
4.	Runs Wilcoxon signed-rank post-hoc tests with Bonferroni correction (α_adj = 0.0042)
5.	Computes effect size r and bootstrap confidence intervals
6.	Generates boxplots, forest plots, and trade-off figures
7.	Outputs a LaTeX-ready results table
Expected output:
PASO 1: FRIEDMAN TEST
  q_{1,3}: χ²=18.0, p<0.001 → post-hoc
  ...
PASO 2: WILCOXON POST-HOC
  q_{1,3} LABC vs PID: W=78, p=0.0005, r=0.88 *
  ...



Step 3 — End-effector MoCap analysis

mocap_analysis.m
This script validates the LABC end-effector precision using the independent OptiTrack measurements:
1.	Loads x_mocap.txt from each valid LABC run
2.	Synchronizes MoCap with controller reference via cross-correlation
3.	Mean-centres signals to remove systematic calibration offset
4.	Applies p95 outlier filter per run to exclude marker occlusions
5.	Computes RMSE and max error over N=10 valid runs
6.	Analyses cycle-to-cycle variability (5 cycles per run, 7.85 s each)
7.	Generates trajectory, error, and boxplot figures
Expected output:
run_ 1: lag=+20.0 ms  RMSE x=3.29 mm  z=2.35 mm
...
x:  RMSE=6.59±1.69 mm   max=14.18±4.38 mm
z:  RMSE=3.17±1.08 mm   max=6.46±1.92 mm
ψ:  RMSE=5.51±0.96 mrad max=11.58±2.49 mrad

________________________________________
# Statistical Methodology
The analysis follows the experimental methodology for robotics described in:
Nehmzow, U. (2006). Scientific Methods in Mobile Robotics: Quantitative Analysis of Agent Behaviour. Springer, London. https://doi.org/10.1007/1-84628-260-8

Key design decisions:
Decision	                  Rationale
Paired design	              Same robot, same trajectory, same session — controls between-session variability
Scalar collapse per run	    Avoids pseudo-replication (Nehmzow 2006, Sec. 2.6.3)
Friedman test	              Non-parametric equivalent of repeated-measures ANOVA for k=3 groups
Wilcoxon signed-rank	      Non-parametric paired test, appropriate for N=12 without normality assumption
Bonferroni correction	      α_adj = 0.05/12 = 0.0042 (3 pairs × 4 joints)
Effect size W	              W=78 = maximum attainable value for N=12, indicating unanimous superiority

________________________________________
# Reference Trajectory
The elliptical rehabilitation trajectory used in all experiments:
x_ref = x₀ + A₁·cos(2πft + φ₀)
z_ref = z₀ + A₂·(1 + sin(2πft + φ₀))
θ_ref = θ₀  (constant)
ψ_ref = ψ₀  (constant)
Parameters: x₀ = −0.0326 m, z₀ = 0.91 m, A₁ = 0.1 m, A₂ = 0.072 m, f = 0.127 Hz, φ₀ = 3π/2, θ₀ = 0.0969 rad, ψ₀ = −0.0407 rad.
The reference is provided as ref_trajectories/ref_cart.txt (5211 rows × 4 columns: x, z, θ, ψ at 100 Hz).

________________________________________
# Hardware Setup
•	Robot: 3UPS+RPU 4-DOF parallel robot (2T2R)
•	Actuators: 4 × Maxon RE40 DC motors with ball-screw
•	Encoders: ENC DEDL 9149 (2000 counts/rev, 0.18°/count)
•	Control PC: Intel Core i7-7700, 8 GB RAM, Linux Ubuntu + ROS 2
•	DAQ: Advantech PCI-1720 (DAC) + PCI-1784 (encoder)
•	Sample time: 10 ms (100 Hz)
•	MoCap: 14 × OptiTrack Flex 13 cameras, accuracy < 0.1 mm, 120 Hz


________________________________________
# Citation
If you use this data or code in your research, please cite:
@article{scaglia2025labc,
  author  = {Scaglia, Gustavo and Vallés, Marina and Díaz-Rodríguez, Miguel
             and Pulloquinga, José L. and Zamora-Ortiz, Pau and Valera, Angel},
  title   = {Linear algebra-based trajectory control of a 4-{DOF} parallel
             robot with stability analysis and experimental validation},
  journal = {},
  year    = {},
  doi     = {}
}

________________________________________
# License
The data and scripts in this repository are released under the Creative Commons Attribution 4.0 International License (CC BY 4.0).

________________________________________
# Contact
Corresponding author: Angel Valera — giuprog@isa.upv.es
Instituto de Automática e Informática Industrial (ai2) Universitat Politècnica de València Camino de Vera s/n, Valencia 46022, Spain

