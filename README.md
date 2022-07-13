# Switching from weak to strong cortical attachment of microtubules accounts for the transition from nuclear centration to spindle elongation in metazoans

## Authors
- Shohei Tada: <>;
- Yoshitaka Yamazaki: <>;
- Takahiro G Yamada: <yamada@fun.bio.keio.ac.jp>;
- Akira Funahashi: <funa@bio.keio.ac.jp>;

Last Modified: Wed, 13 Jul 2022 23:32:52 +0900

Copyright (c) 2013-2022 Funahashi Lab., Keio University.

## Introduction
This program simulates centrosome movement based on a mechanical model for two phenomena, centrosome centralization and spindle elongation, observed in the first cell division of fertilized eggs of C. elegans and sea urchins. In particular, the mechanical model in this program considers the force associated with the contact of microtubules extending from the centrosome with the cell membrane as the central force, and constructs a model in which the angle of microtubules is fixed (MTFixed) due to the weak coupling of these microtubules with the cell membrane, and a model in which the angle is variable (MTVariable) due to the strong coupling of microtubules with the cell membrane.
## All code and how to use
### Requirements
- make (confirmed to work with 3.81)
- ctags (confirmed to work with 5.8)
- gcc  (confirmed to work with 11.3.0)
- X11 (confirmed to work with 1.20.11)
- python3 (confirmed to work with 3.9.12)
  - py39-matplotlib (confirmed to work with 0.1.3_0)
  - py39-pandas (confirmed to work with 1.3.3_0)
- R (confirmed to work with 4.0.4)
- gnuplot (confirmed to work with 5.4)
 

### NCCCentration
This program performs simulation and analysis of the movement and rotation of the centrosome in the nucleus-centrosome complex centralization based on MT Fixed and MT Variable. Analysi\_MovementAngle compares the steady state between MT Fixed and MT Variable for the translation and rotation of the central body, and Analysis_SteadyState analyzes the steady state for a wide range of initial positions and rotation angles for the centrosome.

#### Analysis_MovementAngle
##### How to use
```sh
% cd NCCCentration/Analysis_MovementAngle/
% ./main.sh
```
##### Results
- ``NCCCentration/Analysis_MovementAngle/Result/Simulation/MTFixed/``
	
	Simulation results on centrosome movement and roation in nucleus-centrosome complex centralization based on MTFixed
- ``NCCCentration/Analysis_MovementAngle/Result/Simulation/MTVariable/``

	Simulation results on centrosome movement and roation in nucleus-centrosome complex centralization based on MTVariable	
- ``NCCCentration/Analysis_MovementAngle/Result/Analysis/PositionAngle/position.pdf``

	Comparison results of steady-state movement distances for MTFixed and MTVariable based simulations of central body movement.
- ``NCCCentration/Analysis_MovementAngle/Result/Analysis/PositionAngle/angle.pdf ``

	Comparison results of steady-state rotation angle for MTFixed and MTVariable based simulations of central body movement.
- `` NCCCentration/Analysis_MovementAngle/Result/Analysis/StatisticalTest/wilcoxon_ranksum_test.txt``

	Results of the Wilcoxon rank sum test for steady-state movement distance and rotation angle for simulations of centrosome movement based on MTFixed and MTVariable
	
#### Analysis_SteadyState
##### How to use
```sh
% cd NCCCentration/Analysis_SteadyState
% ./main.sh
```

##### Results
- ``NCCCentration/Analysis_SteadyState/Result/Simulation/``

	Extracted steady-state movement distance and rotation angle of the centrosome from simulation results based on MTFixed for various changes in the initial position and initial angle of the centrosome	
- ``NCCCentration/Analysis_SteadyState/Result/Analysis/rads15.pdf``

	Vector fields and steady-state visualization results for movement and rotation obtained from simulations at various initial positions and angles of the centrosome in an embryo with a long axis of 25 μm and a short axis of 15 μm.
- ``NCCCentration/Analysis_SteadyState/Result/Analysis/rads15_rads10.pdf``

	Vector fields and steady-state visualization results for movement and rotation obtained from simulations at various initial positions and angles of the centrosome in an embryo with a long axis of 25 μm and a short axis of 15 μm and of 10 μm.
- ``NCCCentration/Analysis_SteadyState/Result/Analysis/rads25.pdf``

	Vector fields and steady-state visualization results for movement and rotation obtained from simulations at various initial positions and angles of the centrosome in an embryo with a long axis of 25 μm and a short axis of 25 μm.

### SpindleElongation
This program simulates the movement of the nuclus-centrosome complex in spindle elongation based on MT Fixed and MT Variable and compares their steady states. The steady-state comparison is performed by changing the shape of the embryo in which the nucleus-centrosome complex exists in terms of aspect ratio, which is the ratio of the long axis to the short axis, over a wide range.
#### Analysis_ElongatedSpindle
##### How to use
```sh
% cd cd SpindleElongation/Analysis_ElongatedSpindle/
% ./main.sh
```

##### Results
- ``SpindleElongation/Analysis_ElongatedSpindle/Result/Simulation/``

	Simulation results for centrosome movement in spindle elongation based on MTFixed and MTVariable
- ``SpindleElongation/Analysis_ElongatedSpindle/Result/Analysis/Cel_sim_vs_exp.pdf``

	Steady-state positions obtained from simulations of centrosome movement during spindle elongation based on MTFixed and MTVariable in embryos of various aspect ratios for C. elegans compared to those obtained from experiments.
- ``SpindleElongation/Analysis_ElongatedSpindle/Result/Analysis/Su_sim_vs_exp.pdf``

	Steady-state positions obtained from simulations of centrosome movement during spindle elongation based on MTFixed and MTVariable in embryos of various aspect ratios for sea urchin compared to those obtained from experiments.	

####  Experiment_ElongatedSpindle
The directory stores measurements of the aspect ratio of each embryo and the position of the centrosome at convergence, based on images taken of spindle elongation for C. elegans and sea urchins.	

##### Results
- ``SpindleElongation/Experiment_ElongatedSpindle/exp_cel.csv``
	
	Measurements of the aspect ratio of each embryo and the position of the centrosome at convergence, based on images taken of spindle elongation relative to C. elegans.

- ``SpindleElongation/Experiment_ElongatedSpindle/exp_su.csv``

	Measurements of the aspect ratio of each embryo and the position of the centrosome at convergence, based on images taken of spindle elongation relative to sea urchin.
	
	
## Citation
Coming soon...
