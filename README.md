# MTSim: Microtubule Simulator
This code was used in the research for the paper titled "Switching from weak to strong cortical attachment of microtubules accounts for the transition from nuclear centration to spindle elongation in metazoans".

## Screenshot
![screenshot](https://github.com/funalab/MTSim/blob/images/screenshot.png)

## Authors
- Shohei Tada: <tadashohei.24@gmail.com>;
- Yoshitaka Yamazaki: <yoshi65.ymzk@gmail.com>;
- Takahiro G Yamada: <yamada@fun.bio.keio.ac.jp>;
- Akira Funahashi: <funa@bio.keio.ac.jp>;
- Akatsuki Kimura: <akkimura@nig.ac.jp>;

Last Modified: Thu, 21 Jul 2022 20:09:40 +0900

Copyright (c) 2013-2022 Funahashi Lab., Keio University.

## Introduction
This program simulates centrosome movement based on a mechanical model for two phenomena, centrosome centralization and spindle elongation, observed in the first cell division of fertilized eggs of C. elegans and sea urchins. In particular, the mechanical model in this program considers the force associated with the contact of microtubules extending from the centrosome with the cell membrane as the central force, and constructs a model in which the angle of microtubules is fixed (MTFixed) due to the weak coupling of these microtubules with the cell membrane, and a model in which the angle is variable (MTVariable) due to the strong coupling of microtubules with the cell membrane.
## All code and how to use
### Requirements
- When running on a local machine
  - make (confirmed to work with 3.81)
  - ctags (confirmed to work with 5.8)
  - gcc  (confirmed to work with 10.2.1)
  - X11 (confirmed to work with 1.20.11)
  - xwd (confirmed to work with 1.0.7)
  - ImageMagick (confirmed to work with 6.9.11-60)
  - python3 (confirmed to work with 3.9.2)
    - py39-matplotlib (confirmed to work with 3.3.4)
    - py39-pandas (confirmed to work with 1.1.5)
  - R (confirmed to work with 4.0.4)
  - gnuplot (confirmed to work with 5.4)
  - zsh (confirmed to work with 5.8)
  - cm-super (confirmed to work with 0.3.4-15)
  - dvipng (confirmed to work with 1.15)
  - texlive-latex-base (confirmed to work with 2020.20210202-3)
  - texlive-latex-extra (confirmed to work with 2020.20210202-3)
  - texlive-fonts-recommended (confirmed to work with 2020.20210202-3)
  
- When running on docker
  - Please follow the instructions in this [README.md](./docker/)

### NCCCentration
This program performs simulation and analysis of the movement and rotation of the centrosome in the nucleus-centrosome complex centralization based on MT Fixed and MT Variable. MovementAngle compares the steady state between MT Fixed and MT Variable for the translation and rotation of the central body, and SteadyState analyzes the steady state for a wide range of initial positions and rotation angles for the centrosome.

#### MovementAngle
##### How to use
```sh
% cd NCCCentration/MovementAngle/
% ./main.sh
```
##### Results
- ``NCCCentration/MovementAngle/Result/Simulation/MTFixed/``
	
	Simulation results on centrosome movement and roation in nucleus-centrosome complex centralization based on MTFixed
- ``NCCCentration/MovementAngle/Result/Simulation/MTVariable/``

	Simulation results on centrosome movement and roation in nucleus-centrosome complex centralization based on MTVariable	
- ``NCCCentration/MovementAngle/Result/Analysis/PositionAngle/position.pdf``

	Comparison results of steady-state movement distances for MTFixed and MTVariable based simulations of central body movement.
- ``NCCCentration/MovementAngle/Result/Analysis/PositionAngle/angle.pdf ``

	Comparison results of steady-state rotation angle for MTFixed and MTVariable based simulations of central body movement.
- `` NCCCentration/MovementAngle/Result/Analysis/StatisticalTest/wilcoxon_ranksum_test.txt``

	Results of the Wilcoxon rank sum test for steady-state movement distance and rotation angle for simulations of centrosome movement based on MTFixed and MTVariable
	
#### SteadyState
##### How to use
```sh
% cd NCCCentration/SteadyState
% ./main.sh
```

##### Results
- ``NCCCentration/SteadyState/Result/Simulation/``

	Extracted steady-state movement distance and rotation angle of the centrosome from simulation results based on MTFixed for various changes in the initial position and initial angle of the centrosome	
- ``NCCCentration/SteadyState/Result/Analysis/rads15.pdf``

	Vector fields and steady-state visualization results for movement and rotation obtained from simulations at various initial positions and angles of the centrosome in an embryo with a long axis of 25 μm and a short axis of 15 μm.
- ``NCCCentration/SteadyState/Result/Analysis/rads15_rads10.pdf``

	Vector fields and steady-state visualization results for movement and rotation obtained from simulations at various initial positions and angles of the centrosome in an embryo with a long axis of 25 μm and a short axis of 15 μm and of 10 μm.
- ``NCCCentration/SteadyState/Result/Analysis/rads25.pdf``

	Vector fields and steady-state visualization results for movement and rotation obtained from simulations at various initial positions and angles of the centrosome in an embryo with a long axis of 25 μm and a short axis of 25 μm.

### SpindleElongation
This program simulates the movement of the nuclus-centrosome complex in spindle elongation based on MT Fixed and MT Variable and compares their steady states. The steady-state comparison is performed by changing the shape of the embryo in which the nucleus-centrosome complex exists in terms of aspect ratio, which is the ratio of the long axis to the short axis, over a wide range.
##### How to use
```sh
% cd SpindleElongation/
% ./main.sh
```

##### Results
- ``SpindleElongation/Result/Simulation/``

	Simulation results for centrosome movement in spindle elongation based on MTFixed and MTVariable
- ``SpindleElongation/Result/Analysis/Cel_sim_vs_exp.pdf``

	Steady-state positions obtained from simulations of centrosome movement during spindle elongation based on MTFixed and MTVariable in embryos of various aspect ratios for C. elegans compared to those obtained from experiments.
- ``SpindleElongation/Result/Analysis/Su_sim_vs_exp.pdf``

	Steady-state positions obtained from simulations of centrosome movement during spindle elongation based on MTFixed and MTVariable in embryos of various aspect ratios for sea urchin compared to those obtained from experiments.	

####  Experiment
The directory stores measurements of the aspect ratio of each embryo and the position of the centrosome at convergence, based on images taken of spindle elongation for C. elegans and sea urchins.	

##### Results
- ``SpindleElongation/Experiment/exp_cel.csv``
	
	Measurements of the aspect ratio of each embryo and the position of the centrosome at convergence, based on images taken of spindle elongation relative to C. elegans.

- ``SpindleElongation/Experiment/exp_su.csv``

	Measurements of the aspect ratio of each embryo and the position of the centrosome at convergence, based on images taken of spindle elongation relative to sea urchin.
	
	
## Citation
Coming soon...
