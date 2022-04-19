# Turbine Blades Design - Computer Experiment Designs
Authors: Yuxuan Chen, Xige Huang, Gaojia Xu, Guanqi Zeng

## Background and Introduction

In contemporary society, power is needed everywhere in our daily life. Although not familiar to most of us, to
create power, the generator needs momentum from the engine that produces kinetic energy from some sources.
One type of combustion engine, the gas turbine, is a crucial component of the power generating process
across all industries, such as aircrafts, trains, and electricity plants. Connected to the generator, the gas
turbine drives the generator to function by pressurized gas: the air is drawn to the turbine from its suction
and heated by fuel source combustion; then, as the heated air expands, the motion is made by the turbine
and is connected to the generator to produce electricity. Transforming pressure energy to kinetic energy, the
turbine blades are the key element of a gas turbine. As it extracts energy from the high temperature and
high pressure gas, it is also the limiting part of a gas turbine, as it needs to withstand the severe environment.
To protect the turbine blades, an internal cooling system is incorporated in the turbine, extracting cooling air
from the compressor and passing through the airfoils to cool the blades.

While blade deformation to some extent is acceptable, an extreme deformation can largely reduce the length of
the engine life cycle, and the life of the blades themselves will decrease by half only by a 30 degree celsius off in
the blade temperature prediction. As a result, the design of the turbine blades requires careful consideration
on the choice of materials and the cooling schedule.

In this project, we aim to find an optimal design for the turbine blades to minimize blade stress and
deformation. Since physical experiments are expensive to prototype and run, computer experiments are
used instead to reduce costs. The deterministic black-box simulator makes use of a finite-element model
(FEM) solver which models the physical deformation of the blade during operation. There are six inputs in
the simulator, which will be elaborate on the next section. There are mainly two goals we need to achieve:
firstly, we want to train an surrogate model on stress that yields good prediction performance over the design
space that can be used for a variety of downstream tasks, including model calibration and system control;
secondly, we want to identify good blade material properties and cooling conditions which minimize maximum
blade stress over the turbine blade. More importantly, we also need to consider the fact that a maximum
displacement greater than d  = $1.3 × 10^3$ will result in immediate failure of the engine; therefore, we need to
take the black-box constraint into account while improving emulation accuracy and optimizing the black-box
simulator.

Another main challenge that need to be carefully taken care of in every stage of our design and computer
experiments is that running the simulator is fairly time-consuming. Therefore, given limiting computational
resources and also to make our algorithms and results reproducible, the project is limited to a maximum of
150 computer experiment runs.

## The repository
**Functions used in this project are in `functions.R`. Example of use can be found in `test_runs.rmd`.**
**Note that running the simulator on window and macos systems will yield slightly different results.**
