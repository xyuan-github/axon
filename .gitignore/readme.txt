Algorithm for extraction of axonal arbors

The files included in this folder can be used to test the algorithm for the extraction of axonal arbors from neuronal footprints of electrical activity. An example dataset includes data of three exemplary neurons.

The files included in this folder are:
    - Main.m %  main Matlab file for testing 
    - axon_velocity_auto.m %  main extraction algorithm
    - data.mat %  example dataset, including the footprints of three exemplary neurons
    - ntk.mat %  array information needed for the algorithm


Please simply run Main.m in Matlab to test the algorithm according to the following steps:
1. Load example data set (data.mat % ) and array information (ntk.mat %).
2. Set parameters for extraction, including searching radius, minimum number of electrodes, initial frames, etc. 
3. Use axon_velocity_auto function to extract the axonal arbors.

After running the scripts in Main.m, one should automatically obtain the output plots of the algorithm, including the amplitude map, delay map, extracted axonal segments and propagation velocity (same as can be seen in Fig. 4a).

The software should be run with Matlab R2017a or newer versions. There is no special requirement for the operating system.
