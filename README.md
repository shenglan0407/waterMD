# waterMD

Analysis code for water box MD simulations

Using compute_correlator.py

This is the important script that computes C4. To run it use the following command and pass on the required parameters:
"python compute_correlator.py -i <runname> -o <outputfile> -q <q_inv> -p <nphi> -s <fstart> -e <fend>"]

Parameters:

-i: name of the simulation run to use, e.g. 'run1'

-o: the name output file to store the computed correlator, e.g. "corr_run1.csv". I also try to include other parameters in the name. The outputfile is by default saved in /computed_results. If the file already exist, the program will refuse to overwrite it and propt the user for a new outputfile name.

-q: defines the magnitude of the q vectors. At this point, the program only computes of qs with equal magnitude. The magnitude is defined as 2*pi/(input) and input is assumed to have the unit nm. I find defining this parameter in the unit of nm more intuitive for me. But I am also considering changing it to angstrom-inverse.

-p: number of phi for which to compute. The phi are spaced linearly over the range of 0 to pi.

-s: frame in simulation to start computing

-e: fram in simulation to end computing

