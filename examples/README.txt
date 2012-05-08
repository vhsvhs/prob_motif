###################################################

Files in the example:

1. fordyce.MATA2
The position weight matrix for MATA2 from S. cerevisiae, experimentally derived by Polly Fordyce using the MITOMI system.

2. fordyce.MATA2.mlib
The motif library for MATA2, calculated using the python script "build_mlib.py" with cutoff = 5.2 bits.

3. run_fordyce_example.sh
A shell script to invoke the python program "mlib_to_cdf.py" for MATA2.  This script will run for 1000 mutations, with 100 replicates.

4. fordyce.MATA2.10000gen.pdf
This is the graphic output of running #3 (above) with 10,000 generations.

########################################################
