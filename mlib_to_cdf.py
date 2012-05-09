#################################################
# 
# mlib_to_cdf.py
#
# Victor Hanson-Smith
# victorhansonsmith@gmail.com
# 2012
#
# What is the probability of a cis-regulatory motif occurring by chance?
#
# USAGE:
# $> python mlib_to_cdf.py --urslen L --nmutations M --nsamples N --mlibpath <motif library path> --runid ID
#
#    . . . where L is the length (in nucleotides) of the upstream regulatory sequence.
#    . . . where M is the number of mutations that will be introduced -- one by one -- to the URS.
#    . . . where N is the number of replicates that will be sampled.
#    . . . where <motif library path> is a text file containing a list of motifs, one motif per line.
#            See the script build_mlib.py to generate a motif library from a position weight matrix.
#    . . . where ID is the unique run ID for this execution.
#
# OPTIONAL ARGUMENTS:
# --stride S
#    . . . where S is the stride size in mutations.  The occurance of motifs will be counted
#            only every S mutations.  This is useful for reducing computation time when the number
#            of mutations is large (for example, 10^6 or larger), but if S is too large,
#            then there will be a significant reduction in the accuracy of the cummulative
#            distribution.
# --loadsaved <path>
#    . . . where <path> points to a Python pickled file written by this script.  If --loadsaved
#            is enabled, then this script will skip re-computing the cummulative distribution
#            and will directly plot the data in the pickled file.
# --makeplots True/False
#    . . . by default, this is set to False.  The cummulative distribution will be plotted using
#            matplotlib if this is set to True.
#
# IMPORTANT NOTES:
# 1. This script requires non-standard Python libraries: Mpi4py, and matplotlib
#    You can circumvent the need for matplotlib by specifying "--makeplots False" at the command-line.
#
# 2. The accuracy of your results will suffer if you use a small value for N (nsamples).
#
###########################################################

from argparser import *
import pickle
from pwm_tools import *
import random
from stattools import *
from time_estimate import *

def make_random_urs(len):
    urs = ""
    for i in range(0, len):
        urs += random.sample(["A", "C", "T", "G"], 1)[0]
    return urs

def mutate_urs(urs):
    newurs = ""
    i = random.randint(1, urs.__len__()-1)
    for j in range(0, i):
        newurs += urs[j].__str__() 
    alph = ["A", "C", "T", "G"]
    alph.remove(urs[i])
    newurs += random.sample(alph, 1)[0]
    for j in range(i+1, urs.__len__()):
        newurs += urs[j].__str__()
    return newurs

def summarize_msc(msc, mlib, nmuts, nsamples):
    m_mean_c = {} # m_mean_c[motif] = [] array of mean values for each nmutation
    m_sem_c = {}
        
    for m in mlib:
        m_mean_c[m] = []
        m_sem_c[m] = []
    
    for n in range(0, nsamples):
        for m in mlib:
            m_mean_c[m].append( mean(msc[m][n]) )
            m_sem_c[m].append( sem(msc[m][n]) )
    
    mean_c = []
    sem_c = []
    for n in range(0, nsamples):
        this_n_set = []
        for m in mlib:
            for i in range(0, nmuts):
                this_n_set.append( msc[m][n][i] )
        mean_c.append( mean( this_n_set ) )
        sem_c.append( sem( this_n_set ) )
    
    return [m_mean_c, m_sem_c, mean_c, sem_c]

def sc_to_cdf(sc, nmuts, nsamples, mliblen):
    count_zero = 0
    pvals = [] # key = motif, value = Prob from samples that count @ nmut > 0
    for i in range(0, nmuts):
        count_zero = 0
        for n in range(0, nsamples):
            if sc[n][i] == 0:
                count_zero += 1
        this_pval = 1.0 - (count_zero * 1.0 / (nsamples) )
        pvals.append(this_pval)
        print "mut", i, count_zero
    return pvals


def plot2(msc, mlib, nmuts, nsamples):
    import matplotlib.pyplot as plt
    x = [] # x axis
    for i in range(0, nmuts):
        x.append(i)
    
    """Find max"""
    maxy = 0.0
    for m in mlib:
        for n in range(0, nsamples):
            for i in range(0, nmuts):
                #print "msc", msc[m][n][i]
                if msc[m][n][i] > maxy:
                    maxy = msc[m][n][i]
    
    for m in mlib:
        for n in range(0, nsamples):
            plt.plot(x, msc[m][n] )
    plt.axis([0, nmuts, 0, maxy])
    plt.legend()
    plt.show()
    
def plotcdf(pvals, nmuts):
    import matplotlib.pyplot as plt
    x = [] # x axis
    for i in range(0, nmuts):
        x.append(i)
            
    plt.plot(x, pvals)
    plt.axis([0, nmuts, 0, 1.0])
    plt.legend()
    plt.show()


def cf2(mlib, urslen, nmuts, nsamples, nmut_stride, rank):
    """This method computes and returns sc, which is an array.
    sc[sample number][number of mutations] = the cummulative count of the occurances
        of a satisfactory motif in the URS up to this mutation."""
    #rank = comm.Get_rank()
    
    jobs = []
    sc = [] # key = motif, value = [], one for each sample, value = [] of cumm. count
    ms_locs = {} # key = motif, value = array of last-seen locations

    """Build data structures. . ."""
    count_mem_bytes = 0.0
    for s in range(0, nsamples):
        jobs.append(s)
        sc.append([])
        for n in range(0, nmuts):
            sc[s].append(0.0)
    

    starttime = time.time()
    
    for s in range(0, nsamples):    
        if False == is_my_item(s, jobs, rank):
            continue
        
        sample_start_time = time.time()
        
        if comm.Get_size() == 1:
            eta_str = "calculating"
            if s > 0:
                eta_t = get_time_remainaing(time.time() - starttime, s, nsamples-1 )
                eta_str = eta_t.__str__()
            print ". sample", s, "of", nsamples-1, "from node", rank, " [ est. sec. remaining = ", eta_str, "]"
        
        ms_locs = {}
        for m in mlib:
            ms_locs[m] = []

        urs = make_random_urs(urslen)
        
        for n in range(0, nmuts):
            #print ". mutant", n
            #print ". new urs", urs
            urs = mutate_urs(urs)
            
            if n%nmut_stride != 0:
                if n > 0:
                    sc[s][n] = sc[s][n-1]                    
                continue

            """First, update our info. about previously found copies of the motif."""
            for m in mlib:
                #print "ms_locs", m, ms_locs[m]
                for l in ms_locs[m]:
                    if False == urs[l:l+m.__len__()].__contains__(m):
                        ms_locs[m].remove( l )
                    
            if n > 0:
                sc[s][n] = sc[s][n-1]
            
            """ Next, count unique occurances of motif in urs """
            for i in range(0, urs.__len__()):
                for m in mlib:
                    """If m fits on the sequence..."""
                    if i + m.__len__() < urs.__len__():
                        """Does motif exist in seq?"""
                        if urs[i:i+m.__len__()].__contains__(m):
                            #print "maybe match", m, urs[i:i+m.__len__()]
                            """But we haven't seen this motif in previous mutants?"""
                            if False == ms_locs[m].__contains__(i):
                                    ms_locs[m].append(i)
                                    sc[s][n] += 1
                                    #print "found match", i, m, urs[i:i+m.__len__()]
    
        sample_elapsed_time = time.time() - sample_start_time
        print "\t-> sample", s, "of", nsamples-1, "from node", rank, "took %.3f"%sample_elapsed_time, "seconds."
        
        
    """Return our copy of sc to the master...."""
    if rank > 0:
        sc_p = pickle.dumps(sc)
        comm.send(sc_p, dest=0, tag=11)
    
    if rank == 0:
        """Wait for data from slaves."""
        for slave in range(1, comm.Get_size()):
            sc_p = comm.recv(source=slave, tag=11)
            slave_sc = pickle.loads(sc_p)
            for s in range(0, nsamples):
                if is_my_item(s, jobs, slave):
                    sc[s] = slave_sc[s]
    return sc

def read_mlib(mlibpath):
    mlib = []
    fin = open(mlibpath, "r")
    for l in fin.readlines():
        if l.__len__() > 2:
            mlib.append(l.strip())
    fin.close()
    return mlib

def post_cf2(sc, nmuts, nsamples, mliblen):
    print "\n. Calculating PDF. . ."
    pvals = sc_to_cdf(sc, nmuts, nsamples, mliblen)

    if makeplots == True:
        print "plotting cdf"
        plotcdf(pvals, nmuts)

#################################################
#
# main starts here. . .
#


# Read command-line arguments. . .
ap = ArgParser(sys.argv)

use_mpi = ap.getOptionalArg("--use_mpi")
if use_mpi != False:
    from mpilibs import *
    rank = comm.Get_rank()
    mpi_check()
else:
    rank = 0

urslen = int( ap.getArg("--urslen") )
nmuts = int( ap.getArg("--nmutations") )
nsamples = int( ap.getArg("--nsamples") )
stride = ap.getOptionalArg("--stride")
if stride == False:
    stride = 1
else:
    stride = int(stride)
mlibpath = ap.getArg("--mlibpath")
mlib = read_mlib(mlibpath)
runid = ap.getArg("--runid")
makeplots = ap.getOptionalArg("--makeplots")
if makeplots == "True":
    makeplots = True

loadsaved = ap.getOptionalArg("--loadsaved")
if loadsaved and rank == 0:
    print "\n. Loading saved data from", loadsaved
    fin = open(loadsaved, "r")
    sc = pickle.load(fin)
    fin.close()
elif loadsaved == False:    
    if rank == 0: 
        print "\n. OK, I'm using this motif library:"
        print mlib
        print "\n. Calculating the CDF. . ."
        
    sc = cf2(mlib, urslen, nmuts, nsamples, stride, rank)
    comm.Barrier()
    
    if rank == 0:
        """Save sc."""
        sc_p = pickle.dumps(sc)
        fout = open(runid + ".sc.pickle", "w")
        fout.write(sc_p)
        fout.close()
        
post_cf2(sc, nmuts, nsamples, mlib.__len__())


