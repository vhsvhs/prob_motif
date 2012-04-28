"""
prob_motif

    Victor Hanson-Smith
    victorhansonsmith@gmail.com

What is the probability of a cis-regulatory motif occurring by chance?
"""
from argparser import *
import matplotlib.pyplot as plt
import random
from pwm_tools import *
from stattools import *

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

def msc_to_cdf(msc, mlib, nmuts, nsamples):
    pvals = [] # key = motif, value = Prob from samples that count @ nmut > 0
    for i in range(0, nmuts):
        for m in mlib:
            count_zero = 0
            for n in range(0, nsamples):
                if msc[m][n][i] == 0:
                    count_zero += 1
        this_pval = 1.0 - (count_zero * 1.0 / (nsamples + mlib.__len__() ) )
        pvals.append(this_pval)
    return pvals


def plot2(msc, mlib, nmuts, nsamples):
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
    x = [] # x axis
    for i in range(0, nmuts):
        x.append(i)
            
    plt.plot(x, pvals)
    plt.axis([0, nmuts, 0, 1.0])
    plt.legend()
    plt.show()

def find_occurances(motif, urses):
    count = 0
    locs = []
    ret_counts = []
    mlen = motif.__len__()
    for urs in urses:
        for i in range(0, urs.__len__()):
            if i + mlen < urs.__len__():
                if urs[i:i+mlen].__contains__(motif):
                    if locs.__len__() == 0:
                        locs.append( i )
                        count += 1
                    elif locs[ locs.__len__()-1 ] != i:
                        locs.append(i)
                        #print "locs", locs[ locs.__len__()-2 ], i
                        count += 1
                    #else:
                        #print "already seen it."
        ret_counts.append( count )
    return ret_counts
    

"""cumulative probability of a motif appearing,
given the number of mutations."""
def cf(mlib, urslen, nmutations, nsamples):

    msc = {} #msc[motif] = [] array of arrays
    for m in mlib:
        msc[m] = []
    
    for n in range(0, nsamples):
        print "sample", n
        for m in mlib:
            msc[m].append( [] )
            
        """For each sample, make a random URS"""
        urs = make_random_urs(urslen)
        urses = []
        for i in range(0, nmutations):
            urs = mutate_urs(urs)
            urses.append(urs)
        for m in mlib:
            msc[m][n] = find_occurances( m, urses )
    return msc

""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" Main... """

ap = ArgParser(sys.argv)

urslen = int( ap.getArg("--urslen") )
nmuts = int( ap.getArg("--nmutations") )
nsamples = int( ap.getArg("--nsamples") )
pwm_dir = ap.getOptionalArg("--pwmdir")
pwmpath = ap.getOptionalArg("--pwmpath")
if pwm_dir == False and pwmpath == False:
    print "Ooops, you forgot to specify --pwmdir, or --pwmpath.  Pick one."
cutoffpath = ap.getOptionalArg("--cutoffpath")
cutoffval = float(ap.getOptionalArg("--cutoff"))
if cutoffpath == False and cutoffval == False:
    print "Ooops, you forgot to specify --cutoffpath, or --cutoff.  Pick one."

pwms = {}
pwm_paths = {}

if pwm_dir:
    for f in os.listdir(pwm_dir):
        tokens = f.split(".")
        pwmid = tokens[1]
        pwm_paths[pwmid] = pwm_dir + "/" + f
        print "Reading PWM from", pwm_paths[pwmid]
        pwms[pwmid] = read_pwm_from_file(pwm_paths[pwmid]  )
else:
    tokens = pwmpath.split(".")
    pwmid = tokens[1]
    pwm_paths[pwmid] = pwmpath
    pwms[pwmid] = read_pwm_from_file(pwm_paths[pwmid]  )

mcutoffs = {}
if cutoffpath:
    mcutoffs = read_cutoffs_from_scertf(cutoffpath)
    for m in pwm_paths:
        print "GENE: ", m, "cutoff = ", mcutoffs[m]
else:
    for m in pwm_paths:
        mcutoffs[m] = cutoffval
        print "GENE: ", m, "cutoff = ", mcutoffs[m]

mlibs = {}
for m in pwms:
    print "Calculting mlibs for", m
    print pwms[m]
    mlibs[m] = mlib_from_pwm(pwms[m], mcutoffs[m])
    print "Calculating msc for", m
    msc = cf(mlibs[m], urslen, nmuts, nsamples)
    plot2(msc, mlibs[m], nmuts, nsamples)
    print "Calculating PDF for", m
    
    pvals = msc_to_cdf(msc, mlibs[m], nmuts, nsamples)
    #print pvals
    #print "mlibs for ", m, "=", mlibs[m]
    plotcdf(pvals, nmuts)


