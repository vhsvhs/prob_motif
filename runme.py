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
    count_zero = 0
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

#def find_occurances(motif, urses):
#    count = 0
#    locs = []
#    ret_counts = []
#    mlen = motif.__len__()
#    for urs in urses:
#        for i in range(0, urs.__len__()):
#            if i + mlen < urs.__len__():
#                if urs[i:i+mlen].__contains__(motif):
#                    if locs.__len__() == 0:
#                        locs.append( i )
#                        count += 1
#                    elif locs[ locs.__len__()-1 ] != i:
#                        locs.append(i)
#                        #print "locs", locs[ locs.__len__()-2 ], i
#                        count += 1
#                    #else:
#                        #print "already seen it."
#        ret_counts.append( count )
#    return ret_counts
#    

#"""cumulative probability of a motif appearing,
#given the number of mutations."""
#def cf(mlib, urslen, nmutations, nsamples, stride):
#
#    msc = {} #msc[motif] = [] array of arrays
#    for m in mlib:
#        msc[m] = []
#    
#    for n in range(0, nsamples):
#        print "sample", n
#        for m in mlib:
#            msc[m].append( [] )
#            
#        """For each sample, make a random URS"""
#        urs = make_random_urs(urslen)
#        urses = []
#        for i in range(0, nmutations):
#            urs = mutate_urs(urs)
#            if stride%i == 0:
#                urses.append(urs)
#                for m in mlib:
#                    msc[m].append( = find_occurances( m, urses ) )
#    return msc

def cf2(mlib, urslen, nmutations, nsamples):
    msc = {} # key = motif, value = [] one for each sample, value = [] of cumm. count
    ms_locs = {} # key = motif, value = array of last-seen locations

    for m in mlib:
        msc[m] = []
        for s in range(0, nsamples):
            msc[m].append([])
            for n in range(0, nmuts):
                msc[m][s].append(0.0)
    
    for s in range(0, nsamples):
        print ". sample", s, "of", nsamples-1
        ms_locs = {}
        for m in mlib:
            ms_locs[m] = []

        urs = make_random_urs(urslen)
        
        for n in range(0, nmuts):
            #print ". mutant", n
            #print ". new urs", urs

            """First, update our info. about previously found copies of the motif."""
            for m in mlib:
                #print "ms_locs", m, ms_locs[m]
                for l in ms_locs[m]:
                    if False == urs[l:l+m.__len__()].__contains__(m):
                        ms_locs[m].remove( l )
                    
                if n > 0:
                    msc[m][s][n] = msc[m][s][n-1]
            
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
                                    msc[m][s][n] += 1
                                    #print "found match", i, m, urs[i:i+m.__len__()]
                                    
                                    
            urs = mutate_urs(urs)
    return msc

""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" Main... """

ap = ArgParser(sys.argv)

urslen = int( ap.getArg("--urslen") )
nmuts = int( ap.getArg("--nmutations") )
nsamples = int( ap.getArg("--nsamples") )
mlibpath = ap.getArg("--mlibpath")

pwms = {}
pwm_paths = {}

mlib = []
fin = open(mlibpath, "r")
for l in fin.readlines():
    if l.__len__() > 2:
        mlib.append(l.strip())
fin.close()
  
print "\n. OK, I'm using this motif library:"
print mlib

print "\n. Calculating MSC. . ."
msc = cf2(mlib, urslen, nmuts, nsamples)
#plot2(msc, mlibs[m], nmuts, nsamples)
print "\n. Calculating PDF. . ."
pvals = msc_to_cdf(msc, mlib, nmuts, nsamples)
#print pvals
#print "mlibs for ", m, "=", mlibs[m]
plotcdf(pvals, nmuts)


