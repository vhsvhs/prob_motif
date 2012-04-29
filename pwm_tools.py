import itertools
import os, sys

DNA = ["A", "C", "T", "G"]

def read_pwm_from_file(f):
    """Using the ScerTF file format."""
    pwm = {}
    fin = open(f, "r")
    for l in fin.readlines():
        if l.__len__() > 1:
            tokens = l.split()
            state = tokens[0]
            site = 0
            for val in tokens[2:]:
                if site not in pwm.keys():
                    pwm[site] = {}
                pwm[site][state] = float(val)
                site += 1
    fin.close()
    return pwm

def read_cutoffs_from_scertf(f):
    fin = open(f, "r")
    mcutoff = {}
    for l in fin.readlines():
        if l.__len__() > 2:
            tokens = l.split()
            if False == tokens[0].__contains__("#"):
                mid = tokens[0]
                pwm_cutoff = float(tokens[4])
                mcutoff[mid] = pwm_cutoff
    return mcutoff
    fin.close()

def mlib_from_pwm(pwm, cutoff):
    mlib = []
    
    for mlen in range(0, pwm.__len__()):
        #x = itertools.combinations_with_replacement(["A", "C", "G", "T"], mlen)
        x = itertools.product(["A", "C", "G", "T"], ["A", "C", "G", "T"], repeat = mlen)
        for i in x:
            """Turn i into string."""
            print i
            m = ""
            for l in range(0, mlen):
                m += i[l]
            """Count bits."""
            bits = 0.0
            for j in range(0, m.__len__()):
                bits += pwm[j][ m[j] ]
            if bits > cutoff:
                print m, bits
                #print "*"
                mlib.append(m)
    return mlib