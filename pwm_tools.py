import itertools
import os, sys

DNA = ["A", "C", "G", "T"]

def read_pwm_from_file(f):
    """Using the ScerTF file format."""
    pwm = {}
    try:
        fin = open(f, "r")
        for l in fin.readlines():
            if l.__len__() > 3:
                tokens = l.split()
                if tokens.__len__() <= 1:
                    continue
                if tokens[1] != "|":
                    continue
                state = tokens[0]
                site = 0
                for val in tokens[2:]:
                    if site not in pwm.keys():
                        pwm[site] = {}
                    pwm[site][state] = float(val)
                    site += 1
        fin.close()
        if pwm.keys().__len__() == 0:
            return False
        return pwm
    except IOError:
        return False

def get_information_content( p ):
    sum = 0.0
    for site in p:
        for state in p[site]:
            if p[site][state] > 0.0:
                sum += p[site][state]
    return sum

def get_mochipwms(path):
    """Parse a MochiView-formatted list of PWMs."""
    pwms = {}
    fin = open(path, "r")
    this_pwm = {}
    this_name = None
    count_site = 0
    for l in fin.xreadlines():
        if l.startswith("#"):
            continue
        if l.__contains__("NAME"):
            tokens = l.split()
            #if tokens[0].__contains__("=[JP]"):
            #    continue
            if tokens[2].startswith("TYPE"):
                this_name = tokens[1]
            else:
                this_name = tokens[1] + "." + tokens[2]
                if pwms.keys().__contains__(this_name):
                    print "I found a duplicate for", this_name
            pwms[ this_name ] = {}
            count_site = 0
        elif this_name != None and l.split().__len__() == 4:
            tokens = l.split()
            states = ['A', 'C', 'G', 'T']
            pwms[ this_name ][count_site] = {}
            for i in range(0, 4):
                #print tokens
                #print this_name, count_site, states[i]
                pwms[ this_name ][count_site][ states[i] ] = float( tokens[i] )
            count_site += 1  
    fin.close()
    return pwms

def read_cutoffs_from_scertf(f):
    """Reads the recommended cutoff file that you've pre-downloaded from the ScerTF database."""
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

def pwm_to_tree(pwm, cutoff):
    return pwm_to_tree_recur(pwm, cutoff, 0, 0.0, "")

def pwm_to_tree_recur(pwm, cutoff, site, parent_cost, parent_string):
    """Recursively builds a syntax tree from the PWM, with bit costs at each node."""
    if site == pwm.__len__() - 1:
        ret = {}
        mlibs = []
        for state in pwm[site].keys():
            motif = parent_string + state
            cost = parent_cost + pwm[site][state]

            if cost > cutoff:
                mlibs.append(motif)
                print motif, cost # this explodes runtime!!
            ret[state] = (cost, None)
        return (ret, mlibs)
    else:
        ret = {}
        mlibs = []
        for state in pwm[site].keys():
            my_cost = parent_cost + pwm[site][state]
            my_string = parent_string + state
            if my_cost > cutoff:
                mlibs.append(my_string)
                print my_string, my_cost # this explodes runtime!!
            recur_results = pwm_to_tree_recur(pwm, cutoff, site + 1, my_cost,my_string)
            child_pointer = recur_results[0]
            child_mlibs = recur_results[1]
            mlibs = mlibs + child_mlibs
            ret[state] = (my_cost, child_pointer)
        return (ret, mlibs)

def mlib_from_pwm(pwm, cutoff):
    """Builds a syntax tree, and then finds all subtrees that satisfy the cutoff."""
    mlib = []
    t = pwm_to_tree(pwm)
    for mlen in range(1, pwm.__len__()):
        x = itertools.product(["A", "C", "G", "T"], ["A", "C", "G", "T"], repeat = mlen)
        for i in x:
            """Turn i into string."""
            #print i, i.__len__(), mlen # this is useful for debugging, but it explodes the computation time!
            m = ""
            for l in i:
                m += i[l]
            """Count bits."""
            bits = 0.0
            for j in range(0, m.__len__()):
                bits += pwm[j][ m[j] ]
            if bits > cutoff:
                #print "*"
                if False == mlib.__contains__(m):
                    mlib.append(m)
                    print m, bits
    return mlib