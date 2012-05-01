from argparser import *
import matplotlib.pyplot as plt
import random
from pwm_tools import *
from stattools import *

ap = ArgParser(sys.argv)

pwmpath = ap.getOptionalArg("--pwmpath") 
pwmdir = ap.getOptionalArg("--pwmdir") 
mlib_out_dir = ap.getArg("--mlib_out_dir")
if False == os.path.exists(mlib_out_dir):
    os.system("mkdir " + mlib_out_dir)

cutoffpath = ap.getOptionalArg("--cutoffpath")
cutoff = float( ap.getOptionalArg("--cutoff"))

if pwmpath == False and pwmdir == False:
    print "Ooops. You need to specify either --pwmpath or --pwmdir."
    exit(1)
if cutoffpath == False and cutoff == False:
    print "Oops. You need to specify either --cutoff or --cutoffpath"
    exit(1)

if cutoff == False:
    cutoffs = read_cutoffs_from_scertf(cutoffpath)
print "\n. OK, I found these cutoff values. . .\n"
print cutoffs

def is_this_pwm(path):
    fin = open(path, "r")
    lines = fin.readlines()
    fin.close()
    count_lines_with_content = 0
    for l in lines:
        if False == l.startswith("#"):
            if l.__len__() > 1:
                count_lines_with_content += 1
    if count_lines_with_content != 4:
        return False
    return True

pwms = {}
for f in os.listdir(pwmdir):
    print "\n. Analyzing PWM", f, ". . ."
    if False == is_this_pwm(pwmdir + "/" + f):
        continue
    pwm = read_pwm_from_file(pwmdir + "/" + f)
    mid = f.split(".")[1]
    pwms[mid] = pwm
    if pwm.__len__() > 8:
        print "I'm skipping PWM", f
        continue
    print "\n. I found the following PWM"
    sites = pwm.keys()
    sites.sort()
    for s in sites:
        print pwm[s]
    print "\n. Searching for motifs with cutoff >", cutoffs[mid]
    (t, mlib) = pwm_to_tree(pwm, cutoffs[mid])
    fout = open(mlib_out_dir + "/" + f + ".mlib", "w")
    for m in mlib:
        fout.write(m + "\n")
    fout.close()