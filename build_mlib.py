#############################################################
#
# build_mlib.py
#
# Victor Hanson-Smith
# victorhansonsmith@gmail.com
# 2012
#
# INPUT: a position weight matrix (PWM) and a cutoff (in bits)
# OUTPUT: a motif library, i.e. a list of all nucleotide motifs that satisfy
#    the PWM's cutoff.
#
# USAGE:
# $> python build_mlib.py --pwmpath <filepath> --cutoff <N> --mlib_out_dir <output path>
# 
#    . . . where <filepath> is the path to a PWM.  See the examples folder for formatting criterion.
#    . . . and <N> is the cutoff in bits.
#    . . . and <output path> is a folder into which the motif library will be written.
#
# ALTERNATIVE USAGE, for batch processing many PWMs.
# $> python build_mlib.py --pwmdir <directory> --cutoffpath <path> --mlib_out_dir <output path>
#
#    . . . where <directory> is the path to a folder containing multiple PWM files.
#    . . . and <path> is the path to a file listing cutoff values for those PWMs.
#    . . . and <output path> is a folder into which the motif library will be written.
#
# IMPORTANT NOTES:
# 1. The filenames of PWMs should be formatted as XXXXX.GENE, where XXXXX can be any identifier
#     of your choosing, and GENE should be the name of the transcription factor.
#
# 2. The alternative usage is designed to work with the PWMs from the ScerTF database.
#    I suggest you look at the example directory to learn about file formatting.
#
#####################################################################


from argparser import *
import matplotlib.pyplot as plt
import random
from pwm_tools import *
from stattools import *

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

def pwmpath_to_tfid(path):
    """Get the transcription factor name from the filepath.  See IMPORTANT NOTES in
    the comments at the top of this file."""
    tokens = path.split("/")
    return tokens[ tokens.__len__() - 1].split(".")[1]

############################################
#
# main:
#
ap = ArgParser(sys.argv)
pwmpath = ap.getOptionalArg("--pwmpath") 
pwmdir = ap.getOptionalArg("--pwmdir") 
mlib_out_dir = ap.getArg("--mlib_out_dir")
cutoffpath = ap.getOptionalArg("--cutoffpath")
cutoff = float( ap.getOptionalArg("--cutoff"))

# Sanity check. . .
if pwmpath == False and pwmdir == False:
    print "Ooops. You need to specify either --pwmpath or --pwmdir."
    exit(1)
if pwmpath != False and pwmdir != False:
    print "Ooops. You specified both --pwmpath and --pwmdir.  Please use only one."
    exit(1)
if cutoffpath == False and cutoff == False:
    print "Oops. You need to specify either --cutoff or --cutoffpath"
    exit(1)
if cutoffpath != False and cutoff != False:
    print "Oops. You specified both --cutoff and --cutoffpath.  Please use only one."
    exit(1)
if cutoffpath != False and pwmdir == False:
    print "Oops. You specified --cutoffpath, but you didn't specify --pwmdir."
    exit(1)
if cutoffpath == False and pwmdir != False:
    print "Oops. You specified --pwmdir, but you didn't specify --cutoffpath."
    exit(1)
if cutoff == False and pwmpath != False:
    print "Oops. You specified --pwmpath, but you didn't specify --cutoff."
    exit(1)
if cutoff != False and pwmpath == False:
    print "Oops. You specified --cutoff, but you didn't specify --pwmpath."
    exit(1)


# Setup the output folder. . .
if False == os.path.exists(mlib_out_dir):
    os.system("mkdir " + mlib_out_dir)

# Get the cutoff(s). . .
cutoffs = {}
if cutoff == False:
    cutoffs = read_cutoffs_from_scertf(cutoffpath)
    print "\n. OK, I found these cutoff values. . .\n"
    print cutoffs
else:
    print "\n. OK, I'm using cutoff = ", cutoff
    tfid = pwmpath_to_tfid(pwmpath)
    cutoffs[tfid] = cutoff

# Get the PWM(s). . .
pwm_files = []
if pwmpath != False:
    pwmdir = ""
    pwm_files.append( pwmpath )
elif pwmdir != False:
    pwmdir += "/"
    for f in os.listdir( pwmdir ):
        pwm_files.append( f )
pwms = {} # key = filename, value = PWM hash (see the methods in pwm_tools.py)
for f in pwm_files:
    print "\n. I'm analyzing the file:", f, ". . ."
    if False == is_this_pwm(pwmdir + "" + f):
        q = pwmdir + "" + f
        print "\n. Hmmm, I don't think this file contains a PWM:",  q
        print "---> I'm skipping this file."
        continue
    pwm = read_pwm_from_file(pwmdir + "" + f)
    
    tfid = pwmpath_to_tfid(f)
    
    pwms[tfid] = pwm
    
    if pwm.__len__() > 8:
        print "I'm skipping PWM", f
        continue
    
    print "\n. I found the following PWM:"
    sites = pwm.keys()
    sites.sort()
    for s in sites:
        print pwm[s]
        
    print "\n. Searching for motifs with cutoff >", cutoffs[tfid]
    (t, mlib) = pwm_to_tree(pwm, cutoffs[tfid]) # see pwm_tools.py
    
    # Finally, write the output. . .
    output_path = f + ".mlib"
    tokens = output_path.split("/")
    output_path = tokens[ tokens.__len__() - 1]
    output_path = mlib_out_dir + "/" + output_path
    print "\n. I'm writing the motif library to", output_path, "\n"
    fout = open(output_path, "w")
    for m in mlib:
        fout.write(m + "\n")
    fout.close()