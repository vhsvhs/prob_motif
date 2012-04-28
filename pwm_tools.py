import os, sys

DNA = ["A", "C", "T", "G"]

BASE2 = "01"
BASE4 = "0123"
BASE10 = "0123456789"
BASE16 = "0123456789ABCDEF"
BASE62 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz"

def baseconvert(number,fromdigits,todigits):
    """See http://code.activestate.com/recipes/111286-numeric-base-converter-that-accepts-arbitrary-digi/"""
    """ converts a "number" between two bases of arbitrary digits

    The input number is assumed to be a string of digits from the
    fromdigits string (which is in order of smallest to largest
    digit). The return value is a string of elements from todigits
    (ordered in the same way). The input and output bases are
    determined from the lengths of the digit strings. Negative 
    signs are passed through.

    decimal to binary
    >>> baseconvert(555,BASE10,BASE2)
    '1000101011'

    binary to decimal
    >>> baseconvert('1000101011',BASE2,BASE10)
    '555'

    integer interpreted as binary and converted to decimal (!)
    >>> baseconvert(1000101011,BASE2,BASE10)
    '555'

    base10 to base4
    >>> baseconvert(99,BASE10,"0123")
    '1203'

    base4 to base5 (with alphabetic digits)
    >>> baseconvert(1203,"0123","abcde")
    'dee'

    base5, alpha digits back to base 10
    >>> baseconvert('dee',"abcde",BASE10)
    '99'

    decimal to a base that uses A-Z0-9a-z for its digits
    >>> baseconvert(257938572394L,BASE10,BASE62)
    'E78Lxik'

    ..convert back
    >>> baseconvert('E78Lxik',BASE62,BASE10)
    '257938572394'

    binary to a base with words for digits (the function cannot convert this back)
    >>> baseconvert('1101',BASE2,('Zero','One'))
    'OneOneZeroOne'

    """

    if str(number)[0]=='-':
        number = str(number)[1:]
        neg=1
    else:
        neg=0

    # make an integer out of the number
    x=long(0)
    for digit in str(number):
       x = x*len(fromdigits) + fromdigits.index(digit)
    
    # create the result in base 'len(todigits)'
    res=""
    while x>0:
        digit = x % len(todigits)
        res = todigits[digit] + res
        x /= len(todigits)
    if neg:
        res = "-"+res

    return res

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
    
    nmotifs = 4**pwm.__len__()
    print "There are", nmotifs, "motifs to try."
    #exit(1)
    for i in range(0, nmotifs):
        x = baseconvert(i, BASE10, BASE4)
        xstr = x.__str__()
        m = ""
        bits = 0.0
        for j in range(0, xstr.__len__()):
            state = DNA[ int(xstr[j]) ]
            m += state
            bits += pwm[j][state]
        #if bits > 0:
        #    print m, bits
        if bits > cutoff:
            print m, bits
            mlib.append(m)
    print ". . . of which", mlib.__len__(), "satisfy the cutoff."
    return mlib