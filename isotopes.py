#!/usr/bin/env python
# encoding: utf-8
"""
isotopes.py

adapted from an example by FX Coudert and from
Kubinyi, H.
Calculation of isotope distributions in mass spectrometry. A trivial solution for a non-trivial problem.
Anal Chim Acta 247, 107–119 (1991).

Created by DELSUC Marc-André on 2014-03-28.
Copyright (c) 2014 IGBMC. All rights reserved.
"""

import re
from collections import defaultdict, namedtuple
import copy
global  elem_t, name_t, isotope_t
THRESHOLD = 1E-8

class Isotope(namedtuple("Isotope", "element isotop, mass, abund")):
    "a micro class for isotopes entries"
    pass

class Ion(object):      #as if Ion = namedlist("Ion", "mass, proba", verbose=False)
    """
    a class to hold Ion(mass, proba),
    Distribution is built from it
    """
    def __init__(self, mass=0.0, proba=1.0):
        self.mass = mass
        self.proba = proba
    def __repr__(self):
        return "%f  %.10f"%(self.mass, 100*self.proba)

def load_elements(filename="elements.asc"):
    """
    load the nist element file into dictionnaries
    ...
    5   B   10   10.0129370(4)      0.199(7)      10.811(7)        g,m,r
            11   11.0093054(4)      0.801(7)    
    _________________________________________________________________________
    6   C   12   12.0000000(0)      0.9893(8)     12.0107(8)       g,r
            13   13.0033548378(10)  0.0107(8)   
            14   14.003241989(4)                
    ...
    """
    #-------------------
    def readval(st):
        "0.9893(8) => float(0.9893)"
        v = st.split('(')
        try:
            val = float(v[0])
        except ValueError:  # typically radioactive only elemt
            val = 0.0
        return val
    isotope_t = defaultdict(list)     # a dict by element_number of a list of isotopes data
    name_t = {}                     # dict by name of element_number
    with open(filename,'r') as F:
        elem = 0    # current element
        while F:        # implement a simple parser through a state machine
            l = F.readline()
            f = l.split()
            if re.match("[1-9]",l):     # 6   C   12   12.0000000(0)      0.9893(8)     12.0107(8)       g,r
                elem = int(f[0])        # 6
                if elem>100:
                    break
                name = f[1]             # C
                isotop = int(f[2])      # 12
                mass = readval(f[3])    # 12.0000
                abund = readval(f[4])   # 0.9893
                name_t[name] = elem
                isotope_t[elem].append(Isotope(elem, isotop, mass, abund))  # (6, 12, 12.0000000, 0.9893)
            elif elem>0:
                if re.match("[1-9]",f[0]):  #     13   13.0033548378(10)  0.0107(8)
                    isotop = int(f[0])      # 13
                    mass = readval(f[1])    # 13.0033548378
                    try:
                        abund = readval(f[2])   # 0.0107
                    except IndexError:
                        abund = 0.0
                    isotope_t[elem].append(Isotope(elem, isotop, mass, abund))  # 
            
    elem_t = dict(zip(name_t.values(), name_t.keys()))      # dict by element_number of name
    return (elem_t, name_t, isotope_t )

def print_t():
    " print out the table read by load_elements()"
    for k in isotope_t.keys():
        print elem_t[k], k
        for i in isotope_t[k]:
            print "    ", i


def parse_seq(st):
    """
    parse a raw formula "st" to a list of atom and stoechio
    "CCl4" returns  ["C":1, "Cl":4]
    """
    # to rewritten recursive as a list processing to accept ()
    formula = defaultdict(int)
    for i in range(len(st)):   # through the string
#        print i,st[i:i+2]
        if st[i] == " ":    # skip blanks
            continue
        if st[i:i+2] in name_t:     # double letter elem
            k = st[i:i+2]
            m = re.match("[0-9]+",st[i+2:])
            if m:
                formula[k] += int(m.group(0))
            else:
                formula[k] += 1
        elif st[i] in name_t:       # single letter
            k = st[i]
            m = re.match("[0-9]+",st[i+1:])
            if m:
                formula[k] += int(m.group(0))
            else:
                formula[k] += 1
    return formula

def addformula(f1,f2):
    """add the content of f2 to the content of f1"""
    for i in f2.keys():
        f1[i] += f2[i]
    
def parse_pep(st):
    "compute the formula of a peptide/protein given par one letter code"
    formula = parse_seq("NH2")   # starts with H2N-...
    cterm = parse_seq("COOH")   # end with ..-COOH
    pbound = parse_seq("CO NH")
    H2O = parse_seq("H2O")
    for ires in range(len(st)):
        res = st[ires]
        if res == "A":          f = parse_seq("CH CH3")
        elif res == "C":        f = parse_seq("CH CH2 SH")
        elif res == "D":        f = parse_seq("CH CH2 COOH")
        elif res == "E":        f = parse_seq("CH CH2 CH2 COOH")
        elif res == "F":        f = parse_seq("CH CH2 C6H5")
        elif res == "G":        f = parse_seq("CH2")
        elif res == "H":        f = parse_seq("CH CH2 C3 N2 H3")
        elif res == "I":        f = parse_seq("CH CH CH3 CH2 CH3")
        elif res == "K":        f = parse_seq("CH CH2 CH2 CH2 CH2 NH2")
        elif res == "L":        f = parse_seq("CH CH2 CH CH3 CH3")
        elif res == "M":        f = parse_seq("CH CH2 CH2 S CH3")
        elif res == "N":        f = parse_seq("CH CH2 CONH2")
        elif res == "P":        f = parse_seq("C CH2 CH2 CH2")
        elif res == "Q":        f = parse_seq("CH CH2 CH2 CONH2")
        elif res == "R":        f = parse_seq("CH CH2 CH2 CH2 N C NH2 NH2")
        elif res == "S":        f = parse_seq("CH CH2 OH")
        elif res == "T":        f = parse_seq("CH CHOH CH3")
        elif res == "V":        f = parse_seq("CH CH CH3 CH3")
        elif res == "W":        f = parse_seq("CH CH2 C8 N H6")
        elif res == "Y":        f = parse_seq("CH CH2 C6H4OH")
        elif res == "*":    # star notes phosphate
            f = parse_seq("PO4")
            for i in H2O.keys():
                formula[i] -= H2O[i]    # remove one H2O
        elif res == "a":    # c notes acetate
            f = parse_seq("COOCH3")
            for i in H2O.keys():
                formula[i] -= H2O[i]    # remove one H2O
        elif res == "n":    # n notes amide in Cter
            f = parse_seq("NH2")
            for i in H2O.keys():
                formula[i] -= H2O[i]    # remove one H2O
        else:
            raise(Exception("Unknown residue code"))
        addformula(formula, f)
        if not res in ('*','c','n'):    # these are not a peptide bound
            if ires<len(st)-1:
                addformula(formula, pbound)
    # for i in pbound.keys():
    #     formula[i] -= pbound[i]    # remove last one
    addformula(formula,cterm)
    return formula

def printformula(formula):
    "nice print of a formula"
    st =""
    for k in sorted(formula.keys()):
        st += "%s_%d "%(k, formula[k])
    print st
def monoisotop(formula):
    "returns monoisotopic mass from a formula"
    mass = 0.0
    for el in formula.keys():
        iso =  isotope_t[name_t[el]]
        mono = sorted(iso, key=lambda e: e.abund, reverse=True)[0]     # find most abundant
        mass += mono.mass*formula[el]
    return mass

def average(formula):
    "returns average mass from a formula"
    mass = 0.0
    for el in formula.keys():
        iso =  isotope_t[name_t[el]]
        ave = sum([e.mass*e.abund for e in iso] )
        mass += ave*formula[el]
    return mass


class Distribution(object):
    """
    handle and compute isotopic distribution
    
    """
    def __init__(self, formula=None, isotope=None):    # init either from a formula or a isotope
        self.threshold = THRESHOLD
        self.distrib = [Ion(mass=0,proba=1)]    # starts with empty        
        if formula:
            self.compute(formula)
        if isotope:
            d = []
            for el in isotope:
                d.append(Ion(el.mass,el.abund))
            self.distrib = d
    def len(self):
        """len of the Distribution"""
        return len(self.distrib)
    def __repr__(self, threshold=1E-6):
        "print the distribution, only peaks above threshold"
        st = ""
        for i in self.distrib:
            if i.proba>threshold:
                st += repr(i)+"\n"
        return st
    def sort_by_intens(self):
        " sort distrib by decreasing intensity"
        self.distrib.sort(key=lambda ion: ion.proba, reverse=True)
    def sort_by_mass(self):
        " sort distrib by increasing mass"
        self.distrib.sort(key=lambda ion: ion.mass)

    def combine(self,dist2):
        """combine two Distribution"""
        d = [Ion(mass=0.0, proba=0.0) for i in range((self.len()+dist2.len()-1)) ]
        for i in range(self.len()):
            for j in range(dist2.len()):
                d[i+j].proba += self.distrib[i].proba*dist2.distrib[j].proba
        for i in range(self.len()):
            for j in range(dist2.len()):
#      res->p[i+j].mass += dist1.p[i].prob * dist2.p[j].prob * (dist1.p[i].mass + dist2.p[j].mass);
                d[i+j].mass += (self.distrib[i].mass + dist2.distrib[j].mass) \
                                * self.distrib[i].proba * dist2.distrib[j].proba
        for ion in d:
            if ion.proba > 0.0:
                ion.mass /= ion.proba
        self.distrib = copy.deepcopy(d)
        self.normalize()
        self.prune()

    def compute(self, formula):
        """compute a Distribution from a formula returned by parse_seq()"""
        for el in formula.keys():
            iso =  isotope_t[name_t[el]]
            d = Distribution(isotope=iso)
            dd = Distribution()
            for i in range(formula[el]):
                dd.combine(d)
            self.combine(dd)

    def normalize(self):
        """normalize distrib such that max() == 1"""
        M = max([self.distrib[i].proba for i in range(self.len())])
        for i in range(self.len()):
            self.distrib[i].proba /= M

    def prune(self):
        """prune distrib to remove entries below threshold"""
        while self.distrib[0].proba < self.threshold:      # low
            self.distrib.pop(0)
        while self.distrib[-1].proba < self.threshold:      # high
            self.distrib.pop()

def draw(D, width=0.2, charge=1):
    "quick draw a given distribution"
    import numpy as np
    import matplotlib.pyplot as plt
    D.sort_by_mass()
    x = np.zeros(3*D.len()+2)
    y = np.zeros(3*D.len()+2)
    x[0] = D.distrib[0].mass/charge-5
    x[-1] = D.distrib[-1].mass/charge+5
    for i in range(D.len()):
        x[1+3*i] =  (D.distrib[i].mass)/charge - width
        x[2+3*i] =  (D.distrib[i].mass)/charge
        y[2+3*i] =  100*D.distrib[i].proba
        x[3+3*i] =  (D.distrib[i].mass)/charge + width
    plt.plot(x,y)
    plt.xlabel("$m/z$")
    plt.ylabel("intensity")
    plt.axis(ymin=-5, ymax=105)
    plt.show()
    
if __name__ == '__main__':
    (elem_t, name_t, isotope_t ) = load_elements()
#    print_t()
    test = "CH4 C2H6 CH3COOH PO4"
#    test = "K23 I22"
#    test = "W5"
    test = "C254 H377 N65 O75 S6"   # insuline
#    test = "C1185 H1850 N282 0339 S18"  # cytochrome oxydase
    test2 = "RPKPQQFFGLM"   # substance P
    test2 = "KELCKAVSVSM"
#    test2 = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"
    form = parse_pep(test2)
    addformula(form, form)
#    form = parse_seq(test)
    printformula( form)
    print monoisotop(form), average(form)
    D = Distribution(form)
    print "By mass\n",D,"\n"
    D.sort_by_intens()
    print "By intensities\n",D,"\n"
    draw(D)
