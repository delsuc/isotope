#!/usr/bin/env python
# encoding: utf-8
"""
isotopes.py

adapted from 
Kubinyi, H.
Calculation of isotope distributions in mass spectrometry. A trivial solution for a non-trivial problem.
Anal Chim Acta 247, 107–119 (1991).

First version by FX Coudert,
Python rewrite by DELSUC Marc-André on 2014-03-28.
Copyright (c) 2014 CNRS. All rights reserved.
"""

import re
import os
from collections import defaultdict, namedtuple
import copy
global  elem_t, name_t, isotope_t
THRESHOLD = 1E-8
import matplotlib.pyplot as plt
import numpy as np



class Formula( defaultdict ):
    """
    a micro class for chemical formula entries
    hold a  a dictionnary {"elem":number} (with default value of zero)
    "CH3 CH2 OH" is coded as  ["C":2, "O":1, "H":6]
    """
    def __init__(self,*arg):
        defaultdict.__init__(self,int)
    def __repr__(self):
        return printformula(self)
    def monoisotop(self):
        "return monoisotopique mass"
        return monoisotop(self)
    def average(self):
        "return average mass"
        return average(self)
    def distribution(self):
        "return the mass distribution as a Distribution object"
        return Distribution(self)
    
# class Isotope(namedtuple("Isotope", "element isotop mass abund")):
#     "a micro class for isotopes entries"
#     pass
class Isotope(object):
    """
    a micro class for isotopes entries - 4 attributes
        e.g. 
        self.element = 6
        self.isotop = 12
        self.mass = 12.0
        self.abund = 0.9893

    """
    def __init__(self, element, isotop, mass, abund):
        self.element = element
        self.isotop = isotop
        self.mass = mass
        self.abund = abund
    def __str__(self):
        return "Isotope(element=%s, isotop=%s, mass=%s, abund=%s)"%(self.element, self.isotop, self.mass, self.abund)

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

def load_elements(filename=None):
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
    loads from filename, if empty, loads elements.asc located next to the source prgm.


    returns (elem_t, name_t, isotope_t )
    where
    elem_t : {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C',...   correspondance   Z / name
    name_t : {'H': 1,'He': 2,'Li': 3, ... is the reverse table name /Z 
    isotope_t : {1: (list, of, Isotope_objects), 2: ...    is the isotopic table
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
    #-------------
    # if filename empty, 
    if filename == None:
        filename = os.path.join(os.path.dirname(__file__),"elements.asc")
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

def enrich(element="C",isotop=13, ratio=1.0):
    """
    modifies element abundance by modifying the internal isotopic table

    enrich(element="C",isotop=13, ratio=0.95)  indicates a 95% 13C enrichment
    """
    k = name_t[element]
    for i in isotope_t[k]:  # find it
        if i.isotop == isotop:
            prev = i.abund
    for i in isotope_t[k]:  # modify it
        if i.isotop == isotop:
            i.abund = ratio
        else:
            i.abund = i.abund*(1-ratio)/(1-prev)
        print "    ", i
    
def print_t():
    " print out the table read by load_elements()"
    for k in isotope_t.keys():
        print elem_t[k], k
        for i in isotope_t[k]:
            print "    ", i


def parse_formula(st):
    """
    parse a raw formula "st" to a list of atom and stoechio
    "CCl4" returns  ["C":1, "Cl":4]
    """
    # to rewritten recursive as a list processing to accept ()
    formula = Formula()
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
    
def parse_peptide(st):
    """
    compute the formula of a peptide/protein given par one letter code
    formula = parse_peptide("ACDEY*GH")     # e.g.
    letter code is standard 1 letter code for amino-acides + additional codes for Post Translational Modifications (PTM)
    * posphrylation
    a acetylation
    m methoxylation
    n amidation
    - deamidation
    h hydroxylation
    o oxydation
    does not verify the chemical coherence of the PTM !
    """
    formula = parse_formula("NH2")   # starts with H2N-...
    cterm = parse_formula("COOH")   # end with ..-COOH
    pbound = parse_formula("CO NH")
    AA = "ACDEFGHIKLMNPQRSTVWY"
    PTM = "*amn-ho+"
    for ires in range(len(st)):
        res = st[ires]
        if res == " ":
            continue
        elif res in AA:
            if res == "A":          f = parse_formula("CH CH3")
            elif res == "C":        f = parse_formula("CH CH2 S H")
            elif res == "D":        f = parse_formula("CH CH2 COOH")
            elif res == "E":        f = parse_formula("CH CH2 CH2 COOH")
            elif res == "F":        f = parse_formula("CH CH2 C6H5")
            elif res == "G":        f = parse_formula("CH2")
            elif res == "H":        f = parse_formula("CH CH2 C3 N2 H3")
            elif res == "I":        f = parse_formula("CH CH CH3 CH2 CH3")
            elif res == "K":        f = parse_formula("CH CH2 CH2 CH2 CH2 NH2")
            elif res == "L":        f = parse_formula("CH CH2 CH CH3 CH3")
            elif res == "M":        f = parse_formula("CH CH2 CH2 S CH3")
            elif res == "N":        f = parse_formula("CH CH2 CONH2")
            elif res == "P":        f = parse_formula("C CH2 CH2 CH2")
            elif res == "Q":        f = parse_formula("CH CH2 CH2 CONH2")
            elif res == "R":        f = parse_formula("CH CH2 CH2 CH2 N C NH2 NH2")
            elif res == "S":        f = parse_formula("CH CH2 OH")
            elif res == "T":        f = parse_formula("CH CHOH CH3")
            elif res == "V":        f = parse_formula("CH CH CH3 CH3")
            elif res == "W":        f = parse_formula("CH CH2 C8 N H6")
            elif res == "Y":        f = parse_formula("CH CH2 C6H4OH")
            addformula(formula, f)
        # special codes
        elif res in PTM:
            if res == "*":    # star notes phosphate
                f = parse_formula("PO4")
                to_rem = parse_formula("OH")
            elif res == "a":    # c notes acetate
                f = parse_formula("COOCH3")
                to_rem = parse_formula("H2O")
            elif res == "n":    # amidation (in Cter)
                f = parse_formula("NH2")
                to_rem = parse_formula("OH")
            elif res == "-":    # deamination
                f = parse_formula("OH")
                to_rem = parse_formula('NH2')
            elif res == "h":    # hydroxylation (eg Prolines)
                f = parse_formula("O")
                to_rem = {}
            elif res == "+":    # protonation
                f = parse_formula("H")
                to_rem = {}
            elif res == "o":    # oxydation (eg methionine)
                f = parse_formula("O")
                to_rem = {}
            elif res == "m":    # methylation
                f = parse_formula("CH2")
                to_rem = {}
            # remove 
            for i in to_rem.keys():
                formula[i] -= to_rem[i]
            addformula(formula, f)
        else:
            raise(Exception("Unknown residue code"))
        # then add pbound
        if res in AA: # add pbound only if AA
            if [r for r in st[ires+1:] if r in AA ]:    # and if there is another AA next on st
                addformula(formula, pbound)
    addformula(formula,cterm)
    return formula

def printformula(formula):
    "nice print of a formula"
    st =""
    for k in sorted(formula.keys()):
        if formula[k] >1:
            sto = "_%d"%(formula[k])
        else:
            sto = ""
        st += "%s%s "%(k, sto)
    return st
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
    
    stored in self.distrib as a list of Ion() : [Ion(),...]
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
        """combine two Distributions"""
        d = [Ion(mass=0.0, proba=0.0) for i in range((self.len()+dist2.len()-1)) ]
        for i in range(self.len()):
            for j in range(dist2.len()):
                d[i+j].proba += self.distrib[i].proba*dist2.distrib[j].proba
        for i in range(self.len()):
            for j in range(dist2.len()):
                d[i+j].mass += (self.distrib[i].mass + dist2.distrib[j].mass) \
                                * self.distrib[i].proba * dist2.distrib[j].proba
        for ion in d:
            if ion.proba > 0.0:
                ion.mass /= ion.proba
        self.distrib = copy.deepcopy(d)
        self.normalize()
        self.prune()

    def compute(self, formula):
        """compute a Distribution from a formula returned by parse_formula()"""
        for el in formula.keys():
            iso =  isotope_t[ name_t[el] ]
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
    #---------------- utilities --------------
    def draw(self, title=None, width=0.2, charge=1, R=None):
        "quick draw of distribution"
        import numpy as np
        if R:
            width = 0.25*self.distrib[0].mass/R
        self.sort_by_mass()
        x = np.zeros(3*self.len()+2)
        y = np.zeros(3*self.len()+2)
        x[0] = self.distrib[0].mass/charge-5
        x[-1] = self.distrib[-1].mass/charge+5
        for i in range(self.len()):
            x[1+3*i] =  (self.distrib[i].mass)/charge - width
            x[2+3*i] =  (self.distrib[i].mass)/charge
            y[2+3*i] =  100*self.distrib[i].proba
            x[3+3*i] =  (self.distrib[i].mass)/charge + width
        plt.plot(x,y)
        plt.xlabel("$m/z$")
        plt.ylabel("intensity")
        plt.axis(ymin=-5, ymax=105)
        if title:
            plt.title(title)

    def draw_lowres(self, title=None, charge=1):
        """
        draw the low resolution peak
        """
        spl = self.enveloppe()
        xenv = np.linspace(self.distrib[0].mass-5, self.distrib[-1].mass+5, 1000)
        env = [100*spl(xi) for xi in xenv]
        plt.plot(xenv/charge, env)
        if title:
            plt.title(title)

    def enveloppe(self):
        """
        compute the smoothed enveloppe for a distribution
        f = D.enveloppe()
        returns a function f which computes the enveloppe on any x points : f(x)
        """
        from scipy.interpolate import UnivariateSpline
        from scipy.optimize import newton
        from functools import partial
        m0 = self.distrib[0].mass
        x = [m0-1]
        y = [0]
        x += [dd.mass for dd in self.distrib]
        y += [dd.proba for dd in self.distrib]
        x += [self.distrib[0].mass+5]
        y += [0]
        x = np.array(x)
        y = np.array(y)
        spl = UnivariateSpline(x, y, s=0)
        deriv =  partial( spl, nu=1)
        second = partial( spl, nu=2)
        x0 = sum(x*y)/sum(y)        # average mass
        # summit of curve is :  newton(deriv, x0, fprime=second)
        return lambda x: max(0.0,spl(x)) if x>m0-1 else 0    # 

def test1():
    " example with elemental formula"
#    test = "C6H14"
#    test = "K23 I22 S30"
    #   test = "W5"
    test = "C254 H377 N65 O75 S6"   # insuline
    #    test = "C1185 H1850 N282 O339 S18"  # cytochrome oxydase
    # test = "C769 H1212 N210 O218 S2 H20"  # myoglobine 20+
    form = parse_formula(test)
    print "insuline", printformula( form), monoisotop(form), average(form)

def test2():
    " example with protein formula"
#    test = "RPKPQQFFGCLMn"   # substance P
#    test = "MKVLWAALLV TFLAGCQAKV EQAVETEPEP ELRQQTEWQS GQRWELALGR FWDYLRWVQT LSEQVQEELL SSQVTQELRA LMDETMKELK AYKSELEEQL TPVAEETRAR LSKELQAAQA RLGADMEDVC GRLVQYRGEV QAMLGQSTEE LRVRLASHLR KLRKRLLRDA DDLQKRLAVY QAGAREGAER GLSAIRERLG PLVEQGRVRA ATVGSLAGQP LQERAQAWGE RLRARMEEMG SRTRDRLDEV KEQVAEVRAK LEEQAQQIRL QAEAFQARLK SWFEPLVEDM QRQWAGLVEK VQAAVGTSAA PVPSDNH" # ApoE
    test = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG" # ubiquitine
    form = parse_peptide(test)
#    addformula(form, parse_formula("H12")) # if you want o compute [M+H_12]12+
    
    print "Ubiquitine",printformula( form)
    print monoisotop(form), average(form)
    return Distribution(form)

def insu():
    cc = parse_peptide("GIVEQCCASVCSLYQLENYCN")
    cd = parse_peptide("FVNQHLCGSHLVEALYLVCGERGFFYTPKA")
    print "Chain, C", monoisotop(cc)
    print "Chain, D", monoisotop(cd)
    insuline = copy.deepcopy(cc)
    addformula(insuline,cd)   # add both chains
    insuline["H"] -= 6      # remove 6 H for 3 disulfides
    print "insuline", monoisotop(insuline)
    Ins = Distribution(insuline)
    Ins.draw(width=0.1,title="Insuline")
    
# load table at import
(elem_t, name_t, isotope_t ) = load_elements()
if __name__ == '__main__':
    D = test2()
    print "By mass\n",D,"\n"
    D.draw(charge=12, R=1E5, title="Ubiquitine 12+")
    D.draw_lowres(charge=12)
    plt.figure()
    insu()
    plt.show()
    