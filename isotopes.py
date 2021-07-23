#!/usr/bin/env python
# encoding: utf-8
"""
isotopes.py

This module allows the precise computation of isotopic distributions, average mass, and monoisotopic masses
of arbitrary chemical formula.
Isotopic natural abundance can be changed with the enrich() function.

The handy determination from peptide and protein primary sequences is provided.

The possibility to draw the isotopic profil is also given.

Molecules can be described by chemical formula or peptide primary sequence (one letter code)
so the flollowing entries are legal
"CH3 CH2 OH"  
"(NO3)3 C6H3"    # parenthesis can be nested.
"KEVDFS*GHI"     # for a peptide (* here means phosphorylated (see parse_peptide for the list of PTM)

The main class is Formula from which everything is computed.
Arithmetic on formula is possible

Typical use is :

#parsing
molecule = "K23 I22 S30"
mol_formula = parse_formula(molecule)        # formula object (in fact a dictionary)
# arithmetic
hydrated = mol_formula + 2*parse_formula("H2O")
# masses
print "average mass", hydrated.average()
print "monoisotopic mass", hydrated.monoisotop()
# compute isotopic distrib
distrib = mol_formula.distribution()
distrib.draw()

There is a cascade of data representations

String           either raw atomic composition or one-letter peptide sequence)
    from which you compute the Formula()
Formula()        contains the raw formula with the number of each atomic species
    - you can combine and do arithmetics with them
    from which you can compute
    - monoisotopic()
    - average()  - these are fast
    - Distribution()

Distribution()   contains the list of each exact mass with intensities
    - print
    - draw()
    - sort by intensities or masses

of course the isotopic distribution (and average mass) depend on isotopic abundance, which are either looaded from natural distribution (from NIST db) or modified with the enrich() function.

There are also low level structures ( Ion() and Isotop() that you can usually ignore )

adapted from 
Kubinyi, H.
Calculation of isotope distributions in mass spectrometry. A trivial solution for a non-trivial problem.
Anal Chim Acta 247, 107-119 (1991).

This module depends at import time on a isotope file to be present in the same directory
by default called "elements.asc"

It is a copy of the file found at
NIST | Physical Measurement Laboratory | Physical Reference Data | Atomic Weights and Isotopic Compositions Main Page
http://www.nist.gov/

First version of algo by FX Coudert,
Python rewrite by DELSUC Marc-Andre on 2014-03-28.
madelsuc@unistra.fr
Copyright (c) 2014 CNRS. All rights reserved.

recursive descent chemical parser weakly inspired from
Tim Peters tim_one@email.msn.com
https://mail.python.org/pipermail/tutor/1999-March/000083.html

code for fine isotopic pattern require the external neutronstar prgm - from https://bitbucket.org/orserang/neutronstar.git

The version 1.1 introduces the fast computation of fine isotopic pattern using the neutronstar algo.
"""

from __future__ import division, print_function
import re
import os
from collections import defaultdict, namedtuple
import copy
import unittest
import re
import warnings
from subprocess import check_output
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np

version = "1.1"

#####  Configuration ##############################
# code for fine isotopic pattern
# need to install neutronstar - from https://bitbucket.org/orserang/neutronstar.git
# then enter in NS the absolute location of the neutronstar binary code
NS = "/home/mad/neutronstar/src/neutronstar"

# THRESHOLD is the default value used for pruning when computing isotopic patterns
THRESHOLD = 1E-8

#####  END Configuration ##############################

_lexer = re.compile(r"[A-Z][a-z]*|\d+|\+-|[()]|<EOS>")

global  elem_t, name_t, isotope_t

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
    def __len__(self):
        return len(self.keys())
    def __imul__(self, other):
        if isinstance(other, int):
            multformula(self, other)
            return self
        else:
            return NotImplemented
    def __iadd__(self, other):
        if isinstance(other, Formula):
            addformula(self, other)
            return self
        else:
            return NotImplemented
    def __add__(self, other):
        if isinstance(other, Formula):
            res = copy.deepcopy(self)
            addformula(res, other)
            return res
        else:
            return NotImplemented
    def __isub__(self, other):
        if isinstance(other, Formula):
            rmformula(self, other)
            return self
        else:
            return NotImplemented        
    def __sub__(self, other):
        if isinstance(other, Formula):
            res = copy.deepcopy(self)
            rmformula(res, other)
            return res
        else:
            return NotImplemented        
    def __mul__(self, other):
        if isinstance(other, int):
            res = copy.deepcopy(self)
            multformula(res, other)
            return res
        else:
            return NotImplemented
    def __rmul__(self, other):
        if isinstance(other, int):
            res = copy.deepcopy(self)
            multformula(res, other)
            return res
        else:
            return NotImplemented

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
                    if isotope_t[elem][-1].isotop != isotop-1:
#                        print('*** curring', name)
#   adding 0 abundance isotopes when there is a gap in isotopic scale, 
#   to cure a bug in the algorithm
                        isotope_t[elem].append(Isotope(elem, isotop-1, float(isotop-1), 0.0))
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
        print ("    ", i)

def print_t():
    " print out the table read by load_elements()"
    for k in isotope_t.keys():
        print (elem_t[k], k)
        for i in isotope_t[k]:
            if i.abund != 0.0:
                print ("    ", i)

""
def parse(st):
    """
    generic parser, first try parse_peptide(), then parse_formula(), then give up
    """
    f = None
    try:
        f = parse_peptide(st)
    except:
        try:
            f = parse_formula(st)
        except:
            pass
    if f is None:
        raise Exception("Could not understand the input, not a peptide nor a raw formula")
    return f

def parse_formula(st):
    """
    parse a raw formula "st" to a list of atom and stoechio
    "CCl4" returns  ["C":1, "Cl":4]
    possible formula are
        "CH3 CH2 OH"   "(CH3)3 C Cl"  "HO (CH2 CH2 O)30 H" etc...
        charges are ignored
    """
    def _parse_formula(tklist, match=None):
        """
        recursive part of parse_formula()
        match is either None or '(' and tells what is to be matched on exit
        """
        def countfromtk():
            "utility for parsing counts"
            if tklist:      # if tklist not empty
                nexttoken = tklist.pop(0)   # is next a count
                if re.match("[0-9]+",nexttoken):  # count
                    count = int(nexttoken)
                else:
                    tklist.insert(0,nexttoken)    # put it back
                    count = 1
            else:
                count = 1
            return count
            #----
        formula = Formula()     #initialize
        while tklist:
            token = tklist.pop(0)
            form = Formula()
            if token in name_t:     # element
                elem = token
                count = countfromtk()
                form[elem] = count
            elif token  == "(":
                form = _parse_formula(tklist, token)
                count = countfromtk()
                if count != 1:
                    multformula(form,count)
            elif token == ")":    # end of series
                if match == "(":
                    break
                else:
                    raise Exception("Unmatching parenthesis")
            elif token == "<EOS>":
                if match == "<EOS>":
                    break
                else:
                    raise Exception("Unmatching parenthesis")
            elif token in (" ","\t","+","-"):
                pass
            else:
                raise Exception("Unkown element : ",token)
            addformula(formula,form)
        else:
            raise Exception("Undecipherable chain")
        return formula
        #--------- end of _parse_formula()
    tklist = re.findall(_lexer, st+"<EOS>")     # tokenize
    form = _parse_formula(tklist, match="<EOS>")   # then call recursive function
    return form


def addformula(f1,f2):
    """add inplace the content of f2 to the content of f1"""
    for i in f2.keys():
        f1[i] += f2[i]

def rmformula(f1,f2, check=True):
    """remove inplace the content of f2 to the content of f1
    check enforces a positive atome count"""
    for i in f2.keys():
        f1[i] -= f2[i]
        if check:
            if f1[i]<0 : raise Exception("%s: Negative atom count !"%i)
def multformula(f1, scalar):
    """multiply inplace the content of f1 by scalar"""
    for i in f1.keys():
        f1[i] *= scalar

def parse_peptide(st, extended=False, starts="NH2", ends="COOH"):
    """
    compute the formula of a peptide/protein given by one letter code
    
    formula = parse_peptide("ACDEY*GH")     # e.g.
    letter code is standard 1 letter code for amino-acids
    
    for xyz/abc fragmentation

    - starts is either 
        - "NH2" (default) - for a regular peptide
        - "x"  "y"  "z" - for MS fragments 
        - or any formula
    - ends is either
        - "COOH" (default) - for a regular peptide
        - "a"  "b"  "c" - for MS fragments 
        -  or any formula
    
    additional codes for Post Translational Modifications (PTM)

    * phosphorylation
    a acetylation
    n amidation
    - deamidation
    h hydroxylation
    o oxydation
    + protonation
    m methoxylation
    b beta-mercaptoethanol adduct
    does not verify the chemical coherence of the PTM !
    
    if extended is True, will also interpret U: Seleno-Cysteine and O: Pyrolysine
    
    
    """
    if starts == "x":
        formula = parse_formula("CONH")
    elif starts == "y":
        formula = parse_formula("NH3")  # shouldn't it be "NH" - no it seems, standard y ions are -NH3+ 
    elif starts == "z":
        formula = Formula()
    else:
        formula = parse_formula(starts)    # default is starts with H2N-...
        
    if ends == "a":
        cterm = Formula()
    elif ends == "b":
        cterm = parse_formula("CO")
    elif ends == "c":
        cterm = parse_formula("CONH3")  # idem y ions
    else:
        cterm = parse_formula(ends)  # default is end with ..-COOH
        
    pbound = parse_formula("CO NH")
    AA={}
    AA["A"] = parse_formula("CH CH3")
    AA["C"] = parse_formula("CH CH2 S H")
    AA["D"] = parse_formula("CH CH2 COOH")
    AA["E"] = parse_formula("CH (CH2)2 COOH")
    AA["F"] = parse_formula("CH CH2 C6H5")
    AA["G"] = parse_formula("CH2")
    AA["H"] = parse_formula("CH CH2 C3 N2 H3")
    AA["I"] = parse_formula("CH CH CH3 CH2 CH3")
    AA["K"] = parse_formula("CH (CH2)4 NH2")
    AA["L"] = parse_formula("CH CH2 CH CH3 CH3")
    AA["M"] = parse_formula("CH (CH2)2 S CH3")
    AA["N"] = parse_formula("CH CH2 CONH2")
    AA["P"] = parse_formula("C (CH2)3")
    AA["Q"] = parse_formula("CH (CH2)2 CONH2")
    AA["R"] = parse_formula("CH (CH2)3 N C (NH2)2")
    AA["S"] = parse_formula("CH CH2 OH")
    AA["T"] = parse_formula("CH CHOH CH3")
    AA["V"] = parse_formula("CH CH CH3 CH3")
    AA["W"] = parse_formula("CH CH2 C8 N H6")
    AA["Y"] = parse_formula("CH CH2 C6H4 OH")
    if extended:
        AA["U"] = parse_formula("CH CH2 Se H")  # Selenocysteine
        AA["O"] = parse_formula("CH (CH2)4 NH CO CH CH CH3 CH2 CH N")  # Pyrrolysine
    AAk = AA.keys()
    #PTM coded as a pair of formula [to_add, to_remove]
    PTM = {}    
    PTM["*"] = [parse_formula("PO4H2"), parse_formula("OH")]    # star notes phosphate
    PTM["a"] = [parse_formula("COO CH3"), parse_formula("H2O")]    # c notes acetate
    PTM["b"] = [parse_formula("SH CH2 CH2 OH"), parse_formula("H2")]    # b is for beta mercapto (on Cys)
    PTM["n"] = [parse_formula("NH3"), parse_formula("OH")]    # amidation (in Cter) assumes NH3+
    PTM["-"] = [parse_formula("OH"), parse_formula('NH3')]    # deamination  assumes NH3+
    PTM["h"] = [parse_formula("O"), {}]    # hydroxylation (eg Prolines)
    PTM["+"] = [parse_formula("H"), {}]    # protonation
    PTM["o"] = [parse_formula("O"), {}]    # oxydation (eg methionine)
    PTM["m"] = [parse_formula("CH2"), {}]    # methylation
    PTMk = PTM.keys()
    for ires in range(len(st)):
        res = st[ires]
        if res == " ":
            continue
        elif res in AAk:
            addformula(formula, AA[res])
        # special codes
        elif res in PTMk:
            to_add, to_rem = PTM[res]
            addformula(formula, to_add)
            rmformula(formula, to_rem)
        else:
            raise(Exception("Unknown residue code"))
        # then add pbound
        if res in AAk: # add pbound only if AA
            if [r for r in st[ires+1:] if r in AAk ]:    # and if there is another AA next on st
                addformula(formula, pbound)
    addformula(formula,cterm)
    return formula

def rawformula(formula):
    "normalized print of a formula - needed for neutronstar"
    st =""
    for k in sorted(formula.keys()):
        if formula[k] == 0:  # appears sometimes when doing formula arithmetics !
            continue
        else:
            sto = "%d"%(formula[k])
        st += "%s%s"%(k, sto)
    return st.strip()

def printformula(formula):
    "nice print of a formula"
    st =""
    for k in sorted(formula.keys()):
        if formula[k] == 0:  # appears sometimes when doing formula arithmetics !
            continue
        if formula[k] >1:
            sto = "_%d"%(formula[k])
        else:
            sto = ""
        st += "%s%s "%(k, sto)
    return st.strip()

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

def spectrify(masses, probs, RP=10000, normalize=100.0, point_spacing=None):
    '''
    Compute a gaussian mass spectrum spectrum from a distribution.

    masses is an ordered array of the masses.
    probas is an array of the peak intensities.
    RP is the resolving power of the gaussian
    normalize is the value of the largest point
    point_spacing is the resolution of the mass axis
        computed from RP if missing

    returns mass,spec
    plot(mass,spec) will draw the spectrum

    inspired from a code by P.O'Connor
    '''
    spread = masses[-1]-masses[0]+1
    c = masses[0]/(2.3548*RP)   # 2.3548 is the FWMH of a gaussian
    var = 2*(c**2)

    if point_spacing is None:
        point_spacing = 0.25*np.sqrt(var)
    if RP>10000:
        margin = 0.5
    elif RP>1000:
        margin = 3
    else:
        margin = 8
    mass = np.arange(masses[0]-margin,masses[-1]+margin,point_spacing)
    spec = np.zeros_like(mass)

    for m,p in zip(masses,probs):
        diff = (mass-m)**2
        spec +=  p * np.exp(-diff/var)
    spec *= normalize/spec.max()
    return mass,spec
""
class Distribution(object):
    """
    handle and compute isotopic distribution
    
    stored in self.distrib as a list of Ion() : [Ion(),...]
    """
    def __init__(self, formula=None, isotope=None):    # init either from a formula or a isotope
        self.threshold = THRESHOLD
        self.distrib = [Ion(mass=0,proba=1)]    # starts with empty        
        if formula is not None:
            self.compute(formula)
        if isotope is not None:
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
        dprint = lambda *x: None
        d = [Ion(mass=0.0, proba=0.0) for i in range((self.len()+dist2.len()-1)) ]
        dprint('self',self.distrib)
        dprint(2,d)
        for i in range(self.len()):
            for j in range(dist2.len()):
                d[i+j].proba += self.distrib[i].proba*dist2.distrib[j].proba
                dprint("p",i,j,d)
        for i in range(self.len()):
            for j in range(dist2.len()):
                d[i+j].mass += (self.distrib[i].mass + dist2.distrib[j].mass) \
                                * self.distrib[i].proba * dist2.distrib[j].proba
                dprint("m",i,j,d)
        for iion in d:
            if iion.proba > 0.0:
                iion.mass /= iion.proba
        self.distrib = copy.deepcopy(d)
        dprint("done",self.distrib)
        self.normalize()
        dprint("normed",self.distrib)
        self.prune()
        dprint("pruned",self.distrib)

    def compute(self, formula):
        """compute a Distribution from a formula returned by parse_formula()"""
        for el in formula.keys():
            iso =  isotope_t[ name_t[el] ]
            d = Distribution(isotope=iso)
            dd = Distribution()
            for i in range(formula[el]):
                dd.combine(d)
#                print(i,dd) ###
            self.combine(dd)

    def normalize(self):
        """normalize distrib such that max() == 1"""
        M = max([ion.proba for ion in self.distrib])
        for i in range(self.len()):
            self.distrib[i].proba /= M

    def prune(self):
        """prune distrib to remove entries below threshold"""
        while self.distrib[0].proba < self.threshold:      # low
            self.distrib.pop(0)
        while self.distrib[-1].proba < self.threshold:      # high
            self.distrib.pop()
    #---------------- utilities --------------
    def bar(self, largest=100, charge=1, label=None, color='g'):
        """
        draw the distribution as bars
        """
        masses = [ion.mass for ion in self.distrib ]
        probas = [ion.proba for ion in self.distrib]
        pmax = max(probas)
        ratio = largest/pmax
        for m,p in zip(masses,probas):
            plt.plot([m/charge,m/charge], [0, p*ratio], color, lw=0.8)
        plt.plot([(masses[0]-0.5)/charge,(masses[-1]+0.5)/charge], [0, 0], color, lw=0.8, label=label)
        plt.xlabel("$m/z$")
        plt.ylabel("intensity")
#        plt.axis(ymin=-5, ymax=105)
    def spectrify(self, largest=100, RP=10000 ):
        """
        generate x,y to draw the distribution at a given resolving power RP
        """
        masses = [ion.mass for ion in self.distrib ]
        probas = [ion.proba for ion in self.distrib]
        m,y = spectrify(masses, probas, normalize=largest, RP=RP)
        return (m,y)
    def draw(self, largest=100, charge=1, RP=10000, label=None, color='g'):
        """
        draw the distribution at a given charge and resolving power RP
        """
        self.sort_by_mass()
        m,y = self.spectrify(largest=largest, RP=RP)
        # m and y are nparray, so we can make arithmetics
        plt.plot(m/charge, y, color, label=label)
        plt.xlabel("$m/z$")
        plt.ylabel("intensity")
#        plt.axis(ymin=-5, ymax=105)

    def _draw(self, title=None, width=0.2, charge=1, R=None, label=None):
        "obsolete version of the draw() method"
        "quick draw of distribution"
        warnings.warn("Obsolete version", UserWarning, stacklevel=1)

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
        plt.plot(x, y, label=label)
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
        probably useless - maybe for teaching
        """
        from scipy.interpolate import UnivariateSpline
        from scipy.optimize import newton
        from functools import partial
        self.sort_by_mass()
        m0 = self.distrib[0].mass
        x = [m0-5]
        y = [0]
        x += [dd.mass for dd in self.distrib]
        y += [dd.proba for dd in self.distrib]
        x += [self.distrib[-1].mass+5]
        y += [0]
        x = np.array(x)
        y = np.array(y)
        spl = UnivariateSpline(x, y, s=0)
        deriv =  partial( spl, nu=1)
        second = partial( spl, nu=2)
        x0 = sum(x*y)/sum(y)        # average mass
        # summit of curve is :  newton(deriv, x0, fprime=second)
        return lambda x: max(0.0,spl(x)) if x>m0-1 else 0    # 

    def tolist(self):
        """
        return the distribution as a plain list [ (mass,intensity), () .. ]
        """
        return [ (ion.mass, ion.proba) for ion in self.distrib ]

def fineisotopicdistrib(formula, length=1000, full_output=False):
    """
    compute a peak a list of all the peaks in a fine isotopic pattern.
    if not full_output
        returns D a normalized iso.Distribution, with the largest isotopes (truncated to length), sorted by abundance
    else
        returns (D, coverage, processing_time_in_sec)
    where coverage is ratio of isotopes covered by D
    
    Do D.prune() on the returned distribution if you want to reduce the number of isotopes
       to only those higher than D.threshold (default to isotope.THRESHOLD - usually 1E-8) 

    It uses the super fast "neutronstar" algorithm by Oliver Serang et al,
        https://dx.doi.org/10.1021/acs.analchem.0c01670
        Anal. Chem. 2020, 92, 10613−10619
    which has to be installed before hand  from https://bitbucket.org/orserang/neutronstar.git

    code probably not valid on windows !
    """
    if os.path.exists(NS) and os.access(NS, os.X_OK):   # check if exists and runable
        retcode = check_output([NS, rawformula(formula), '-k', str(length), '1.05', 'sort_by_abundance'])
    else:
        raise Exception("fine isotopic distribution is not available - please install neutronstar program")

    # parse results
    result = retcode.split(b'\n')     # split lines
    timline = result.pop(0).decode()  # get time
    tim = float(timline.split()[1])
#    print(tim, 'seconds')
    result.pop(0)                     # skip head
    # parse distribution into an isotopes Distribution()
    D = Distribution()
    D.distrib = []
    #D.threshold = 1E-12
    for r in result:
        line = r.decode().split()
        if line[0] == 'TOTAL':        # last line !  
            ltotal = float(line[3])   # contains log(coverage)
            break
        labund = float(line[1])
        mass = float(line[3])
        D.distrib.append(Ion(mass, np.exp(labund)))
#    D.sort_by_mass()
    D.normalize()
    if full_output:
        return (D, np.exp(ltotal), tim)
    else:
        return D
    #pprint( retcode.split(b'\n')[0:10])

class Test(unittest.TestCase):
    """tests """
    def test_formula(self):
        " test elemental formula, masses and operations"
        molecule = "K23 I22 S30 W5"     # does not exist !
        form = parse_formula(molecule)
        self.assertAlmostEqual(form.average(), 5572.31183942, 8)    # test down to the 8th digit
        self.assertAlmostEqual(form.monoisotop(), 5566.98044564, 8)
        self.assertEqual(printformula(form), "I_22 K_23 S_30 W_5")
        form += parse_formula("K2")
        form -= parse_formula("S2")
        self.assertEqual(form["S"],28)
        self.assertEqual(form["K"],25)
        form = parse_formula( " ( (HO (CH2 CH2 O)10 )3 N)2 Cl" )    # Check parenthesis
        form = form + 2*parse_formula("H+")
        self.assertEqual(form["C"],2*3*10*2)
        self.assertEqual(form["H"],2+2*3*(1+10*4))
        self.assertEqual(form["Cl"],1)

    def test_residu(self):
        "test sum of all aa in extended list"
        AAlist = "ACDEFGHIKLMNPQRSTVWYOU"
        masstot = sum( [parse_peptide(aa, extended=True).monoisotop()  for aa in AAlist] )
        self.assertAlmostEqual(masstot, 3160.44812713)

    def test_prot(self):
        " test with protein formula"
        test = "MKVLWAALLV TFLAGCQAKV EQAVETEPEP ELRQQTEWQS GQRWELALGR FWDYLRWVQT LSEQVQEELL SSQVTQELRA LMDETMKELK AYKSELEEQL TPVAEETRAR LSKELQAAQA RLGADMEDVC GRLVQYRGEV QAMLGQSTEE LRVRLASHLR KLRKRLLRDA DDLQKRLAVY QAGAREGAER GLSAIRERLG PLVEQGRVRA ATVGSLAGQP LQERAQAWGE RLRARMEEMG SRTRDRLDEV KEQVAEVRAK LEEQAQQIRL QAEAFQARLK SWFEPLVEDM QRQWAGLVEK VQAAVGTSAA PVPSDNH" # ApoE
        form = parse_peptide(test)
        monotest = monoisotop(form)
        self.assertEqual(printformula(form), "C_1569 H_2559 N_477 O_483 S_10")
        pepPTM = "A*CaDmEnF-GhHoIK+"    #that's an heavy PTM !
        formPTM = parse_peptide(pepPTM)
        self.assertEqual(printformula(formPTM), "C_47 H_69 N_12 O_20 P S")
        mH2 = parse_formula('H2').monoisotop()
        test1 = test[:103]
        test2 = test[103:]
        mt1 = monoisotop(parse_peptide(test1,ends="b"))
        mt2 = monoisotop(parse_peptide(test2,starts="y"))
        self.assertAlmostEqual(mt1+mt2, monotest+mH2)
        mt1 = monoisotop(parse_peptide(test1,ends="a"))
        mt2 = monoisotop(parse_peptide(test2,starts="x"))
        self.assertAlmostEqual(mt1+mt2, monotest)
        mt1 = monoisotop(parse_peptide(test1,ends="c"))
        mt2 = monoisotop(parse_peptide(test2,starts="z"))
        self.assertAlmostEqual(mt1+mt2, monotest+mH2)

    def test_insu(self):
        " test distribution on insuline convalent dimer "
        cc = parse_peptide("GIVEQCCASVCSLYQLENYCN")
        cd = parse_peptide("FVNQHLCGSHLVEALYLVCGERGFFYTPKA")
        #print "Chain, C", monoisotop(cc)
        #print "Chain, D", monoisotop(cd)
        insuline = cc + cd - parse_formula("H6")    # remove 6 H for 3 disulfides
        #insuline 5729.60086987
        self.assertAlmostEqual(monoisotop(insuline), 5729.60086987, 8)
        Ins = Distribution(insuline)
        self.assertEqual(Ins.len(), 23)
        ion10 = Ins.distrib[10]
        #print "Distribution"
        #print Ins
        self.assertAlmostEqual(ion10.mass, 5739.618326, 2)
        self.assertAlmostEqual(100*ion10.proba, 3.86529, 4)

# load table at import
def initialize():
    global elem_t, name_t, isotope_t 
    (elem_t, name_t, isotope_t ) = load_elements()

initialize()
import scipy.constants as cts
m_e = cts.physical_constants['electron mass in u'][0]    # electron mass in u ~0.00054

def demo0():
    print("Handling molecular formula")
    mol = "SH (CH2)11 (OCH2CH2)3 OCH2 COOH"
    form = parse_formula( mol )    # Check
    print("""
Mol: {0}
Formula: {1}
Monoisotopic mass: {2:.6f}
Average mass: {3:.3f}
""".format(mol, form, form.monoisotop(), form.average()))

def demo1():
    print("Handling proteins")
    prot = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG" # ubiquitine
    form = parse_peptide(prot)
    print("Ubiquitin:", prot)
    print("parsed as:",form)
    print ("average mass: %.3f  monoisotopic mass: %.3f"%(form.average(), form.monoisotop()))
    D = form.distribution()
    print ("Simplified Isotopic pattern:\n",D,"\n")
    D.draw(charge=6, RP=70000)
    D.draw_lowres(charge=6)
    plt.title("Ubiquitine 6+  RP=70k")

    plt.show()

def demo2():
    print("molecular arithmetics")
    form1 = parse_peptide("CCCC")
    form2 = parse_peptide("VVAVG")
    print("peptide P1: CCCC", form1)
    print("peptide P2: VVAVG", form2)
    # you can do formula arithmetics
    somme = form1 + 2*form2 - 2*parse_formula("H2O")
    print("ligation P2-P1-P2 - 2 H2O:", somme)
    print ("monoisotopic masses P1: %.3f  P2: %.3f  ligation:: %.3f"%(form1.monoisotop(), form2.monoisotop(), somme.monoisotop()))

def demo3():
    print("fine isotopic pattern")
    " test distribution on insuline covalent dimer "
    cc = parse_peptide("GIVEQCCASVCSLYQLENYCN")
    cd = parse_peptide("FVNQHLCGSHLVEALYLVCGERGFFYTPKA")
    #print "Chain, C", monoisotop(cc)
    #print "Chain, D", monoisotop(cd)
    insuline = cc + cd - parse_formula("H6")    # remove 6 H for 3 disulfides

    D1 = insuline.distribution()
    try:
        DD1 = fineisotopicdistrib(insuline)
    except:
        print("neutronstar utilities is not installed - check instructions")
        exit()
    # computing with neutronstar is very fast - drawing at high RP is slow !
    plt.figure(1,figsize=(8,3))
    D1.draw(RP=1E6, color='b', label="simplified pattern - RP=1E6")
    DD1.draw(RP=1E6, color='r', label="fine pattern - RP=1E6")
    DD1.bar(color='r')
    plt.legend()
    plt.title("comparing simple and detailled isotopic patterns for insuline")
    plt.figure(2,figsize=(8,3))
    for i,ion in enumerate([0,1,4,8]):
        plt.subplot(1,4,i+1)
        D1.draw(RP=1E6, color='b', label="ion %d"%ion)
        DD1.draw(RP=1E6, color='r')
        DD1.bar(color='r')
        plt.legend()
        plt.xlim(xmin=D1.distrib[ion].mass-0.03, xmax=D1.distrib[ion].mass+0.03)
    plt.title("zooming on peaks")
    plt.show()

def demo4():
    print("Not only organic chemistry !")
    print("Isotopic distribution for Palladium:")
    print(parse_formula('Pd').distribution())
    print("monoisotopic mass of atoms is defined as the most naturally abundant, hence 106Pd here")
    Pdalloy = "Au2 Ca10 Ga10 Pd76"
    print("Palladium alloy:", Pdalloy)
    alloy = parse(Pdalloy)
    print("parsed formula",alloy)
    print ("average mass: %.3f  monoisotopic mass: %.3f"%(alloy.average(), alloy.monoisotop()))

    alloy.distribution().bar()

if __name__ == '__main__':
    unittest.main()
    # demo0()
    # demo1()
    # demo2()
    # demo3()
    # demo4()
