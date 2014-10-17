#!/usr/bin/env python
# encoding: utf-8
"""
proteins.py

Tools for manipulating peptides and proteins in MS



Created by DELSUC Marc-Andre on 2014-04-06.
Copyright (c) 2014 IGBMC. All rights reserved.
"""

import numpy as np
import matplotlib.pyplot as plt
import isotopes as iso

def averagine(mass):
    """
    given an average mass, returns the averagine formula closest to this mass
    
    from - Senko 1995 JASMS
    C4.9384 H7.7583 N1.3577 O1.4773 S0.0417 
    
    i.e.
    >>> print averagine(20000), iso.average(averagine(20000))
    C_889 H_1395 N_243 O_266 S_8  19999.6104529
    
    """
    averagine = {"C":4.9384, "H":7.7583, "N":1.3577, "O":1.4773, "S":0.0417 }
    if mass < 1000:
        raise(Exception("mass should be larger than 1000"))
    ratio = mass/111.123648381
    f = iso.Formula()
    for e,m in averagine.items():
        f[e] = int(round(ratio*m))
    def correc(add, elem, val):
        """
        used locally to add o substract a few atom to bring mass within 1 dalton
        adds or substract "add" "elem" to running formula if mass is away by more than "val"
        """
        mf = iso.average(f)
        delta = mf-mass
        if delta<-val:
            f[elem] += add
            mf += val
        elif delta>val:
            f[elem] -= add
            mf -= val
    # the following loop allow some tuning of the formula to get a mass within less than 1 dalton to the target.
    # this will try to change by +/- 1 each heavy atom and by +/- 1 to 11 for hydrogens
    for (a,e,v) in ((1,"S",34),(1,"O",16),(1,"N",14),(1,"C",12),(8,"H",7.5),(4,"H",3.5),(2,"H",1.5),(1,"H",0.5)):
        correc(a,e,v)   # the ".5" allow to be within less than 0.5 in nearly all cases (max is ~0.51)
    return f

def test():
    N = 10000
    mm = np.linspace(1000,50000,N)
    ave = np.zeros(N)
    mono = np.zeros(N)
    for i in range(N):
        a = averagine(mm[i])
        ave[i] = iso.average(a)
        mono[i] = iso.monoisotop(a)
#    plt.plot(mm,ave-mm)
    plt.plot (mm,ave-mono)
#    plt.plot(mm,mono)
    print max(abs(ave-mm)), (ave-mm).mean(), (ave-mm).std()
    plt.show()
    
def test2():
    fb = averagine(13000.0)
    iso.Distribution(fb).draw()
    iso.enrich("N", 15, 0.50)
    iso.Distribution(fb).draw()
    iso.enrich("N", 15, 0.65)
    iso.Distribution(fb).draw()
    iso.Distribution(fb).enveloppe()
    plt.show()
if __name__ == '__main__':
    test2()