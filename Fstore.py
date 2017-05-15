"""
Utilities to handle distribution of isotopic deconvolution

"""
from __future__ import division, print_function
import copy
import isotopes as iso
from proteins import averagine

class Fstore(object):
    """
    a store for isotopic distributions from formulas
    can be speeded up by implmenting a sorted list, (would be in O(log n) rather that in O(n))
    see
    https://code.activestate.com/recipes/577197-sortedcollection/
    """
    def __init__(self):
        "initialize with a base formula ~around 990 dalton"
        self.Memo = {}
#        self.base = iso.parse_formula("C40 H60 N10 O10 S") - iso.parse_formula("S")
        # with this trick, base contains S_0 !
#        self.set( self.base, self.base.distribution() )
    def hash(self,form):
        "a hash for storing formulas "
        return tuple(form.items())
    def hashtoform(self,h):
        " back from hash to formula"
        f = iso.Formula()
        for i,j in h:
            f[i] = j
        return f
    def content(self):
        return [self.hashtoform(f) for f in self.Memo.keys()]
    def set(self, form, dist):
        "setting in the store"
        self.Memo[self.hash(form)] = dist
    def get(self, form):
        "getting from the store"
        return copy.copy(self.Memo[self.hash(form)])
    def closest(self,form):
        "given a formula, returns the closest (smaller) formula in the store"
        bestsofar = sum( form.values() )  # number of atoms
        sol = None
        for h in self.Memo.keys():
            f = self.hashtoform(h)
            diff = [form[a]-f[a] for a in ('C','H','N','O','S')]  # diff of nb of atoms
            if all( [i>=0 for i in diff] ):  # if all pos
                sdiff = sum(diff)         # distance in nb of atoms
                if sdiff<bestsofar:       # a new best
                    bestsofar = sdiff
                    sol = f
                    if bestsofar<2:
                        break
        return sol
    def distrib(self, form):
        "compute a distribution for form"
        try:
            b = self.get(form)
        except KeyError:
            b = None
        if b is None:
            close = self.closest(form)
            if close is None:
                print(form)
                b = form.distribution()
            else:
                print(close)
                mdif = form-close   # formula arithmetics
                try:
                    b = self.get(close)  # stored distribution
                except KeyError:    # should never happen !!!!  A VERIFIER   MAD
                    b = close.distribution()
                    print("***")
                mddist = mdif.distribution()    # distrib to add
                b.combine(mddist)   # combine() adds mddist to b
            self.set(form,b)
        return b
    def pattern(self, mass):
        """
        isotopic pattern definition, estimated from mass, using averagine hypothesis

        returns a isotopic pattern for a given protein mass
        return (deltaMass, intensities)
        """
        molecule = averagine(mass)
        # distribution de la molecule
        U = self.distrib(molecule)
        U.sort_by_mass()
        # determiner le signal de distribution
        deltaM = [ion.mass-U.distrib[0].mass for ion in U.distrib]
        S = sum([I.proba for I in U.distrib])
        intens = [ion.proba/S for ion in U.distrib]
        return (deltaM, intens)

if __name__ == '__main__':
    FF = Fstore()
    mm = iso.parse_formula('C46 H70 O12 N12 S')
    print (FF.distrib(mm))
    print ("\n",mm.distribution())