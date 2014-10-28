# isotopes

This module allows the precise computation of isotopic distributions, average mass, and monoisotopic masses
of arbitrary chemical formula.

The handy determination from peptide and protein primary sequences is provided.

The possibility to draw the isotopic profil is also given.

Isotopic natural abundance can be changed with the enrich() function.

Typical use is :
molecule = "K23 I22 S30"
formula = parse_formula(molecule)   # formula object (in fact a dictionary)
print "average mass", formula.average()
print "monoisotopic mass", formula.monoisotop()
distrib = formula.distribution()    # distribution object
print distrib
distrib.draw()


adapted from 
Kubinyi, H.
Calculation of isotope distributions in mass spectrometry. A trivial solution for a non-trivial problem.
Anal Chim Acta 247, 107-119 (1991).

This module depends at import time on a isotope file to be present in the same directory
by default called "elements.asc"

It is a copy of the file found at
NIST | Physical Measurement Laboratory | Physical Reference Data | Atomic Weights and Isotopic Compositions Main Page
http://www.nist.gov/

First version of algo by FX Coudert, Python rewrite by DELSUC Marc-Andre on 2014-03-28.
Copyright (c) 2014 CNRS. All rights reserved.