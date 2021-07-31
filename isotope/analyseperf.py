#!/usr/bin/env python
# encoding: utf-8
"""Analyser of performance of some main functions of isotopes """
import time
import random
#~ import numpy as np
import matplotlib.pyplot as plt

import isotopes as iso

def get_xy(couple):
    """from a list of couple return a couple of list """
    x = []
    y = [] 
    for (z, k) in couple:
        x.append(z)
        y.append(k)
    return (x, y)


class Analyser(object):
    """ This classe is used to analyse the performances of some main
    function of isotopes """


    def __init__(self, lsizes=[2 ** i for i in range(1, 10)]):
        self.listofsize = lsizes
        self.seq = ""
        self.forms = []
        self.average_time = []
        self.monoisotop_time = []
        self.distribution_time = []
        self.ppeptide_time = []
        

    def generate_sequence(self, size):
        """ Generate a random sequence of size AAs """
        aas = "ACDEFGHIKLMNPQRSTVWY"
        self.seq = ""
        random.seed()
        for i in range(size):
            self.seq += aas[random.randrange(0, len(aas))]
        return self.seq


    def analyse_parse_peptide(self):
        """ Return a list of couple of (size, time) for the function
        parse_peptide """
        listoftime = []
        for s in self.listofsize:
            self.generate_sequence(s)
            start = time.time()
            self.forms.append(iso.parse_peptide(self.seq))
            realtime = time.time()-start
            listoftime.append(realtime)
        self.ppeptide_time = list(zip(self.listofsize, listoftime))

    def analyse_standard(self, func):
        """ Return a list of couple of (size, time) for a generic function
        with formula in params"""
        listoftime = []
        for formula in self.forms:
            start = time.time()
            func(formula)
            realtime = time.time()-start
            listoftime.append(realtime)
        res = list(zip(self.listofsize, listoftime))
        return res


    def analyse_monoisotopic(self):
        """ Return a list of couple of (size, time) for the monoisotopic
        function"""
        self.monoisotop_time = self.analyse_standard(iso.monoisotop)

    def analyse_average(self):
        """ Return a list of couple of (size, time) for the average
        function"""
        self.average_time = self.analyse_standard(iso.average)

    def analyse_distribution(self):
        """ Return a list of couple of (size, time) for the init of
        distribution"""
        self.distribution_time = self.analyse_standard(iso.Distribution)

    def analyse_function(self, name):
        """ analyse a function if it can be analyse """
        if name == "parse_peptide":
            self.analyse_parse_peptide()
        elif name == "monoisotop":
            if len(self.forms) == 0:
                raise Exception\
                ("Formulas need to be compute before monoisotopic function")
            self.analyse_monoisotopic()
        elif name == "average":
            if len(self.forms) == 0:
                raise Exception\
                ("Formulas need to be compute before average function")
            self.analyse_average()
        elif name == "distribution":
            if len(self.forms) == 0:
                raise Exception\
                ("Formulas need to be compute before Distribution function")
            self.analyse_distribution()

    def analyse_all(self):
        """ analyse all functions """
        self.analyse_parse_peptide()
        self.analyse_monoisotopic()
        self.analyse_average()
        self.analyse_distribution()
        #~ print self.ppeptide_time
        #~ print self.monoisotop_time
        #~ print self.average_time
        #~ print self.distribution_time

    def get_all_xy(self):
        """ return the list of couple of list """
        (x, y) = get_xy(self.ppeptide_time)
        ppep, = plt.loglog(x, y, 'co-', basex=2, basey=2)
        (x, y) = get_xy(self.monoisotop_time)
        monoi, = plt.loglog(x, y, 'go-', basex=2, basey=2)
        (x, y) = get_xy(self.average_time)
        aver, = plt.loglog(x, y, 'bo-', basex=2, basey=2)
        (x, y) = get_xy(self.distribution_time)
        dist, = plt.loglog(x, y, 'ro-', basex=2, basey=2)
        return [ppep, monoi, aver, dist]

    def create_graph(self):
        """ create the loglog graph with a curve for each function """
        plt.title("loglog")
        plt.xlabel("sizesequence")
        plt.ylabel("time")
        curvelist = self.get_all_xy()
        plt.legend(curvelist, ["parse_peptide", "monoisotop",
        "average", "distribution"], bbox_to_anchor=(0., 1.1, 0.2, 0.))
        plt.show()



if __name__ == '__main__':
    TEST = Analyser()
    TEST.analyse_all()
    TEST.create_graph()
