#!/usr/bin/env python
# encoding: utf-8
import time
import random
import numpy as np 
import matplotlib.pyplot as plt

import isotopes as iso

"""Analyser of performance of some main functions of isotopes """

class Analyser(object):
    """ This classe is used to analyse the performances of some main
     function of isotopes """
    
    
    def __init__(self, lsizes = [2 ** i for i in range(1,10)]):
        self.listofsize = lsizes
        self.seq = ""
        self.forms = []
    
    """ Generate a random sequence of size AAs """ 
    def generateSequence(self, size):
        aas = "ACDEFGHIKLMNPQRSTVWY"
        self.seq = ""
        random.seed()
        for i in range(size):
            self.seq += aas[random.randrange(0,len(aas))]
        return self.seq
    
    """ Return a list of couple of (size, time) for the function parse_peptide """
    def analyseParsePeptide(self):
        listoftime = []
        for s in self.listofsize:
            self.generateSequence(s)
            start = time.clock()
            self.forms.append(iso.parse_peptide(self.seq))
            realtime = time.clock()-start
            listoftime.append(realtime)
        self.ppeptide_time = zip(self.listofsize, listoftime)
    
    """ Return a list of couple of (size, time) for a generic function with formula in params"""
    def analyseStandard(self, func):
        listoftime = []
        for formula in self.forms:
            start = time.clock()
            func(formula)
            realtime = time.clock()-start
            listoftime.append(realtime)
        res = zip(self.listofsize, listoftime)
        return res
        
    """ Return a list of couple of (size, time) for the monoisotopic function"""
    def analyseMonoisotopic(self):
        self.monoisotop_time = self.analyseStandard(iso.monoisotop)
    
    """ Return a list of couple of (size, time) for the average function"""    
    def analyseAverage(self):
        self.average_time =self.analyseStandard(iso.average)
    
    """ Return a list of couple of (size, time) for the init of distribution"""
    def analyseDistribution(self):
        self.distribution_time = self.analyseStandard(iso.Distribution)
    
    """ analyse a function if it can be analyse """
    def analyseFunction(self, name):
        if name == "parse_peptide":
            self.analyseParsePeptide()
        elif name == "monoisotop":
            if len(self.forms) == 0:
                raise(Exception\
                ("Formulas need to be compute before monoisotopic function"))
            self.analyseMonoisotopic()
        elif name == "average":
            if len(self.forms) == 0:
                raise(Exception\
                ("Formulas need to be compute before average function"))
            self.analyseAverage()
        elif name == "distribution" :
            if len(self.forms) == 0:
                raise(Exception\
                ("Formulas need to be compute before Distribution function"))
            self.analyseDistribution()
    
    """ analyse all functions """
    def analyseAll(self):
        self.analyseParsePeptide()
        self.analyseMonoisotopic()
        self.analyseAverage()
        self.analyseDistribution()
        #~ print self.ppeptide_time
        #~ print self.monoisotop_time
        #~ print self.average_time
        #~ print self.distribution_time
    
    """from a list of couple return a couple of list """
    def getXY(self, couple):
        x = []
        y = []
        for (seq , time) in couple:
            x.append(seq)
            y.append(time)
        return (x, y)
    
    """ return the list of couple of list """
    def getAllXY(self):
        (x, y) = self.getXY(self.ppeptide_time)
        ppep, = plt.loglog(x,y,'co-',basex=2,basey=2)
        (x, y) = self.getXY(self.monoisotop_time)
        monoi, = plt.loglog(x,y,'go-',basex=2,basey=2)
        (x, y) = self.getXY(self.average_time)
        aver, = plt.loglog(x,y,'bo-',basex=2,basey=2)
        (x, y) = self.getXY(self.distribution_time)
        dist, = plt.loglog(x,y,'ro-',basex=2,basey=2)
        return [ppep, monoi, aver, dist]
        
    """ create the loglog graph with a curve for each function """
    def createGraph(self):
        plt.title("loglog")
        plt.xlabel("sizesequence")
        plt.ylabel("time")
        curvelist = self.getAllXY()
        plt.legend(curvelist, ["parse_peptide", "monoisotop", 
        "average", "distribution"], bbox_to_anchor=(0., 1.1, 0.2, 0.))
        plt.show()
        
        
        
if __name__ == '__main__':
    truc = Analyser()
    truc.analyseAll()
    truc.createGraph()
    
    
   
