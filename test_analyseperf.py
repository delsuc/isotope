import unittest
import numpy.testing as nt

import analyseperf

#~ lauch all test "python -m unittest discover" 

class TestAnalysePerf(unittest.TestCase):
    
    def setUp(self):
        self.analyse = analyseperf.Analyser()
    
    def tearDown(self): #do the setup to each test
        pass
    
    def test_init(self):
        self.assertGreater(len(self.analyse.listofsize), 0)
        
    def test_01_generateSequence(self):
        self.analyse.generateSequence(20)
        nt.assert_equal(len(self.analyse.seq), 20)
    
    def test_02_analyseParsePeptide(self):
        self.analyse.analyseParsePeptide()
        nt.assert_equal(len(self.analyse.ppeptide_time),len(self.analyse.listofsize))
        for (x, y) in self.analyse.ppeptide_time:
            self.assertIsNotNone(x)
            self.assertIsNotNone(y)
        nt.assert_equal(len(self.analyse.forms), len(self.analyse.listofsize))
    
    def test_03_analyseMonoisotopic(self):
        self.analyse.analyseParsePeptide()
        self.analyse.analyseMonoisotopic()
        nt.assert_equal(len(self.analyse.monoisotop_time), len(self.analyse.listofsize))
        for value in self.analyse.monoisotop_time:
            self.assertIsNotNone(value)
    
    def test_03_analyseAverage(self):
        self.analyse.analyseParsePeptide()
        self.analyse.analyseAverage()
        nt.assert_equal(len(self.analyse.average_time), len(self.analyse.listofsize))
        for value in self.analyse.average_time:
            self.assertIsNotNone(value)
    
    def test_03_analyseDistribution(self):
        self.analyse.analyseParsePeptide()
        self.analyse.analyseDistribution()
        nt.assert_equal(len(self.analyse.distribution_time), len(self.analyse.listofsize))
        for value in self.analyse.distribution_time:
            self.assertIsNotNone(value)
    
    def test_03_getXY(self):
        couple = [ (1, 1), (1, 1), (1, 1)]
        (x,y) = self.analyse.getXY(couple)
        nt.assert_equal(len(x), len(couple))
        nt.assert_equal(len(y), len(couple))
        

if __name__ == '__main__':
    unittest.main()
