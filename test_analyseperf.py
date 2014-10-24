import unittest

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
        self.assertEqual(len(self.analyse.seq), 20)
    
    def test_02_analyseParsePeptide(self):
        self.analyse.analyseParsePeptide()
        self.assertEqual(len(self.analyse.ppeptide_time),len(self.analyse.listofsize))
        for (x, y) in self.analyse.ppeptide_time:
            self.assertIsNotNone(x)
            self.assertIsNotNone(y)
        self.assertEqual(len(self.analyse.listofsize), len(self.analyse.forms))
    
    def test_03_analyseMonoisotopic(self):
        self.analyse.analyseParsePeptide()
        self.analyse.analyseMonoisotopic()
        self.assertEqual(len(self.analyse.listofsize), len(self.analyse.monoisotop_time))
        for value in self.analyse.monoisotop_time:
            self.assertIsNotNone(value)
    
    def test_03_analyseAverage(self):
        self.analyse.analyseParsePeptide()
        self.analyse.analyseAverage()
        self.assertEqual(len(self.analyse.listofsize), len(self.analyse.average_time))
        for value in self.analyse.average_time:
            self.assertIsNotNone(value)
    
    def test_03_analyseDistribution(self):
        self.analyse.analyseParsePeptide()
        self.analyse.analyseDistribution()
        self.assertEqual(len(self.analyse.listofsize), len(self.analyse.distribution_time))
        for value in self.analyse.distribution_time:
            self.assertIsNotNone(value)
    
    def test_03_getXY(self):
        couple = [ (1, 1), (1, 1), (1, 1)]
        (x,y) = self.analyse.getXY(couple)
        self.assertEqual(len(x),len(y))
        self.assertEqual(len(couple), len(x))

if __name__ == '__main__':
    unittest.main()
