#!/usr/bin/env python
# encoding: utf-8
"""Tests with unittest"""
import unittest
#~ import numpy.testing as nt
from numpy.testing import assert_equal

import analyseperf

#~ lauch all test "python -m unittest discover" 

class TestAnalysePerf(unittest.TestCase):
    """This class is used to do unit tests on the analyser
    test has done in the alphabetical order"""
    
    def setUp(self):
        """To do before each tests"""
        self.analyse = analyseperf.Analyser()
    
    def tearDown(self):
        """Clear after each tests"""
        pass
    
    def test_00_init(self):
        """test the init"""
        self.assertGreater(len(self.analyse.listofsize), 0)
        
    def test_01_generate_sequence(self):
        """test the sequence"""
        self.analyse.generate_sequence(20)
        assert_equal(len(self.analyse.seq), 20)
    
    def test_02_analyse_parse_peptide(self):
        """test the analyse_parse_peptide"""
        self.analyse.analyse_parse_peptide()
        assert_equal(len(self.analyse.ppeptide_time),\
        len(self.analyse.listofsize))
        for (x, y) in self.analyse.ppeptide_time:
            self.assertIsNotNone(x)
            self.assertIsNotNone(y)
        assert_equal(len(self.analyse.forms),\
        len(self.analyse.listofsize))
    
    def test_03_analyse_monoisotopic(self):
        """test the analyse_monoisotopic"""
        self.analyse.analyse_parse_peptide()
        self.analyse.analyse_monoisotopic()
        assert_equal(len(self.analyse.monoisotop_time),\
         len(self.analyse.listofsize))
        for value in self.analyse.monoisotop_time:
            self.assertIsNotNone(value)
    
    def test_03_analyse_average(self):
        """test the analyse_average"""
        self.analyse.analyse_parse_peptide()
        self.analyse.analyse_average()
        assert_equal(len(self.analyse.average_time),\
         len(self.analyse.listofsize))
        for value in self.analyse.average_time:
            self.assertIsNotNone(value)
    
    def test_03_analyse_distribution(self):
        """test the analyse_distribution"""
        self.analyse.analyse_parse_peptide()
        self.analyse.analyse_distribution()
        assert_equal(len(self.analyse.distribution_time),\
         len(self.analyse.listofsize))
        for value in self.analyse.distribution_time:
            self.assertIsNotNone(value)
    
    def test_03_get_xy(self):
        """test the get_xy"""
        couple = [(1, 1), (1, 1), (1, 1)]
        (x, y) = analyseperf.get_xy(couple)
        assert_equal(len(x), len(couple))
        assert_equal(len(y), len(couple))
        

if __name__ == '__main__':
    unittest.main()
