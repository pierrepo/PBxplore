#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for PBxplore.

Tests functions from different programs.

2014 - P. Poulain
"""

#===============================================================================
# load modules
#===============================================================================
import unittest
import collections

import PBlib as PB


#===============================================================================
# classes for tests
#===============================================================================
class TestPBlib(unittest.TestCase):
    """
    Tests for PBlib
    """

    def test_get_dihedral(self):
        """
        Test for get_dihedral()
        """
        Result = collections.namedtuple('Result', ['A', 'B', 'C', 'D', 'torsion'])
        results = (Result((-7.28, -9.262, 5.077), 
                          (-7.526, -10.643, 5.529),
                          (-6.221, -11.438, 5.555),
                          (-6.289, -12.685, 5.931),
                           -179.663656153),
                   Result((-1.373, -8.817, -4.389), 
                          (-1.203, -8.335, -5.792),
                          (-1.891, -6.977, -5.927),
                          (-1.918, -6.429, -7.107),
                           -176.048770127),
                   Result((-0.533, -8.42, -3.47  ), 
                          (-1.373, -8.817, -4.389),
                          (-1.203, -8.335, -5.792),
                          (-1.891, -6.977, -5.927),
                           -84.8356057692),
                   Result((-1.918, -6.429, -7.107), 
                          (-2.609, -5.125, -7.305),
                          (-4.108, -5.392, -7.331),
                          (-4.469, -6.494, -7.911),
                           -36.8942888266),
                   Result((-11.285, 6.472, -7.44 ), 
                          (-12.62, 5.829, -7.425 ),
                          (-13.585, 6.626, -6.544),
                          (-13.098, 7.621, -5.858),
                           -6.58786169376),
                   Result((-11.284, -0.971, -2.679), 
                          (-12.65, -0.794, -3.226),
                          (-13.665, -1.664, -2.479),
                          (-13.262, -2.363, -1.452),
                           3.91626706556),
                   Result((-2.004, -10.892, -2.611), 
                          (-1.87, -9.835, -1.853),
                          (-0.726, -8.877, -2.011),
                          (-0.533, -8.42, -3.47),
                           50.065196067),
                   Result((11.174, -6.725, 0.458), 
                          (10.732, -7.258, -0.86),
                          (9.27, -6.869, -1.096),
                          (8.741, -7.185, -2.245),
                           175.872397707))
        
        for res in results:
            torsion = PB.get_dihedral(res.A, res.B, res.C, res.D)
            self.assertAlmostEqual(torsion, res.torsion)

    def test_read_fasta(self):
        headers, sequences = PB.read_fasta("test_data/1BTA.PB.fasta")
        self.assertEqual(headers, ['test_data/1BTA.pdb | chain A'])
        self.assertEqual(sequences, ['ZZdddfklonbfklmmmmmmmmnopafklnoiakl'
                                     'mmmmmnoopacddddddehkllmmmmngoilmmmm'
                                     'mmmmmmmmnopacdcddZZ'])

if __name__ == '__main__':
    unittest.main()


