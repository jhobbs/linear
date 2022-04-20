from unittest import TestCase
import numpy as np

from gaussian import is_reduced_row_echelon


class TestIsReducedRowEchelon(TestCase):
    def test_empty_array(self):
        array = np.ndarray((0,))
        self.assertTrue(is_reduced_row_echelon(array))

    def test_one_row(self):
        array = np.asarray([[0,1,2]])
        self.assertTrue(is_reduced_row_echelon(array))

    def test_two_rows(self):
        array = np.asarray([[0,1,2], [0,0,1]])
        self.assertTrue(is_reduced_row_echelon(array))
    
    def test_one_row_not_reduced(self):
        array = np.asarray([[0,2,2]])
        self.assertFalse(is_reduced_row_echelon(array))
    
    def test_two_rows_not_reduced(self):
        array = np.asarray([[1,0,0],[0,2,2]])
        self.assertFalse(is_reduced_row_echelon(array))
    
    def test_two_rows_wrong_order(self):
        array = np.asarray([[0,1,0],[1,0,0]])
        self.assertFalse(is_reduced_row_echelon(array))