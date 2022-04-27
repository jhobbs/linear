import sys
import numpy as np
from sympy import Rational, Matrix, zeros, oo
from typing import List, Dict, Tuple


class LinearSystemOfEquations:
    def __init__(self, raw_rows):
        self._matrix = self.rows_to_matrix(raw_rows)

    def rows_to_matrix(self, input_rows):
        rows = []
        constants = []
        for row in input_rows:
            row_variables, constant = self.parse_raw_row(row)
            rows.append(row_variables)
            constants.append(constant)
        var_col = self.variables_to_cols(rows)
        matrix = zeros(len(rows), len(var_col) + 1)
        for i, row in enumerate(rows):
            for variable, coefficient in row.items():
                matrix[i, var_col[variable]] = coefficient
            matrix[i, -1] = constants[i]
        return matrix

    @staticmethod
    def parse_term(term: str):
        if not term[-1].isalpha():
            raise Exception(f"Invalid term: {term}")

        variable = term[-1]
        sign = term[0]
        if term[0].isdigit():
            coefficient = Rational(term[0:-1])
        elif len(term[1:-1]) == 0:
            coefficient = Rational(1)
        else:
            coefficient = Rational(term[1:-1])

        return (
            coefficient,
            variable,
        )

    def parse_variable_terms(self, variable_part: str):
        split_row = variable_part.split()
        terms = {}
        for raw_term in split_row:
            coefficient, variable = self.parse_term(raw_term)
            if variable in terms:
                terms[variable] += coefficient
            else:
                terms[variable] = coefficient
        return terms

    def parse_raw_row(self, raw_row: str):
        variable_part, constant_part = raw_row.split("=")
        terms = self.parse_variable_terms(variable_part)
        constant = Rational(int(constant_part.strip()))
        return terms, constant

    @staticmethod
    def variables_to_cols(rows: List[dict]):
        # collect a set of variables
        # sort them and number them by order
        variable_set = set()
        for row in rows:
            for variable in row.keys():
                variable_set.add(variable)

        return {
            variable: position for position, variable in enumerate(sorted(variable_set))
        }

    def is_reduced_row_echelon(self):
        return self._inner_is_reduced_row_echelon(self._matrix)

    def _inner_is_reduced_row_echelon(self, matrix):
        if len(matrix) == 0:
            return True

        for i in range(len(matrix[0])):
            if matrix[0, i] == 0:
                if not np.all(matrix[1:, i] == 0):
                    return False
            elif matrix[0, i] == 1:
                if not np.all(matrix[1:, i] == 0):
                    return False
                return self._inner_is_reduced_row_echelon(matrix[1:, i:])
            else:
                return False
        return True

    @staticmethod
    def compare_rows(row_a, row_b):
        for i in range(len(row_a)):
            if row_a[i] == row_b[i]:
                continue
            if row_a[i] == 0 or (abs(row_a[i]) < abs(row_b[i])):
                return -1
            return 1
        return 0

    @staticmethod
    def swap_rows(matrix, row1, row2):
        return matrix.elementary_row_op("n<->m", row1=row1, row2=row2)

    def sort_rows(self, matrix):
        # bubble sort, but 0's are less than negatives
        while True:
            swapped = False
            for i in range(matrix.rows):
                if i == matrix.rows - 1:
                    break
                if self.compare_rows(matrix.row(i), matrix.row(i + 1)) == -1:
                    matrix = self.swap_rows(matrix, i, i + 1)
                    swapped = True

            if swapped == False:
                break
        return matrix

    def reduce_rows(self):
        self._matrix = self._inner_reduce_rows(self.sort_rows(self._matrix))
        self.back_substitute()

    def _inner_reduce_rows(self, matrix):
        if matrix.is_zero_matrix:
            return

        leading_coeff = matrix[0, 0]
        if leading_coeff == 0:
            raise Exception(f"Unexpected 0 leading coefficient: {i} {matrix}")
        matrix[0, :] = matrix[0, :] * 1 / leading_coeff
        for i in range(matrix.rows - 1):
            row_leading_coeff = matrix[i + 1, 0]
            if row_leading_coeff == 0:
                continue
            matrix[i + 1, :] = matrix[i + 1, :] + (-row_leading_coeff * matrix[0, :])
        reduced_submatrix = self._inner_reduce_rows(matrix[1:, 1:])
        if reduced_submatrix is not None:
            matrix[1:, 1:] = reduced_submatrix
        return matrix

    def number_of_solutions(self):
        non_zero_rows = 0
        for i in range(self._matrix.rows):
            # check for contradiction
            if self._matrix[i, :-1].is_zero_matrix and self._matrix[i, -1] != 0:
                return 0

            if not self._matrix[i, :].is_zero_matrix:
                non_zero_rows += 1
        if non_zero_rows == len(self._matrix[i, :]) - 1:
            return 1
        return oo

    def back_substitute(self):
        for i in range(self._matrix.rows - 1, 0, -1):
            for j in range(len(self._matrix[i, :-1])):
                if self._matrix[i, j] != 1:
                    continue
                for i2 in range(i):
                    if self._matrix[i2, j] != 0:
                        self._matrix[i2, :] = (
                            self._matrix[i2, :]
                            + -self._matrix[i2, j] * self._matrix[i, :]
                        )

    @classmethod
    def from_file(cls, path):
        with open(path) as matrix_file:
            lines = matrix_file.readlines()
        return cls(lines)


def main():
    system = LinearSystemOfEquations.from_file(sys.argv[1])
    print(system._matrix)
    system.reduce_rows()
    print(system._matrix)
    print(system.number_of_solutions())


if __name__ == "__main__":
    main()
