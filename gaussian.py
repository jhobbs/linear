import sys
import numpy as np
from sympy import Rational, Matrix, symbols, zeros, oo, FiniteSet, init_printing
from typing import List, Dict, Tuple


class LinearSystemOfEquationsParser:
    @classmethod
    def rows_to_matrix(cls, input_rows):
        rows = []
        constants = []
        if input_rows[0].startswith("order:"):
            operation_rows = input_rows[1:]
            col_var = input_rows[0][6:].split()
            var_col = {variable: col for col, variable in enumerate(col_var)}
        else:
            operation_rows = input_rows
            col_var = None
            var_col = None

        for row in operation_rows:
            row_variables, constant = cls.parse_raw_row(row)
            rows.append(row_variables)
            constants.append(constant)

        if var_col is None:
            var_col = cls.variables_to_cols(rows)
            col_var = cls.cols_to_vars(var_col)

        matrix = zeros(len(rows), len(var_col) + 1)
        for i, row in enumerate(rows):
            for variable, coefficient in row.items():
                matrix[i, var_col[variable]] = coefficient
            matrix[i, -1] = constants[i]
        return matrix, symbols(col_var)

    @classmethod
    def parse_raw_row(cls, raw_row: str):
        variable_part, constant_part = raw_row.split("=")
        terms = cls.parse_variable_terms(variable_part)
        try:
            constant = Rational(int(constant_part))
        except ValueError:
            constant = symbols(constant_part)

        return terms, constant

    @classmethod
    def parse_variable_terms(cls, variable_part: str):
        split_row = variable_part.split()
        terms = {}
        for raw_term in split_row:
            coefficient, variable = cls.parse_term(raw_term)
            if variable in terms:
                terms[variable] += coefficient
            else:
                terms[variable] = coefficient
        return terms

    @staticmethod
    def parse_term(term: str):
        if not term[-1].isalpha():
            raise Exception(f"Invalid term: {term}")

        variable = term[-1]
        if term[0].isdigit():
            coefficient = Rational(term[0:-1])
        else:
            if len(term[1:-1]) == 0:
                coefficient = Rational(1)
            else:
                coefficient = Rational(term[1:-1])
            sign = term[0]
            if sign == "-":
                coefficient = coefficient * -1

        return (
            coefficient,
            variable,
        )

    @staticmethod
    def variables_to_cols(rows: List[dict]):
        # collect a set of variables
        # sort them and number them by order
        variable_set = set()
        for row in rows:
            for variable in row.keys():
                variable_set.add(variable)
        mapping = {
            variable: position for position, variable in enumerate(sorted(variable_set))
        }
        return mapping

    @staticmethod
    def cols_to_vars(vars_to_cols):
        result = [
            variable
            for variable, _ in list(sorted(vars_to_cols.items(), key=lambda x: x[1]))
        ]
        return result


class LinearSystemOfEquations:
    def __init__(self, matrix, column_variables):
        self._matrix = matrix
        self._column_variables = column_variables
        self._echelon_matrix = self._gaussian_elimination(
            self._sort_rows(self._matrix.copy())
        )
        self._normalized_echelon_matrix = self._normalize_rows(
            self._echelon_matrix.copy()
        )
        self._reduced_matrix = self._reduce_rows(self._normalized_echelon_matrix.copy())

    @staticmethod
    def _compare_rows(row_a, row_b):
        for i in range(len(row_a)):
            if row_a[i] == row_b[i]:
                continue
            if row_a[i] == 0 or (abs(row_a[i]) < abs(row_b[i])):
                return -1
            return 1
        return 0

    @staticmethod
    def _swap_rows(matrix, row1, row2):
        return matrix.elementary_row_op("n<->m", row1=row1, row2=row2)

    def _sort_rows(self, matrix):
        # bubble sort, but 0's are less than negatives
        while True:
            swapped = False
            for i in range(matrix.rows):
                if i == matrix.rows - 1:
                    break
                if self._compare_rows(matrix.row(i), matrix.row(i + 1)) == -1:
                    matrix = self._swap_rows(matrix, i, i + 1)
                    swapped = True

            if swapped == False:
                break
        return matrix

    def _reduce_rows(self, matrix):
        """Convert matrix to reduced row echelon form."""
        self._back_substitute(matrix)
        return matrix

    def _normalize_rows(self, matrix):
        for i in range(matrix.rows):
            for j in range(matrix.cols):
                if matrix[i, j] != 0:
                    matrix[i, :] = (1 / matrix[i, j]) * matrix[i, :]
                    break
        return matrix

    def _gaussian_elimination(self, matrix):
        """Convert matrix to echelon form."""
        if matrix[:, :-1].is_zero_matrix:
            return matrix

        leading_coeff = matrix[0, 0]
        if leading_coeff == 0:
            raise Exception(f"Unexpected 0 leading coefficient: {i} {matrix}")
        for i in range(matrix.rows - 1):
            row_leading_coeff = matrix[i + 1, 0]
            if row_leading_coeff == 0:
                continue
            matrix[i + 1, :] = matrix[i + 1, :] + (
                -(row_leading_coeff / leading_coeff) * matrix[0, :]
            )
        reduced_submatrix = self._gaussian_elimination(matrix[1:, 1:])
        if reduced_submatrix is not None:
            matrix[1:, 1:] = reduced_submatrix
        return matrix

    def _back_substitute(self, matrix):
        for i in range(matrix.rows - 1, 0, -1):
            for j in range(len(matrix[i, :-1])):
                if matrix[i, j] == 1:
                    for i2 in range(i):
                        if matrix[i2, j] != 0:
                            matrix[i2, :] = (
                                matrix[i2, :] + -matrix[i2, j] * matrix[i, :]
                            )
                    break
        return matrix

    def number_of_solutions(self):
        """Return the number of solutions to the system."""
        non_zero_rows = 0
        for i in range(self._reduced_matrix.rows):
            # check for contradiction
            if (
                self._reduced_matrix[i, :-1].is_zero_matrix
                and self._reduced_matrix[i, -1] != 0
            ):
                return 0

            if not self._reduced_matrix[i, :].is_zero_matrix:
                non_zero_rows += 1
        if non_zero_rows == len(self._reduced_matrix[i, :]) - 1:
            return 1
        return oo

    def partition_variables(self):
        """Partition the variables into leading and free variables."""
        matrix = self._reduced_matrix
        leading = []
        for i in range(matrix.rows - 1):
            for j in range(matrix.cols - 1):
                if matrix[i, j] != 0:
                    leading.append(self._column_variables[j])
                    break
        free = [var for var in self._column_variables if var not in leading]
        return leading, free

    def _multiple_solutions_set(self, matrix):
        parameters = []
        for i in range(matrix.rows):
            parameter = matrix[i, -1]
            for j in range(matrix.cols - 1):
                if matrix[i, j] == 1:
                    break
            else:
                continue
            for k in range(j + 1, matrix.cols - 1):
                parameter += self._column_variables[k] * -matrix[i, k]
            parameters.append(parameter)
        _, free_variables = self.partition_variables()
        parameters.extend(list(free_variables))
        return FiniteSet(tuple(parameters))

    def _multiple_solutions_set_as_vectors(self, matrix):
        leading, free_variables = self.partition_variables()
        column_vectors = {
            variable: Matrix.zeros(rows=len(self._column_variables)).col(0)
            for variable in free_variables
        }
        for free_variable in free_variables:
            col = self._column_variables.index(free_variable)
            for row in range(len(self._column_variables)):
                if self._column_variables[row] in leading:
                    column_vectors[free_variable][row] = -matrix[row, col]
                elif self._column_variables[row] == free_variable:
                    column_vectors[free_variable][row] = 1
        particular = Matrix.zeros(rows=len(self._column_variables)).col(0)
        for row in range(matrix.rows):
            particular[row] = matrix[row, -1]
        display = [f"{repr(particular)}"]
        for variable in free_variables:
            display.append(f"{repr(column_vectors[variable])} * {variable}")
        print(" + ".join(display))
        vector_sum = particular
        for variable in free_variables:
            vector_sum += variable * column_vectors[variable]
        return FiniteSet(vector_sum)

    def solution_set(self):
        matrix = self._reduced_matrix
        if self.number_of_solutions() == 0:
            return FiniteSet()
        if self.number_of_solutions() == 1:
            for i in range(matrix.rows):
                if matrix[i, :].is_zero_matrix:
                    boundary_row = i
                    break
            else:
                boundary_row = i + 1
            return FiniteSet(tuple(matrix[:boundary_row, -1]))
        return self._multiple_solutions_set_as_vectors(matrix)

    def display_solution(self):
        vars = ", ".join(map(str, self._column_variables))
        print(f"({vars}) = {self.solution_set()}")

    @classmethod
    def from_file(cls, path):
        with open(path) as matrix_file:
            return cls.from_string(matrix_file.read(), sep="\n")

    @classmethod
    def from_string(cls, string, sep=","):
        lines = string.split(sep)
        matrix, column_variables = LinearSystemOfEquationsParser.rows_to_matrix(lines)
        return cls(matrix, column_variables)


def main():
    init_printing(use_unicode=True)
    system = LinearSystemOfEquations.from_file(sys.argv[1])
    print(repr(system._matrix))
    print(repr(system._echelon_matrix))
    print(repr(system._normalized_echelon_matrix))
    print(repr(system._reduced_matrix))
    leading, free = system.partition_variables()
    system.display_solution()


if __name__ == "__main__":
    main()
