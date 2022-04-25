import numpy as np
from sympy import Rational
from typing import List


def parse_term(term: str):
    if not term[-1].isalpha():
        raise Exception(f"Invalid term: {term}")

    variable = term[-1]
    coefficient = 0
    for i, c in enumerate(reversed(term[:-1])):
        if not c.isdigit():
            if c == "+":
                break
            if c == "-":
                coefficient = coefficient * -1
            else:
                raise Exception(f"Invalid term: {term}")
        else:
            coefficient += int(c) * 10 ** (i)

    coefficient = Rational(coefficient)

    print(f"coeff: {coefficient} var: {variable}")
    return (
        coefficient,
        variable,
    )


def parse_variable_terms(variable_part: str):
    split_row = variable_part.split()
    terms = {}
    for raw_term in split_row:
        coefficient, variable = parse_term(raw_term)
        if variable in terms:
            terms[variable] += coefficient
        else:
            terms[variable] = coefficient
    return terms


def parse_raw_row(raw_row: str):
    variable_part, constant_part = raw_row.split("=")
    terms = parse_variable_terms(variable_part)
    constant = Rational(int(constant_part.strip()))
    return terms, constant


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


def rows_to_matrix(input_rows):
    rows = []
    constants = []
    for row in input_rows:
        row_variables, constant = parse_raw_row(row)
        rows.append(row_variables)
        constants.append(constant)
    var_col = variables_to_cols(rows)
    matrix = np.ndarray((len(rows), len(var_col) + 1))
    for i, row in enumerate(rows):
        for variable, coefficient in row.items():
            matrix[i][var_col[variable]] = coefficient
        matrix[i][-1] = constants[i]
    return matrix


def get_rows_from_stdin():
    rows = []
    while True:
        raw_row = input("Next row: ")
        if len(raw_row) == 0:
            break
        rows.append(raw_row)
    return rows


def is_reduced_row_echelon(matrix):
    if len(matrix) == 0:
        return True

    for i in range(len(matrix[0])):
        if matrix[0][i] == 0:
            if not np.all(matrix[1:, i] == 0):
                return False
        elif matrix[0][i] == 1:
            if not np.all(matrix[1:, i] == 0):
                return False
            return is_reduced_row_echelon(matrix[1:, i:])
        else:
            return False
    return True


def gauss_jordan(matrix):
    while not is_reduced_row_echelon(matrix):
        pass


def compare_rows(row_a, row_b):
    for i in range(len(row_a)):
        if row_a[i] == row_b[i]:
            continue
        if row_a[i] == 0 or row_a[i] < row_b[i]:
            return -1
        return 1
    return 0


def sort_rows(matrix):
    # bubble sort
    while True:
        swapped = False
        for i in range(len(matrix)):
            if i == len(matrix) - 1:
                break
            if compare_rows(matrix[i], matrix[i + 1]) == -1:
                matrix[[i, i + 1]] = matrix[[i + 1, i]]
                swapped = True

        if swapped == False:
            break
    return matrix


def reduce_rows(matrix):
    if len(matrix) == 0:
        return

    leading_coeff = matrix[0][0]
    if leading_coeff == 0:
        raise Exception(f"Unexpected 0 leading coefficient: {i} {matrix}")
    matrix[0] = matrix[0] * 1 / leading_coeff
    for i in range(len(matrix) - 1):
        row_leading_coeff = matrix[i + 1][0]
        if row_leading_coeff == 0:
            continue
        matrix[i + 1] = matrix[i + 1] + (-row_leading_coeff * matrix[0])
    reduce_rows(matrix[1:, 1:])


def get_matrix_from_file():
    with open("matrix.txt") as matrix_file:
        lines = matrix_file.readlines()
    matrix = rows_to_matrix(lines)
    return matrix


def main():
    matrix = get_matrix_from_file()
    sort_rows(matrix)
    reduce_rows(matrix)
    print(matrix)


if __name__ == "__main__":
    main()
