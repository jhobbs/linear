import pytest
from sympy import FiniteSet, Matrix, Rational, EmptySet, symbols

from gaussian import LinearSystemOfEquations


@pytest.mark.parametrize(
    "system, solution",
    (
        (
            # 1.19a
            LinearSystemOfEquations.from_string("2x +2y = 5,x -4y = 0"),
            FiniteSet(
                (
                    Rational("2"),
                    Rational("1/2"),
                )
            ),
        ),
        (
            # 1.19b
            LinearSystemOfEquations.from_string("-x +y = 1, x +y = 2"),
            FiniteSet(
                (
                    Rational("1/2"),
                    Rational("3/2"),
                )
            ),
        ),
        (
            # 1.19d
            LinearSystemOfEquations.from_string("-x -y = 1, -3x -3y = 2"),
            EmptySet,
        ),
        (
            # 1.19e
            LinearSystemOfEquations.from_string(
                "4y +z = 20, 2x -2y +z = 0, x +z = 5, x +y -z = 10"
            ),
            FiniteSet(
                (
                    5,
                    5,
                    0,
                )
            ),
        ),
        (
            # 1.20a
            LinearSystemOfEquations.from_string("x +y +z = 5, x -y = 0, y +2z = 7"),
            FiniteSet(
                (
                    1,
                    1,
                    3,
                )
            ),
        ),
        (
            # 1.20c
            LinearSystemOfEquations.from_string(
                "x +3y +z = 0, -x -y = 2, -x +y +2z = 8"
            ),
            FiniteSet(
                (
                    -1,
                    -1,
                    4,
                )
            ),
        ),
    ),
)
def test_finite_solutions(system, solution):
    assert system.solution_set() == solution


@pytest.mark.parametrize(
    "system",
    (
        # 1.19c
        LinearSystemOfEquations.from_string("x -3y +z =1, x +y +2z =14"),
        # 1.19f
        LinearSystemOfEquations.from_string(
            "2x +z +w = 5, y -w = -1, 3x -z -w = 0, 4x +y +2z +w = 9"
        ),
        # 1.20b
        LinearSystemOfEquations.from_string("3x +z = 7, x -y +3z = 4, x +2y -5z = -1"),
    ),
)
def test_infinite_solutions(system):
    solution = system.solution_set()


a, b, c, d, e, x, y, z, w = symbols("a b c d e x y z w")


@pytest.mark.parametrize(
    "system, solution",
    (
        (
            # One.1.2.19a
            LinearSystemOfEquations.from_string("2x +y -z =1, 4x -y =3"),
            FiniteSet(
                (
                    Matrix(
                        [
                            [Rational("2/3")],
                            [Rational("-1/3")],
                            [0],
                        ]
                    )
                    + Matrix(
                        [
                            [Rational("1/6")],
                            [Rational("2/3")],
                            [1],
                        ]
                    )
                    * z
                )
            ),
        ),
        (
            # One.1.2.19b
            LinearSystemOfEquations.from_string(
                "order: x y z w, x -z = 1, y +2z -w = 3, x +2y +3z -w = 7"
            ),
            FiniteSet(
                (
                    Matrix(
                        [
                            [1],
                            [3],
                            [0],
                            [0],
                        ]
                    )
                    + Matrix(
                        [
                            [1],
                            [-2],
                            [1],
                            [0],
                        ]
                    )
                    * z
                )
            ),
        ),
        (
            # One.1.2.19c
            LinearSystemOfEquations.from_string(
                "order: x y z w, x -y +z = 0, y +w =0, 3x -2y +3z +w = 0, -y -w = 0"
            ),
            FiniteSet(
                (
                    Matrix([[0], [0], [0], [0]])
                    + Matrix([[-1], [0], [1], [0]]) * z
                    + Matrix([[-1], [-1], [0], [1]]) * w
                )
            ),
        ),
        (
            # One.1.2.19d
            LinearSystemOfEquations.from_string(
                "a +2b +3c +d -e = 1, 3a -b +c +d +e = 3"
            ),
            FiniteSet(
                (
                    Matrix([[1], [0], [0], [0], [0]])
                    + Matrix(
                        [[Rational("-5 / 7")], [Rational("-8 / 7")], [1], [0], [0]]
                    )
                    * c
                    + Matrix(
                        [[Rational("-3 / 7")], [Rational("-2 / 7")], [0], [1], [0]]
                    )
                    * d
                    + Matrix([[Rational("-1 / 7")], [Rational("4 / 7")], [0], [0], [1]])
                    * e
                )
            ),
        ),
        (
            # One.1.2.20b
            LinearSystemOfEquations.from_string(
                "x +y -2z = 0, x -y = -3, 3x -y -2z = -6, 2y -2z = 3"
            ),
            FiniteSet(
                (
                    Matrix([[Rational("-3 / 2")], [Rational("3 / 2")], [0]])
                    + Matrix([[1], [1], [1]]) * z
                )
            ),
        ),
    ),
)
def test_infinite_solutions_with_solution_set(system, solution):
    assert system.solution_set() == solution
