import pytest
from sympy import FiniteSet, Rational, EmptySet

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
