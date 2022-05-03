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
    ),
)
def test_finite_solutions(system, solution):
    assert system.solution_set() == solution
