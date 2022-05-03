import pytest
from sympy import FiniteSet, Rational

from gaussian import LinearSystemOfEquations


@pytest.mark.parametrize(
    "system, solution",
    (
        (
            LinearSystemOfEquations.from_string("2x +2y = 5,x -4y = 0"),
            FiniteSet(
                (
                    Rational("2"),
                    Rational("1/2"),
                )
            ),
        ),
        (
            LinearSystemOfEquations.from_string("-x +y = 1, x +y = 2"),
            FiniteSet(
                (
                    Rational("1/2"),
                    Rational("3/2"),
                )
            ),
        ),
    ),
)
def test_finite_solutions(system, solution):
    assert system.solution_set() == solution
