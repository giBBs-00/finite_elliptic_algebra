from functools import reduce
from itertools import zip_longest
from math import gcd
from numbers import Integral
from typing import Iterator, List, Tuple, Union


class Polynomial:
    def __init__(self, coefficients: List[Integral]) -> None:
        if not isinstance(coefficients, list) or not all(isinstance(c, Integral) for c in coefficients):
            raise TypeError("coefficients must be a list of integers")

        # normalize the coefficients
        while coefficients and coefficients[-1] == 0:
            coefficients.pop()
        if not coefficients:
            coefficients = [0]  # The zero polynomial

        self._coefficients = coefficients

    @property
    def coefficients(self) -> List[Integral]:
        """
        Get the polynomial's coefficients.

        :return: A copy of the coefficient list to prevent modification.
        """
        return self._coefficients.copy()

    @property
    def degree(self) -> Integral:
        """
        Return the degree of the polynomial.

        The degree of the zero polynomial is 0.

        :return: The degree of the polynomial

        :examples:
        >>> Polynomial([0]).degree
        0
        >>> Polynomial([0, 0, 1]).degree
        2
        """
        for i in range(len(self.coefficients) - 1, -1, -1):
            if self.coefficients[i] != 0:
                return i
        return 0  # The zero polynomial has degree 0

    def __str__(self) -> str:
        """Return a string representation of the polynomial in descending order of degree."""
        if not self.coefficients:
            return "0"

        # Start with the highest degree term (reverse order)
        terms = []
        for i in range(len(self.coefficients) - 1, -1, -1):
            c = self.coefficients[i]
            if c == 0:
                continue

            # Handle coefficient formatting
            if i == 0:  # Constant term
                terms.append(str(c))
            elif i == 1:  # Linear term
                if c == 1:
                    terms.append("x")
                elif c == -1:
                    terms.append("-x")
                else:
                    terms.append(f"{c}x")
            else:  # Higher power terms
                if c == 1:
                    terms.append(f"x^{i}")
                elif c == -1:
                    terms.append(f"-x^{i}")
                else:
                    terms.append(f"{c}x^{i}")

        # Join terms with proper operators
        if not terms:
            return "0"

        result = terms[0]
        for term in terms[1:]:
            if term[0] == '-':
                result += f" {term}"  # Already has the minus sign
            else:
                result += f" + {term}"

        return result

    def __repr__(self) -> str:
        return f'Polynomial({self.coefficients})'

    def __bool__(self) -> bool:
        return any(c != 0 for c in self.coefficients)

    def __eq__(self, other):
        if not isinstance(other, Polynomial):
            return NotImplemented
        return self.coefficients == other.coefficients

    def __ne__(self, other):
        if not isinstance(other, Polynomial):
            return NotImplemented
        return self.coefficients != other.coefficients

    def __len__(self) -> int:
        return len(self.coefficients)

    def __iter__(self) -> Iterator[int]:
        return iter(self.coefficients)

    def __reversed__(self) -> Iterator[int]:
        return reversed(self.coefficients)

    def __contains__(self, item: int) -> bool:
        return item in self.coefficients

    def __call__(self, x: float) -> float:
        """
        Evaluate the polynomial at a given value.

        :param x: The value at which to evaluate the polynomial.
        :return: The result of evaluating the polynomial at x.
        """
        return sum(c * x ** i for i, c in enumerate(self.coefficients))

    def __add__(self, other: Union["Polynomial", int]) -> "Polynomial":
        if isinstance(other, int):
            return self.__radd__(other)
        if isinstance(other, Polynomial):  # Fixed type check
            result = [
                x + y for x, y in zip_longest(self.coefficients, other.coefficients, fillvalue=0)]
            return Polynomial(result)
        return NotImplemented

    def __radd__(self, other: int) -> "Polynomial":
        if isinstance(other, int):
            coefficients = self.coefficients.copy()
            coefficients[0] += other
            return Polynomial(coefficients)
        return NotImplemented

    def __sub__(self, other: Union["Polynomial", int]) -> "Polynomial":
        if isinstance(other, int):
            return self.__rsub__(other)
        if isinstance(other, type(Polynomial)):
            result = [
                x - y for x, y in zip_longest(self.coefficients, other.coefficients, fillvalue=0)]
            return Polynomial(result)
        return NotImplemented

    def __rsub__(self, other: int) -> "Polynomial":
        if isinstance(other, int):
            coefficients = [-c for c in self.coefficients]
            coefficients[0] += other
            return Polynomial(coefficients)
        return NotImplemented

    def __mul__(self, other: Union["Polynomial", int]) -> "Polynomial":
        if isinstance(other, int):
            return self.__rmul__(other)
        if isinstance(other, Polynomial):
            result = [0] * (len(self.coefficients) +
                            len(other.coefficients) - 1)
            for i, c1 in enumerate(self.coefficients):
                for j, c2 in enumerate(other.coefficients):
                    result[i + j] += c1 * c2
            return Polynomial(result)
        return NotImplemented

    def __rmul__(self, other: int) -> "Polynomial":
        if isinstance(other, int):
            result = [other * c for c in self.coefficients]
            return Polynomial(result)
        return NotImplemented

    def __neg__(self) -> "Polynomial":
        return Polynomial([-c for c in self.coefficients])

    def __pos__(self) -> "Polynomial":
        return Polynomial([+c for c in self.coefficients])

    def __abs__(self) -> "Polynomial":
        return Polynomial([abs(c) for c in self.coefficients])

    def __pow__(self, n: int) -> "Polynomial":
        if not isinstance(n, int):
            raise TypeError("exponent must be an integer")
        if n < 0:
            raise ValueError("negative powers are not supported")
        if n == 0:
            return Polynomial([1])
        result = Polynomial([1])
        for _ in range(n):
            result *= self
        return result

    def _divmod_helper(self, dividend, divisor):
        """Helper method for polynomial division."""
        if not divisor:
            raise ZeroDivisionError("division by zero polynomial")

        # Normalize divisor to have leading coefficient 1
        divisor_degree = max(
            (i for i, c in enumerate(divisor) if c != 0), default=0)
        divisor_lc = divisor[divisor_degree]

        # Create a copy of the dividend
        remainder = dividend.copy()

        # Calculate the degree of the quotient
        dividend_degree = max(
            (i for i, c in enumerate(remainder) if c != 0), default=0)
        quotient_degree = max(0, dividend_degree - divisor_degree)
        quotient = [0] * (quotient_degree + 1)

        # Perform polynomial long division
        for i in range(quotient_degree, -1, -1):
            if dividend_degree < divisor_degree:
                break

            # Calculate the coefficient of the quotient
            coef = remainder[dividend_degree] // divisor_lc
            quotient[i] = coef

            # Subtract coef * divisor * x^i from the remainder
            for j in range(divisor_degree + 1):
                if j < len(divisor) and divisor[j] != 0:
                    idx = i + j
                    if idx < len(remainder):
                        remainder[idx] -= coef * divisor[j]

            # Recalculate the degree of the remainder
            while remainder and remainder[-1] == 0:
                remainder.pop()

            dividend_degree = len(remainder) - 1

        return quotient, remainder

    def divmod(self, other: "Polynomial") -> Tuple["Polynomial", "Polynomial"]:
        """
        Perform polynomial long division and return quotient and remainder.

        :param other: The divisor polynomial
        :return: A tuple (quotient, remainder) of Polynomial objects
        """
        if not isinstance(other, Polynomial):
            raise TypeError("other must be a Polynomial object")

        quotient, remainder = self._divmod_helper(
            self.coefficients.copy(), other.coefficients.copy())
        return Polynomial(quotient), Polynomial(remainder)

    def __floordiv__(self, other: Union["Polynomial", int]) -> "Polynomial":
        """Integer division of polynomials."""
        if isinstance(other, int):
            return Polynomial([c // other for c in self.coefficients])

        quotient, _ = self.divmod(other)
        return quotient

    def __mod__(self, other: Union["Polynomial", int]) -> "Polynomial":
        """Modulo operation for polynomials."""
        if isinstance(other, int):
            return Polynomial([c % other for c in self.coefficients])

        _, remainder = self.divmod(other)
        return remainder

    def _content(self) -> int:
        """
        Calculate the content of the polynomial (GCD of all coefficients).

        :return: The content (an integer)
        """
        non_zero = [abs(c) for c in self.coefficients if c != 0]
        if not non_zero:
            return 0

        return reduce(gcd, non_zero)

    def gcd(self, other: "Polynomial") -> "Polynomial":
        """
        Calculate the greatest common divisor of this polynomial and another.

        Uses the Euclidean algorithm.

        :param other: Another polynomial
        :return: The GCD polynomial
        """
        if not isinstance(other, Polynomial):
            raise TypeError("other must be a Polynomial object")

        # Handle special cases
        if not self:
            return other
        if not other:
            return self

        a, b = self, other

        # Ensure a has higher or equal degree
        if a.degree() < b.degree():
            a, b = b, a

        # Euclidean algorithm
        while b:
            _, r = a.divmod(b)
            a, b = b, r

        # Normalize the result to have a positive leading coefficient
        if a.coefficients and a.coefficients[-1] < 0:
            a = -a

        return a
