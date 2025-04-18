from functools import cached_property
from typing import List, Union
from utils.helpers import extended_gcd
from prime import is_prime
from cyclic import Cyclic
from itertools import zip_longest


class Fp:
    def __init__(self, p: int, a: int) -> None:

        # Ensure that the field modulus p is a prime number.
        if not is_prime(p):
            raise ValueError("Field modulus p must be a prime number.")

        # Check that the element a is an integer.
        if not isinstance(a, int):
            raise TypeError(f"The element a of F{p} must be an integer.")

        # Store the modulus and reduce a modulo p to keep it in the field.
        self.p = p
        self.a = a % p

    @cached_property
    def elements(self) -> List["Fp"]:

        # Returns a list of all elements in the finite field F_p.
        # Computed once and cached for efficiency.

        return [Fp(self.p, i) for i in range(self.p)]

    @property
    def inverse(self) -> "Fp":
        """
        The multiplicative inverse of an element in Fp.
        """
        # The modular inverse of 0 does not exist in any field.
        if self.a == 0:
            raise ValueError(f"No inverse exists for 0 in F{self.p}.")

        # Use the extended Euclidean algorithm to find the inverse of a mod p.
        _, x, _ = extended_gcd(self.a, self.p)

        # Ensure the inverse is in the correct field by reducing mod p.
        return Fp(self.p, x)

    # Internal method: check if `other` is an element of the same field F_p.
    def __in_field(self, other: "Fp") -> bool:
        return isinstance(other, Fp) and self.p == other.p

    # Internal method: raise an error if `other` is not in the same field.
    def __assert_in_field(self, other: "Fp") -> None:
        if not isinstance(other, Fp):
            raise TypeError("Operand must be an instance of Fp.")
        if self.p != other.p:
            raise ValueError(
                f"Operands are from different fields: F{self.p} and F{other.p}.")

    # Field addition: (a + b) mod p
    def __add__(self, other: "Fp"):
        self.__assert_in_field(other)
        return Fp(self.p, self.a + other.a)

    # Field subtraction: (a - b) mod p
    def __sub__(self, other: "Fp"):
        self.__assert_in_field(other)
        return Fp(self.p, self.a - other.a)

    # Field multiplication: (a * b) mod p
    def __mul__(self, other: "Fp"):
        self.__assert_in_field(other)
        return Fp(self.p, self.a * other.a)

    # Field equality: check same field and same value
    def __eq__(self, other: "Fp"):
        return self.__in_field(other) and self.a == other.a

    # Field division: (a / b) mod p using b⁻¹
    def __truediv__(self, other: "Fp"):
        self.__assert_in_field(other)
        return Fp(self.p, self.a * other.inverse.a)

    # Field exponentiation: a^n mod p
    def __pow__(self, other: int):
        if not isinstance(other, int):
            raise TypeError("Exponent must be an integer.")

        # Exponentiation by repeated multiplication (manual implemntation)
        result = self.a
        for _ in range(other - 1):
            result = self.a * result
        return Fp(self.p, result)

    def __hash__(self):
        return hash((self.p, self.a))

    # Human-readable string representation (e.g. just prints the value)
    def __str__(self) -> str:
        return f"{self.a}"

    # Official debug representation (includes field info)
    def __repr__(self) -> str:
        return f"Fp({self.p}, {self.a})"


class GF:
    def __init__(self, p, irred: List[Union[Fp, int]]):
        # Ensure p is a prime
        if not is_prime(p):
            raise TypeError("Field characteristic p must be prime.")
        self.p = p

        # Convert integer coefficients to Fp instances if needed
        irred = [Fp(self.p, x) if isinstance(x, int) else x for x in irred]

        # Strip leading zeros from the polynomial
        while irred and irred[0].a == 0:
            irred.pop(0)

        # If the irreducible polynomial is all zero, use [0] to avoid empty list
        self.irred = irred if irred else [Fp(self.p, 0)]

        # Set degree of the extension
        self.n = len(self.irred) - 1

        # Check if the polynomial is primitive
        if not self.is_primitive:
            raise ValueError(
                f"The provided polynomial is not primitive over GF({p}^{self.n}).")

    def a(self, a: int) -> "GF.GFe":
        # Return the non-zero element corresponding to α^a from the multiplicative group
        if not isinstance(a, int):
            raise TypeError("Exponent must be an integer.")
        if not (0 <= a < len(self.alpha_rep)):
            raise ValueError(
                f"Exponent must be in range 0 to {len(self.alpha_rep) - 1}.")

        # This representation is only for the multiplicative group.
        # α^a is at index a + 1
        elem = self.elements[a + 1].poly
        return GF.GFe(self, elem)

    class GFe:
        def __init__(self, gf: "GF", poly: List[Union[Fp, int]]):
            # Ensure the field definition is valid
            if not isinstance(gf, GF):
                raise TypeError("gf must be an instance of GF.")

            self.gf = gf

            # Convert all integer coefficients to Fp elements
            # self.p is defined outside __init__
            poly = [Fp(self.p, x) if isinstance(x, int) else x for x in poly]

            # Remove leading zeros to keep polynomial in canonical form
            while poly and poly[0].a == 0:
                poly.pop(0)

            # If all coefficients were zero, keep as zero element
            if not poly:
                poly = [Fp(self.p, 0)]

            self.poly = poly

        @property
        def p(self) -> int:
            # Characteristic of the base field Fp
            return self.gf.p

        @property
        def n(self) -> int:
            # Degree of the polynomial representing this element
            # Not the extension degree
            return len(self.poly) - 1

        @property
        def a(self) -> "Cyclic":

            # Return the α-exponent representation of this nonzero element
            if self.poly == [Fp(self.p, 0)]:
                return GF.GFe(self.gf, [Fp(self.p, 0)])

            # Get rid of 0 elements in GF(p^n)
            elements = self.gf.elements.pop(0)
            alpha = self.gf.alpha_rep

            try:
                index = elements.index(self)
            except ValueError:
                raise ValueError(
                    "Element not found in the field’s element list.")

            return alpha[index]

        def __add__(self, other: "GF.GFe") -> "GF.GFe":
            # Ensure elements are from the same field
            if not isinstance(other, GF.GFe):
                raise TypeError(
                    "Both operands must be elements of the same GF.")

            # Reverse to add lowest-degree terms first
            a = list(reversed(self.poly))
            b = list(reversed(other.poly))

            # Pad with zeros and add elementwise
            result = [x + y for x,
                      y in zip_longest(a, b, fillvalue=Fp(self.p, 0))]

            # Reverse back to normal order
            return self.gf(list(reversed(result)))

        # Same idea as __add__ method
        def __sub__(self, other: "GF.GFe") -> "GF.GFe":
            if not isinstance(other, GF.GFe):
                raise TypeError(
                    "Both operands must be elements of the same GF.")

            a = list(reversed(self.poly))
            b = list(reversed(other.poly))

            result = [x - y for x,
                      y in zip_longest(a, b, fillvalue=Fp(self.p, 0))]

            return self.gf(list(reversed(result)))

        def __mul__(self, other: "GF.GFe") -> "GF.GFe":
            # Handle zero case early
            if self == GF.GFe(self.gf, [0]) or other == GF.GFe(self.gf, [0]):
                return GF.GFe(self.gf, [Fp(self.p, 0)])

            # Ensure operand is a valid field element
            if not isinstance(other, GF.GFe):
                raise TypeError("Operand must be an instance of GF.GFe.")

            elements = self.gf.elements
            alpha = self.gf.alpha_rep

            # Find indices of self and other
            i = elements.index(self)
            j = elements.index(other)

            # Multiply as exponents: α^i * α^j = α^(i + j)
            a = alpha[i - 1] * alpha[j - 1]

            # Map result back to polynomial form
            return elements[a.n + 1]

        def __pow__(self, other: int) -> "GF.GFe":
            # Exponentiation for elements in GF(p^n)
            if not isinstance(other, int):
                raise TypeError("Exponent must be an integer.")

            if self == GF.GFe(self.gf, [0]):
                if other == 0:
                    # 0^0 = 1 by convention
                    return GF.GFe(self.gf, [Fp(self.p, 1)])
                # 0^e = 0 for e > 0
                return GF.GFe(self.gf, [Fp(self.p, 0)])

            elements = self.gf.elements
            alpha = self.gf.alpha_rep

            # Get index of self
            i = elements.index(self)

            # Raise α^i to power, α^(i * other)
            result_alpha = alpha[i - 1] ** other

            # Return α^(i * other) as element from lookup table
            return elements[result_alpha.n + 1]

        def __truediv__(self, other: "GF.GFe") -> "GF.GFe":
            # Type and field consistency check
            if not isinstance(other, GF.GFe):
                raise TypeError("Operand must be an instance of GF.GFe.")

            # 0 / x = 0 for any x ≠ 0
            if self == GF.GFe(self.gf, [0]):
                return GF.GFe(self.gf, [Fp(self.p, 0)])

            # x / 0 is undefined
            if other == GF.GFe(self.gf, [0]):
                raise ZeroDivisionError(
                    "Cannot divide by zero in a finite field.")

            elements = self.gf.elements
            alpha = self.gf.alpha_rep

            # Get indices: self = α^i, other = α^j → result = α^(i - j)
            i = elements.index(self)
            j = elements.index(other)

            result_alpha = alpha[i - 1] / alpha[j - 1]

            return elements[result_alpha.n + 1]

        def __iter__(self):
            return iter(self.poly)

        def __hash__(self):
            return hash(tuple(self.poly))

        def __eq__(self, other):
            if not isinstance(other, GF.GFe):
                return False

            return self.poly == other.poly

        def __repr__(self):
            return str(self)

        # Return a human-readable string representation
        # of the element as a polynomial over Fp
        def __str__(self):
            terms = []
            for i, coeff in enumerate(self.poly):
                if coeff.a == 0:
                    continue
                power = self.n - i
                if power == 0:
                    terms.append(f"{coeff}")
                elif power == 1:
                    terms.append(f"{coeff if coeff.a > 1 else ''}x")
                else:
                    superscript = "".join(
                        "⁰¹²³⁴⁵⁶⁷⁸⁹"[int(d)] for d in str(power))
                    terms.append(
                        f"{coeff if coeff.a > 1 else ''}x{superscript}")

            return " + ".join(terms) if terms else "0"

    @cached_property
    def alpha_rep(self) -> List["Cyclic"]:
        # Precompute the multiplicative group representation: α^1 to α^{p^n - 1}
        # Skips α^0 = 1, which is not needed for exponent indexing
        return Cyclic(1, len(self.elements) - 1).elements

    @cached_property
    def elements(self) -> List["GF.GFe"]:
        # Generate all elements of GF(p^n), starting with 0 followed by powers of α
        return self.__generate()

    @cached_property
    def is_primitive(self) -> bool:
        # Check if the irreducible polynomial generates the full field GF(p^n)
        # verifies the number of generated elements matches the expected size
        return len(self.elements) == self.p**self.n

    def __generate(self) -> List["GF.GFe"]:
        # Start with zero element
        zero = [Fp(self.p, 0)]

        # Generate monomials x, x^2, ..., x^{n-1}
        pows = [
            [Fp(self.p, 1)] + [Fp(self.p, 0) for _ in range(j)]
            for j in range(self.n)
        ]

        # Initial elements: [0, x^0, x^1, ..., x^{n-1}]
        elements = [GF.GFe(self, zero)] + [GF.GFe(self, p) for p in pows]

        # Remove leading coefficient from irreducible (assumed to be 1)
        irred = self.irred[1:]

        # Setup for modular reduction: -f(x) with f = irreducible
        reduct = [Fp(self.p, -1) * coeff for coeff in irred]
        mult = reduct.copy()

        # Check if reduct is already in elements
        reduct_check = GF.GFe(self, reduct)
        for element in elements:
            if reduct_check == element:
                return elements

        elements.append(GF.GFe(self, reduct))

        # Iteratively build new elements using modular multiplication
        while len(elements) < self.p**self.n:
            # Multiply current polynomial (mult) by x, reducing mod irreducible
            next_poly = [
                a + b
                for a, b in zip(
                    mult[1:] + [Fp(self.p, 0)],
                    [mult[0] * c for c in reduct]
                )
            ]

            # Check for cycle completion by seeing if next element is already in elements
            # Not just checking against 1 incase something unintended happens
            next_check = GF.GFe(self, next_poly)
            for element in elements:
                if next_check == element:
                    return elements

            elements.append(GF.GFe(self, next_poly))

            mult = next_poly

        return elements

    def fromList(self, coefficients: List[Union[Fp, int]]) -> "GF.GFe":
        # Create an instance of GFe from a list of Fp or int coefficients
        return GF.GFe(self, coefficients)

    def fromString(self, input: str) -> "GF.GFe":
        # Parse a polynomial string like "x^2 + 3x + 1" into a GFe element
        terms = [term.strip() for term in input.split("+")]
        coeffs = []

        for term in terms:
            if "x" in term:
                coeff, power = term.split("x")
                coeff = coeff.strip()
                power = power.strip() if power else ""

                # Handle cases like "x", "-x"
                if coeff == "":
                    coeff = 1
                elif coeff == "-":
                    coeff = -1
                else:
                    coeff = int(coeff)

                if "^" in power:
                    power = int(power.split("^")[1])
                elif power == "":
                    power = 1
                else:
                    power = int(power)
            else:
                coeff = int(term) % self.p
                power = 0

            # Ensure coeffs list is long enough
            while len(coeffs) <= power:
                coeffs.append(Fp(self.p, 0))

            coeffs[power] = Fp(self.p, coeff)

        return GF.GFe(self, coeffs)

    def __call__(self, input: Union[str, List[Union[Fp, int]]]) -> "GF.GFe":
        # Create a GFe element from either a string (like "x^2 + 1")
        # or a list of coefficients (like [1, 0, 1])
        if isinstance(input, str):
            return self.fromString(input)
        elif isinstance(input, list):
            return self.fromList(input)
        else:
            raise TypeError(
                "Input must be a string or a list of integers/Fp elements.")
