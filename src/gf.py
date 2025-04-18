from collections import defaultdict
from typing import List, Tuple, Union
from functools import cached_property
import itertools
from itertools import zip_longest
import math
from math import isqrt
import statistics
import pandas as pd
from pandas import DataFrame
import cmath
from cmath import tau
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from prime import is_prime

# Need to find divisors sometimes
def divisors(n):
    """Return all positive divisors of n."""
    divs = set()
    # only need to check up sqrt(n) + 1 as we add to divs in factor pairs
    for i in range(1, isqrt(n) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(divs)

# Extended Euclidean Algorithm:
# Returns gcd(a, b), x, y such that: ax + by = gcd(a, b)
def extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """Return gcd(a, b) and integers x, y such that ax + by = gcd(a, b)."""
    # Base case: gcd(a, 0) = a
    if b == 0:
        return a, 1, 0

    # Recursive case: compute gcd and coefficients
    gcd, x1, y1 = extended_gcd(b, a % b)

    # Update x and y
    x = y1
    y = x1 - (a // b) * y1

    return gcd, x, y

# Need for clean character representation globally 
def clean_complex(z: complex, dp: int = 4, epsilon: float = 1e-10) -> complex:
    """Round real and imaginary parts of a complex number, and snap to integer n, if abs(n - val) < epsilon."""
    def snap(x):
            nearest = round(x)
            return nearest if abs(x - nearest) < epsilon else round(x, dp)
        
    return complex(snap(z.real), snap(z.imag))


class AG:
    def __init__(self, *orders):
        """Initialize an abelian group of given orders (i.e., Z/mZ × Z/nZ)."""
        self.orders = tuple(sorted(orders))
        self.rank = len(orders)
        self.order = math.prod(orders)

    def generate(self):
        """Generate all elements of the abelian group AG."""
        if self.rank == 1:
            return [self.AGE(self, i) for i in range(self.order)]

        if self.rank == 2:
            a = self.AGE(self, 1, 0)
            b = self.AGE(self, 0, 1)
            elements = []
            for i in range(self.orders[0]):
                for j in range(self.orders[1]):
                    elements.append(a * i + b * j)
            return elements

        raise NotImplementedError("Rank 3 and higher not supported")

    def subgroups(self):
        """Return a list of all distinct cyclic subgroups of the group."""
        elems = self.generate()
        subgroups = []
        seen = set()

        for g in elems:
            subgroup = [i * g for i in range(g.order)]
            check = frozenset(subgroup)

            if check not in seen:
                seen.add(check)
                subgroups.append(subgroup)

        return subgroups

    def torsion(self, n: int):
        """Return elements of the group with exact order n."""
        torsion = []
        for a in self.generate():
            if a.order == n:
                torsion.append(a)
        return torsion

    def characters(self):
        """Return all characters χ: G → ℂ* as functions, each corresponding to a group element."""
        group_elements = self.generate()
        character_functions = []

        # For each character, represented by an element a in G
        for a in group_elements:
            def make_char(a):
                def χ(g):
                    # Pairwise dot product mod order for each coordinate
                    val = 0
                    for ai, gi, ni in zip(a.coords, g.coords, self.orders):
                        val += (ai * gi) / ni
                    return cmath.exp(tau * 1j * val)
                return χ
            char_func = make_char(a)
            character_functions.append((a, char_func))

        return character_functions

    def sum_char(self):
        """Compute the sum of each character χ over all group elements."""
        characters = self.characters()
        elements = self.generate()
        char_sums = []
        for char in characters:
            char_sum = 0
            for g in elements:
                char_sum += char(g)
            char_sums.append(char_sum)
        return char_sums

    class AGE:
        def __init__(self, group: "AG", *coords):
            """Create an abelian group element from coordinates, within a specific group."""
            self.group = group
            self.coords = tuple(coords)

            if not isinstance(group, AG):
                raise TypeError(
                    "An abelian group element must live in an abelian group")
            if len(coords) != len(group.orders):
                raise ValueError("Coordinate count does not match group rank.")


        @property
        def order(self):
            """Return the order of this element (smallest n such that n⋅g = 0)."""
            divs = divisors(self.group.order)
            for div in divs:
                if div * self == AG.AGE(self.group, *([0] * len(self.coords))):
                    return div

        @property
        # The signature of an element in a group in the sense of characters
        def signature(self):
            """Return the character signature: values χ(g) for all characters χ."""
            characters = self.group.characters()
            return [char(self) for char in characters]

        def __add__(self, other):

            if not isinstance(other, AG.AGE):
                raise TypeError("Operands must be in the same Abelian Group")

            return AG.AGE(self.group, *((a + b) % n for a, b, n in zip(self.coords, other.coords, self.group.orders)))

        def __mul__(self, scalar):

            if not isinstance(scalar, int):
                raise TypeError("Scalar must be an integer")

            return AG.AGE(self.group, *(((scalar % n) * a) % n for a, n in zip(self.coords, self.group.orders)))

        def __rmul__(self, scalar):

            if not isinstance(scalar, int):
                raise TypeError("Scalar must be an integer")

            return AG.AGE(self.group, *(((scalar % n) * a) % n for a, n in zip(self.coords, self.group.orders)))

        def __eq__(self, other):
            return isinstance(other, AG.AGE) and self.coords == other.coords and self.group.orders == other.group.orders

        def __hash__(self):
            return hash(self.coords)

        def __repr__(self):
            if len(self.coords) == 1:
                return f"({self.coords[0]} mod {self.group.orders[0]})"
            else:
                return f"({self.coords[0]} mod {self.group.orders[0]}, {self.coords[1]} mod {self.group.orders[1]})"


class Fp:
    def __init__(self, p: int, a: int) -> None:
        """Create an element of the prime field Fp (mod p)."""
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
        """List all elements in the field Fp(0) to Fp(p - 1)."""
        # Returns a list of all elements in the finite field F_p.
        # Computed once and cached for efficiency.

        return [Fp(self.p, i) for i in range(self.p)]

    @property
    def inverse(self) -> "Fp":
        """Return the multiplicative inverse of the element."""
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
        """Initialize a finite field GF(p^n) using a primitive irreducible polynomial over F(P), of degree n."""
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
        if self.n <= 0:
            raise ValueError("Irreducible polynomial must be of degree ≥ 1.")


    def a(self, a: int) -> "GF.GFe":
        """Return α^a as a field element (based on multiplicative group representation)."""
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
    
    def map_to_abelian(self) -> Tuple[dict, dict, AG]:
        nonzero_elements = self.elements[1:]  # skip 0
        order = len(nonzero_elements)  # = p^n - 1
        group = AG(order)
        ag_elements = group.generate()
        
        gf_to_ag = {gf_elem: ag_elem for gf_elem, ag_elem in zip(nonzero_elements, ag_elements)}
        ag_to_gf = {v: k for k, v in gf_to_ag.items()}
        return gf_to_ag, ag_to_gf
    
    def all_signatures(self):
        """Return the character signature for each elements of GF(p^n)"""
        # Compute the signature of all elements of GF(p^n)
        signatures = []
        elements = self.elements[1:]
        print(elements)
        for g in elements:
            signatures.append(g.signature)
        return signatures

    class GFe:
        def __init__(self, gf: "GF", poly: List[Union[Fp, int]]):
            """Create a GF(p^n) element from a polynomial over Fp."""
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

            if len(poly) > self.gf.n:
                raise ValueError(f"Element polynomial degree exceeds field extension degree {self.gf.n}.")


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
            """Return the α-power representation of this nonzero field element."""
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
        
        @property
        def signature(self):
            """Return the AG character signature of this point.""" 
            if self == self.gf([0]):
                return "Undefined"
            gf_ag, ag_gf = self.gf.map_to_abelian()
            ag = next(iter(ag_gf)).group
            chars = ag.characters()
            signature = []
            for name, func in chars:
                signature.append((ag_gf[name], func(gf_ag[self])))
            return signature


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
        """Return the powers of a generator α of GF(p^n)^*.""" 
        # Precompute the multiplicative group representation: α^1 to α^{p^n - 1}
        # Skips α^0 = 1, which is not needed for exponent indexing
        return Cyclic(1, len(self.elements) - 1).elements

    @cached_property
    def elements(self) -> List["GF.GFe"]:
        """Return all elements in the field."""  
        # Generate all elements of GF(p^n), starting with 0 followed by powers of α
        return self.__generate()

    @cached_property
    def is_primitive(self) -> bool:
        """Check whether the field was constructed using a primitive polynomial."""
        # Check if the irreducible polynomial generates the full field GF(p^n)
        # verifies the number of generated elements matches the expected size
        return len(self.elements) == self.p**self.n

    def __generate(self) -> List["GF.GFe"]:
        """Generate all elements in GF(p^n) via modular polynomial arithmetic."""
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
        """Create a GFe from list of coefficients. Entries can be integers or instances of Fp."""
        # Create an instance of GFe from a list of Fp or int coefficients
        return GF.GFe(self, coefficients)

    def fromString(self, input: str) -> "GF.GFe":
        """Parse a polynomial string and return a GFe element."""  
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
        """Convert from string or list to GFe element (dispatch function)."""
        # Create a GFe element from either a string (like "x^2 + 1")
        # or a list of coefficients (like [1, 0, 1] or [Fp(p, 1), Fp(p, 0), Fp(p, 1)])
        if isinstance(input, str):
            return self.fromString(input)
        elif isinstance(input, list):
            return self.fromList(input)
        else:
            raise TypeError(
                "Input must be a string or a list of integers/Fp elements.")


class Cyclic:
    def __init__(self, n: int, k: int):
        """Return the α-power representation of this nonzero field element."""
        # k is the order of the cyclic group
        # n is the exponent (reduced mod k)
        if not isinstance(n, int) or not isinstance(k, int):
            raise TypeError("Both n and k must be integers.")
        if k <= 0:
            raise ValueError("Group order k must be a positive integer.")

        self.k = k
        self.n = n % k  # Ensure n is always reduced mod k

    @cached_property
    def elements(self):
        """Return a list of all k elements in this cyclic group."""
        # Return all elements a(0) to a(k - 1)
        return [Cyclic(i, self.k) for i in range(self.k)]

    @property
    def inverse(self):
        # The inverse of a(n) is a(-n) ≡ a(k - n) mod k
        return Cyclic(self.k - self.n, self.k)

    def __mul__(self, other: "Cyclic") -> "Cyclic":
        # Multiplication in the cyclic group: add exponents mod k
        if self.k != other.k:
            raise ValueError("Cyclic elements must belong to the same group.")
        return Cyclic(self.n + other.n, self.k)

    def __pow__(self, other: int) -> "Cyclic":
        # Exponentiation: multiply exponent n by scalar power
        if not isinstance(other, int):
            raise TypeError("Exponent must be an integer.")
        return Cyclic(self.n * other, self.k)

    def __truediv__(self, other: "Cyclic") -> "Cyclic":
        # Division: multiply by the inverse
        if self.k != other.k:
            raise ValueError("Cyclic elements must belong to the same group.")
        return self * other.inverse

    def __str__(self) -> str:
        return f"a({self.n})"

    def __repr__(self) -> str:
        return f"a({self.n})"


class EllipticCurve:
    def __init__(self, a: GF.GFe, b: GF.GFe):
        """Create an elliptic curve over GF(p^n) defined by parameters a, b."""
        # Curve defined by y² = x³ + ax + b over GF(p^n), p > 2
        # Curve defined by y² + xy = x³ + ax² + b over GF(p^n), p = 2
        if not isinstance(a, GF.GFe) or not isinstance(b, GF.GFe):
            raise TypeError("Curve parameters must be instances of GF.GFe.")

        self.gf = a.gf
        self.a = a
        self.b = b

        if a.gf != b.gf:
            raise ValueError("a and b must be in the same field")
        if self.is_singular:
            raise ValueError(
                f"The chosen curve is singular over GF({self.gf.p}^{self.gf.n})")

        if self.a.p == 3:
            raise NotImplementedError(
                "Characteristic 3 curves not yet supported")

    @property
    def is_singular(self) -> bool:
        """Check whether the curve is singular.""" 
        # Check singularity based on discriminant
        if self.a.p == 2:
            return self.b == self.gf([0])
        Δ = self.gf([4]) * self.a**3 + self.gf([27]) * self.b**2
        return Δ == self.gf([0])

    @property
    def infinity(self):
        """Return the point at infinity (identity element)."""
        # Identity of the additive group - point at infinity
        return (math.inf, math.inf)

    @property
    def points_poly(self) -> List[Tuple[GF.GFe, GF.GFe]]:
        """Return all affine (x, y) points satisfying the EC equation.""" 
        # Return all affine points (x, y) on the curve over GF(p^n),
        # where x, y ∈ GF(p^n), including the point at infinity

        # Point at infinity at start so .pop(0) can be used if needed
        points = [self.infinity]

        # Cartesian product of all field elements (x, y)
        all_pts = itertools.product(self.gf.elements, repeat=2)

        for x, y in all_pts:
            if self.a.p == 2:
                lhs = y**2 + x * y
                rhs = x**3 + self.a * x**2 + self.b
            else:
                lhs = y**2
                rhs = x**3 + self.a * x + self.b

            if lhs == rhs:
                points.append((x, y))

        return points

    @property
    def points_alpha(self) -> List[Tuple[int, int]]:
        """Return the affine points as powers of alpha."""
        # Return points as powers of alpha: coordinates (x.a, y.a)
        # Mostly for easy reading or exporting
        # Same functionality as points_poly method

        points = [self.infinity]
        all_pts = itertools.product(self.gf.elements, repeat=2)

        for x, y in all_pts:
            if self.a.p == 2:
                lhs = y**2 + x * y
                rhs = x**3 + self.a * x**2 + self.b
            else:
                lhs = y**2
                rhs = x**3 + self.a * x + self.b

            if lhs == rhs:
                points.append((x.a, y.a))

        return points

    @property
    def group_structure(self):
        """Determine and return the abelian group structure of the EC (e.g., Z/nZ or Z/mZ × Z/nZ)."""
        # Take the amount of points on the curve
        pts = [self.ECpoint(self, x, y) for x, y in self.points_poly]
        N = len(pts)
        structure = []

        # If N is prime then an abelian group is cyclic of order N
        if is_prime(N):
            for pt in pts:
                if pt.order == N:
                    point = pt
                    break
            string_value = f"Z/{N}Z"
            for i in range(len(pts)):
                structure.append(point * i)
            return string_value, structure, N

        # Keep track of max order computed
        max_order = 0
        point1 = None

        # Find the order of each point
        for P in pts:
            o = P.order
            if o == N:
                # If a point has order N exit early
                string_value = f"Z/{N}Z (cyclic)"
                point = P
                for i in range(len(pts)):
                    structure.append(point * i)
                return string_value, structure, N
            if o > max_order:
                max_order = o
                point1 = P

        # Return the group structure
        m = N // max_order

        point1_sub = point1.subgroup

        for pt in pts:
            if pt.order == m and pt not in point1_sub:
                point2 = pt
                break

        string_value = f"Z/{m}Z × Z/{max_order}Z"
        for i in range(m):
            for j in range(max_order):
                structure.append(point1 * j + point2 * i)

        return string_value, structure, m, max_order

    def map_to_abelian(self) -> Tuple[dict["ECpoint", AG.AGE], dict[AG.AGE, "ECpoint"]]:
        """Return dictionaries mapping EC points to AG elements and vice versa."""
        group_structure = self.group_structure
        if len(group_structure) == 3:
            _, ec_elements, n = group_structure
            abelian_elements = AG(n).generate()
        else:
            _, ec_elements, m, n = group_structure
            sort = sorted((m, n))
            abelian_elements = AG(sort[0], sort[1]).generate()

        mapping = [i for i in zip(ec_elements, abelian_elements)]
        ec_to_ag = {ec: ag for ec, ag in mapping}
        ag_to_ec = {ag: ec for ec, ag in mapping}
        return ec_to_ag, ag_to_ec

    def torsion(self, n: int):
        """Return all torsion points of order n on the curve.""" 
        _, ag_ec = self.map_to_abelian()
        ag = next(iter(ag_ec)).group
        return [ag_ec[ag] for ag in ag.torsion(n)]

    def subgroups(self):
        """Return all subgroups of the EC based on its AG isomorphism."""
        _, ag_ec = self.map_to_abelian()
        ag = next(iter(ag_ec)).group
        ag_subgroups = ag.subgroups()
        ec_subgroups = []
        for subgroup in ag_subgroups:
            ec_subgroup = []
            for element in subgroup:
                ec_subgroup.append(ag_ec[element])
            ec_subgroups.append(ec_subgroup)
        return ec_subgroups

    @property
    def j_invariant(self):
        """Return the j-invariant of the elliptic curve."""
        if self.a.p == 2:
            if self.b == self.gf([0]):
                raise ValueError("j-invariant undefined for singular curve")
            # Return the j-invariant for char 2 curves
            return self.a**3 / self.b

        num = self.gf([1728]) * (self.gf([4]) * self.a**3)
        denom = self.gf([4]) * self.a**3 + self.gf([27]) * self.b**2
        if denom == self.gf([0]):
            raise ValueError("j-invariant undefined for singular curve")
        return num / denom

    def char_sum(self):
        """Sum all characters over the EC group (via AG isomorphism)."""
        # Compute the sum of a character over all points, for all characters
        _, ag_ec = self.map_to_abelian()
        ag = next(iter(ag_ec)).group
        return ag.sum_char()

    def all_signatures(self):
        """Return the character signature for each point on the elliptic curve."""
        # Compute the signature of all points on the curve
        signatures = []
        for x, y in self.points_poly:
            point = EllipticCurve.ECpoint(self, x, y)
            signatures.append(point.signature)
        return signatures

    def curve_char(self):
        a = self.a
        b = self.b
        N = len(self.points_poly)
        def χ(P):
            val = (a * P.x + b * P.y) / self.gf([N])
            return cmath.exp(tau * 1j * val)
        return χ

    class ECpoint:
        def __init__(self, EC: "EllipticCurve", x, y):
            """Create a point on the given elliptic curve from coordinates (x, y)."""
            # Point on an elliptic curve over GF(p^n)
            self.EC = EC
            self.a = EC.a
            self.b = EC.b
            self.x = x
            self.y = y

            # Handle point at infinity first
            if x == math.inf and y == math.inf:
                return
            
            if (x == math.inf or y == math.inf) and not (x == y == math.inf):
                raise ValueError("Invalid representation of point at infinity.")

            # Check coordinate types
            if not isinstance(x, GF.GFe) or not isinstance(y, GF.GFe):
                raise TypeError(
                    "x and y must be elements of GF(p^n), or math.inf")

            # Check ordinates are in the same field
            if x.gf != y.gf:
                raise ValueError(
                    "x and y must belong to the same finite field.")

            # Check that point satisfies the correct curve equation
            if self.EC.a.p == 2:
                lhs = y**2 + x * y
                rhs = x**3 + self.a * x**2 + self.b
            else:
                lhs = y**2
                rhs = x**3 + self.a * x + self.b

            if lhs != rhs:
                raise ValueError(f"Point ({x}, {y}) is not on the curve.")

        @cached_property
        def order(self) -> int:
            """Return the order of this point in the EC group."""  
            # Compute the order of the point (i.e., the smallest n such that nP = 𝒪)
            # by testing all divisors of the curve’s group order.
            points = self.EC.points_poly
            N = len(points)

            for d in divisors(N):
                if (self * d).is_inf():
                    return d

        @property
        def subgroup(self):
            """Return the cyclic subgroup generated by this point."""
            sg = []
            for i in range(self.order):
                sg.append(self * i)
            return sg

        @property
        def signature(self):
            """Return the AG character signature of this point.""" 
            ec_ag, ag_ec = self.EC.map_to_abelian()
            ag = next(iter(ag_ec)).group
            chars = ag.characters()
            signature = []
            for name, func in chars:
                signature.append((ag_ec[name], func(ec_ag[self])))
            return signature

        def __add__(self, other):

            # Point at infinity + Point at infinity = ∞
            if self.is_inf() and other.is_inf():
                return EllipticCurve.ECpoint(self.EC, math.inf, math.inf)

            # ∞ + P = P
            if self.is_inf():
                return EllipticCurve.ECpoint(self.EC, other.x, other.y)

            # P + ∞ = P
            if other.is_inf():
                return EllipticCurve.ECpoint(self.EC, self.x, self.y)

            # Special case: characteristic 2
            if self.EC.a.p == 2:
                if self == other:
                    # Doubling at x = 0 gives ∞
                    if self.x == self.EC.a.gf([0]):
                        return EllipticCurve.ECpoint(self.EC, math.inf, math.inf)

                    # Point doubling (char 2 form): y² + xy = x³ + ax² + b
                    m = self.x + self.y / self.x
                    x = m**2 + m + self.a
                    y = self.x**2 + (m + self.EC.a.gf([1])) * x
                    return EllipticCurve.ECpoint(self.EC, x, y)

                # Same x but y1 ≠ y2 → vertical line → ∞
                if self.x == other.x:
                    return EllipticCurve.ECpoint(self.EC, math.inf, math.inf)

                # Regular addition (char 2): P + Q
                m = (self.y + other.y) / (self.x + other.x)
                x = m**2 + m + self.a + self.x + other.x
                y = m * (self.x + x) + x + self.y
                return EllipticCurve.ECpoint(self.EC, x, y)

            # General case: x1 ≠ x2 → regular addition
            if self.x != other.x:
                m = (other.y - self.y) / (other.x - self.x)
                x = m**2 - self.x - other.x
                y = m * (self.x - x) - self.y
                return EllipticCurve.ECpoint(self.EC, x, y)

            # Same x, but y1 ≠ y2 → vertical line → ∞
            if self.y != other.y:
                return EllipticCurve.ECpoint(self.EC, math.inf, math.inf)

            # y == 0 → vertical tangent → ∞
            if self.y == self.EC.a.gf([0]):
                return EllipticCurve.ECpoint(self.EC, math.inf, math.inf)

            # Point doubling for p > 2
            m = ((self.EC.a.gf([3]) * self.x**2 + self.a) /
                 (self.EC.a.gf([2]) * self.y))

            x = m**2 - self.EC.a.gf([2]) * self.x
            y = m * (self.x - x) - self.y
            return EllipticCurve.ECpoint(self.EC, x, y)

        def __mul__(self, scalar: int) -> "EllipticCurve.ECpoint":
            # Scalar multiplication using double-and-add algorithm.
            # Computes scalar * point, where scalar ∈ ℕ and point ∈ EC.

            if not isinstance(scalar, int):
                raise TypeError("Can only multiply ECpoint by an integer.")

            # Start from the identity point (point at infinity)
            result = EllipticCurve.ECpoint(self.EC, math.inf, math.inf)
            addend = self

            while scalar > 0:
                if scalar & 1:
                    result = result + addend
                addend = addend + addend
                # shift right: floor of division by 2
                scalar >>= 1

            return result

        def __rmul__(self, scalar: int) -> "EllipticCurve.ECpoint":
            # Scalar multiplication using double-and-add algorithm.
            # Computes scalar * point, where scalar ∈ ℕ and point ∈ EC.

            if not isinstance(scalar, int):
                raise TypeError("Can only multiply ECpoint by an integer.")

            # Start from the identity point (point at infinity)
            result = EllipticCurve.ECpoint(self.EC, math.inf, math.inf)
            addend = self

            while scalar > 0:
                if scalar & 1:
                    result = result + addend
                addend = addend + addend
                # shift right: floor of division by 2
                scalar >>= 1

            return result

        def is_inf(self) -> bool:
            """Check if the point is the point at infinity."""
            # Return True if this is the point at infinity, indicated by x being math.inf.
            # Checking if float first to get round comparing GFe instance to math.inf
            return isinstance(self.x, float) and self.x == math.inf

        def __eq__(self, other: object) -> bool:
            # Check for equality by comparing x and y coordinates.
            # If other is not an ECpoint, return False.
            if isinstance(other, EllipticCurve.ECpoint):
                return self.x == other.x and self.y == other.y
            return False

        def __hash__(self):
            return hash((self.x, self.y))

        def __str__(self) -> str:
            # Return a human-readable string representation of the point.
            # This shows the coordinates, for example: "(x, y)".
            return f"({self.x}, {self.y})"

        def __repr__(self) -> str:
            return f"({self.x}, {self.y})"

    class FreeECPoint:
        def __init__(self, x: GF.GFe, y: GF.GFe):
            """Create a coordinate pair (x, y) not yet bound to a specific curve."""
            # A point (x, y) in GF(p^n) that is not yet bound to any curve.
            self.x = x
            self.y = y

            # Ensure x and y are in the same field
            if x.gf != y.gf:
                raise ValueError("x and y must be in the same finite field")

        def candidate_curves(self) -> List[List[GF.GFe]]:
            """Return all non-singular ECs this point lies on."""
            # Return all valid elliptic curves (a, b) that this point lies on
            # and are non-singular over GF(p^n)

            F = self.x.gf
            curves = []

            if self.x.p == 2:
                # Check y² + xy = x³ + ax² + b
                for a in F.elements:
                    b = self.y**2 + self.x * self.y + self.x**3 + a * self.x**2
                    # Must be non-singular
                    if b != self.x.gf([0]):
                        curves.append([a, b])
            else:
                # Check y² = x³ + ax + b
                for a in F.elements:
                    b = self.y**2 - self.x**3 - a * self.x
                    Δ = self.x.gf([4]) * a**3 + self.x.gf([27]) * b**2
                    # Must be non-singular
                    if Δ != self.x.gf([0]):
                        curves.append([a, b])

            return curves

        def singular_curves(self) -> List[List[GF.GFe]]:
            """Return all singular ECs this point lies on."""
            # Return all singular elliptic curves (a, b) that this point lies on

            F = self.x.gf
            curves = []

            if self.x.p == 2:
                for a in F.elements:
                    b = self.y**2 + self.x * self.y + self.x**3 + a * self.x**2
                    # Must be singular
                    if b == self.x.gf([0]):
                        curves.append([a, b])
            else:
                for a in F.elements:
                    b = self.y**2 - self.x**3 - a * self.x
                    Δ = self.x.gf([4]) * a**3 + self.x.gf([27]) * b**2
                    # Must be singular
                    if Δ == self.x.gf([0]):
                        curves.append([a, b])

            return curves


gf = GF(5, [1, 1, 2])
EC = EllipticCurve(gf([1]), gf([2]))
points = [EllipticCurve.ECpoint(EC, x, y) for x, y in EC.points_poly]
a = EC.a
b = EC.b
print(a, b)
print(a.signature)
print(b.signature)
for x, y in EC.points_poly[1:]:
    point = EllipticCurve.ECpoint(EC, x, y)
    
    print(x,y)

    print(x.signature)

    print(y.signature)
    
    print(point.signature)
    break

