import math
from cmath import exp, tau
from functools import cached_property

from elliptic_algebra.utils import divisors


class AG:
    def __init__(self, *orders):
        """Initialize an abelian group of given orders (i.e., Z/mZ × Z/nZ)."""
        self.orders = tuple(sorted(orders))
        self.rank = len(orders)
        self.order = math.prod(orders)

    @cached_property
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

    @cached_property
    def characters(self):
        """Return all characters χ_g: G → ℂ* as functions, for all g in G."""
        group_elements = self.generate()
        character_functions = []

        for a in group_elements:
            def make_char(a):
                def χ(g):
                    val = 0
                    for ai, gi, ni in zip(a.coords, g.coords, self.orders):
                        val += (ai * gi) / ni
                    return exp(tau * 1j * val)
                return χ
            char_func = make_char(a)
            character_functions.append((a, char_func))

        return character_functions

    def sum_char(self):
        """Compute the sum of each character χ_g(h) over all elements g,h in the group."""
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
            """Create an element of the abelian group Z\nZ x Z\mZ, given by its coordinates. Can be 1 or 2 dimensional"""
            self.group = group
            self.coords = tuple(coords)

            if not isinstance(group, AG):
                raise TypeError(
                    "An abelian group element must live in an abelian group")
            if len(coords) != len(group.orders):
                raise ValueError("Coordinate count does not match group rank.")

        @property
        def order(self):
            """Return the order of this element (smallest n such that n*g = 0)."""
            divs = divisors(self.group.order)
            for div in divs:
                if div * self == AG.AGE(self.group, *([0] * len(self.coords))):
                    return div

        @property
        def signature(self):
            """Return the character signature: values χ_g(h) for all g in G and a given h in G."""
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


class Cyclic:
    def __init__(self, n: int, k: int):
        """Deal with multiplicaive cyclic groups i.e. modular reduction is done by the exponent."""
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
