import math
import itertools
from typing import List, Tuple
from finite_field.field import GF
from finite_field.abelian import AG
from utils.helpers import divisors
from prime import is_prime


class EllipticCurve:
    def __init__(self, a: GF.GFe, b: GF.GFe):
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
        # Check singularity based on discriminant
        if self.a.p == 2:
            return self.b == self.gf([0])
        Δ = self.gf([4]) * self.a**3 + self.gf([27]) * self.b**2
        return Δ == self.gf([0])

    @property
    def infinity(self):
        # Identity of the additive group - point at infinity
        return (math.inf, math.inf)

    @property
    def points_poly(self) -> List[Tuple[GF.GFe, GF.GFe]]:
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
        # Return points as raw integer coordinates (x.a, y.a)
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
        for pt in pts:
            if pt.order == m:
                point2 = pt
                break
        string_value = f"Z/{m}Z × Z/{max_order}Z"
        for i in range(m):
            for j in range(max_order):
                structure.append(point1 * j + point2 * i)
        return string_value, structure, m, max_order

    def map_to_abelian(self):
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
        _, ag_ec = self.map_to_abelian()
        ag = next(iter(ag_ec)).group
        return [ag_ec[ag] for ag in ag.torsion(n)]

    def subgroups(self):
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
        # Compute the sum of a character over all points - sanity check as should be 0 for all other than char_0
        _, ag_ec = self.map_to_abelian()
        ag = next(iter(ag_ec)).group
        ag_elemets = ag.generate()
        ag_characters = ag.characters()
        char_sums = []
        count = 0
        for char in ag_characters:
            char_sum = 0
            for i in ag_elemets:
                char_sum += char(i)
            char_sums.append(
                f"char_{count}={char_sum.real:.4f} + {char_sum.imag:.4f}i")
            count += 1
        return char_sums

    class ECpoint:
        def __init__(self, EC: "EllipticCurve", x, y):
            # Point on an elliptic curve over GF(p^n)
            self.EC = EC
            self.a = EC.a
            self.b = EC.b
            self.x = x
            self.y = y

            # Handle point at infinity first
            if x == math.inf and y == math.inf:
                return

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

        @property
        def order(self) -> int:
            # Compute the order of the point (i.e., the smallest n such that nP = 𝒪)
            # by testing all divisors of the curve’s group order.
            points = self.EC.points_poly
            N = len(points)

            for d in divisors(N):
                if (self * d).is_inf():
                    return d

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
            # A point (x, y) in GF(p^n) that is not yet bound to any curve.
            self.x = x
            self.y = y

            # Ensure x and y are in the same field
            if x.gf != y.gf:
                raise ValueError("x and y must be in the same finite field")

        def candidate_curves(self) -> List[List[GF.GFe]]:
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
