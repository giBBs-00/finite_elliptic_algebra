
from elliptic_algebra import EllipticCurve, GF

# ---------------------------------------------------
# ELLIPTIC CURVE OVER GF(p^n): FULL WALKTHROUGH DEMO
# ---------------------------------------------------
# This is a full demo of how to make an elliptic curve over GF(5²),
# explore its structure, print all its points, analyze torsion, and list subgroups.
# It's minimal. It's meant to be remixed.

def main():
    # Define GF(5^2) using aa primitive irreducible polynomial x^2 + x + 2
    F = GF(5, [1, 1, 2])  # x² + x + 1

    # Define the curve: y² = x³ + (z + 1)x + z, for z being the root of our primitive polynomial
    a = F([1, 1])
    b = F([1, 0])
    E = EllipticCurve(a, b)

    print(f"Elliptic curve over GF(5^2) with a = {a}, b = {b}")
    print("Points on the curve:")
    for point in E.points_poly:
        print(point)

    print("\nGroup structure:")
    structure = E.group_structure
    print(structure[0])  # structure name like "Z/36Z" or "Z/6Z x Z/6Z"

    print("\nTorsion points of order 3:")
    torsion_pts = E.torsion(3)
    for pt in torsion_pts:
        print(pt)

    print("\nCyclic subgroups:")
    subgroups = E.subgroups()
    for subgroup in subgroups:
        print(subgroup)


if __name__ == "__main__":
    main()
