from math import isqrt
from typing import Tuple


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
    """Round real and imaginary parts of a complex number, and snap to integer n, if abs(n - val) < epsilon. 
    Sorts out floating point errors in calculations."""
    def snap(x):
        nearest = round(x)
        return nearest if abs(x - nearest) < epsilon else round(x, dp)

    return complex(snap(z.real), snap(z.imag))


def is_prime(n: int) -> bool:
    """Check if the given value is a prime number."""
    if not isinstance(n, int):
        raise TypeError(f"Expected an integer, got {type(n).__name__}")

    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6

    return True
