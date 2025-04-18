from typing import Tuple
from math import isqrt

def is_prime(n: int) -> bool:
    """
    Check if the given value is a prime number.

    :param n: The value to check. Must be an integer.
    :return: True if the value is prime, False otherwise.
    :raises TypeError: If the input is not an integer.
    """
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


def divisors(n):
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
    # Base case: gcd(a, 0) = a
    if b == 0:
        return a, 1, 0

    # Recursive case: compute gcd and coefficients
    gcd, x1, y1 = extended_gcd(b, a % b)

    # Update x and y
    x = y1
    y = x1 - (a // b) * y1

    return gcd, x, y
