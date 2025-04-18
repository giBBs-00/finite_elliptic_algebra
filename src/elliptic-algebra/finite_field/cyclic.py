from functools import cached_property


class Cyclic:
    def __init__(self, n: int, k: int):
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
