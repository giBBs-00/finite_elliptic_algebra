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


class Prime(int):
    """
    A class to represent a prime number.

    This class ensures that the number is a prime integer.
    It inherits from the built-in `int` class.

    :example:
        >>> a = Prime(7)
        >>> print(a)
        7
        >>> b = Prime(4)  # Raises ValueError
        ValueError: 4 is not a prime number
        >>> c = Prime(3.5)  # Raises TypeError
        TypeError: Expected an integer, got float
    """

    @staticmethod
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

    def __new__(cls, value: int):
        """
        Create a new Prime instance.

        :param value: The value to be checked and stored.
        :raises TypeError: If the value is not an integer.
        :raises ValueError: If the value is not a prime number.
        """
        if not cls.is_prime(value):
            raise ValueError(f"{value} is not a prime number")
        return super().__new__(cls, value)
