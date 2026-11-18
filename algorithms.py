import random

# Euclid's extended algorithm
# Implementation of Cohen's Algorithm 1.3.6
def gcdExtended(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return 1, 0, a

    # Step 1: Initialise
    u = 1
    d = a
    v1 = 0
    v3 = b

    # Step 2: Finished?
    while v3 != 0:

        # Step 3: Euclidean step
        q = d // v3
        t3 = d % v3
        t1 = u - q * v1
        u = v1
        d = v3
        v1 = t1
        v3 = t3

    return u, (d - a * u) // b, d

# Tonelli-Shanks square root mod prime p algorithm
# Implementation of Cohen's Algorithm 1.5.1
def modSqrt(a: int, p: int):
    a %= p
    if a == 0:
        return 0

    # Writing p - 1 as 2^e * q
    q = p - 1
    e = 0
    while q % 2 == 0:
        q //= 2
        e += 1

    # Step 1: Find a quadratic non-residue n
    n: int
    while True:
        n = random.randint(1, p)
        if pow(n, (p - 1) // 2, p) == p - 1:
            break
    z = pow(n, q, p)

    # Step 2: Initialize
    y = z
    r = e
    x = pow(a, (q - 1) // 2, p)
    b = (a * x * x) % p
    x = (a * x) % p

    # Step 3: Find exponent
    while b % p != 1:
        m = 0
        temp = b
        while temp != 1 and m < r:
            temp = (temp * temp) % p
            m += 1
        if m == r:
            return None

        # Step 4: Reduce exponent
        t = pow(y, 1 << (r - m - 1), p)
        y = (t * t) % p
        r = m % p
        x = (x * t) % p
        b = (b * y) % p

    return x

# Miller-Rabin primality test
# Implementation of Cohen's Algorithm 8.2.2
def isPrime(p: int):
    if p < 1:
        raise ValueError("p must be a positive integer")
    if p in [2, 3]:
        return True
    if p % 2 == 0 or p == 1:
        return False
    
    # Step 1: Initialise
    q = p - 1
    t = 0
    while q % 2 == 0:
        q //= 2
        t += 1
    c = 20

    # Step 2: Choose a new a
    while c > 0:
        a = random.randint(2, p - 1)
        e = 0
        b = pow(a, q, p)

        # Step 3: Squarings
        if b != 1:
            while b != 1 and b != p - 1 and e <= t - 2:
                b = pow(b, 2, p)
                e += 1
            if b != p - 1:
                return False

        # Step 4: Decrease counter
        c -= 1
    return True