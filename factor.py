import math
import random
from algorithms import isPrime, modSqrt
from PositiveDefiniteForm import PositiveDefiniteForm

def factor(N: int) -> tuple[int, int]:
    
    # ensure N is positive
    N = abs(N)

    # ensure N is odd
    if N % 2 == 0:
        return 2, 0

    # check if N is a perfect square
    sqrt = int(round(math.sqrt(N)))
    if sqrt * sqrt == N:
        return sqrt, 0

    # check that N is composite
    if isPrime(N):
        return N, 0

    # choose the appropriate discriminant
    D = -N if -N % 4 == 1 else -4 * N

    # store square roots of forms in the 2-Sylow subgroup
    S: dict[str, PositiveDefiniteForm] = {}

    # for testing purposes, stores the number of BSGS computations made, which is the most expensive part of the computation
    BSGS_count = 0

    # sample primes until a factorisation is found
    while True:
        p: int
        while True:
            p = random.randint(3, abs(D))
            if isPrime(p):
                break

        # if p divides D or the primality test gave a false positive, we may get a shared factor between p and D, which gives a non-primitive form
        x = math.gcd(abs(D), p)
        if x != 1:
            if not x in [2, 4, D, 2 * D, 4 * D]:
                while x % 2 == 0:
                    x //= 2
                return x, BSGS_count
            p //= x

        # find b which squares to D mod p (if it exists)
        b = modSqrt(D, p)
        if b == None:
            continue
        if (b*b - D) % 4 != 0:
            b += p
        if ((b*b - D) % (4 * p)) != 0:
            continue

        # compute the corresponding form and its order
        f = PositiveDefiniteForm(p, b, (b * b - D) // (4 * p)).reduce()
        d = f.order()
        BSGS_count += 1

        # if d is even compute the corresponding ambiguous form
        if d % 2 != 0:
            continue
        g = pow(f, d // 2)

        # extract a factor from the ambiguous form
        x = (2 * g.a) - g.b if g.a == g.c else g.a
        if not x in [1, 2, 4, D, 2 * D, 4 * D]:
            while x % 2 == 0:
                x //= 2
            return x, BSGS_count

        # if the factorisation is trivial, compute d = 2^A * B
        A = 0
        B = d
        while B % 2 == 0:
            B //= 2
            A += 1

        # h is in the 2-sylow subgroup with order 2^A
        h = pow(f, B)
            
        for _ in range(A):
            h2 = h * h
            if str(h2) in S.keys():
                g = h.inverse() * S[str(h2)]
                x = (2 * g.a) - g.b if g.a == g.c else g.a
                if not x in [1, 2, 4, D, 2 * D, 4 * D]:
                    while x % 2 == 0:
                        x //= 2
                    return x, BSGS_count
            else:
                S[str(h2)] = h
            h = h2

# Simple CLI
if __name__ == "__main__":
    while True:
        print("\nEnter an integer to factor")
        print("N = ", end="")
        N = input()

        try:
            N = int(N)
        except:
            print("N must be in integer")
            continue

        primeFactors = []
        compositeFactors = [N]
        while len(compositeFactors) > 0:
            current = compositeFactors[0]
            compositeFactors = compositeFactors[1:]
            if isPrime(current):
                primeFactors.append(current)
            else:
                x, _ = factor(current)
                compositeFactors.append(x)
                if x != current:
                    compositeFactors.append(current // x)
        primeFactors.sort()
        condensed = []
        while len(primeFactors) > 0:
            currentPrime = primeFactors[0]
            currentPower = primeFactors.count(currentPrime)
            primeFactors = primeFactors[currentPower:]
            condensed.append(str(currentPrime) if currentPower == 1 else f"{currentPrime}^{currentPower}")
        print(f"N = {" * ".join(condensed)}")