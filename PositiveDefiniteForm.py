import math
from algorithms import gcdExtended

class PositiveDefiniteForm:

    @staticmethod
    def Identity(discriminant):
        if discriminant % 4 not in (0, 1):
            raise ValueError("Invalid quadratic form discriminant")
        if discriminant % 4 == 0:
            return PositiveDefiniteForm(1, 0, -discriminant // 4)
        else:
            return PositiveDefiniteForm(1, 1, (1 - discriminant) // 4)

    def __init__(self, a, b, c):
        self.a = int(a)
        self.b = int(b)
        self.c = int(c)

        self.discriminant = self.b**2 - 4*self.a*self.c

        self._validate()

    def _validate(self):
        if self.a <= 0:
            raise ValueError("a must be positive")

        if self.discriminant >= 0:
            raise ValueError("Discriminant must be negative")

        if self.discriminant % 4 not in (0, 1):
            raise ValueError("Invalid quadratic form discriminant")

        if math.gcd(math.gcd(abs(self.a), abs(self.b)), abs(self.c)) != 1:
            raise ValueError("Form must be primitive")

    def is_reduced(self):
        return (
            abs(self.b) <= self.a <= self.c
            and not (
                (abs(self.b) == self.a or self.a == self.c)
                and self.b < 0
            )
        )

    def is_ambiguous(self):
        f = self.reduce()
        return (f.b == 0) or (f.a == f.b) or (f.a == f.c)

    def reduce(self):
        f = PositiveDefiniteForm(self.a, self.b, self.c)
        while f.is_reduced() == False:

            # T-reduction
            k = (f.a - f.b) // (2 * f.a)
            f.c = (k * k * f.a) + (k * f.b) + f.c
            f.b = (2 * k * f.a) + f.b

            # S-reduction
            if (f.a > f.c) or (f.b < 0 and f.a == f.c):
                f.a, f.c = f.c, f.a
                f.b = -f.b

        return f

    # Multiplication is composition of forms
    def __mul__(self, other) -> PositiveDefiniteForm:

        if not isinstance(other, PositiveDefiniteForm):
            return NotImplemented

        if self.discriminant != other.discriminant:
            raise ValueError("Forms must have the same discriminant")

        D = self.discriminant

        f1 = self.reduce()
        f2 = other.reduce()

        a1, b1 = f1.a, f1.b
        a2, b2 = f2.a, f2.b

        B = (b1 + b2) // 2

        u1, v1, n1 = gcdExtended(a1, a2)
        z, w, n = gcdExtended(n1, B)
        u = z * u1
        v = z * v1

        m = (a1 * b2 * u) + (a2 * b1 * v) + (((b1 * b2) + D) * w // 2)

        a = a1 * a2 // (n**2)
        b = m // n
        c = (m**2 - ((n**2) * D)) // (4 * a1 * a2)

        result = PositiveDefiniteForm(a, b, c)
        return result.reduce()

    def __pow__(self, value: int):
        g = self
        h = PositiveDefiniteForm.Identity(self.discriminant)
        for i in range(math.ceil(math.log2(value)) + 1):
            if (value >> i) % 2 == 1:
                h *= g
            g = g * g
        return h

    def inverse(self) -> PositiveDefiniteForm:
        f = PositiveDefiniteForm(self.a, -self.b, self.c)
        return f.reduce()

    # Baby Step Giant Step algorithm to compute order
    def order(self, bound = None) -> int:

        # bound on the class number
        H = bound if bound != None else math.ceil(math.sqrt(abs(self.discriminant)) * (math.log(abs(self.discriminant)) + 2) / math.pi)

        # baby steps
        b = self
        babySteps = { str(b) : 1 }
        for i in range(2, math.ceil(math.sqrt(H)) + 1):
            b *= self
            if b.a == 1:
                return i
            babySteps[str(b)] = i

        # giant steps
        G = PositiveDefiniteForm(b.a, -b.b, b.c)
        g = G
        for i in range(1, math.ceil(math.sqrt(H)) + 1):
            if str(g) in babySteps.keys():
                return babySteps[str(g)] + math.ceil(math.sqrt(H)) * i
            g *= G

    # Equality is defined up to equivalence
    def __eq__(self, other):
        if not isinstance(other, PositiveDefiniteForm):
            return False

        self_reduced = self.reduce()
        other_reduced = other.reduce()

        return (
            self_reduced.a == other_reduced.a
            and self_reduced.b == other_reduced.b
            and self_reduced.c == other_reduced.c
        )

    def __repr__(self):
        return f"PositiveDefiniteForm(a={self.a}, b={self.b}, c={self.c}, D={self.discriminant})"

    def __str__(self):
        return f"({self.a}, {self.b}, {self.c})"