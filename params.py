#!/usr/bin/env python
from estimator import *

from copy import copy
from math import ceil, inf, log2, prod, sqrt
from sage.all import is_prime


def chunck(q, w):
    chuncks = []

    cnt = ceil(len(q) / w)
    for i in range(w):
        chuncks.append(q[i * cnt:(i + 1) * cnt])

    return chuncks


def genprime(init, n, reverse=False):
    p = init - (init % (2 * n)) + 2 * n + 1
    assert p > init and p % (2 * n) == 1

    while not is_prime(p):
        if reverse:
            p -= 2 * n
        else:
            p += 2 * n

    return p


def level0(Binit, q, P, drop, check=False):
    # copy args
    Binit = Binit.copy()

    Bsum = Binit.copy()
    for _ in range(1, N):
        Bsum.add(Binit)
    Bsum.ms(q, drop, check)

    return Bsum,


def level1(Bsum, q, P, drop, check=False):
    # copy args
    Bsum = Bsum.copy()

    # input MAC
    Bmac = Bsum.copy()
    Bmac.mul(Bsum).flood()

    # matrix multiplication
    BA1 = Bsum.copy()
    BA2 = Bsum.copy()
    BB = Bsum.copy()

    # initial rotation
    BA2.ks(q, P, check)

    # first multiplication
    Bmat = BA2.copy()
    Bmat.mul(BB)

    # other multiplications
    for _ in range(1, d):
        BA1.ks(q, P, check)
        BA2.ks(q, P, check)
        BB.ks(q, P, check)

        Bphi0 = BA1.copy().const()
        Bphi1 = BA2.copy().const()
        Bphi = Bphi0.add(Bphi1)
        Bphi.ms(q, drop, check)

        Bpsi = BB.copy()
        Bpsi.ms(q, drop, check)

        Bmul = Bphi.mul(Bpsi)
        Bmat.add(Bmul)

    return Bsum, Bmac, Bmat


def level2(Bsum, Bmac, Bmat, q, P, drop, check=False):
    # copy args
    Bsum = Bsum.copy()
    Bmac = Bmac.copy()
    Bmat = Bmat.copy()

    Bmat.ks(q, P, check)
    Bmat.ms(q, drop, check)

    Bmac2 = Bmat.copy()
    Bmac2.mul(Bsum).flood()

    Bmat.flood()

    return Bmac, Bmat, Bmac2


N = 5
d = 128
n = 2**15
t = genprime(2**128, n)
eflood = 2**80

bits = 61
mods = 3
levels = [level0, level1, level2]
drops = [0] * len(levels)
w = 3
msmany = False
ksmany = True
verbose = False


class Bound:
    # builtin
    def __init__(self, n, t, w, eflood):
        self.n = n
        self.t = t
        self.w = w
        self.eflood = eflood
        self.value = 0
        self.enc()

    def __gt__(self, other):
        if type(self) == type(other):
            return self.value > other.value
        else:
            return self.value > other

    def __lt__(self, other):
        if type(self) == type(other):
            return self.value < other.value
        else:
            return self.value < other

    def __repr__(self):
        return str(self.value)

    # util
    def copy(self):
        return copy(self)

    def log2(self):
        return log2(self.value)

    def mods(self, bits):
        if self.log2() == inf:
            return inf
        return ceil(self.log2() / bits)

    # noise bounds
    def enc(self):
        self.value = self.t * sqrt(3 * self.n * (126 * self.n + 127))
        return self

    def dec(self, q):
        return self.value < prod(q) // 2

    def add(self, other):
        self.value += other.value
        return self

    def mul(self, other):
        self.value *= other.value
        return self

    def const(self):
        self.value *= self.t * sqrt(3 * self.n)
        return self

    def ms(self, q, drop, check):
        if drop == 0:
            return self

        delta = self.t / 2 * sqrt(3 * self.n * (2 * self.n + 4))
        if msmany:
            qdrop = [prod(q[:drop])]
            delta *= sqrt(drop)
        else:
            qdrop = q[:drop]

        new = self.value
        for i, mod in enumerate(qdrop):
            if new + mod * delta >= prod(q[i:]) // 2:
                if check:
                    assert False
                return self

            new /= mod
            new += delta

        self.value = new
        return self

    def ks(self, q, P, check):
        qmax = max(map(lambda x: prod(map(float, x)), chunck(q, self.w)))
        bv = lambda cnt: qmax * self.t * self.n * sqrt(63 * cnt / 2)

        addTrue = bv(len(q) + len(P) + self.w - 1)
        addFalse = bv(len(q) + self.w)

        if len(P) == 0:
            self.value += addFalse
            return self

        delta = self.t / 2 * sqrt(3 * self.n * (2 * self.n + 4))
        if ksmany:
            Pdrop = [prod(P)]
            delta *= sqrt(len(P))
        else:
            Pdrop = P

        new = addTrue
        for i, mod in enumerate(Pdrop):
            if self.value + new + mod * delta >= prod(q + P[i:]) // 2:
                if check:
                    assert False
                self.value += addFalse
                return self

            new /= mod
            new += delta

        self.value += new
        return self

    def flood(self):
        self.value += self.t * eflood
        return self


P = []
depth = 0
while depth < len(levels):
    Binit = Bound(n, t, w, eflood)
    B = Binit.copy(),

    # generate q
    q, p = [], 2**(bits - 1)
    for _ in range(mods):
        p = genprime(p, n)
        q.append(p)

    Pmods = len(P)
    P, p = [], q[-1]
    for _ in range(Pmods):
        p = genprime(p, n)
        P.append(p)
    p = genprime(p, n)

    i = 0
    while i <= depth:
        # improve current level values
        #   1) can we improve P for hybrid key switching?
        #   2) can we improve drop?
        drop = drops[i]
        dropped = sum(drops[:i])
        qlvl = q[dropped:]
        Bcur = levels[i](*B, qlvl, P, drop)

        if verbose:
            print(f"[>] i:      {i}/{depth}/{len(levels)}")
            print(f"[>] q:      {len(qlvl)}/{len(q)}")
            print(f"[>] drops:  {drops}")
            print(f"[>] P:      {len(P)}")
            print(f"[>] Bcur:   {Bcur}")
            print()

        # improve P
        Bnew = levels[i](*B, qlvl, P + [p], drop)
        if max(Bnew).mods(bits) < max(Bcur).mods(bits):
            P.append(p)
            p = genprime(p, n)
            B, i = (Binit,), 0
            continue

        # improve drops
        Bnew = levels[i](*B, qlvl, P, drop + 1)
        if max(Bnew).mods(bits) + drop < max(Bcur).mods(bits):
            drops[i] = drop + 1
            B, i = (Binit,), 0
            continue

        # can we decrypt?
        if not max(Bcur).dec(qlvl[drop:]):
            mods += 1
            B, i = (Binit,), 0
            break

        B = Bcur
        i += 1
    else:
        depth += 1

if verbose:
    print(f"[>] i:      {i}/{depth}/{len(levels)}")
    print(f"[>] q:      {len(qlvl)}/{len(q)}")
    print(f"[>] drops:  {drops}")
    print(f"[>] P:      {len(P)}")
    print(f"[>] Bcur:   {Bcur}")
    print()

CBD = nd.NoiseDistribution.CenteredBinomial
DG = nd.NoiseDistribution.DiscreteGaussian
Xs, Xe = CBD(1), CBD(21)
params = LWE.Parameters(n, prod(q + P), Xs, Xe)
est = log2(LWE.primal_usvp(params).rop)
# est  = log2(min([x.rop for x in LWE.estimate(params).values()]))

print(f"qbits:  {mods * bits}")
print(f"mods:   {mods}")
print(f"bits:   {bits}")
print(f"w:      {w}")
print(f"drops:  {drops}")
print(f"chunck: {list(map(lambda x: len(x) * bits, chunck(q, w)))}")
print(f"est:    {est:.2f}")
print()
print(f"n:      {n}")
print(f"t:      {t}")
print(f"q:      {q}")
print(f"P:      {P}")

# final check
B = Bound(n, t, w, eflood),
if verbose:
    print()
    print(f"[>] B: {B}")

for i, level in enumerate(levels):
    dropped = sum(drops[:i])
    qlvl = q[dropped:]
    B = level(*B, q, P, drops[i], check=True)
    if verbose:
        print(f"[>] B: {B}")
assert max(B).dec(q[sum(drops):])
