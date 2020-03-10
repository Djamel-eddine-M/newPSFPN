# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import secrets
import math
import hashlib
import random
import time


class EllipticCurve3(object):
    """Represents a single elliptic curve defined over a finite field.

    See here:
        http://en.wikipedia.org/wiki/Elliptic_curve
        

    p must be prime, since we use the modular inverse to compute point
    addition.

    """
    def __init__(self, a, b, p):
        self.a = a%p
        self.b = b%p
        self.p = p

    def __eq__(self, C):
        return (self.a, self.b) == (C.a, C.b)

    def has_point(self, x, y):
        return (y ** 2) % self.p == (x ** 3 + self.a * x + self.b) % self.p

    def __str__(self):
        return 'y^2 = x^3 + {}x + {}'.format(self.a, self.b)



class Point(object):
    """A point on a specific curve. all points are associated with a curve
        so we have to test if this point are in the curve 
        and we have another class for the infinite point
    """
    def __init__(self, curve, x, y):
        self.curve = curve
        self.x = x % curve.p
        self.y = y % curve.p
        #we test if the point are on the curve
        if not self.curve.has_point(x, y):
            raise ValueError('{} is not on curve {}'.format(self, self.curve))

    def __str__(self):
        return '({}, {})'.format(self.x, self.y)

    def __getitem__(self, index):
        return [self.x, self.y][index]

    def __eq__(self, Q):
        return (self.curve, self.x, self.y) == (Q.curve, Q.x, Q.y)

    def __neg__(self):
        return Point(self.curve, self.x, -self.y)
    
    def __add__(self, Q):
        """Add two points together.

        for that we need to take care of special cases:
         * Q is the infinity point (0)
         * P == Q
         * The line crossing P and Q is vertical.

        """
        assert self.curve == Q.curve
        # 0 + P = P
        if isinstance(Q, Inf):
            return self

        xp, yp, xq, yq = self.x, self.y, Q.x, Q.y
        m = None

        # P == Q
        if self == Q:
            if self.y == 0:
                R = Inf(self.curve)
            else:
                m = ((3 * xp * xp + self.curve.a) * mod_inverse(2 * yp, self.curve.p)) % self.curve.p

        # Vertical line
        elif xp == xq:
            R = Inf(self.curve)

        # Common case
        else:
            m = ((yq - yp) * mod_inverse(xq - xp, self.curve.p)) % self.curve.p
        #if it's not the point at the infinite
        if m is not None:
            xr = (m ** 2 - xp - xq) % self.curve.p
            yr = (m * (xp - xr) - yp) % self.curve.p
            R = Point(self.curve, xr, yr)


        return R
    
    
    def __mul__(self, n):
        
        """ Multiplication of a point by a integer
           
        """
        assert isinstance(n, int)
        assert n >= 0

    
        if n == 0:
            return Inf(self.curve)
    
        else:
            Q = self
            R = Inf(self.curve)
    
            i = 1
            while i <= n:
                if n & i == i:
                    R = R + Q
    
                Q = Q + Q
    
                i = i << 1
    
        return R

    def __rmul__(self, n):
        return self * n
    
    
class Inf(Point):
    """The custom infinity point."""
    def __init__(self, curve):
        self.curve = curve

    def __eq__(self, Q):
        return isinstance(Q, Inf)

    def __neg__(self):
        """-0 = 0"""
        return self
    def __add__(self, Q):
        """P + 0 = P"""
        return Q
    
    def __str__(self):
        return 'Inf'

def mod_inverse(a, n):
    """Return the inverse of a mod n.

    n must be prime.

    >>> mod_inverse(42, 2017)    


    """
    b = n
    if abs(b) == 0:
        return (1, 0, a)

    x1, x2, y1, y2 = 0, 1, 1, 0
    while abs(b) > 0:
        q, r = divmod(a, b)
        x = x2 - q * x1
        y = y2 - q * y1
        a, b, x2, x1, y2, y1 = b, r, x1, x, y1, y

    return x2 % n

                    

def solution(n,p):
    a=pow(n, (p - 1) // 2, p)
    #print("critere deuler=")
    #print(a)
    if(a==1):
        return True
    else :
        return False
    



def access_bit(data, num):
    base = int(num // 8)
    shift = int(num % 8)
    return (data[base] & (1<<shift)) >> shift

def random_elliptiqueV4(p):
    #while we didnt find a curve
    log2 = math.log(p, 2.0)
    t=math.ceil(log2)
    s=math.floor((t-1)/160)
    h=t-(160*s)
    print("start")
    while(True):
        m=hashlib.sha1()
        seed=secrets.token_bytes(20) #arbitary bit string of 160bits
        m.update(seed)
        hash1=m.digest() #sha1 of seedE
        tmp=[access_bit(hash1,i) for i in range(len(hash1)*8)]#conversion en tableau de bits
        c0=tmp[0:h] #the h rightmost bit of hash1
        w0=c0
        w0[len(w0)-1]=0
        w0.reverse()
        #Forming the w it's the concatenation of the Wi 
        #with Wi=sha1((seed+1)mod 2^160)
        W=[]
        for i in range (s,0,-1):
            m=hashlib.sha1()
            seedE=((seed+i )%( 2 **160))
            m.update(seedE)
            hash2=m.digest()
            tmp=[access_bit(hash2,i) for i in range(len(hash2)*8)]#conversion en tableau de bits
            tmp.reverse()
            W=tmp+W
        W=w0+W
        
        r=0
        for i in range(t):
            r=r+(W[i]*(2**((t-1)-i)))
        r=r%p
        if(r==0):
            a=0
            b=random.randint(1,p-1)#car si comme a  vaut zero (car r vaut 0) =>b=0=>courbe super singulière
        else:
            a=random.randint(1,p-1)#car si a vaut zero =>b=0=>courbe super singulière
            rinv=mod_inverse(r,p)
            res=(pow(a,3,p)*rinv)%p
            #print("r vaut:")
            #print(r)
            #print("a vaut:")
            #print(a)
            #print("rinv vaut")
            #print(rinv)
            #print("res_carre ")
            #print(res)


            while(solution(res,p)==False):#au plus une itration en moyenne
                #print("non")
                a=random.randint(0,p-1)
                res=(pow(a,3,p)*rinv)%p
                #print("a vaut:")
                #print(a)
                #print("rinv vaut")
                #print(rinv)
                #print("res vaut")
                #print(res)

            b=tonelli(res,p)
        #if we pass this test it mean we have a good curve        
        if(((4* pow(a,3,p))+(27*pow(b,2,p))) % p != 0):
            break
        print("next")
    return seed,a,b

def random_point(a,seedE,b,p):
    log2 = math.log(p, 2.0)
    t=math.ceil(log2)
    s=math.floor((t-1)/160)
    h=t-(160*s)
    while(True):
        m=hashlib.sha1()
        seedp=secrets.token_bytes(20)
        m.update(seedp)
        hash1=m.digest()
        tmp=[access_bit(hash1,i) for i in range(len(hash1)*8)]#conversion en tableau de bits
        c0=tmp[0:h] #the h rightmost bit of hash1
        x0=c0
        x0[len(x0)-1]=0
        #print("h vaut:")
        #print(h)
        x0.reverse()
        X=[]
        for i in range (s,0,-1):
            m=hashlib.sha1()
            seed=((seedE+i )%( 2 **160))
            m.update(seed)
            hash2=m.digest()
            tmp=[access_bit(hash2,i) for i in range(len(hash2)*8)]#conversion en tableau de bits        Xbis=concatenation(X)
            tmp.reverse()
            X=tmp+X
        X=x0+X
        xu=0
        for i in range(t):
            xu=xu+(X[i]*(2**((t-1)-i)))
        
        temp=((xu**3)+(a*xu)+b)%p
        if(solution(temp,p)==True):
            break
    yu=tonelli(temp,p)
    u=Point(EllipticCurve3(a,b,p),xu,yu)   
    P=u*h
    return seedp,yu,P

#Algorithme de Tonelli-Shanks

##########################################
'''
#Complexite :O(ln(p)⁴) c.f.  
A course in computational algebraic number theory



https://rosettacode.org/wiki/Tonelli-Shanks_algorithm
Contenu du site est sous licence de documentation libre GNU 1.2
'''

def legendre3(a, p):
    return pow(a, (p - 1) // 2, p)

def tonelli(n, p):
    #if(n==0):
     #return 0
    assert legendre3(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    for z in range(2, p):
        if p - 1 == legendre3(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r
########################################


    
                        #Algrithme de schoof
'''
Complexite O(ln(p)⁸)
'''                                                                         
#######################################################################   
# -*- coding: utf-8 -*-
# Copyright (c) 2010--2012  Peter Dinges <pdinges@acm.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
An implementation of Schoof's algorithm for counting the points of
ptptic curves over finite fields; the implementation avoids unnecessary
computations.

The implementation is faster than the naive version. Yet, it still uses naive
arithmetic and is badly suited for actual counting tasks. It demonstrates the
performance improvement from a smarter choice of moduli and reordered
computations.

@see Schoof, R.
"Elliptic Curves Over Finite Fields and the Computation of Square Roots mod p"
in Mathematics of Computation, Vol. 44, No. 170 (Apr. 1985), pp. 483--494.

@package   schoof.reduced_computation
@author    Peter Dinges <pdinges@acm.org>
"""

from elliptic_curves.l_torsion_group.naive import LTorsionGroup
from support.primes import inverse_primorial, primes_range
from support.quotients import solve_congruence_equations, representative_in_range

def frobenius_trace(curve):
    """
    Compute the trace of the Frobenius endomorphism for the given EllpiticCurve
    @p curve.
    
    This is an implementation of Schoof's original algorithm for counting the
    points of an elliptic curve over a finite field.
    
    @return    The trace @f$ t @f$ of the Frobenius endomorphism. The number of
               points on the curve then is @f$ q + 1 - t @f$, where @f$ q @f$
               is the size of the finite field over which the curve was defined.
    """
    trace_congruences = []
    search_range = hasse_frobenius_trace_range( curve.field() )
    torsion_primes = greedy_prime_factors(
                                 len(search_range),
                                 curve.field().characteristic()
                             )
    
    # To avoid multivariate polynomial arithmetic, make l=2 a special case.
    if 2 in torsion_primes:
        trace_congruences.append( frobenius_trace_mod_2( curve ) )
        torsion_primes.remove( 2 )

    torsion_group = LTorsionGroup( curve )
    for prime in torsion_primes:
        trace_congruences.append(
                frobenius_trace_mod_l( torsion_group( prime ) )
             )
    
    trace_congruence = solve_congruence_equations( trace_congruences )
    return representative_in_range( trace_congruence, search_range )


from rings.integers.naive import Integers
from rings.polynomials.naive import Polynomials
from rings.quotients.naive import QuotientRing
from support.rings import gcd

def frobenius_trace_mod_2(curve):
    """
    Compute the trace of the Frobenius endomorphism modulo 2.
    
    We cannot use the implicit torsion group representation of
    frobenius_trace_mod_l() in the case @f$ l=2 @f$ because the second division
    polynomial is @f$ y @f$, and multivariate polynomial arithmetic is
    unavailable. Implementing it for this single case would be overkill.
    
    Instead, we test whether there are any rational 2-torsion points. If so,
    then @f$ E[2] @f$ is a subgroup of the @p curve and its order divides the
    order of the curve (which is the number of rational points). Now, the
    number of rational points in @f$ E[2] @f$ can only be 2 or 4 in this case, so
    we know that the number of rational points is even. Hence, the Frobenius
    trace must then be 0 modulo 2.
    
    @return    0 if the number of points on the curve is even, 1 otherwise.
               The result is a congruence class modulo 2, that is, an element
               of @c QuotientRing( Integers, 2 ).
    """
    R = Polynomials( curve.field() )
    
    x = R(0, 1)
    A, B = curve.parameters()

    defining_polynomial = x**3 + A*x + B
    rational_characteristic = x**curve.field().size() - x
    
    # gcd() returns a gcd, which may have any unit as leading coefficient.
    # For relatively prime polynomials, the gcd is constant. 
    if gcd( rational_characteristic, defining_polynomial ).degree() == 0:
        # The rational characteristic and the defining polynomial are
        # relatively prime. Thus there is no rational point of order 2.
        return QuotientRing( Integers, 2 )(1)
    else:
        return QuotientRing( Integers, 2 )(0)


from rings.integers.naive import Integers
from rings.quotients.naive import QuotientRing

def frobenius_trace_mod_l(torsion_group):
    """
    Compute the trace of the Frobenius endomorphism modulo @f$ l @f$, where
    @f$ l @f$ is the torsion of @p torsion_group.
    
    The function guesses candidates and verifies whether the function that
    results from applying the characteristic polynomial @f$ \chi_\phi @f$
    to @f$ \phi @f$ maps every point in the @p torsion_group onto the point
    at infinity.
    
    @note      A torsion of 2 requires multivariate polynomial arithmetic,
               which is unavailable.  Therefore @f$ l @f$ must be greater than 2.
               Use frobenius_trace_mod_2() to get the trace modulo 2.
    
    @return    The congruence class of the trace of the Frobenius endomorphism.
               This is an element of @c QuotientRing( Integers, l ). 
    """
    assert torsion_group.torsion() > 2, \
        "torsion 2 requires multivariate polynomial arithmetic"
        
    torsion_quotient_ring = QuotientRing( Integers, torsion_group.torsion() )
    field_size = torsion_group.curve().field().size()

    # Note: Technically, there could be several points so we would have to
    #       filter the one candidate that worked for all points in the end.
    #       Luckily, there is only one point.
    for point in torsion_group.elements():
        frobenius_point = frobenius( point, field_size )
        frobenius2_point = frobenius( frobenius_point, field_size )
        determinant_point = ( field_size % torsion_group.torsion() ) * point
        
        point_sum = frobenius2_point + determinant_point
        if point_sum.is_infinite():
            return torsion_quotient_ring( 0 )
        
        trace_point = frobenius_point
        for trace_candidate in range( 1, (torsion_group.torsion()+1) // 2 ):
            if point_sum.x() == trace_point.x():
                if point_sum.y() == trace_point.y():
                    return torsion_quotient_ring( trace_candidate )
                else:
                    return torsion_quotient_ring( -trace_candidate )
            else:
                trace_point += frobenius_point

    message = "Frobenius equation held for no trace candidate"
    raise ArithmeticError( message )


def frobenius(point, q):
    """
    The Frobenius endomorphism @f$ \phi @f$.
    
    @return    The point @f$ (x^q, y^q) @f$ if @p point is @f$ (x, y) @f$.
    """
    return point.__class__( point.x() ** q, point.y() ** q )


from math import ceil, sqrt

def hasse_frobenius_trace_range(field):
    """
    Return the interval in which the trace of the Frobenius endomorphism
    must reside (using Hasse's theorem).
    
    Hasse's theorem gives a limit for the trace @f$ t @f$ of the Frobenius
    endomorphism on an elliptic curve over a field with  @f$ q @f$ elements:
    @f$ \left|t\right| \leq 2\sqrt{q} @f$.
    
    @return    The interval that must contain the trace of the Frobenius
               endomorphism according to Hasse's theorem.
    """
    # This only depends on the field and holds for any curve over it. 
    l = 2 * ceil( sqrt( field.size() ) )
    return range( -l, l+1 )


from support.primes import primes_range, inverse_primorial

def greedy_prime_factors(n, shunned=0):
    """
    Return a list of the first primes whose product is greater than, or equal
    to @p n, but do not use @p shunned.
    
    For example, if @p n is 14, then the returned list will consist of 3 and
    5, but not 2, because 3 times 5 is greater than 14. The function behaves
    like inverse_primorial() except that it removes unnecessary smaller primes.
    
    @note      Canceling of unnecessary primes follows a greedy algorithm.
               Therefore the choice of primes might be suboptimal; perfect
               choice, however, is an NP-complete problem (KNAPSACK).
    
    @note      This function uses primes_range() to obtain a list of primes.
               See the notes there for use case limitations.
    """
    primes = primes_range( 2, n+1 )
    
    # Find the smallest product of primes that is at least n, but don't use
    # the shunned prime.
    product = 1
    for index, prime in enumerate( primes ):
        if prime != shunned:
            product *= prime
            if product >= n:
                break
    
    # Throw away excess primes
    primes = primes[ : index+1 ]
    if shunned in primes:
        primes.remove( shunned )
    
    # Try to cancel unnecessary primes, largest first.
    # (This greedy search is not optimal; however, we did not set out to solve
    # the KNAPSACK problem, did we?)
    for index, prime in enumerate( reversed( primes ) ):
        canceled_product = product / prime
        if canceled_product >= n:
            product = canceled_product
            primes[ -(index+1) ] = 0
    
    return list( filter( None, primes ) )
    

#------------------------------------------------------------------------------

from fields.finite.naive import FiniteField
from elliptic_curves.naive import EllipticCurve

import sys
from support.running import AlgorithmRunner

def reduced_computation_schoof_algorithm( p, A, B, output=sys.stdout ):
    p, A, B = int(p), int(A), int(B)
    
    message = "Counting points of y^2 = x^3 + {A}x + {B} over GF<{p}>: "
    print( message.format( p=p, A=A, B=B ), end="", file=output )
    output.flush()
    
    order = p + 1 - frobenius_trace( EllipticCurve( FiniteField(p), A, B ) )
    print( order, file=output )
    return order

####################################################################
#reduced_computation_schoof_algorithm(23,1,13)





                    #  Self-Initializing Quadratic Sieve 

####O(exp(sqrt(log n log log n)))  c.f. thèse de Scott Patrick 
### Factoring Integers withthe Self-Initializing Quadratic Sieve 
##########################################################

from math import sqrt, log2, ceil, floor
import random
from fractions import gcd
import sys
from builtins import ValueError


"""This script factorises a natural number given as a command line
parameter into its prime factors. It first attempts to use trial
division to find very small factors, then uses Brent's version of the
Pollard rho algorithm [1] to find slightly larger factors. If any large
factors remain, it uses the Self-Initializing Quadratic Sieve (SIQS) [2]
to factorise those.
[1] Brent, Richard P. 'An improved Monte Carlo factorization algorithm.'
    BIT Numerical Mathematics 20.2 (1980): 176-184.
[2] Contini, Scott Patrick. 'Factoring integers with the self-
    initializing quadratic sieve.' (1997).
"""

# Some tuning parameters
MAX_DIGITS_POLLARD = 30
POLLARD_QUICK_ITERATIONS = 20
MIN_DIGITS_POLLARD_QUICK2 = 45
POLLARD_QUICK2_ITERATIONS = 25
SIQS_TRIAL_DIVISION_EPS = 25
SIQS_MIN_PRIME_POLYNOMIAL = 400
SIQS_MAX_PRIME_POLYNOMIAL = 4000

# Number of iterations for the Miller-Rabin primality test
MILLER_RABIN_ITERATIONS = 50


class Polynomial:
    """A polynomial used for the Self-Initializing Quadratic Sieve."""

    def __init__(self, coeff=[], a=None, b=None):
        self.coeff = coeff
        self.a = a
        self.b = b

    def eval(self, x):
        res = 0
        for a in self.coeff[::-1]:
            res *= x
            res += a
        return res


class FactorBasePrime:
    """A factor base prime for the Self-Initializing Quadratic Sieve."""

    def __init__(self, p, tmem, lp):
        self.p = p
        self.soln1 = None
        self.soln2 = None
        self.tmem = tmem
        self.lp = lp
        self.ainv = None


def lowest_set_bit(a):
    b = (a & -a)
    low_bit = -1
    while (b):
        b >>= 1
        low_bit += 1
    return low_bit


def to_bits(k):
    """Return a generator that returns the bits of k, starting from the
    least significant bit, using True for 1s, and False for 0s.
    """
    k_binary = bin(k)[2:]
    return (bit == '1' for bit in k_binary[::-1])


def pow_mod(a, k, m):
    """Return a^k mod m."""
    r = 1
    b = a
    for bit in to_bits(k):
        if bit:
            r = (r * b) % m
        b = (b * b) % m
    return r


def is_quadratic_residue(a, p):
    """Return whether a is a quadratic residue modulo a prime p."""
    return legendre(a, (p - 1) // 2, 1, p) == 1


def legendre(a, q, l, n):
    x = q ** l
    if x == 0:
        return 1

    z = 1
    a %= n

    while x != 0:
        if x % 2 == 0:
            a = (a ** 2) % n
            x //= 2
        else:
            x -= 1
            z = (z * a) % n
    return z


def sqrt_mod_prime(a, p):
    """Return the square root of a modulo the prime p. Behaviour is
    undefined if a is not a quadratic residue mod p."""
    # Algorithm from http://www.mersennewiki.org/index.php/Modular_Square_Root
    assert a < p
    assert is_probable_prime(p)
    if a == 0:
        return 0
    if p == 2:
        return a
    if p % 2 == 0:
        return None
    p_mod_8 = p % 8
    if p_mod_8 == 1:
        # Shanks method
        q = p // 8
        e = 3
        while q % 2 == 0:
            q //= 2
            e += 1
        while True:
            x = random.randint(2, p - 1)
            z = pow_mod(x, q, p)
            if pow_mod(z, 2 ** (e - 1), p) != 1:
                break
        y = z
        r = e
        x = pow_mod(a, (q - 1) // 2, p)
        v = (a * x) % p
        w = (v * x) % p
        while True:
            if w == 1:
                return v
            k = 1
            while pow_mod(w, 2 ** k, p) != 1:
                k += 1
            d = pow_mod(y, 2 ** (r - k - 1), p)
            y = (d ** 2) % p
            r = k
            v = (d * v) % p
            w = (w * y) % p
    elif p_mod_8 == 5:
        v = pow_mod(2 * a, (p - 5) // 8, p)
        i = (2 * a * v * v) % p
        return (a * v * (i - 1)) % p
    else:
        return pow_mod(a, (p + 1) // 4, p)


def inv_mod(a, m):
    """Return the modular inverse of a mod m."""
    return eea(a, m)[0] % m


def eea(a, b):
    """Solve the equation a*x + b*y = gcd(a,b).
    Return (x, y, +/-gcd(a,b)).
    """
    if a == 0:
        return (0, 1, b)
    x = eea(b % a, a)
    return (x[1] - b // a * x[0], x[0], x[2])


def is_probable_prime(a):
    """Perform the Miller-Rabin primality test to determine whether the
    given number a is a prime. Return True if the number is a prime
    with very high probability, and False if it is definitely composite.
    """
    if a == 2:
        return True
    if a == 1 or a % 2 == 0:
        return False
    return primality_test_miller_rabin(a, MILLER_RABIN_ITERATIONS)


def primality_test_miller_rabin(a, iterations):
    m = a - 1
    lb = lowest_set_bit(m)
    m >>= lb

    for _ in range(iterations):
        b = random.randint(2, a - 1)
        j = 0
        z = pow_mod(b, m, a)
        while not ((j == 0 and z == 1) or z == a - 1):
            if (j > 0 and z == 1 or j + 1 == lb):
                return False
            j += 1
            z = (z * z) % a

    return True


def siqs_factor_base_primes(n, nf):
    """Compute and return nf factor base primes suitable for a Quadratic
    Sieve on the number n.
    """
    global small_primes
    factor_base = []
    for p in small_primes:
        if is_quadratic_residue(n, p):
            t = sqrt_mod_prime(n % p, p)
            lp = round(log2(p))
            factor_base.append(FactorBasePrime(p, t, lp))
            if len(factor_base) >= nf:
                break
    return factor_base


def siqs_find_first_poly(n, m, factor_base):
    """Compute the first of a set of polynomials for the Self-
    Initialising Quadratic Sieve.
    """
    p_min_i = None
    p_max_i = None
    for i, fb in enumerate(factor_base):
        if p_min_i is None and fb.p >= SIQS_MIN_PRIME_POLYNOMIAL:
            p_min_i = i
        if p_max_i is None and fb.p > SIQS_MAX_PRIME_POLYNOMIAL:
            p_max_i = i - 1
            break

    # The following may happen if the factor base is small, make sure
    # that we have enough primes.
    if p_max_i is None:
        p_max_i = len(factor_base) - 1
    if p_min_i is None or p_max_i - p_min_i < 20:
        if( p_min_i is not None ):
            p_min_i=5

    target = sqrt(2 * float(n)) / m
    target1 = target / ((factor_base[p_min_i].p + 
                         factor_base[p_max_i].p) / 2) ** 0.5

    # find q such that the product of factor_base[q_i] is approximately
    # sqrt(2 * n) / m; try a few different sets to find a good one
    best_q, best_a, best_ratio = None, None, None
    for _ in range(30): 
        a = 1
        q = []

        while a < target1:
            p_i = 0
            while p_i == 0 or p_i in q:
                p_i = random.randint(p_min_i, p_max_i)
            p = factor_base[p_i].p
            a *= p
            q.append(p_i)

        ratio = a / target

        # ratio too small seems to be not good
        if (best_ratio is None or (ratio >= 0.9 and ratio < best_ratio) or
                    best_ratio < 0.9 and ratio > best_ratio):
            best_q = q
            best_a = a
            best_ratio = ratio
    a = best_a
    q = best_q

    s = len(q)
    B = []
    for l in range(s):
        fb_l = factor_base[q[l]]
        q_l = fb_l.p
        assert a % q_l == 0
        gamma = (fb_l.tmem * inv_mod(a // q_l, q_l)) % q_l
        if gamma > q_l // 2:
            gamma = q_l - gamma
        B.append(a // q_l * gamma)

    b = sum(B) % a
    b_orig = b
    if (2 * b > a):
        b = a - b

    assert 0 < b
    assert 2 * b <= a
    assert ((b * b - n) % a == 0)

    g = Polynomial([b * b - n, 2 * a * b, a * a], a, b_orig)
    h = Polynomial([b, a])
    for fb in factor_base:
        if a % fb.p != 0:
            fb.ainv = inv_mod(a, fb.p)
            fb.soln1 = (fb.ainv * (fb.tmem - b)) % fb.p
            fb.soln2 = (fb.ainv * (-fb.tmem - b)) % fb.p

    return g, h, B


def siqs_find_next_poly(n, factor_base, i, g, B):
    """Compute the (i+1)-th polynomials for the Self-Initialising
    Quadratic Sieve, given that g is the i-th polynomial.
    """
    v = lowest_set_bit(i) + 1
    z = -1 if ceil(i / (2 ** v)) % 2 == 1 else 1
    b = (g.b + 2 * z * B[v - 1]) % g.a
    a = g.a
    b_orig = b
    if (2 * b > a):
        b = a - b
    assert ((b * b - n) % a == 0)

    g = Polynomial([b * b - n, 2 * a * b, a * a], a, b_orig)
    h = Polynomial([b, a])
    for fb in factor_base:
        if a % fb.p != 0:
            fb.soln1 = (fb.ainv * (fb.tmem - b)) % fb.p
            fb.soln2 = (fb.ainv * (-fb.tmem - b)) % fb.p

    return g, h


def siqs_sieve(factor_base, m):
    """Perform the sieving step of the SIQS. Return the sieve array."""
    sieve_array = [0] * (2 * m + 1)
    for fb in factor_base:
        if fb.soln1 is None:
            continue
        p = fb.p
        i_start_1 = -((m + fb.soln1) // p)
        a_start_1 = fb.soln1 + i_start_1 * p
        lp = fb.lp
        if p > 20:
            for a in range(a_start_1 + m, 2 * m + 1, p):
                sieve_array[a] += lp

            i_start_2 = -((m + fb.soln2) // p)
            a_start_2 = fb.soln2 + i_start_2 * p
            for a in range(a_start_2 + m, 2 * m + 1, p):
                sieve_array[a] += lp
    return sieve_array


def siqs_trial_divide(a, factor_base):
    """Determine whether the given number a can be fully factorised into
    primes from the factors base. If so, return the indices of the
    factors from the factor base. If not, return None.
    """
    divisors_idx = []
    for i, fb in enumerate(factor_base):
        if a % fb.p == 0:
            exp = 0
            while a % fb.p == 0:
                a //= fb.p
                exp += 1
            divisors_idx.append((i, exp))
        if a == 1:
            return divisors_idx
    return None


def siqs_trial_division(n, sieve_array, factor_base, smooth_relations, g, h, m,
                        req_relations):
    """Perform the trial division step of the Self-Initializing
    Quadratic Sieve.
    """
    sqrt_n = sqrt(float(n))
    limit = log2(m * sqrt_n) - SIQS_TRIAL_DIVISION_EPS
    for (i, sa) in enumerate(sieve_array):
        if sa >= limit:
            x = i - m
            gx = g.eval(x)
            divisors_idx = siqs_trial_divide(gx, factor_base)
            if divisors_idx is not None:
                u = h.eval(x)
                v = gx
                assert (u * u) % n == v % n
                smooth_relations.append((u, v, divisors_idx))
                if (len(smooth_relations) >= req_relations):
                    return True
    return False


def siqs_build_matrix(factor_base, smooth_relations):
    """Build the matrix for the linear algebra step of the Quadratic Sieve."""
    fb = len(factor_base)
    M = []
    for sr in smooth_relations:
        mi = [0] * fb
        for j, exp in sr[2]:
            mi[j] = exp % 2
        M.append(mi)
    return M


def siqs_build_matrix_opt(M):
    """Convert the given matrix M of 0s and 1s into a list of numbers m
    that correspond to the columns of the matrix.
    The j-th number encodes the j-th column of matrix M in binary:
    The i-th bit of m[i] is equal to M[i][j].
    """
    m = len(M[0])
    cols_binary = [""] * m
    for mi in M:
        for j, mij in enumerate(mi):
            cols_binary[j] += "1" if mij else "0"
    return [int(cols_bin[::-1], 2) for cols_bin in cols_binary], len(M), m


def add_column_opt(M_opt, tgt, src):
    """For a matrix produced by siqs_build_matrix_opt, add the column
    src to the column target (mod 2).
    """
    M_opt[tgt] ^= M_opt[src]


def find_pivot_column_opt(M_opt, j):
    """For a matrix produced by siqs_build_matrix_opt, return the row of
    the first non-zero entry in column j, or None if no such row exists.
    """
    if M_opt[j] == 0:
        return None
    return lowest_set_bit(M_opt[j])


def siqs_solve_matrix_opt(M_opt, n, m):
    """
    Perform the linear algebra step of the SIQS. Perform fast
    Gaussian elimination to determine pairs of perfect squares mod n.
    Use the optimisations described in [1].
    [1] Koç, Çetin K., and Sarath N. Arachchige. 'A Fast Algorithm for
        Gaussian Elimination over GF (2) and its Implementation on the
        GAPP.' Journal of Parallel and Distributed Computing 13.1
        (1991): 118-122.
    """
    row_is_marked = [False] * n
    pivots = [-1] * m
    for j in range(m):
        i = find_pivot_column_opt(M_opt, j)
        if i is not None:
            pivots[j] = i
            row_is_marked[i] = True
            for k in range(m):
                if k != j and (M_opt[k] >> i) & 1:  # test M[i][k] == 1
                    add_column_opt(M_opt, k, j)
    perf_squares = []
    for i in range(n):
        if not row_is_marked[i]:
            perfect_sq_indices = [i]
            for j in range(m):
                if (M_opt[j] >> i) & 1:  # test M[i][j] == 1
                    perfect_sq_indices.append(pivots[j])
            perf_squares.append(perfect_sq_indices)
    return perf_squares


def siqs_calc_sqrts(square_indices, smooth_relations):
    """Given on of the solutions returned by siqs_solve_matrix_opt and
    the corresponding smooth relations, calculate the pair [a, b], such
    that a^2 = b^2 (mod n).
    """
    res = [1, 1]
    for idx in square_indices:
        res[0] *= smooth_relations[idx][0]
        res[1] *= smooth_relations[idx][1]
    res[1] = sqrt_int(res[1])
    return res


def sqrt_int(n):
    """Return the square root of the given integer, rounded down to the
    nearest integer.
    """
    a = n
    s = 0
    o = 1 << (floor(log2(n)) & ~1)
    while o != 0:
        t = s + o
        if a >= t:
            a -= t
            s = (s >> 1) + o
        else:
            s >>= 1
        o >>= 2
    assert s * s == n
    return s


def kth_root_int(n, k):
    """Return the k-th root of the given integer n, rounded down to the
    nearest integer.
    """
    u = n
    s = n + 1
    while u < s:
        s = u
        t = (k - 1) * s + n // pow(s, k - 1)
        u = t // k
    return s


def siqs_factor_from_square(n, square_indices, smooth_relations):
    """Given one of the solutions returned by siqs_solve_matrix_opt,
    return the factor f determined by f = gcd(a - b, n), where
    a, b are calculated from the solution such that a*a = b*b (mod n).
    Return f, a factor of n (possibly a trivial one).
    """
    sqrt1, sqrt2 = siqs_calc_sqrts(square_indices, smooth_relations)
    assert (sqrt1 * sqrt1) % n == (sqrt2 * sqrt2) % n
    return gcd(abs(sqrt1 - sqrt2), n)


def siqs_find_factors(n, perfect_squares, smooth_relations):
    """Perform the last step of the Self-Initialising Quadratic Field.
    Given the solutions returned by siqs_solve_matrix_opt, attempt to
    identify a number of (not necessarily prime) factors of n, and
    return them.
    """
    factors = []
    rem = n
    non_prime_factors = set()
    prime_factors = set()
    for square_indices in perfect_squares:
        fact = siqs_factor_from_square(n, square_indices, smooth_relations)
        if fact != 1 and fact != rem:
            if is_probable_prime(fact):
                if fact not in prime_factors:
                    print ("SIQS: Prime factor found: %d" % fact)
                    prime_factors.add(fact)

                while rem % fact == 0:
                    factors.append(fact)
                    rem //= fact

                if rem == 1:
                    break
                if is_probable_prime(rem):
                    factors.append(rem)
                    rem = 1
                    break
            else:
                if fact not in non_prime_factors:
                    print ("SIQS: Non-prime factor found: %d" % fact)
                    non_prime_factors.add(fact)

    if rem != 1 and non_prime_factors:
        non_prime_factors.add(rem)
        for fact in sorted(siqs_find_more_factors_gcd(non_prime_factors)):
            while fact != 1 and rem % fact == 0:
                print ("SIQS: Prime factor found: %d" % fact)
                factors.append(fact)
                rem //= fact
            if rem == 1 or is_probable_prime(rem):
                break

    if rem != 1:
        factors.append(rem)
    return factors


def siqs_find_more_factors_gcd(numbers):
    res = set()
    for n in numbers:
        res.add(n)
        for m in numbers:
            if n != m:
                fact = gcd(n, m)
                if fact != 1 and fact != n and fact != m:
                    if fact not in res:
                        print("SIQS: GCD found non-trivial factor: %d" % fact)
                        res.add(fact)
                    res.add(n // fact)
                    res.add(m // fact)
    return res


def siqs_choose_nf_m(d):
    """Choose parameters nf (sieve of factor base) and m (for sieving
    in [-m,m].
    """
    # Using similar parameters as msieve-1.52
    if d <= 34:
        return 200, 65536
    if d <= 36:
        return 300, 65536
    if d <= 38:
        return 400, 65536
    if d <= 40:
        return 500, 65536
    if d <= 42:
        return 600, 65536
    if d <= 44:
        return 700, 65536
    if d <= 48:
        return 1000, 65536
    if d <= 52:
        return 1200, 65536
    if d <= 56:
        return 2000, 65536 * 3
    if d <= 60:
        return 4000, 65536 * 3
    if d <= 66:
        return 6000, 65536 * 3
    if d <= 74:
        return 10000, 65536 * 3
    if d <= 80:
        return 30000, 65536 * 3
    if d <= 88:
        return 50000, 65536 * 3
    if d <= 94:
        return 60000, 65536 * 9
    return 100000, 65536 * 9


def siqs_factorise(n):
    """Use the Self-Initializing Quadratic Sieve algorithm to identify
    one or more non-trivial factors of the given number n. Return the
    factors as a list.
    """
    dig = len(str(n))
    nf, m = siqs_choose_nf_m(dig)

    factor_base = siqs_factor_base_primes(n, nf)

    required_relations_ratio = 1.05
    success = False
    smooth_relations = []
    prev_cnt = 0
    i_poly = 0
    while not success:
        print("*** Step 1/2: Finding smooth relations ***")
        required_relations = round(len(factor_base) * required_relations_ratio)
        print("Target: %d relations" % required_relations)
        enough_relations = False
        while not enough_relations:
            if i_poly == 0:
                g, h, B = siqs_find_first_poly(n, m, factor_base)
            else:
                g, h = siqs_find_next_poly(n, factor_base, i_poly, g, B)
            i_poly += 1
            if i_poly >= 2 ** (len(B) - 1):
                i_poly = 0
            sieve_array = siqs_sieve(factor_base, m)
            enough_relations = siqs_trial_division(
                n, sieve_array, factor_base, smooth_relations,
                g, h, m, required_relations)

            if (len(smooth_relations) >= required_relations or
                i_poly % 8 == 0 and len(smooth_relations) > prev_cnt):
                print("Total %d/%d relations." %
                      (len(smooth_relations), required_relations))
                prev_cnt = len(smooth_relations)

        print("*** Step 2/2: Linear Algebra ***")
        print("Building matrix for linear algebra step...")
        M = siqs_build_matrix(factor_base, smooth_relations)
        M_opt, M_n, M_m = siqs_build_matrix_opt(M)

        print("Finding perfect squares using matrix...")
        perfect_squares = siqs_solve_matrix_opt(M_opt, M_n, M_m)

        print("Finding factors from perfect squares...")
        factors = siqs_find_factors(n, perfect_squares, smooth_relations)
        if len(factors) > 1:
            success = True
        else:
            print("Failed to find a solution. Finding more relations...")
            required_relations_ratio += 0.05

    return factors


def check_factor(n, i, factors):
    while n % i == 0:
        n //= i
        factors.append(i)
        if is_probable_prime(n):
            factors.append(n)
            n = 1
    return n


def trial_div_init_primes(n, upper_bound):
    """Perform trial division on the given number n using all primes up
    to upper_bound. Initialise the global variable small_primes with a
    list of all primes <= upper_bound. Return (factors, rem), where
    factors is the list of identified prime factors of n, and rem is the
    remaining factor. If rem = 1, the function terminates early, without
    fully initialising small_primes.
    """
    print("Trial division and initialising small primes...")
    global small_primes
    is_prime = [True] * (upper_bound + 1)
    is_prime[0:2] = [False] * 2
    factors = []
    small_primes = []
    max_i = sqrt_int(upper_bound)
    rem = n
    for i in range(2, max_i + 1):
        if is_prime[i]:
            small_primes.append(i)
            rem = check_factor(rem, i, factors)
            if rem == 1:
                return factors, 1

            for j in (range(i ** 2, upper_bound + 1, i)):
                is_prime[j] = False

    for i in range(max_i + 1, upper_bound + 1):
        if is_prime[i]:
            small_primes.append(i)
            rem = check_factor(rem, i, factors)
            if rem == 1:
                return factors, 1

    print("Primes initialised.")
    return factors, rem


def pollard_brent_f(c, n, x):
    """Return f(x) = (x^2 + c)%n. Assume c < n.
    """
    x1 = (x * x) % n + c
    if x1 >= n:
        x1 -= n
    assert x1 >= 0 and x1 < n
    return x1


def pollard_brent_find_factor(n, max_iter=None):
    """Perform Brent's variant of the Pollard rho factorisation
    algorithm to attempt to a non-trivial factor of the given number n.
    If max_iter > 0, return None if no factors were found within
    max_iter iterations.
    """
    y, c, m = (random.randint(1, n - 1) for _ in range(3))
    r, q, g = 1, 1, 1
    i = 0
    while g == 1:
        x = y
        for _ in range(r):
            y = pollard_brent_f(c, n, y)
        k = 0
        while k < r and g == 1:
            ys = y
            for _ in range(min(m, r - k)):
                y = pollard_brent_f(c, n, y)
                q = (q * abs(x - y)) % n
            g = gcd(q, n)
            k += m
        r *= 2
        if max_iter:
            i += 1
            if (i == max_iter):
                return None

    if g == n:
        while True:
            ys = pollard_brent_f(c, n, ys)
            g = gcd(abs(x - ys), n)
            if g > 1:
                break
    return g


def pollard_brent_quick(n, factors):
    """Perform up to max_iter iterations of Brent's variant of the
    Pollard rho factorisation algorithm to attempt to find small
    prime factors. Restart the algorithm each time a factor was found.
    Add all identified prime factors to factors, and return 1 if all
    prime factors were found, or otherwise the remaining factor.
    """
    rem = n
    while True:
        if is_probable_prime(rem):
            factors.append(rem)
            rem = 1
            break

        digits = len(str(n))
        if digits < MIN_DIGITS_POLLARD_QUICK2:
            max_iter = POLLARD_QUICK_ITERATIONS
        else:
            max_iter = POLLARD_QUICK2_ITERATIONS

        f = pollard_brent_find_factor(rem, max_iter)
        if f and f < rem:
            if is_probable_prime(f):
                print("Pollard rho (Brent): Prime factor found: %s" % f)
                factors.append(f)
                assert rem % f == 0
                rem //= f
            else:
                print("Pollard rho (Brent): Non-prime factor found: %s" % f)
                rem_f = pollard_brent_quick(f, factors)
                rem = (rem // f) * rem_f
        else:
            print("No (more) small factors found.")
            break
    return rem


def check_perfect_power(n):
    """Check if the given integer is a perfect power. If yes, return
    (r, b) such that r^b == n. If no, return None. Assume that
    global small_primes has already been initialised and that n does
    not have any prime factors from small_primes.
    """
    largest_checked_prime = small_primes[-1]
    for b in small_primes:
        bth_root = kth_root_int(n, b)
        if bth_root < largest_checked_prime:
            break
        if (bth_root ** b == n):
            return (bth_root, b)
    return None


def find_prime_factors(n):
    """Return one or more prime factors of the given number n. Assume
    that n is not a prime and does not have very small factors, and that
    the global small_primes has already been initialised. Do not return
    duplicate factors.
    """

    print("Checking whether %d is a perfect power..." % n)
    perfect_power = check_perfect_power(n)
    if perfect_power:
        print("%d is %d^%d" % (n, perfect_power[0], perfect_power[1]))
        factors = [perfect_power[0]]
    else:
        print("Not a perfect power.")
        digits = len(str(n))
        if digits <= MAX_DIGITS_POLLARD:
            print("Using Pollard rho (Brent's variant) to factorise %d (%d digits)..."
                  % (n, digits))
            factors = [pollard_brent_find_factor(n)]
        else:
            print("Using Self-Initializing Quadratic Sieve to factorise" +
                  " %d (%d digits)..." % (n, digits))
            factors = siqs_factorise(n)

    prime_factors = []
    for f in set(factors):
        for pf in find_all_prime_factors(f):
            prime_factors.append(pf)

    return prime_factors


def find_all_prime_factors(n):
    """Return all prime factors of the given number n. Assume that n
    does not have very small factors and that the global small_primes
    has already been initialised.
    """
    rem = n
    factors = []

    while rem > 1:
        if is_probable_prime(rem):
            factors.append(rem)
            break

        for f in find_prime_factors(rem):
            print("Prime factor found: %d" % f)
            assert is_probable_prime(f)
            assert rem % f == 0
            while rem % f == 0:
                rem //= f
                factors.append(f)

    return factors


def product(factors):
    """Return the product of all numbers in the given list."""
    prod = 1
    for f in factors:
        prod *= f
    return prod


def factorise(n):
    """Factorise the given integer n >= 1 into its prime factors."""

    if type(n) != int or n < 1:
        raise ValueError("Number needs to be an integer >= 1")

    print("Factorising %d (%d digits)..." % (n, len(str(n))))
    if n == 1:
        return []

    if is_probable_prime(n):
        return [n]

    factors, rem = trial_div_init_primes(n, 1000000)

    if factors:
        print("Prime factors found so far: %s" % factors)
    else:
        print("No small factors found.")

    if rem != 1:
        digits = len(str(rem))
        if digits > MAX_DIGITS_POLLARD:
            print("Attempting quick Pollard rho (Brent's variant) to find slightly " +
                  "larger factors...")
            rem = pollard_brent_quick(rem, factors)
        if rem > 1:
            for fr in find_all_prime_factors(rem):
                factors.append(fr)

    factors.sort()
    assert product(factors) == n
    for p in factors:
        assert is_probable_prime(p)
    return factors
######################################################################################
def listeCouple(l):
    lc=[]
    stable=l[0]
    k=0
    tmp=stable
    cpt=0
    while k!=len(l):
        while (tmp==stable):
            cpt=cpt+1
            k=k+1
            if k<len(l):
                tmp=l[k]
            else:
                break
        lc.append((stable,cpt))
        stable=tmp
        cpt=0
    return lc



def OrderPoint(P,OderCurve,facto):
    m=OderCurve
    for i in range(len(facto)):
        p,e=facto[i]
        m=m//(p**e)
        print("m vaut")
        print(m)
        Q=m*P
        while isinstance(Q, Inf)==False:
            Q=p*Q
            m=m*p
    return m

P=1009
seed,a,b=random_elliptiqueV4(P)
courbe=EllipticCurve3(a,b,P)  
print(courbe)
sed,yu,point1=random_point(a,seed,b,P)
print(point1)
Ordercurve=reduced_computation_schoof_algorithm(P,a,b)
print("le nombre de point de la courbe vautvaut")

print(Ordercurve)


ll=factorise(Ordercurve)
l=listeCouple(ll)
print("l vaut")
print(l)
k=OrderPoint(point1,Ordercurve,l)
if isinstance(k*point1, Inf):
    print("OrderPoint fonctionne!!!!")
print("k vaut")
print(k)
'''
print("")
#version naive
def order(P,p):
    if(isinstance(P,Inf)):
        return 1
    else :
        h=P
        i=2
    while(True):
        h=h+P
        print(h)
        if(isinstance(h,Inf)):
            return i
        i=i+1

print(order(point1,P))
print(k)
print(k*point1)

#algorithme rho pour resoudre ECDLP
def rho_probability(P,A,B,Q,p):
    prob=random.randint(1,9999)%3
    if(prob %3 == 0):
        A=2*A % p
        B=2*B % p
        if(A==0):
            A=A+1
        if(B==0):
            B=B+1
        Z=(P*A)+(Q*B)
        return Z,A,B
    elif(prob % 3 == 1):
        A=A+1 % p
        if(A==0):
            A=A+1
        Z=(P*A)+(Q*B)
        return Z,A,B
    else :
        B=B+1 % p
        if(B==0):
            B=B+1
        Z=(P*A)+(Q*B)
        return Z,A,B
    
  

def rho_man(P,Q,p):
    ordre=order(P,p)
    Ai=random.randint(1,ordre-1)
    Bi=random.randint(1,ordre-1)
    Zi=(P*Ai)+(Q*Bi)
    Zibis=Zi
    Aibis=Ai
    Bibis=Bi
    i=0
    while(True):
        Zi,Ai,Bi=rho_probability(P,Ai,Bi,Q,ordre)
        Zibis,Aibis,Bibis=rho_probability(P,Aibis,Bibis,Q,ordre)
        Zibis,Aibis,Bibis=rho_probability(P,Aibis,Bibis,Q,ordre)
        print('Zi,ZIBIS,Ai,Bi,Aibis,Bibis')
        print(Zi,Zibis,Ai,Bi,Aibis,Bibis)
        print(i)
        i=i+1
        if((Bi % ordre) != (Bibis % ordre)):
            if(Zi==Zibis) :
                break
    return ((Ai-Aibis) *mod_inverse(-( Bi-Bibis),ordre)) % ordre 



'''        
    



