
import secrets
import math
import hashlib
import random
import time
import shlex, subprocess



'''
This code is inspired by:
    
Title: elliptic-curves-finite-fields
Author: Jeremy Kun, Software Engineer at Google
Date: 2017
Availability: https://github.com/j2kun/elliptic-curves-finite-fields
'''
#####################################################################
class EllipticCurve3(object):
   
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
#####################################################################
'''
Code personelle:
'''                  

def solution(n,p):
    """
    Return True if n is square in Fp.
    p must be prime.
    """
    
    a=pow(n, (p - 1) // 2, p)
    if(a==1):
        return True
    else :
        return False
    



def access_bit(data, num):
    """
    Return the n-th bit of data.
    """
    
    base = int(num // 8)
    shift = int(num % 8)
    return (data[base] & (1<<shift)) >> shift
 
    
def random_elliptiqueV4(p):
    '''
    Generate a random elliptic curve on Fp
    by using Certicom Algorithm Page 15
    see: https://www.certicom.com/content/dam/certicom/images/pdfs/challenge-2009.pdf
    
    '''
    log2 = math.log(p, 2.0)
    t=math.ceil(log2)
    s=math.floor((t-1)/160)
    h=t-(160*s)
    #print("start")
    
    
    while(True):
        
        
        m=hashlib.sha1()
        seed=secrets.token_bytes(20) #arbitary bit string of 160bits
        m.update(seed)
        hash1=m.digest() 
        tmp=[access_bit(hash1,i) for i in range(len(hash1)*8)]
        c0=tmp[0:h] #the h rightmost bit of hash1
        w0=c0
        w0[len(w0)-1]=0
        w0.reverse()
        W=[]
        
        
        for i in range (s,0,-1):
            m=hashlib.sha1()
            seedE=((seed+i )%( 2 **160))
            m.update(seedE)
            hash2=m.digest()
            tmp=[access_bit(hash2,i) for i in range(len(hash2)*8)]
            tmp.reverse()
            W=tmp+W
        W=w0+W
        
        
        
        r=0
        for i in range(t):
            r=r+(W[i]*(2**((t-1)-i)))
        r=r%p
        
        
        
        if(r==0):
            a=0
            b=random.randint(1,p-1)#because r=0=>a=0 =>[b=0=>courbe super singulière]
        else:
            a=random.randint(1,p-1)#because [a=0 and r!=0 ]=>b=0=>courbe super singulière
            rinv=mod_inverse(r,p)
            res=(pow(a,3,p)*rinv)%p
            


            while(solution(res,p)==False):#at most one iteration on average
                a=random.randint(0,p-1)
                res=(pow(a,3,p)*rinv)%p
                

            b=tonelli(res,p)
        
        
        #if we pass this test it mean we have a good curve        
        if(((4* pow(a,3,p))+(27*pow(b,2,p))) % p != 0):
            break
        #print("next")
    return seed,a,b

def random_point(a,seedE,b,p):
    '''
    Generate a random point on a elleptic curve on Fp
    by using Certicom Algorithm Page 16
    see: https://www.certicom.com/content/dam/certicom/images/pdfs/challenge-2009.pdf
    '''
    
    
    log2 = math.log(p, 2.0)
    t=math.ceil(log2)
    s=math.floor((t-1)/160)
    h=t-(160*s)
    while(True):
        
        
        
        m=hashlib.sha1()
        seedp=secrets.token_bytes(20)
        m.update(seedp)
        hash1=m.digest()
        tmp=[access_bit(hash1,i) for i in range(len(hash1)*8)]
        c0=tmp[0:h] #the h rightmost bit of hash1
        x0=c0
        x0[len(x0)-1]=0
        x0.reverse()
        X=[]
        
        
        
        
        for i in range (s,0,-1):
            m=hashlib.sha1()
            seed=((seedE+i )%( 2 **160))
            m.update(seed)
            hash2=m.digest()
            tmp=[access_bit(hash2,i) for i in range(len(hash2)*8)]
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
    '''
    Find a square root of n over Fp
    if n is a square in Fp
    
    '''
    
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








def Prime(N):
    '''
    Return True if N is a prime number
    
    >>> Prime(23)
    True
    
    >>> Prime(216545618759)
    True
    
    >>> Prime(2222222)
    False
    '''
    z=str(N)
    #print("Prime in python")
    command_line="sage Prime.sage "+z
    args = shlex.split(command_line)
    now = subprocess.Popen(args, stdout=subprocess.PIPE)
    result = now.communicate()
    a=result[0].decode("utf-8")
    a=a[0:len(a)-1]
    if(a=='True'):
        return True
    else:
        return False






def OrderCourbe(a,b,p):
    '''
    Return the number of points
    on the elliptic curve 
    y²=x³+ax+b  over Fp
    
    >>> OrderCourbe(1,1,23)
    28
    
    '''
    z=str(p)+" "+str(a)+" "+str(b)
    command_line="sage OrdreCourbe.sage "+z 
    args = shlex.split(command_line)
    now = subprocess.Popen(args, stdout=subprocess.PIPE)
    result = now.communicate()
    
    a=int(result[0])
    return a






def listeCouple(l):
    '''
    Convert a list
    to a list of ordered pair
    '''
    
    lc=[]
    for i in range(0,len(l),2):
        lc.append((int(l[i]),int(l[i+1])))
    return lc

def Listeprime(interval,nb_primes):
    '''
    Return a list of nb_primes prime numbers
    where the i-th element of the list
    is the (i*interval)-th prime number
    '''
    
    z=str(interval)+" "+str(nb_primes)
    command_line="sage Listprime.sage "+z
    args = shlex.split(command_line)
    now = subprocess.Popen(args, stdout=subprocess.PIPE)
    result = now.communicate()
    a=result[0].decode("utf-8")
    a=a.split(" ")
    a=a[0:len(a)-1]
    return a



def factorisation(N):
    '''
    Return a list of ordered pairs (p,e) 
    which give the integer factorization of N
    (for all (p,e) p^e is a factor of N
    
    
    >>> factorisation(215525392304898663817)
    [(697499261, 1), (308997305597, 1)]
    
    
    '''
    
    command_line="sage Facorisation.sage "+str(N)
    args = shlex.split(command_line)
    now = subprocess.Popen(args, stdout=subprocess.PIPE)
    result = now.communicate()
    a=result[0].decode("utf-8")
    a=a.split(" ")
    a=a[0:len(a)-1]
    return listeCouple(a)









'''
    Algorithm from:
        Cours de Cryptographie Avancée
        Courbes Elliptiques
        Application à la Cryptographie
        Page 23
        
        by
        Stéphane Ballet et Alexis Bonecaze
        Ecole Polytech de Marseille
    Personal implementation:
'''


def OrderPoint(P,OderCurve,facto):
    '''
    Given  a point P on a elliptic curve, 
    the order of the curve 
    and the integer factorisation of
    curve's order
    Return the order of the point P
    
    
    

    >>> p=9678136687
    >>> Prime(p)
    True
    >>> seed,a,b=random_elliptiqueV4(p)
    >>> OrderCurve=OrderCourbe(a,b,p)
    >>> l=factorisation(OrderCurve)
    >>> sed,yu,P=random_point(a,seed,b,p)
    >>> OrderP=OrderPoint(P,OrderCurve,l)
    >>> isinstance(OrderP*P,Inf)
    True
    >>> OrderCurve%OrderP ==0 
    True
    
    
    
    
    
    >>> sed,yu,P=random_point(a,seed,b,p)
    >>> OrderP=OrderPoint(P,OrderCurve,l)
    >>> isinstance(OrderP*P,Inf)
    True
    >>> OrderCurve%OrderP ==0 
    True
    
    >>> sed,yu,P=random_point(a,seed,b,p)
    >>> OrderP=OrderPoint(P,OrderCurve,l)
    >>> isinstance(OrderP*P,Inf)
    True
    >>> OrderCurve%OrderP ==0 
    True

    

    '''
    
    m=OderCurve
    for i in range(len(facto)):
        p,e=facto[i]
        m=m//(p**e)
        
        Q=m*P
        while isinstance(Q, Inf)==False:
            Q=p*Q
            m=m*p
    return m

def P_order_prime(a,seedE,b,p,Ordercurve,facto):
    """
    Return a point P
    with prime order
    
    """
    x,y,P=random_point(a,seedE,b,p)
    OrderP=OrderPoint(P,Ordercurve,facto)
    pr=Prime(OrderP)
    while(not pr):
        #print("Work Point")
        x,y,P=random_point(a,seedE,b,p)
        OrderP=OrderPoint(P,Ordercurve,facto)
        pr=Prime(OrderP)
    return P,OrderP
def Random_Curve_orderprime(p):
    """
    Return a curve
    with prime order
    
    """
    seed,a,b=random_elliptiqueV4(p)
    Ordre=OrderCourbe(a,b,p)
    pr=Prime(Ordre)
    while(not pr):
        #print("Work Curve")
        seed,a,b=random_elliptiqueV4(p)
        Ordre=OrderCourbe(a,b,p)
        pr=Prime(Ordre)
    return seed,a,b,Ordre

        

        
def rho_probability(P,A,B,Q,p):
    '''
    Fonction aléatoire
    pour la méthode rho
    
    
    '''
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
    
  
def rho_man(P,Q,p,ordreP):
    '''
    Rho algorithme
    Return k that make
    Q=kP
    
    We know that y²=x³+1304x+186
    define an equation of a curve 
    with prime order equal to 2411
    Then (by Lagrange theorem)
    the order of a point P
    different of O have got the same order
    
    >>> p=2437
    >>> a,b= 1304,186
    >>> C=EllipticCurve3(a,b,p)
    >>> P,ordreP=Point(C,568,2048),2411
    >>> k=random.randint(1,ordreP-1)
    >>> Q=k*P
    >>> res=rho_man(P,Q,p,ordreP)
    >>> res==k
    True

    '''
    Ai=random.randint(1,ordreP-1)
    Bi=random.randint(1,ordreP-1)
    Zi=(P*Ai)+(Q*Bi)
    Zibis=Zi
    Aibis=Ai
    Bibis=Bi
    i=0
    while(True):
        Zi,Ai,Bi=rho_probability(P,Ai,Bi,Q,ordreP)
        Zibis,Aibis,Bibis=rho_probability(P,Aibis,Bibis,Q,ordreP)
        Zibis,Aibis,Bibis=rho_probability(P,Aibis,Bibis,Q,ordreP)
        #print('Zi,ZIBIS,Ai,Bi,Aibis,Bibis')
        #print(Zi,Zibis,Ai,Bi,Aibis,Bibis)
        #print(i)
        i=i+1
        if((Bi % ordreP) != (Bibis % ordreP)):
            if(Zi==Zibis) :
                break
    return ((Ai-Aibis) *mod_inverse(-( Bi-Bibis),ordreP)) % ordreP




if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)

"""




'''









def Procedure(P,feinte):
    seed,a,b=random_elliptiqueV4(P)
    courbe=EllipticCurve3(a,b,P)  
    print(courbe)
    Ordercurve=OrderCourbe(a,b,P)
    print("le nombre de point de la courbe vaut")
    print(Ordercurve)
    pr=Prime(Ordercurve)
    while(feinte==1 and not pr ):
        seed,a,b=random_elliptiqueV4(P)
        courbe=EllipticCurve3(a,b,P)  
        print(courbe)
        Ordercurve=OrderCourbe(a,b,P)
        print("le nombre de point de la courbe vaut")
        print(Ordercurve)
        pr=Prime(Ordercurve)
        
        if(pr):
            sed,yu,point1=random_point(a,seed,b,P)
            print(point1)
            k=Ordercurve
            if isinstance(k*point1, Inf):
                print("OrderPoint fonctionne!!!!")
            print("k vaut")
            print(k)
            return a,b,P,point1, Ordercurve,k
    
    sed,yu,point1=random_point(a,seed,b,P)
    print(point1)
    
    l=factorisation(Ordercurve)
    print("l vaut")
    print(l)

    k=OrderPoint(point1,Ordercurve,l)
    if isinstance(k*point1, Inf):
        print("OrderPoint fonctionne!!!!")
    print("k vaut")
    print(k)
    return a,b,P,point1,Ordercurve,k

'''



def Procedure_complexity(p,feinte):
    '''
    Return a random elliptic curve over Fp (with a,b,p)
    a point P with prime order,
    the order of the curve,
    the order of P
    
    if feinte=1
    the order of the curve is a prime
    
    '''
    
    if(feinte==1):
        seed,a,b,Ordercurve=Random_Curve_orderprime(p)
        l=factorisation(Ordercurve)
        P,OrderP=P_order_prime(a,seed,b,p,Ordercurve,l)
        return a,b,p,P,Ordercurve,OrderP
    
    
    
    seed,a,b=random_elliptiqueV4(p)
    Ordercurve=OrderCourbe(a,b,p)
    l=factorisation(Ordercurve)
    P,OrderP=P_order_prime(a,seed,b,p,Ordercurve,l)
    return a,b,p,P,Ordercurve,OrderP
    
'''
def Work(interval,nb_primes,feinte):
    #fichier = open("data", "w")
    listeprimes=Listeprime(interval,nb_primes)
    print(listeprimes)
    print("La liste des premiers est générer")
    for i in range(len(listeprimes)):
        a,b,p,P,Ordercurve,ordreP=Procedure_complexity(int(listeprimes[i]),feinte)
        print("L'odre de la courbe  est :")
        print(Ordercurve)
        print("L'odre de P est :")
        print(ordreP)
        k=random.randint(1,ordreP-1)
        Q=k*P
        print("Q=")
        print(k)
        print("fois P")
        res=rho_man(P,Q,p,ordreP)
        print("rho renvoie :")
        print(res)
        if(k==res):
            print("RHO OK p=")
            print(p)
            #fichier.write("Bonjour monde")
            
            
        else:
            print("error p=")
            print(p)
            return
        

        
        
   
Work(10,100,1)
'''
def Stock(interval,nb_primes):
    '''
    Write in a file for every line:
        -parameters of elliptic curve with prime Order (p,a,b) 
        - a point P on the elliptic curve (Xp,Yp)
        -the order of the curve= order of P
    '''
    fichier = open("data", "w")
    listeprimes=Listeprime(interval,nb_primes)
    print(listeprimes)
    print("La liste des premiers est générer")
    for i in range(len(listeprimes)):
        a,b,p,P,Ordercurve,ordreP=Procedure_complexity(int(listeprimes[i]),1)
        zprint=str(p)+" "+str(a)+" "+str(b)+" "+str(P[0])+" "+str(P[1])+" "+str(Ordercurve)
        print(zprint)
        z=zprint+"\n"
        fichier.write(z)
    fichier.close()
   
        

        
Stock(10,100)
        
'''
a=65
b=6
z=str(a)+" "+str(b)+"\n\n"
Stock(10,100,1)
fichier = open("data", "w")
fichier.write(z)
fichier.write("aa")
fichier.close()
C=EllipticCurve3(1,1,23)
P=Point(C,0,22)
print(P[1])
'''
"""