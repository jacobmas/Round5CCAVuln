def testautomorphisms(n,h,q,p,t):
    (a,s,b)=r5_pke_gen(n,h,q,p,t)
    K.<x>=PolynomialRing(ZZ)
    a_K=K(list(a))
    s_K=K(list(s))
    b_K=K(list(b))
    b1=b_K(1)
    s1=s_K(1)
    a1=a_K(1)
    b1mod=(b1 % (q))%(n+1)
    print("b_K(1)={0},{1},{2}".format(b_K(1),b1mod,b1%(n+1)))
    print("a_K(1)={0},s_K(1)={1}".format(a1,s1))
    return (a,s,b)

def round_x(x,q):
    while x>=q/2:
        x=x-q
    while x<-q/2:
        x=x+q
    return x
    

def r5_pke_gen(n,h,q,p,t):
    s=create_secret_vector(n,h)
    a=create_a_vector(n,h,q)
    R=CyclotomicField(n+1)
    sprime=R(s)
    aprime=R(a)
    t=map(lambda x:x%q,list(sprime*aprime))
    # round
    u=map(lambda x:(q/p)*Integer(floor(x*p*1./(q*1.)+1./2.)),t)
    v=[x - y for x,y in zip(t,u)]
    bprime=R(u)
    e=list(bprime-aprime*sprime)
    f=reduce(lambda x,y:x+y,e)
    print("f={0}".format(f%(n+1)))
    return (aprime,sprime,bprime)

def create_a_vector(n,h,q):
    return [randrange(0,q) for i in range(n)]
    
    
def create_secret_vector(n,h):
    vector=[0 for i in range(n)]
    for i in range(h):
        idx=randrange(0,n)
        while vector[idx]!=0:
            idx=randrange(0,n)
        vector[idx] = -1 if (i & 1) else 1
    return vector

s=testautomorphisms(1170,222,(1<<13),(1<<9),(1<<5))
#print("s={0}".format(s))
