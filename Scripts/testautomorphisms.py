def testautomorphisms(n,h):
    s=create_secret_vector(n,h)
    return 0

def create_secret_vector(n,h):
    vector=[0 for i in range(n)]
    for i in range(h):
        idx=randrange(0,n)
        while vector[idx]!=0:
            idx=randrange(0,n)
        vector[idx] = -1 if (i & 1) else 1



