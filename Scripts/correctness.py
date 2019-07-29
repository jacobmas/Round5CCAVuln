# NOTE: The below code needs to be run in the SAGE Jupyter Notebook


####################################################################
##       ESTIMATION OF FAILURE PROBABILITIES IN ROUND5            ##
####################################################################
## This SAGE code estimates the probability of decryption failures
## in Round5 parameters, both for ones that use error correction
## (and thus use polynomial multiplication modulo x^(n+1)-1 in 
## combination with balanced sparse-ternary secrets), and for ones
## that do not use it (thus using polynomial multiplication modulo
## \Phi_{n+1}(x) = x^n + x^(n-1) + ... + 1.)
## 
## In order to obtain accurate estimates (as much as possible), 
## explicit high-precision convolutions are used.
##
## INPUT: Round5 parameters d, n, h, 
##                          log_2(q), log_2(p), log_2(t),
##                          logb, \overline{n}, \overline{m},
##                          f, mu
##        A number of exemplary "high failure" Round5 parameters
##        are provided as part of this code to demonstrate how it works.
##
## OUTPUT: The final output is a plot of how the following two 
##         probabilities behave with increasing Hamming weight (h) 
##         (a) the probability of at least one error occurring,
##         (b) the probability of at least two errors occurring.
##
## The functions below provide further intermediate output that are
## useful in themselves for analyzing the system. These include 
## the overall error distribution for a given parameter set, i.e.,
## the probability of at least 'i' errors occurring, for 1<=i<=mu.
#####################################################################
    


RRR = RealField(150) # 150 bits precision, increase if necessary...

#
# Convolve two given probability distributions A and B, 
# represented as lists, modulo q
#
def myconv(A,B,q):
    C = [0 for _ in range(q) ]
    for i in range(q):
    	if i%(1 << 10)==0:
	   print "myconv,i={0}".format(i)
        for j in range(q):
            # restrict support to Z_q
            k = (i+j)%q
            C[k] = C[k] + A[i] * B[j]
    return C

#
# Iteratively convolve a given probability distribution A,
# that has support in Z_q, with itself, 'p' times.
#
def myconvpower(A,p,q):
    nb = p.nbits()
    #print "myconvpower,A={0},p={1},q={2},nb={3}".format(A,p,q,nb)
    result = [0 for _ in range(q)]
    result[0]=1
    for bit in range(nb-1,-1,-1):
    	print "myconvpower,bit={0}".format(bit)
        if bit < nb - 1:
            result = myconv(result,result,q)
        if (p//2^bit)%2 == 1:
            result = myconv(result,A,q)
    return result

#
# Compute the error distribution for
# for a given set of Round5 parameters. 
# ASSUMPTION: Rounding errors are treated as if sampled from a 
# uniform, bounded distribution. Any dependence between the secret
# and the rounding error is ignored.
#   
def errprobsconv(d, n, h, q, p, t, logb, nbar, mbar, f, mu):
    print "errorprobsconv, n={0},h={1}".format(n,h)
    # Sanity check
    if (not (q/2/p) in ZZ) or (not (q/2/t) in ZZ):
        print 'q/2p or q/2t not in ZZ. Abort.'
        return
    
    
    # Uniform error distribution over -q/2p...q/2p
    # represents error introduced when rounding from Z_q to Z_p
    # 
    Ep = [0 for _ in range(q)]
    for i in range(q/2/p):
        Ep[i]=RRR(p/q)
    Ep[q/2/p]=RRR(p/2/q)
    Ep[q-q/2/p]=RRR(p/2/q)
    for i in range(q-q/2/p+1,q):
        Ep[i]=RRR(p/q)

    # Not-so-Uniform error distribution over -q/2p...q/2p
    # represents error introduced when rounding from Z_q to Z_p
    x=(q/2/p)*(q/2/p)
    Ep2 = [0 for _ in range(q)]
    for i in range(q/2/p):
        Ep2[i]=RRR(i/x)
    Ep2[q/2/p]=RRR((q/2/p)/(2*x))
    Ep2[q-q/2/p]=RRR((q/2/p)/(2*x))
    for i in range(q-q/2/p+1,q):
        Ep2[i]=RRR(i/x)
    
    
    # Uniform error distribution over -q/2t...q/2t
    # represents error introduced when rounding from Z_q to Z_p, then to Z_t
    Et = [0 for _ in range(q)]
    for i in range(q/2/t):
        Et[i]=RRR(t/q)
    Et[q/2/t]=RRR(t/2/q)
    Et[q-q/2/t]=RRR(t/2/q)
    for i in range(q-q/2/t+1,q):
        Et[i]=RRR(t/q)
    

    # Compute the distribution of the error term (j_B*r - s*j_U), 
    # where s, r are the secret-keys, and j_B, j_U are the errors 
    # introduced while computing the public-keys of Alice and Bob
    # respectively; see Section 2.8, 2.8.3 of the specification.
    # ASSUMPTION: Any dependencies between j_B and s, j_U and r ignored.
    conv_fold1 = h
    print "Gadzooks"
    EpH1 = myconvpower(Ep,conv_fold1,q)

    print "Madzooks"
    conv_fold2 = h
    EpH2 = myconvpower(Ep2,conv_fold2,q)

    EpH = myconv(EpH1,EpH2,q)
    
    # Compute the distribution of (j_B*r - s*j_U) + j_v
    EpHt = myconv(EpH,Et,q)
    
    
    # Return the error distribution
    return (EpH,EpHt)

# (bit) failure probabilties for Round5 parameters,
# (a) ring parameters with error correction,
# (b) ring parameters without error correction,
# (c) non-ring parameters.
#
def errprobs(d, n, h, logq, logp, logt, logb, nbar, mbar, f, mu):
    print "In errorprobs, n={0}".format(n)
    # Round5 moduli
    q = 2**logq; p = 2**logp; t = 2**logt;
     # Decryption failure/success thresholds
    # Decryption fails if error terms \in [up_fail, down_fail)
    up_fail = (q + (1<<logb)) // (1<<(logb+1))
    down_fail = (q*( (1<<(logb+1)) - 1) + (1<<logb)) // (1<<(logb+1))

    print "up_fail={0},down_fail={1},logq={2},logp={3},logt={4},logb={5}".format(up_fail,down_fail,logq,logp,logt,logb)
    
    # Get the (rounding) error distribution for given parameters,
    # both the intermediate (before adding j_v) and the final
    # distribution (after adding j_v)
    (EpH,EpHt) = errprobsconv(d, n, h, q, p, t, logb, nbar, mbar, f, mu)
    print "Blop"
    
   
    
    
    # Sanity check
    if logb==1:
        assert up_fail==(q//4)
        assert down_fail==(3*q//4)
    
    
    ##########################################################################
    # Error probability for ring parameters with and without error correction.
    # 
    # For both cases, polynomial multiplication modulo x^(n+1)-1 must be done,
    # hence computing the error distribution in that case, i.e., distribution
    # of (j_B*r - s*j_U) + j_v is always necessary, irrespective of whether 
    # f=0 or f>0.
    #
    # ASSUMING that dependencies between coefficients of the overall error term
    # (j_B*r - s*j_U) + j_v are ignored, the distribution can be approximated 
    # by a binomial distribution with some bit error probability 
    # bfp = Pr[ up_fail <= XX < down_fail ], 
    # where Pr[ XX = x ] = EpHt[x].
    #
    # In more generic terms,
    # let Pcond[a] = Pr[ up_fail <= XX - a mod q < down_fail ],
    # where Pr[ XX = x ] = EpHt[x].
    ##########################################################################
    
    # Init Pcond
    Pcond = [0 for _ in range(q)]
    
    # Compute Pcond[a], differentiating between the two types of ring parameters.
    for a in range(q):
    	print "a={0}".format(a)
        for i in range(q):
	    #print "a={0},i={1}".format(a,i)
            if f==0 and n==d:
                ##########################################################################
                # Error probability for ring parameters without error correction, 
                # i.e., polynomial multiplication modulo x^n+x^(n-1)+...+1).
                # The above multiplication involves one modulo x^(n+1)-1, followed by a 
                # subtraction of the highest order coefficient a = c_n(s,e) from all the
                # other coefficients, see Sec. 2.8.4 of the specification.
                # Thus, compute Pcond[a] = Pr[ up_fail <= XX - a mod q < down_fail ],
                # where Pr[ XX = x ] = EpHt[x]
                ##########################################################################
        
                coeff = (i+q-a)%q
            
            elif f>0 or n==1:
                ##############################################################################
                # Error probability for ring parameters with error correction; i.e., f>0
                # (polynomial multiplication modulo x^(n+1)-1), or for non-ring parameters;
                # i.e., n==1, where a similar independence assumption is made.
                # 
                # Note that here the highest order coefficient a = c_n(s,e) 
                # does not influence the other coefficients.
                #
                # Therefore, as mentioned above, bfp = Pr[ up_fail <= XX < down_fail ], 
                # where Pr[ XX = x ] = EpHt[x].
                ##############################################################################
                coeff = i%q               # same as coeff=i
                assert coeff==i
            # Finally, compute Pcond[a]
            if coeff>=up_fail and coeff<down_fail:
                Pcond[a] = Pcond[a] + EpHt[i]

    # Compute the overall error distribution of (j_B*r - s*j_U) + j_v that can be approximated
    # by a Binomial distribution, with additional influence by the a=c_n(s,e) term in the case
    # f==0, i.e., in parameters without error correction.
    
    # Number of errors that occur during decryption
    ee = [ 0 for i in range(mu+1) ]
    for k in range(mu+1):
        # Note that below, for the case f>0, Pcond[a] always = Pcond[0]
        ee[k]=sum( binomial(mu,k)*(Pcond[a]^k)*((1-Pcond[a])^(mu-k))*EpH[a] for a in range(q) )
    # Overall error distribution
    eetail = [0 for _ in range(mu+1)]
    eetail[-1]=ee[-1]
    for i in range(mu-1,-1,-1):
        eetail[i]=eetail[i+1]+ee[i]
    
    
    # Compute the bit failure probability
    if f==0 and n==d:
        # additional influence due to subtraction of a=c_n(s,e) must be accounted for.
        bfp = eetail[1]/mu
    else:
        # no influence from a=c_n(s,e)
        bfp = sum(EpHt[up_fail:down_fail])
    if n==d and f>0:  assert bfp==Pcond[0]         # Sanity check        
        
    
    # Compute the final error probability after correcting 'f' errors
    ffp = 0.
    if f>0:
        ffp = sum( binomial(mu, j) * (bfp^j) * ((1 - bfp)^(mu - j)) for j in range(f + 1, mu))
    else:
        ffp = None
        
    # Return the error distribution, the bit failure probability, 
    # and the post-error-correction failure probability
    return (eetail,bfp,ffp)


#
# For a given set of Round5 parameters, compute the probability of 
# (a) at least 1 error, (b) at least two errors,
# (c) overall error distribution.
#
def getfail(param_name, d, n, h, logq, logp, logt, logb, nbar, mbar, f, mu):
    
    (etail, bfp, ffp) = errprobs(d, n, h, logq, logp, logt, logb, nbar, mbar, f, mu)
    
    # Probability of at least one error
    log_1fp = log(etail[1],2.0)
    
    # Probability of at least 2 errors
    log_2fp = log(etail[2],2.0)
    
    return (log_1fp,log_2fp,etail,ffp)


def showfail(param_name, d, n, h, logq, logp, logt, logb, nbar, mbar, f, mu):
    (log_1fp,log_2fp,etail,ffp) = getfail(param_name, d, n, h, logq, logp, logt, logb, nbar, mbar, f, mu)
    print "\n\nParameter",param_name,"has failure probability (before error correction):\t", log_1fp
    if ffp:               print "Post-XEf failure probability:\t",log(ffp,2.0),"\n"
    for i in range(6):    print "Parameter",param_name,"-- Probability of at least",i,"errors:\t", log(etail[i],2.0)
    return


#############################################################################################################


print "\n\nWarning: Please make sure this script is being run in the SAGE Jupyter notebook.\nIf not, then please exit (Ctrl+C) and try again.\n\n"


# Print failure rates of Round5 parameters
# Uncomment parameters to run them.

# Be warned that parameters with q>2^15 take *a long time* to compute

# Round5.KEM

# Round5 ring parameters without error correction
#showfail("R5ND_1KEM_0c",618,618,104,11,8,4,1,1,1,0,128)
#showfail("R5ND_3KEM_0c",786,786,384,13,9,4,1,1,1,0,192)
#showfail("R5ND_5KEM_0c",1018,1018,428,14,9,4,1,1,1,0,256)

# Round5 ring parameters with error correction
#showfail("R5ND_1KEM_5c",490,490,162,10,7,3,1,1,1,5,318)
#showfail("R5ND_3KEM_5c",756,756,242,12,8,2,1,1,1,5,410)
#showfail("R5ND_5KEM_5c",940,940,414,12,8,2,1,1,1,5,490)

# Round5 non ring parameters, without error correction
#showfail("R5N1_1KEM_0c",594,1,238,13,10,7,3,7,7,0,43)
#showfail("R5N1_3KEM_0c",881,1,238,13,10,7,3,8,8,0,64)
#showfail("R5N1_5KEM_0c",1186,1,712,15,12,7,4,8,8,0,64)

# Round5 PKE

# Round5 ring parameters without error correction
#showfail("R5ND_1PKE_0c",586,586,182,13,9,4,1,1,1,0,128)
#showfail("R5ND_3PKE_0c",852,852,212,12,9,5,1,1,1,0,192)
showfail("R5ND_5PKE_0c",1170,1170,222,13,9,5,1,1,1,0,256)

# Round5 ring parameters with error correction
#showfail("R5ND_1PKE_5c",508,508,136,10,7,4,1,1,1,5,318)
#showfail("R5ND_3PKE_5c",756,756,242,12,8,3,1,1,1,5,410)
#showfail("R5ND_5PKE_5c",940,940,414,12,8,3,1,1,1,5,490)

# Round5 non ring parameters, without error correction
#showfail("R5N1_1PKE_0c",636,1,114,12,9,6,2,8,8,0,64)
#showfail("R5N1_3PKE_0c",876,1,446,15,11,7,3,8,8,0,64)
#showfail("R5N1_5PKE_0c",1217,1,462,15,12,9,4,8,8,0,64)


