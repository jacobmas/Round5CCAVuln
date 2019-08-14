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
        if i%(1<<11)==0:
            print("i={0}".format(i))
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
    print("p={0},typeof(p)={1}".format(p,type(p)))
    nb = p.nbits()
    print("Fuck")
    result = [0 for _ in range(q)]
    result[0]=1
    for bit in range(nb-1,-1,-1):
        print("bit={0}".format(bit))
        if bit < nb - 1:
            result = myconv(result,result,q)
        if (p//2^bit)%2 == 1:
            result = myconv(result,A,q)
    return result

#
# Compute the error distribution for
# for a given set of Round5 parameters. 
# ASSUMPTION: Special depending on two-boys and one-boys and biased-boys
#
# min_abs denotes the minimum absolute value to be considered "big"
# min_big denotes the minimum number of coefficients out of all n that need to be at least absolute value min_abs for a
# query to be accepted and sent on
#   
def errprobsconv(d, n, h, twos,ones,min_abs,min_big,q, p, t, logb, nbar, mbar, f, mu):
    (exptwos,expones)=getexpectedtwosones(n,h)
    print("twos={0},ones={1},extwos={2},exones={3},min_abs={4},min_big={5}".format(twos,ones,exptwos,expones,min_abs,min_big))
    # Sanity check
    if (not (q/2/p) in ZZ) or (not (q/2/t) in ZZ):
        print 'q/2p or q/2t not in ZZ. Abort.'
        return
    
    
    # Uniform one-boys error distribution over -q/2p...q/2p
    # represents error introduced when rounding from Z_q to Z_p
    Ep = [0 for _ in range(q)]
    for i in range(q/2/p):
        Ep[i]=RRR(p/q)
    Ep[q/2/p]=RRR(p/2/q)
    Ep[q-q/2/p]=RRR(p/2/q)
    for i in range(q-q/2/p+1,q):
        Ep[i]=RRR(p/q)

     # Uniform two-boys error distribution over 2*(-q/2p)...2*(q/2p)
    # represents error introduced when rounding from Z_q to Z_p
    Ep2 = [0 for _ in range(q)]
    for i in range(0,q/p,2):
        Ep2[i]=RRR(p/q)
    Ep2[q/p]=RRR(p/2/q)
    Ep2[q-q/p]=RRR(p/2/q)
    for i in range(q-q/p+2,q,2):
        Ep2[i]=RRR(p/q)

    # compute number of elements in biased big range, counting +-q/2/p each as half
    num_biased_elems=((q/p)/2-min_abs)*2+1

    # Biased one-boys error distribution over -q/2p...q/2p
    # represents error introduced when rounding from Z_q to Z_p
    EpB = [0 for _ in range(q)]
    for i in range(min_abs,q/2/p):
        EpB[i]=RRR(1/num_biased_elems)
    EpB[q/2/p]=RRR(1/(2*num_biased_elems))
    EpB[q-q/2/p]=RRR(1/(2*num_biased_elems))
    for i in range(q-q/2/p+1,q-min_abs+1):
        EpB[i]=RRR(1/num_biased_elems)

    check_EpB=reduce(lambda x,y:x+y,EpB)

    # Biased two-boys error distribution over 2*(-q/2p)...2*(q/2p)
    # represents error introduced when rounding from Z_q to Z_p
    EpB2 = [0 for _ in range(q)]
    for i in range(2*min_abs,q/p,2):
        EpB2[i]=RRR(1/num_biased_elems)
    EpB2[q/p]=RRR(1/(2*num_biased_elems))
    EpB2[q-q/p]=RRR(1/(2*num_biased_elems))
    for i in range(q-q/p+2,q-2*min_abs+2,2):
        EpB2[i]=RRR(1/num_biased_elems)

    check_EpB2=reduce(lambda x,y:x+y,EpB2)

        
    print("check_EpB={0},check_EpB2={1}".format(check_EpB,check_EpB2))
    
    
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
    conv_fold = 2*h

    # EpH1 and EpH2 are biased distribution, expected number of one boys and two boys
    print("cyka+blyat=doing ones now")
    EpH1 = myconvpower(EpB,expones,q)
    print("cyka+blyat=doing twos now")

    EpH2 = myconvpower(EpB2,exptwos,q)
    
    EpH = myconv(EpH1,EpH2,q)

    print("cyka blyat=doing EpHB1")

    # EpHB1 and EpHB2 are uniform distribution, biased number of one boys and two boys 
    EpHB1 = myconvpower(Ep,expones,q)

    print("cyka blyat=doing EpHB2")

    EpHB2 = myconvpower(Ep,exptwos,q)

    print("cyka blyat doing EpHB1,EpHB2")
    EpHB = myconv(EpHB1,EpHB2,q)

    print("cyka blyat doing EpHStar from EpH,EpHB")
    EpHStar=myconv(EpH,EpHB,q)
    
    print("cyka blyat doing EpHt")

    # Compute the distribution of (j_B*r - s*j_U) + j_v
    EpHt = myconv(EpHStar,Et,q)
    
    
    # Return the error distribution
    return (EpHStar,EpHt)


#
# Compute (bit) failure probabilities for Round5 parameters,
# (a) ring parameters with error correction,
# (b) ring parameters without error correction,
# (c) non-ring parameters.
#
def errprobs(d, n, h, twos,ones,min_abs,min_big, logq, logp, logt, logb, nbar, mbar, f, mu):
    
    # Round5 moduli
    q = 2**logq; p = 2**logp; t = 2**logt;
    
    

    # Get the (rounding) error distribution for given parameters,
    # both the intermediate (before adding j_v) and the final
    # distribution (after adding j_v)
    (EpH,EpHt) = errprobsconv(d, n, h,twos,ones,min_abs,min_big, q, p, t, logb, nbar, mbar, f, mu)
        
    
    # Decryption failure/success thresholds
    # Decryption fails if error terms \in [up_fail, down_fail)
    up_fail = (q + (1<<logb)) // (1<<(logb+1))
    down_fail = (q*( (1<<(logb+1)) - 1) + (1<<logb)) // (1<<(logb+1))
    
    
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
        for i in range(q):
            if true:
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


# Compute expected number of twos and ones for a given n,h, assuming n-1 possibilities (=p-2 possibilities) per the
# Lemma in the paper version
# This is going to be ever so slightly optimistic for # of twos by like 1 due to the longer tail but not much
def getexpectedtwosones(n,h):
    twos=h^2*(n-1)/((2*n)*(n-1))*1.
    ones=2*h*(n-h)/(n)*1.
    return (Integer(round(twos)),Integer(round(ones)))
    

#
# For a given set of Round5 parameters, compute the probability of 
# (a) at least 1 error, (b) at least two errors,
# (c) overall error distribution.
#
def getfail(param_name, d, n, h,twos,ones,min_abs,min_big,logq, logp, logt, logb, nbar, mbar, f, mu):

    
    (etail, bfp, ffp) = errprobs(d, n, h, twos,ones,min_abs,min_big, logq, logp, logt, logb, nbar, mbar, f, mu)
    
    # Probability of at least one error
    log_1fp = log(etail[1],2.0)
    
    # Probability of at least 2 errors
    log_2fp = log(etail[2],2.0)
    
    return (log_1fp,log_2fp,etail,ffp)


def showfail(param_name, d, n, h, twos, ones,min_abs,min_big,logq, logp, logt, logb, nbar, mbar, f, mu):
    (log_1fp,log_2fp,etail,ffp) = getfail(param_name, d, n, h, twos,ones,min_abs,min_big, logq, logp, logt, logb, nbar, mbar, f, mu)
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
#          param_name    d     n   h  twos ones min_abs=4 min_big=1152
showfail("R5ND_5PKE_0c",1170,1170,222,88,268,4,1152,13,9,5,1,1,1,0,256)


