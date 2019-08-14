RRR = RealField(150)
def countbadt(h):
    maxpos=h/2
    maxneg=h/2
    maxt=maxpos+maxneg-1
    y=0
    for p in range(0,maxpos+1):
        for n in range(0,maxneg+1):
            for t in range(0,maxpos+maxneg):
                if t<=2*min(p,n):
                    y+=1
    return y


# count all the absolute value 2 using dynamic programming
def countallabs2(n,h):
    full_d={}
    d_2={}
    d_2[(0,2,0,-1,-1)]=1
    d_2[(0,1,0,-1,0)]=1
    d_2[(1,1,1,-1,1)]=1
    d_2[(0,1,0,0,-1)]=1
    d_2[(0,0,0,0,0)]=1
    d_2[(1,0,0,0,1)]=1
    d_2[(1,1,1,1,-1)]=1
    d_2[(1,0,0,1,0)]=1
    d_2[(2,0,0,1,1)]=1
    full_d[2]=d_2
    for s in range(3,n+1):
        compute_dp(s,full_d,n,h)
        if s-1 in full_d:
            del full_d[s-1]
    summarize_result(full_d,n,h)

def get_dp_value(full_d,s,pos,neg,t,x,y):
    tp=(pos,neg,t,x,y)
    if tp not in full_d[s]:
        return 0
    return full_d[s][tp]

def compute_dp(s,full_d,n,h):
    full_d[s]={}
    print("s={0}".format(s))
    for pos in range(0,min(s+1,h/2+1)):
        for neg in range(0,min(s+1-pos,h/2+1)):
            for t in range(0,2*min(pos,neg)+1):
                if neg>=2:
                    tp=(pos,neg,t,-1,-1)
                    full_d[s][tp]=0
                    for i in range(-1,2):
                       full_d[s][tp]+=get_dp_value(full_d,s-1,pos,neg-1,t,i,-1)
                if neg>=1:            
                    tp=(pos,neg,t,-1,0)
                    full_d[s][tp]=0
                    for i in range(-1,2):
                       full_d[s][tp]+=get_dp_value(full_d,s-1,pos,neg,t,i,-1)
                if neg>=1 and pos>=1 and t>=1:            
                    tp=(pos,neg,t,-1,1)
                    full_d[s][tp]=0
                    for i in range(-1,2):
                       full_d[s][tp]+=get_dp_value(full_d,s-1,pos-1,neg,t-1,i,-1)
                if neg>=1:       
                    tp=(pos,neg,t,0,-1)
                    full_d[s][tp]=0
                    for i in range(-1,2):
                       full_d[s][tp]+=get_dp_value(full_d,s-1,pos,neg-1,t,i,0)
                tp=(pos,neg,t,0,0)
                full_d[s][tp]=0
                for i in range(-1,2):
                    full_d[s][tp]+=get_dp_value(full_d,s-1,pos,neg,t,i,0)
                if pos>=1:       
                    tp=(pos,neg,t,0,1)
                    full_d[s][tp]=0
                    for i in range(-1,2):
                       full_d[s][tp]+=get_dp_value(full_d,s-1,pos-1,neg,t,i,0)
                if neg>=1 and pos>=1 and t>=1:
                    tp=(pos,neg,t,1,-1)
                    full_d[s][tp]=0
                    for i in range(-1,2):
                       full_d[s][tp]+=get_dp_value(full_d,s-1,pos,neg-1,t-1,i,1)
                if pos>=1:
                    tp=(pos,neg,t,1,0)
                    full_d[s][tp]=0
                    for i in range(-1,2):
                       full_d[s][tp]+=get_dp_value(full_d,s-1,pos,neg,t,i,1)
                if pos>=2:
                    tp=(pos,neg,t,1,1)
                    full_d[s][tp]=0
                    for i in range(-1,2):
                       full_d[s][tp]+=get_dp_value(full_d,s-1,pos-1,neg,t,i,1)
def summarize_result(full_d,n,h):
    result_dict={}
    log_result_dict={}
    tot_settings=binomial(n,h)*binomial(h,h/2)
    log_tot_settings=log(tot_settings*1.,2.)*1.
    print("Printing counts for t, warning will be LARGE!")
    print("t   ,log_2(Pr[# 2s=t]),count of #2s=t")
    for t in range(0,h):
        result_dict[t]=0
        for x in range(-1,2):
            for y in range(-1,2):
                tp=(h/2,h/2,t,x,y)
                if tp in full_d[n]:
                    result_dict[t]=result_dict[t]+full_d[n][tp]
        print("{0:4},{1}".format(t,result_dict[t]))
    # Sum up


def new_summarize_results(filename,n,h):
    import csv
    tot_settings=binomial(n,h)*binomial(h,h/2)
    extra_tot=500
    tot_settings=tot_settings>>500
    log_tot_settings=log(tot_settings*1.,2.)*1.+500.
    with open(filename, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            extra=0
            x=Integer(row[1])
            if x>(1<<750):
                x=x>>400
                extra+=400
            y=log(x*1.,2.)+extra
            print("{0},{1}".format(row[0],y-log_tot_settings))

# compute the probabilities of getting some extreme parameters for rounding

def compute_rounding_extreme_probs(n,h,min_abs,min_big):
    ret=0
    if min_abs<1:
        print("BAD min_abs, must be >0")
        return -1
    p=RRR(((8.-min_abs)*2+1.)/16.0)
    print("p={0}".format(p))
        
    for i in range(min_big,n+1):
        curr=binomial(n,i)*power(p,i)*power((1-p),(n-i))
 #       print("i={0},curr={1}".format(i,curr,pow(p,i),pow(1-p,n-i)))
        ret=ret+curr
    return ret
    

