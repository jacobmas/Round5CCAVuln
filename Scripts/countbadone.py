# count all the absolute value 2 using dynamic programming
def countallabs1(n,h):
    full_d={}
    d_2={}
    d_2[(0,2,0,0,-1,-1)]=1
    d_2[(0,1,0,1,-1,0)]=1
    d_2[(1,1,1,0,-1,1)]=1
    d_2[(0,1,0,1,0,-1)]=1
    d_2[(0,0,0,0,0,0)]=1
    d_2[(1,0,0,1,0,1)]=1
    d_2[(1,1,1,0,1,-1)]=1
    d_2[(1,0,0,1,1,0)]=1
    d_2[(2,0,0,0,1,1)]=1
    full_d[2]=d_2
    for s in range(3,n+1):
        compute_dp1(s,full_d,n,h)
        if s-1 in full_d:
            del full_d[s-1]
    summarize_result1(full_d,n,h)

def get_dp_value1(full_d,s,pos,neg,t,ones,x,y):
    tp=(pos,neg,t,ones,x,y)
    if tp not in full_d[s]:
        return 0
    return full_d[s][tp]

def compute_dp1(s,full_d,n,h):
    full_d[s]={}
    print("s={0}".format(s))
    for pos in range(0,min(s+1,h/2+1)):
        for neg in range(0,min(s+1-pos,h/2+1)):
            for t in range(0,2*min(pos,neg)+1):
                for ones in range(0,2*(pos+neg)+1):
                    do_stuff_dp1(s,full_d,n,h,pos,neg,t,ones)
    #print("full_d[s]={0}".format(full_d[s]))

def do_stuff_dp1(s,full_d,n,h,pos,neg,t,ones):
    #print("do_stuff_dp1,s={0},pos={1},neg={2},t={3},ones={4}".format(s,pos,neg,t,ones))
    if neg>=2:
        tp=(pos,neg,t,ones,-1,-1)
        full_d[s][tp]=0
        for i in range(-1,2):
            full_d[s][tp]+=get_dp_value1(full_d,s-1,pos,neg-1,t,ones,i,-1)
    if neg>=1:            
        tp=(pos,neg,t,ones,-1,0)
        full_d[s][tp]=0
        for i in range(-1,2):
            full_d[s][tp]+=get_dp_value1(full_d,s-1,pos,neg,t,ones-1,i,-1)
    if neg>=1 and pos>=1 and t>=1:            
        tp=(pos,neg,t,ones,-1,1)
        full_d[s][tp]=0
        for i in range(-1,2):
            full_d[s][tp]+=get_dp_value1(full_d,s-1,pos-1,neg,t-1,ones,i,-1)
    if neg>=1:       
        tp=(pos,neg,t,ones,0,-1)
        full_d[s][tp]=0
        for i in range(-1,2):
            full_d[s][tp]+=get_dp_value1(full_d,s-1,pos,neg-1,t,ones-1,i,0)
    tp=(pos,neg,t,ones,0,0)
    full_d[s][tp]=0
    for i in range(-1,2):
        full_d[s][tp]+=get_dp_value1(full_d,s-1,pos,neg,t,ones,i,0)
    if pos>=1:       
        tp=(pos,neg,t,ones,0,1)
        full_d[s][tp]=0
        for i in range(-1,2):
            full_d[s][tp]+=get_dp_value1(full_d,s-1,pos-1,neg,t,ones-1,i,0)
    if neg>=1 and pos>=1 and t>=1:
        tp=(pos,neg,t,ones,1,-1)
        full_d[s][tp]=0
        for i in range(-1,2):
            full_d[s][tp]+=get_dp_value1(full_d,s-1,pos,neg-1,t-1,ones,i,1)
    if pos>=1:
        tp=(pos,neg,t,ones,1,0)
        full_d[s][tp]=0
        for i in range(-1,2):
            full_d[s][tp]+=get_dp_value1(full_d,s-1,pos,neg,t,ones-1,i,1)
    if pos>=2:
        tp=(pos,neg,t,ones,1,1)
        full_d[s][tp]=0
        for i in range(-1,2):
            full_d[s][tp]+=get_dp_value1(full_d,s-1,pos-1,neg,t,ones,i,1)
def summarize_result1(full_d,n,h):
    result_dict={}
    log_result_dict={}
    tot_settings=binomial(n,h)*binomial(h,h/2)
    log_tot_settings=log(tot_settings*1.,2.)*1.
    print("Printing counts for t, warning will be LARGE!")
    print("t   ,log_2(Pr[# 2s=t]),count of #2s=t")
    for t in range(0,h):
        for ones in range(0,h+1):
            result_dict[(t,ones)]=0
                    
            for x in range(-1,2):
                for y in range(-1,2):
                    tp=(h/2,h/2,t,ones,x,y)
                    if tp in full_d[n]:
                        result_dict[(t,ones)]=result_dict[(t,ones)]+full_d[n][tp]
            print("{0},{1},{2}".format(t,ones,result_dict[(t,ones)]))
    # Sum up
