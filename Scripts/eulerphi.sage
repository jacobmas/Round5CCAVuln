
def compute_euler_list(max_val):
    my_lst=[[] for y in range(max_val)]
    for y in range(2,max_val):
        z=euler_phi(y)
        my_lst[z].append(y)
    for y in range(2,max_val/20,2):
        if len(my_lst[y])==0:
            print("{0},{1},{2}".format(y,len(my_lst[y]),my_lst[y]))
