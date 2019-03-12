import numpy as np

N=20
a=np.arange(N)
ind=np.copy(a)
np.random.shuffle(ind)
# a=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],dtype=np.int)
# ind=np.array([2, 5, 0, 4, 3, 9, 7, 6, 1, 8],dtype=np.int)
frag_pos=ind*10.0

print "inizio:"
print "ordine: ",a
print "indici: ",ind
print "dati:   ",frag_pos

counter=0
for i in range(N):
    if ind[i]!=i:
        oldf=frag_pos[i]
        next=ind[i]
        print "Nella casellina: [ind=%d,d=%d]"%(next,oldf)
        while (next != i):
            print "scambio %d e casellina"%next
            tmp=frag_pos[next]
            frag_pos[next]=oldf
            oldf=tmp
            tmp=ind[next]
            ind[next]=next
            next=tmp
            #print "ordine: ",a
            #print "indici: ",ind
            #print "dati:   ",frag_pos
            #print "Nella casellina: [ind=%d,d=%d]"%(next,oldf)
            counter+=1
            if counter>N:
                print "ALLUPPATO!"
                break
        frag_pos[i]=oldf
        ind[i]=next
    else:
        print "%d a posto"%i
    print "Adesso %d e` a posto"%i
print "ordine: ",a
print "indici: ",ind
print "dati:   ",frag_pos


print "finito"
