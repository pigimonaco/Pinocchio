import numpy as np

safe=3
Lx=10
Ly=15
Lz=25
Lxwbl=Lx+2*safe
Lywbl=Ly+2*safe
Lzwbl=Lz+2*safe
Np=Lxwbl*Lywbl*Lzwbl
Len=(Np/32)+(Np%32!=0)

print "Np=%d, Len=%d, Len*32=%d"%(Np,Len,Len*32)

frag_map=np.zeros(Len,dtype=np.uint32)

def set_bit(i,j,k):
    pos=i+(j+k*Lywbl)*Lxwbl
    frag_map[pos/32] |= np.uint32(1<<pos%32)

def get_bit(i,j,k):
    pos=i+(j+k*Lywbl)*Lxwbl
    bit=(frag_map[pos/32] & np.uint32(1<<pos%32))>>pos%32
    return bit

for i in range(safe-1,Lx+safe+1):
    for j in range(safe-1,Ly+safe+1):
        for k in range(safe-1,Lz+safe+1):
            set_bit(i,j,k)


def visualize(k):
    for j in range(Lywbl):
        string=''
        for i in range(Lxwbl):
            string+=str(get_bit(i,j,k))
        print string


