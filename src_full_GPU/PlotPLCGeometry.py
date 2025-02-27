import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import mpl_toolkits.mplot3d.art3d as art3d

def plot_cube(C,L,fate,ax):


    if (fate>0):
        color={ 0: "gray",
                1: "blue",
                2: "cyan",
                3: "palegreen",
                4: "turquoise" }
        alpha={ 0: 0.01,
                1: 0.05,
                2: 0.05,
                3: 0.05,
                4: 0.05 }

        r = [0.0,1.0]
        vertices=np.array(list(product(r,r,r)))
        for i in range(3):
            vertices[:,i]*=L[i]
            vertices[:,i]+=C[i]
        X,Y=np.meshgrid(r,r)
        ax.scatter3D(vertices[:, 0], vertices[:, 1], vertices[:, 2])
        ax.plot_surface(X*L[0]+C[0], Y*L[1]+C[1],        C[2], alpha=alpha[fate], color=color[fate])
        ax.plot_surface(X*L[0]+C[0], Y*L[1]+C[1],   L[2]+C[2], alpha=alpha[fate], color=color[fate])
        ax.plot_surface(X*L[0]+C[0],        C[1], Y*L[2]+C[2], alpha=alpha[fate], color=color[fate])
        ax.plot_surface(X*L[0]+C[0],   L[1]+C[1], Y*L[2]+C[2], alpha=alpha[fate], color=color[fate])
        ax.plot_surface(       C[0], X*L[1]+C[1], Y*L[2]+C[2], alpha=alpha[fate], color=color[fate])
        ax.plot_surface(  L[0]+C[0], X*L[1]+C[1], Y*L[2]+C[2], alpha=alpha[fate], color=color[fate])

def plot_dist(C,L,V,r,ax):

    vers=(C+L/2-V)/np.linalg.norm(C+L/2-V)
    ax.plot(np.array([V[0],V[0]+vers[0]*r]),
            np.array([V[1],V[1]+vers[1]*r]),
            np.array([V[2],V[2]+vers[2]*r]))


def plot_cone(V,D1,A,ax,rstart=0.0,rstop=1.0,Nr=5,Nphi=10):

    deltar=(rstop-rstart)/np.float(Nr)

    D = D1/np.linalg.norm(D1)

    ivers=np.array([1.,0.,0.])
    p = np.cross(D,ivers)
    if np.linalg.norm(p)==0:
        ivers=np.array([0.,0.,1.])
        p = np.cross(D,ivers)
    p /= np.linalg.norm(p)
    q = np.cross(D,p)

    AA=min(A,np.pi)
    if A<180.:
        for iphi in range(Nphi):

            phi1 = iphi * 2.0*np.pi / np.float(Nphi)
            phi2 = (iphi+1) * 2.0*np.pi / np.float(Nphi)

            for ir in range(Nr):

                r1=ir*deltar+rstart
                r2=(ir+1)*deltar+rstart

                rad1=r1*np.sin(AA)
                rad2=r2*np.sin(AA)

                X=np.array([[V[0]+r1*np.cos(AA)*D[0]+rad1*(p[0]*np.cos(phi1)+q[0]*np.sin(phi1)),
                             V[0]+r2*np.cos(AA)*D[0]+rad2*(p[0]*np.cos(phi1)+q[0]*np.sin(phi1))],
                            [V[0]+r1*np.cos(AA)*D[0]+rad1*(p[0]*np.cos(phi2)+q[0]*np.sin(phi2)),
                             V[0]+r2*np.cos(AA)*D[0]+rad2*(p[0]*np.cos(phi2)+q[0]*np.sin(phi2))]])
                Y=np.array([[V[1]+r1*np.cos(AA)*D[1]+rad1*(p[1]*np.cos(phi1)+q[1]*np.sin(phi1)),
                             V[1]+r2*np.cos(AA)*D[1]+rad2*(p[1]*np.cos(phi1)+q[1]*np.sin(phi1))],
                            [V[1]+r1*np.cos(AA)*D[1]+rad1*(p[1]*np.cos(phi2)+q[1]*np.sin(phi2)),
                             V[1]+r2*np.cos(AA)*D[1]+rad2*(p[1]*np.cos(phi2)+q[1]*np.sin(phi2))]])
                Z=np.array([[V[2]+r1*np.cos(AA)*D[2]+rad1*(p[2]*np.cos(phi1)+q[2]*np.sin(phi1)),
                             V[2]+r2*np.cos(AA)*D[2]+rad2*(p[2]*np.cos(phi1)+q[2]*np.sin(phi1))],
                            [V[2]+r1*np.cos(AA)*D[2]+rad1*(p[2]*np.cos(phi2)+q[2]*np.sin(phi2)),
                             V[2]+r2*np.cos(AA)*D[2]+rad2*(p[2]*np.cos(phi2)+q[2]*np.sin(phi2))]])

                ax.plot_surface(X,Y,Z, color='r', alpha=0.2)

    ax.plot([V[0],V[0]+rstop*D[0]],[V[1],V[1]+rstop*D[1]],[V[2],V[2]+rstop*D[2]],lw=3,color='k')

    Ncap=max(int(1.5*Nphi*AA/np.pi)+1,5)
    a=np.complex(0,Nphi+1)
    b=np.complex(0,Ncap)
    u, v = np.mgrid[0:2*np.pi:a, 0:AA:b]

    x=V[0]+rstop*(np.cos(v)*D[0]+np.sin(v)*(p[0]*np.cos(u)+q[0]*np.sin(u)))
    y=V[1]+rstop*(np.cos(v)*D[1]+np.sin(v)*(p[1]*np.cos(u)+q[1]*np.sin(u)))
    z=V[2]+rstop*(np.cos(v)*D[2]+np.sin(v)*(p[2]*np.cos(u)+q[2]*np.sin(u)))
    ax.plot_wireframe(x, y, z, color="r")
    x=V[0]+rstart*(np.cos(v)*D[0]+np.sin(v)*(p[0]*np.cos(u)+q[0]*np.sin(u)))
    y=V[1]+rstart*(np.cos(v)*D[1]+np.sin(v)*(p[1]*np.cos(u)+q[1]*np.sin(u)))
    z=V[2]+rstart*(np.cos(v)*D[2]+np.sin(v)*(p[2]*np.cos(u)+q[2]*np.sin(u)))
    ax.plot_wireframe(x, y, z, color="b")



file=open("pinocchio.cdm.geometry.out","r")
line=file.readline()
Ncubes=np.int(line.split()[3])
line=file.readline()
rstart=np.float(line.split()[3])
rstop=np.float(line.split()[4])
line=file.readline()
V=np.array(line.split()[3:]).astype(np.float)
line=file.readline()
D=np.array(line.split()[3:]).astype(np.float)
line=file.readline()
L=np.array(line.split()[3:]).astype(np.float)
line=file.readline()
A=np.float(line.split()[3])
line=file.readline()
IPD=np.float(line.split()[3])

line=file.readline()

replication=[]
rmin=[]
rmax=[]
fate=[]
for i in range(Ncubes):
    line=file.readline()
    replication.append(np.array(line.split()[1:4]).astype(np.int))
    rmin.append(np.float(line.split()[4]))
    rmax.append(np.float(line.split()[5]))
    fate.append(np.int(line.split()[6]))

replication=np.asarray(replication)
rmin=np.asarray(rmin)
rmax=np.asarray(rmax)
fate=np.asarray(fate)

A *= np.pi/180.
V *= IPD
L *= IPD
rstart *= IPD
rstop *= IPD
rmin *= IPD
rmax *= IPD

file.close()

rlargest=np.max(rmax)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")
ax.set_xlim([-rlargest,rlargest])
ax.set_ylim([-rlargest,rlargest])
ax.set_zlim([-rlargest,rlargest])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# draw cubes
for i in range(Ncubes):
    plot_cube(replication[i]*L,L,fate[i],ax)
    # plot_dist(replication[i]*L,L,V,rmin[i],ax)
    # plot_dist(replication[i]*L,L,V,rmax[i],ax)

# draw cone
plot_cone(V,D,A,ax,rstart=rstart,rstop=rstop,Nr=10,Nphi=50)

plt.show()

