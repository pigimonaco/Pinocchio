# *****************************************************************
# *                        PINOCCHIO  V4.1                        *
# *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
# *****************************************************************

# This code was written by
# Pierluigi Monaco
# Dipartimento di Fisica, Universita` di Trieste
# Copyright (C) 2016

# web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

# The original code was developed with:
# Tom Theuns           Institute for Computational Cosmology, University of Durham 
# Giuliano Taffoni     INAF - Osservatorio Astronomico di Trieste

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



# This python script gives three plots of the example run provided in the distribution 
# and of the same run performed as a test


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


(m,nm,fit)=np.loadtxt('pinocchio.0.0000.example.mf.out',unpack=True,usecols=(0,1,5))

plt.figure()

plt.plot(m,m*nm,label='pinocchio MF',ls='-',lw=3,c='green')
plt.plot(m,m*fit,label='Watson fit',c='blue')


plt.xlim([0.8e13,1.e16])
plt.ylim([1.e-8,1.e-3])
plt.xscale('log')
plt.yscale('log')

plt.legend(frameon=True)
plt.title('Mass function at z=0')
plt.xlabel(r'M (M$_\odot$)',fontsize=16)
plt.ylabel(r'M n(M) (Mpc$^{-3}$)',fontsize=16)

plt.savefig('mf.png')
plt.show()

(x,y,z)=np.loadtxt('pinocchio.0.0000.example.catalog.out',unpack=True,usecols=(5,6,7))

plt.figure()

index=(z<100)
plt.plot(x[index],y[index],'o',c='green',label='pinocchio halos')

plt.xlim([0,500])
plt.ylim([0,500])
plt.xscale('linear')
plt.yscale('linear')

plt.legend(frameon=True)
plt.title('Large-scale structure at z=0')
plt.xlabel(r'x (Mpc/h)',fontsize=16)
plt.ylabel(r'y (Mpc/h)',fontsize=16)

plt.savefig('lss.png')
plt.show()

(x,y,z,m)=np.loadtxt('pinocchio.example.plc.out',unpack=True,usecols=(2,3,4,8))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

index=(m>1.e14)
ax.scatter(x[index],y[index],z[index],marker='o',c='green')

ax.set_xlabel('x (Mpc/h)')
ax.set_ylabel('y (Mpc/h)')
ax.set_zlabel('z (Mpc/h)')

plt.title('Past-light cone in comoving coordinates')

plt.savefig('plc.png')
plt.show()
