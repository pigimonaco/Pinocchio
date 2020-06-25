import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpmath import *
from scipy.interpolate import InterpolatedUnivariateSpline

omegam = 0.3
omegal = 1. - omegam
H0     = 70.
x0     = 1.e-3

B         = omegam**(-1.) - 1.
fR0       = 1.e-4
H0_over_c = 0.00023333333333333333
R         = 10.

##########################
# Cosmological functions #
##########################

def E(a):
    return np.sqrt(omegal + omegam / a**3.)

def Om(a):
    return omegam / (E(a)**2. * a**3.)

def Ol(a):
    return omegal / E(a)**2.

def delta(la):
    return ((1. - la[0]) * (1. - la[1]) * (1. - la[2]))**(-1.) - 1.

# MG functions
def scr(x):
    ff = x**3 - 3. * x**2 + 3. * x
    return (1./3.) * min(1., ff)

def xi(a, la1, la2, la3, dd, Rin):
    x1 = a * (1. - la1)
    x2 = a * (1. - la2)
    x3 = a * (1. - la3)
    aa = x1*x2*x3
    bb = 4. * omegal / omegam
    coef = fR0 / (omegam * H0_over_c * H0_over_c * Rin * Rin)
    ff = coef * a**7 * (dd + 1.)**(1./3.) * (((1. + bb)/(1. + bb * a**3.))**2. - ((1. + bb)/(dd + 1. + bb * a**3.))**2.)
    return ff


# Eigenvalues of the deformation tensor;
# MUST BE lam1 > lam2 > lam3
lam1 = 1.5 
lam2 = 1.3
lam3 = 1.

# Initial conditions
f0 = [  lam1 * x0, \
        lam2 * x0, \
        lam3 * x0, \
       -lam1 * x0, \
       -lam2 * x0, \
       -lam3 * x0, \
        lam1 * x0, \
        lam2 * x0, \
        lam3 * x0]

# Step between X values
dx = 1.e-3

# First derivatives of the lambda parameters
def LA(i, a, d, la, lv, ld):
    return lv[i] * (la[i] - 1.) / a

def LV(i, a, d, la, lv, ld, F):
    return (Ol(a) - 0.5 * Om(a) - 1.5 * Om(a) * (1. + F) * ld[i] + 1.5 * Om(a) * (1. + lv[i]) - (1. + lv[i])**2.) / a

def LD(i, a, d, la, lv, ld):
    summ = 0.
    for j in range(3):
        if i == j:
            continue
        if la[i] == la[j]:
            summ = 0
        else:    
            summ += (ld[j] - ld[i]) * ((1. - la[i])**2. * (1. + lv[i]) - (1. - la[j])**2. * (1. + lv[j])) / ((1. - la[i])**2. - (1. - la[j])**2.)
    return (- (1. + d) * (d + 2.5)**(-1.) * (ld[i] + 5. / 6.) * sum(lv) + (ld[i] + 5. / 6.) * (3. + sum(lv)) - (d + 2.5) * (1. + lv[i]) + summ) / a


# Function returning a vector with the rhs of the ODE
# p = [l_a1, l_a2, l_a3, l_v1, l_v2, l_v3, l_d1, l_d2, l_d3]
# returns [l_a1', l_a2', l_a3', l_v1', l_v2', l_v3', l_d1', l_d2', l_d3']
def rhs(a, p):
    la = p[0:3]
    lv = p[3:6]
    ld = p[6:9]
    d  = sum(ld)
    F = scr(xi(a, p[0], p[1], p[2], d, R))
    #print d
    return  [LA(0, a, d, la, lv, ld), \
             LA(1, a, d, la, lv, ld), \
             LA(2, a, d, la, lv, ld), \
             LV(0, a, d, la, lv, ld, F), \
             LV(1, a, d, la, lv, ld, F), \
             LV(2, a, d, la, lv, ld, F), \
             LD(0, a, d, la, lv, ld), \
             LD(1, a, d, la, lv, ld), \
             LD(2, a, d, la, lv, ld)  ]          

# Choose solver:Runge-Kutta of order 4
solver = ode(rhs).set_integrator('dopri5')

# Set initial conditions (NOTE: the set_initial_value method reset the t value
# to zero, hence we must set t AFTER the set_initial_value call).
solver.set_initial_value(f0)
solver.t = x0

# Output arrays: X and F values
A   = np.empty(0)
LA1 = np.empty(0)
LA2 = np.empty(0)
LA3 = np.empty(0)
LV1 = np.empty(0)
LV2 = np.empty(0)
LV3 = np.empty(0)
LD1 = np.empty(0)
LD2 = np.empty(0)
LD3 = np.empty(0)

# Integration
while solver.successful():
    solver.integrate(solver.t + dx)
    if solver.y[0] > 1.:
        break
    A   = np.append(A, solver.t)
    LA1 = np.append(LA1, solver.y[0])
    LA2 = np.append(LA2, solver.y[1])
    LA3 = np.append(LA3, solver.y[2])
    LV1 = np.append(LV1, solver.y[3])
    LV2 = np.append(LV2, solver.y[4])
    LV3 = np.append(LV3, solver.y[5])
    LD1 = np.append(LD1, solver.y[6])
    LD2 = np.append(LD2, solver.y[7])
    LD3 = np.append(LD3, solver.y[8])

#Interpolation
s = InterpolatedUnivariateSpline(LA1, A, k=3) # Extrapolation with cubic spline
print s(1.)

#Plot
fig = plt.figure(figsize = (10,8))
ax  = plt.axes()
ax.plot(A, A             , color = 'grey', ls = '--' , label = '$a$')
ax.plot(A, A * (1. - LA1), color = 'r',    ls = '-', label = '$a_1$')
ax.plot(A, A * (1. - LA2), color = 'g',    ls = '-', label = '$a_2$')
ax.plot(A, A * (1. - LA3), color = 'b',    ls = '-', label = '$a_3$')
ax.set_ylabel(r'$a_i(a)$')
ax.set_xlabel(r'$a$')
plt.ylim((0.,0.5))
plt.xlim((0.,0.5))
ax.legend(loc = 2)
plt.show()


fig = plt.figure(figsize = (10,8))
ax  = plt.axes()
ax.set_yscale('log')
ax.plot(A, DELTA, color='r', label='NGS')
ax.plot(X, DELTA2, color='b', label='num')
ax.plot(A, DELTA3, color='g', label='sum ld')
plt.show()
