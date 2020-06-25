import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpmath import *
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from math import pi
from scipy.interpolate import InterpolatedUnivariateSpline
import time

# Cosmological parameters
omegam = 1.
omegal = 1. - omegam
omegar = 0.
omegak = 0.
H0     = 67.8
B      = omegam**(-1.) - 1.

# Cosmo functions
def H(a):
    return H0 * np.sqrt(omegar / a**4. + omegam / a**3. + omegak / a**2. + omegal)

def Om(a):
    return omegam * H0**2. / (H(a)**2. * a**3.)

def Ol(a):
    return omegal * H0**2. / H(a)**2.

def delta(la):
    return ((1. - la[0]) * (1. - la[1]) * (1. - la[2]))**(-1.) - 1.

# Function that returns value of scale factor at collapse of the first axis
def coll_time(lam1, lam2, lam3):

    # Initial conditions
    x0 = 1.e-3
    f0 = [ lam1 * x0, \
           lam2 * x0, \
           lam3 * x0, \
           lam1 * x0 / (lam1 * x0 - 1.), \
           lam2 * x0 / (lam2 * x0 - 1.), \
           lam3 * x0 / (lam3 * x0 - 1.), \
           lam1 * x0, \
           lam2 * x0, \
           lam3 * x0]

    # Step between X values
    dx = 1.e-3

    # First derivatives of the lambda parameters
    def LA(i, a, d, la, lv, ld):
        return lv[i] * (la[i] - 1.) / a


    def LV(i, a, d, la, lv, ld):
        return (-1.5 * Om(a) * ld[i] + lv[i] * (-1. + Om(a)/2. - Ol(a) - lv[i])) / a


    def LD(i, a, d, la, lv, ld):
        summ = 0.
        for j in range(3):
            if i == j:
                continue
            if la[i] == la[j]:
                summ += 0.
            else:    
                summ += (ld[j] - ld[i]) * ((1. - la[i])**2. * (1. + lv[i]) - (1. - la[j])**2. * (1. + lv[j])) / ((1. - la[i])**2. - (1. - la[j])**2.)
        return (- (1. + d) * (d + 2.5)**(-1.) * (ld[i] + 5. / 6.) * sum(lv) + (ld[i] + 5. / 6.) * (3. + sum(lv)) - (d + 2.5) * (1. + lv[i]) + summ) / a


    # Function returning a vector with the rhs of the ODE. NOTE: the integrator
    # works only on first order ODE; to integrate a second order ODE I will use a
    # temporary variable equal to F', and integrate both functions simultaneously.

    # p = [l_a1, l_a2, l_a3, l_v1, l_v2, l_v3, l_d1, l_d2, l_d3]
    # returns [l_a1', l_a2', l_a3', l_v1', l_v2', l_v3', l_d1', l_d2', l_d3']
    def rhs(a, p):
        la = p[0:3]
        lv = p[3:6]
        ld = p[6:9]
        d  = sum(ld)
        return  [LA(0, a, d, la, lv, ld), \
                 LA(1, a, d, la, lv, ld), \
                 LA(2, a, d, la, lv, ld), \
                 LV(0, a, d, la, lv, ld), \
                 LV(1, a, d, la, lv, ld), \
                 LV(2, a, d, la, lv, ld), \
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
    LD1 = np.empty(0)
    LD2 = np.empty(0)
    LD3 = np.empty(0)

    # Integration. NOTE: the overdensity might never collapse, it is
    # necessary to manually set a time at which the integration is stopped
    while solver.successful() and solver.t <= 5.5:
        solver.integrate(solver.t + dx)
        if solver.y[0] > 1.:
            break
        A   = np.append(A, solver.t)
        LA1 = np.append(LA1, solver.y[0])
        LD1 = np.append(LD1, solver.y[0])
        LD2 = np.append(LD2, solver.y[0])
        LD3 = np.append(LD3, solver.y[0])

    # Interpolate to find a for which la_1(a) = 1.
    if A[-1] > 5.4:
        c_time = -0.1
    else:
        s = InterpolatedUnivariateSpline(LA1, A, k=3) # Extrapolation with cubic spline
        c_time = s(1.)

    return [c_time, LA1, A, LD1]
 
# Timing one integration:
start = time.time()
coll_time(1./3., 1./3., 1./3.)
end = time.time()
print "Evaluating one CT takes", end-start, "s"

# For a fixed value of d_l = delta, find collapse time for different lambdas 
AC = np.empty(0)
L1 = np.empty(0)
L2 = np.empty(0)
L3 = np.empty(0)

x = np.linspace(0.,1.,10.) # l1 - l2
y = np.linspace(0.,2.,10.) # l2 - l3
d_l = 1.

for i in range(len(x)):
    for j in range(len(y)):
        l1 = (d_l + 2. * x[i] + 1. * y[j]) / 3.
        l2 = (d_l - 1. * x[i] + 1. * y[j]) / 3.
        l3 = (d_l - 1. * x[i] - 2. * y[j]) / 3.
        elle = np.sort([l1, l2, l3])
        a_c = coll_time(elle[2], elle[1], elle[0])[0]
        AC = np.append(AC, a_c)
        L1 = np.append(L1, elle[2])
        L2 = np.append(L2, elle[1])
        L3 = np.append(L3, elle[0])
        print x[i], y[j], AC[-1]

# Save data
np.savetxt('/data/my/phd/ell_coll/cosmodep/eds_sng_NL_d5.dat', \
           np.column_stack([L1, L2, L3, AC]))
