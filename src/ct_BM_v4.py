import numpy as np
import time
import sys
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpmath import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import InterpolatedUnivariateSpline

# Cosmological parameters
omegam = 1.
omegal = 1. - omegam
omegar = 0.
omegak = 0.
H0     = 67.8
B      = omegam**(-1.) - 1.

# Cosmo functions
def E(a):
    return np.sqrt(omegar / a**4. + omegam / a**3. + omegak / a**2. + omegal)

def Eprime(a):
    return (-4. * omegar / a**5. - 3. * omegam / a**4. - 2. * omegak / a**3.) / (2. * E(a))

def b(a1, a2, a3):    # Elliptic integral
    return (2./3.) * (a1 * a2 * a3 * elliprd(a1**2., a2**2., a3**2.) - 1.)
        
def c(d,b,l):
    return (1. + d) / 3. + (d * b / 2.) + l
                
def delta(a, a1, a2, a3):
    return a**3. / (a1 * a2 * a3) - 1.
                
# Function that returns value of scale factor at collapse of the first axis
def coll_time(lam1,lam2,lam3):

	# Initial conditions
	x0   = 1.e-3
	f0 = [ \
	       (1. -      x0 * lam1) * x0,  \
	        1. - 2. * x0 * lam1      ,  \
	       (1. -      x0 * lam2) * x0,  \
	        1. - 2. * x0 * lam2      ,  \
	       (1. -      x0 * lam3) * x0,  \
	        1. - 2. * x0 * lam3 ]
	
	# Step between X values
	dx = 1.e-3

	# Function to integrate
        # x = first derivative of a_i,
        # y = a_i, e = E(a), f = Eprime(a), c = C_i
        def d2a(x, y, a, e, f, c):
            return - x * (1./a + f / e) - y * (3. * omegam * c / (2. * a**5. * e**2.) - omegal / (a**2. * e**2.))

	# Function returning the right hand side of the ODE. NOTE: the integrator
	# works only on first order ODE; to integrate a second order ODE I will use a
	# temporary variable equal to F', and integrate both functions simultaneously.
	def rhs(x, p):
	    d = delta(x, p[0], p[2], p[4])
	    b1 = b(p[2], p[4], p[0])
	    b2 = b(p[4], p[0], p[2])
	    b3 = b(p[0], p[2], p[4])
            l1 = 5. * b1 / 4.    # NL approx
            l2 = 5. * b2 / 4.
            l3 = 5. * b3 / 4.
            #l1 = (bc(x)/ x0) * (x0 * lam1 - delta(x0, f0[0], f0[2], f0[4]) / 3.)    # lin approx
            #l2 = (bc(x)/ x0) * (x0 * lam2 - delta(x0, f0[0], f0[2], f0[4]) / 3.)
            #l3 = (bc(x)/ x0) * (x0 * lam3 - delta(x0, f0[0], f0[2], f0[4]) / 3.)
	    C1    = c(d, b1, l1)
	    C2    = c(d, b2, l2)
	    C3    = c(d, b3, l3)
	    # The p array contains the last computed F value, and its first
	    # derivative.  This function must return the first and second order
	    # derivatives, hence the first value must always be p[1].
	    return  [p[1], d2a(p[1], p[0], x, E(x), Eprime(x), C1), \
	             p[3], d2a(p[3], p[2], x, E(x), Eprime(x), C2), \
	             p[5], d2a(p[5], p[4], x, E(x), Eprime(x), C3)]
	
	# Choose solver:Runge-Kutta of order 4
	solver = ode(rhs).set_integrator('dopri5')
	
	# Set initial conditions (NOTE: the set_initial_value method reset the t value
	# to zero, hence we must set t AFTER the set_initial_value call).
	solver.set_initial_value(f0)
	solver.t = x0
	
	# Output arrays: X and F values
	X  = np.empty(0)
	F1 = np.empty(0)
        F2 = np.empty(0)
	
        # Integration. NOTE: the overdensity might never collapse, it is
        # necessary to manually set a time at which the integration is stopped
	while solver.successful() and solver.t <= 5.5:
            solver.integrate(solver.t + dx)
            if solver.y[0] < 0.:
                break
            X  = np.append(X, solver.t)
            F1 = np.append(F1, solver.y[0])
            F2 = np.append(F2, solver.y[1])

        # Interpolate to find a for which a_1(a) = 0.
        if F2[-1] > 0.:
            c_time = -0.1
        else:
            Xcut  = X[np.argmax(F1):]
            F1cut = F1[np.argmax(F1):]
            order = F1cut.argsort()
            x = F1cut[order]
            y = Xcut[order]
            s = InterpolatedUnivariateSpline(x, y, k=3) # Extrapolation with cubic spline
            c_time = s(0.)
        
        return c_time

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
                a_c = coll_time(elle[2], elle[1], elle[0])
                AC = np.append(AC, a_c)
                L1 = np.append(L1, elle[2])
                L2 = np.append(L2, elle[1])
                L3 = np.append(L3, elle[0])
                print l1, l2, l3, a_c


## Plot
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot_trisurf(X, Y, BC, color = 'b', linewidth=0.5, antialiased=True, alpha=0.3, label='$a_c$')
#ax.plot_trisurf(X, Y, AC, color = 'r', linewidth=0.5, antialiased=True, alpha=0.3, label='$b_c$')
#ax.legend(loc=2)
#ax.set_xlim(0.,3.)
#ax.set_ylim(0.,6.)
#ax.set_zlim(0.,1.)
#ax.set_xlabel('$\lambda_1 - \lambda_2$')
#ax.set_ylabel('$\lambda_2 - \lambda_3$')
#ax.set_zlabel('$b_c$')
#plt.show()

# Save data
np.savetxt('/data/my/phd/ell_coll/cosmodep/eds_numerical_NL_d-1.dat', \
           np.column_stack([L1, L2, L3, AC]))
