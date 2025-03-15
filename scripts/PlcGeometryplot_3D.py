
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
import argparse


########################################## Functions definition ##########################################

# Plot 3D cubes 

def plot_cube( shifted_origin , side_lenght , fate , ax ):

    # Dictionary with a different colors for a different fate of the cube face

    if ( fate > 0 ):

        color={ 0: "gray", # Perché stiamo usando anche zero se c'è un if che dice > 0 ?
                1: "blue", # Starting cube 
                2: "cyan", # 
                3: "palegreen", #
                4: "turquoise" } #
        
    # Control Transparency

        alpha={ 0: 0.1,
                1: 0.3,
                2: 0.3,
                3: 0.3,
                4: 0.3 }

    # Generating vertices of a unitary cube 

        base_vector = [ 0.0 , 1.0 ]
        vertices = np.array(list(product( base_vector, base_vector , base_vector ))) # Cartesian product of the base vector 
        
    # Translation of generic vertices on the vertices of box replications that are used to sample the light cone 
     
        for i in range(3):

            vertices[ : , i ] *= side_lenght[i]     # L is the side length of the cube                                      
            vertices[ : , i ] += shifted_origin[i]  # C is the starting vertex of the cube (starting and replications)
        
    # Scatter plot of all vertices generated before

        ax.scatter3D( vertices[: , 0], vertices[: , 1] , vertices[: , 2] , s=5 )  

    # Creating X , Y , Z base data point for plotting the cubes using plot_surface
    # X and Y are a 2D array of points x and y while Z is used to indicate the 2D array of heights for x and y points
 
        X , Y = np.meshgrid( base_vector , base_vector ) 
        Z = np.ones(( 2 , 2))

    # Plotting the faces of the cubes 
    # We have to plot each face separately because it could potentially have a different color from the others according to the fate value
                                                                                                                  
        ax.plot_surface(X*side_lenght[0] + shifted_origin[0] , Y*side_lenght[1] + shifted_origin[1] , Z*shifted_origin[2]                    ,  alpha = alpha[fate], color = color[fate])   # Bottom Z face             
        ax.plot_surface(X*side_lenght[0] + shifted_origin[0] , Y*side_lenght[1] + shifted_origin[1] , Z*(side_lenght[2] + shifted_origin[2]) ,  alpha = alpha[fate], color = color[fate])   # Top Z face    
        ax.plot_surface(X*side_lenght[0] + shifted_origin[0] ,                    shifted_origin[1] , Y*side_lenght[2] +  shifted_origin[2]  ,  alpha = alpha[fate], color = color[fate])   # Right Y face        
        ax.plot_surface(X*side_lenght[0] + shifted_origin[0] ,   side_lenght[1] + shifted_origin[1] , Y*side_lenght[2] +  shifted_origin[2]  ,  alpha = alpha[fate], color = color[fate])   # Left Y face
        ax.plot_surface(                   shifted_origin[0] , X*side_lenght[1] + shifted_origin[1] , Y*side_lenght[2] +  shifted_origin[2]  ,  alpha = alpha[fate], color = color[fate])   # Left X face
        ax.plot_surface(  side_lenght[0] + shifted_origin[0] , X*side_lenght[1] + shifted_origin[1] , Y*side_lenght[2] +  shifted_origin[2]  ,  alpha = alpha[fate], color = color[fate])   # Right X face

# Plot distance

def plot_dist( C , L , V , r , ax ):

    vers  = ( C + L/2 - V ) /np.linalg.norm( C + L/2 - V )
    
    ax.plot( np.array([ V[0] , V[0] + vers[0] * r ]), 
             np.array([ V[1] , V[1] + vers[1] * r ]),
             np.array([ V[2] , V[2] + vers[2] * r ]))
    
# Plot the cone lateral surface and the dome 

def plot_cone( starting_vertex , cone_direction , cone_aperture , ax, rstart = 0.0 , rstop = 1.0 , N_r = 5.0 , N_phi = 10.0 ):
    
    # Step between radius of concentric circumferences  

    deltar = ( rstop - rstart) / N_r  
    
    # Normalizing the cone axis of symmetry

    norm_direction = cone_direction/np.linalg.norm( cone_direction ) 
   
    # Find the perpendicular direction to the axis of symmetry 
     
    x_versor = np.array([ 1. , 0. , 0. ])    
    p = np.cross( norm_direction , x_versor ) # Cross product of norm_direction and x_versor in R^3 is a vector perpendicular to both norm_direction and x_versor 
    
    # Check perpendicularity
    
    if np.linalg.norm(p) == 0: # i.e. p vector is parallel to norm_direction

        z_versor = np.array([ 0. , 0. , 1.])
        p = np.cross( norm_direction , z_versor ) 
    
    # Normalizing p 

    p /= np.linalg.norm(p)

    # Find the perpendicular direction to p 

    q = np.cross( norm_direction , p ) 
    
    # p and q form a basis for the plane perpendicular to the cone axis

    # Plot lateral surface of the cone
    
    AA = min( cone_aperture , np.pi ) # The aperture angle must be less than or equal to pi. If the aperture angle is less than pi, the cone is truncated at a certain height
    
    if cone_aperture < np.pi : 

        # Loop over polar angle 
        
        for i_phi in range( N_phi ):  # N_phi is the number of polar angles at which the generatrixes of the cone will be sampled
            
            phi_1 = i_phi * 2.0 * np.pi / (N_phi)   # Start value
            phi_2 = ( i_phi + 1 ) * 2.0 * np.pi / (N_phi)  # Increased value
            
            # Loop over the radius

            for i_r in range( N_r ): # N_r is the number of radii

                r_1 = i_r * deltar + rstart # Start value
                r_2= ( i_r + 1 ) * deltar + rstart # Increased value
                
                # For each combination of angles and radii, it calculates the position of two points on the surface of the cone using trigonometry and the basis vectors calculated earlier 
                # It then plots a line segment connecting these two points on the given matplotlib axis object. The result is a series of line segments that together form the lateral
                # surface of the cone. The aperture angle is taken into account when determining the positions of the points, and the number of samples along the surface of the cone is
                # determined by the input parameters N_r and N_phi

                X = np.array([[ starting_vertex[0] + r_1 * np.cos(AA) * norm_direction[0] + r_1 * np.sin(AA) * ( p[0] * np.cos(phi_1) + q[0] * np.sin(phi_1) ),  
                                starting_vertex[0] + r_2 * np.cos(AA) * norm_direction[0] + r_2 * np.sin(AA) * ( p[0] * np.cos(phi_1) + q[0] * np.sin(phi_1) )],  

                              [ starting_vertex[0] + r_1 * np.cos(AA) * norm_direction[0] + r_1 * np.sin(AA) * ( p[0] * np.cos(phi_2) + q[0] * np.sin(phi_2) ),
                                starting_vertex[0] + r_2 * np.cos(AA) * norm_direction[0] + r_2 * np.sin(AA) * ( p[0] * np.cos(phi_2) + q[0] * np.sin(phi_2) )]])
                
                Y = np.array([[ starting_vertex[1] + r_1 * np.cos(AA) * norm_direction[1] + r_1 * np.sin(AA) * ( p[1] * np.cos(phi_1) + q[1] * np.sin(phi_1) ),
                                starting_vertex[1] + r_2 * np.cos(AA) * norm_direction[1] + r_2 * np.sin(AA) * ( p[1] * np.cos(phi_1) + q[1] * np.sin(phi_1) )],

                              [ starting_vertex[1] + r_1 * np.cos(AA) * norm_direction[1] + r_1 * np.sin(AA) * ( p[1] * np.cos(phi_2) + q[1] * np.sin(phi_2) ),
                                starting_vertex[1] + r_2 * np.cos(AA) * norm_direction[1] + r_2 * np.sin(AA) * ( p[1] * np.cos(phi_2) + q[1] * np.sin(phi_2) )]])
                
                Z = np.array([[ starting_vertex[2] + r_1 * np.cos(AA) * norm_direction[2] + r_1 * np.sin(AA) * ( p[2] * np.cos(phi_1) + q[2] * np.sin(phi_1) ),
                                starting_vertex[2] + r_2 * np.cos(AA) * norm_direction[2] + r_2 * np.sin(AA) * ( p[2] * np.cos(phi_1) + q[2] * np.sin(phi_1) )],

                              [ starting_vertex[2] + r_1 * np.cos(AA) * norm_direction[2] + r_1 * np.sin(AA) * ( p[2] * np.cos(phi_2) + q[2] * np.sin(phi_2) ),
                                starting_vertex[2] + r_2 * np.cos(AA) * norm_direction[2] + r_2 * np.sin(AA) * ( p[2] * np.cos(phi_2) + q[2] * np.sin(phi_2) )]])
                
                # Plot of each line segment of the lateral surface

                ax.plot_wireframe(X,Y,Z, rcount=1, ccount=1, color='yellow' , alpha=0.3)  

    # Plot the dome of the cone as a 3D wireframe surface with axis along the vector norm_direction, centered at the point starting_vertex
    # The wireframe is generated by sampling points on the surface using a polar coordinate system and then converting them to Cartesian coordinates 
    
    N_cap = max( int( 1.5 * N_phi * AA / np.pi) + 1, 5 ) # Determines the number of points along the circular cross-section of the wireframe that is being plotted
  
    # Specify the step sizes of the meshgrid in a single variable
    
    a = complex( 0 , N_phi+1 ) # real part 0 and imaginary part N_phi+1, which specifies the step size along the phi dimension
    b = complex( 0 , N_cap ) # real part 0 and imaginary part N_cap, which specifies the step size along the r dimension
    
    # Creates a dense grid of points in the polar coordinate system
    
    u , v = np.mgrid[ 0 : 2*np.pi : a, 0 : AA : b ] # u and v represent the polar coordinates of the sampled points. Meshgrid_shape = ( N_phi+1 , Ncap )
                                                    # u is the angle [0:2*np.pi]
                                                    # v is the radius (in rad)
    
    # Cartesian coordinates 

    x = starting_vertex[0] + rstop * ( np.cos(v) * norm_direction[0] + np.sin(v) * ( p[0] * np.cos(u) + q[0] * np.sin(u) ) )
    y = starting_vertex[1] + rstop * ( np.cos(v) * norm_direction[1] + np.sin(v) * ( p[1] * np.cos(u) + q[1] * np.sin(u) ) )
    z = starting_vertex[2] + rstop * ( np.cos(v) * norm_direction[2] + np.sin(v) * ( p[2] * np.cos(u) + q[2] * np.sin(u) ) )
    
    # Plot the dome surface 

    ax.plot_wireframe(x,y,z, rstride=5, cstride=5, color='yellow', alpha=0.5) 


# Input #
 
parser = argparse.ArgumentParser()
parser.add_argument('--Geometry_cat', type=str, help='Input plc geometry catalog')
args = parser.parse_args()


################################# Read Input ###################################################

catalog = open( args.Geometry_cat , "r" )  # Open catalog
lines = catalog.readlines() # read lines

# Gathering the cube and PLC geometric infromation 

# Cube

N_cubes = int( lines[0].split()[3] )                                 # Number of replicated cubes
Cube_side_lenght = np.array( lines[4].split()[3:] ).astype( float )  # Side lenght of the cube
IPD = float( lines[6].split()[3] )                                   # Intra particle distance

# PLC

rstart = float( lines[1].split()[3] )                                # Starting distance i.e. cone vertex
rstop = float( lines[1].split()[4] )                                 # Max distance from cone vertex  

Plc_vertex = np.array( lines[2].split()[3:] ).astype( float )        # Plc vertex
Plc_direction = np.array( lines[3].split()[3:] ).astype( float )     # Plc simmetry axis 

Aperture_angle = float( lines[5].split()[3] )                        # Aperture angle

# Gathering simulation output 

replication_origin, rmin, rmax, fate = [], [], [], []

for i in range(8 , 8 + N_cubes) :
    
    # Origin of the replicated cube

    replication_origin.append( np.array(lines[i].split()[1:4]).astype(int) ) 
    
    # Min and max distance from the origin

    rmin.append( float( lines[i].split()[4] ) )                    
    rmax.append( float( lines[i].split()[5] ) )                       
    
    # Fate of the cube. See plot_cube for more information

    fate.append( float( lines[i].split()[6] ) )                       

# Closing geometry catalog

catalog.close()

# Formatting simulation output as array for plotting aims 
    
replication_origin = np.asarray( replication_origin )
rmin = np.asarray( rmin )
rmax = np.asarray( rmax )
fate = np.asarray( fate )

# Convert the aperture angle in rad

Aperture_angle *= np.pi/180.

# Rescaling all the geometric parameters using the IPD

Plc_vertex *= IPD
Cube_side_lenght *= IPD
rstart *= IPD
rstop *= IPD
rmin *= IPD
rmax *= IPD

################################# Make Plot ###################################################

# General Plot style

plt.style.use('dark_background') 
fig = plt.figure(figsize=(12, 9), dpi=800)
fig.tight_layout()

##################### First Plot ############################

# 3D view

# General settings

ax = fig.add_subplot(221, projection='3d')
ax.set_aspect("equal")
ax.grid(True)

# X-axis

ax.xaxis.pane.fill = False
ax.set_xlabel('X')
ax.xaxis.set_tick_params(labelsize=6)

# Y-axis

ax.yaxis.pane.fill = False
ax.set_ylabel('Y')
ax.yaxis.set_tick_params(labelsize=6)

# Z-axis

ax.zaxis.pane.fill = False
ax.set_zlabel('Z')
ax.zaxis.set_tick_params(labelsize=6)

# Point of view

ax.view_init(25, -155) # First value elevation, second azimuth

# Plot cubes

for i in range( N_cubes ) :
    
    plot_cube( replication_origin[i] * Cube_side_lenght , Cube_side_lenght , fate[i] , ax )

# Plot cone surface

plot_cone( Plc_vertex , Plc_direction , Aperture_angle , ax , rstart=rstart , rstop = rstop , N_r = 10 , N_phi = 50 )

# Optimized scale on the axis

plt.autoscale()

##################### Second Plot ############################

# X-Y plane projection

# General settings

ax = fig.add_subplot(222, projection='3d')
ax.set_aspect("equal")
ax.grid(True)

# X-axis

ax.xaxis.pane.fill = False
ax.set_xlabel('X')
ax.xaxis.set_tick_params(labelsize=6)

# Y-axis

ax.yaxis.pane.fill = False
ax.set_ylabel('Y')
ax.yaxis.set_tick_params(labelsize=6)

# Z-axis

ax.zaxis.pane.fill = False
ax.set_zticks([])  # For graphical reason we set zticks to void

# Point of view == X-Y projection

ax.view_init(90, -90)

# Plot cubes

for i in range( N_cubes ) :
    
    plot_cube( replication_origin[i] * Cube_side_lenght , Cube_side_lenght , fate[i] , ax )

# Plot cone surface

plot_cone( Plc_vertex , Plc_direction , Aperture_angle , ax , rstart=rstart , rstop = rstop , N_r = 10 , N_phi = 50 )

# Optimized scale on the axis

plt.autoscale()

##################### Third Plot ############################

# Y-Z plane projection

# General settings

ax = fig.add_subplot(223, projection='3d')
ax.set_aspect("equal")
ax.grid(True) 

# X-axis

ax.xaxis.pane.fill = False
ax.set_xticks([]) # For graphical reason we set xticks to void

# Y-axis

ax.yaxis.pane.fill = False
ax.set_ylabel('Y')
ax.yaxis.set_tick_params(labelsize=6)

# Z-axis

ax.zaxis.pane.fill = False
ax.set_zlabel('Z')
ax.zaxis.set_tick_params(labelsize=6)

# Point of view == Y-Z projection

ax.view_init(0, 180)

# Plot cubes

for i in range( N_cubes ) :
    
    plot_cube( replication_origin[i] * Cube_side_lenght , Cube_side_lenght , fate[i] , ax )

# Plot cone surface

plot_cone( Plc_vertex , Plc_direction , Aperture_angle , ax , rstart=rstart , rstop = rstop , N_r = 10 , N_phi = 50 )

# Optimized scale on the axis

plt.autoscale()

##################### Fourth Plot ############################

# X-Z plane projection

# General settings

ax = fig.add_subplot(224, projection='3d')
ax.set_aspect("equal")
ax.grid(True) 

# X-axis

ax.xaxis.pane.fill = False
ax.set_xlabel('X')
ax.xaxis.set_tick_params(labelsize=6)

# Y-axis

ax.yaxis.pane.fill = False
ax.set_yticks([]) 

# Z-axis

ax.zaxis.pane.fill = False
ax.set_zlabel('Z')
ax.zaxis.set_tick_params(labelsize=6)

# Point of view == X-Z projection

ax.view_init(0, - 90)

# Plot cubes

for i in range( N_cubes ) :
    
    plot_cube( replication_origin[i] * Cube_side_lenght , Cube_side_lenght , fate[i] , ax )

# Plot cone surface

plot_cone( Plc_vertex , Plc_direction , Aperture_angle , ax , rstart=rstart , rstop = rstop , N_r = 10 , N_phi = 50 )

# Optimized scale on the axis

plt.autoscale()

############################# Save Plot ########################

plt.subplots_adjust( wspace=0.005, hspace=0.005) # Settings space between the four subplots 

# Set name for the output plot

name = args.Geometry_cat.split('.')[1]

plt.savefig(f'{name}_PLC.png')
