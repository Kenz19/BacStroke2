'''
This script is the main simulation script for BacStroke. It runs a simulation
of a single bacterium within a clinostat based on parameters found within the
configuration file. 

The bacteriums position within the clinostat is updated at each timestep
within the simulation. Where it moves to depends on its velocity, this depends
on four components:
    
    Gravity
    Diffusion
    Bacterial swimming
    Bacterial Tumbling
    Centripetal force
    Clinostat rotation
    
Any one of these parameters can be turned off via the configuation file, with the
exception of swimming. This is turned off by setting the swimming
velocity in the initial conditions file to 0 (last element in the line).
'''

# Imports #####################################################################

# modules
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
from matplotlib.patches import Rectangle
import sys
import random
import matplotlib as mpl

# external files
import Functions as f 
from Bacteria3D import Bacteria3D as bac # class

# style sheet

plt.rcParams.update({
    'text.usetex': False,
    'font.family': 'roman',
})

# #handle graph formatting and style
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed
MYLW=1.0 #line width

random.seed(66)

###############################################################################

def main(config_file, output_file, tumble_file, time_file, swimming_file, figure_output_file):
    
    # 1. read input files #####################################################
    
    # read config_file, obtain constants describing clinostat system
    
    # opening config file containing all file paths needed to execute code
    file = open(config_file, "r")
    
    # reading every entry from configuration file (containing constants etc)
    config = file.readlines()[1::3]
    
    # removing trailing new line (\n)
    for i in range(len(config)):
        config[i] = config[i].rstrip()
    
    # obtaining file names from config file
    inputfile_name = config[0]# input data file (containing becterium initial conditions)
    
    # Time parmeters of simulation
    time = 0.0 # [s]
    dt = float(config[1]) # length of each timestep [s]
    total_time = float(config[2]) # total length of similation [s]
    numstep = round(total_time/dt) # number of steps that simulation will take, rounded to int as used as np sizing
    
    # clinostat parameters 
    clino_rotation_rate = float(config[3]) # rotation rate of clinostat [RPM]
    omega = clino_rotation_rate*2*np.pi/60 # converting rotation rate to angular velocity [rad/s]
    
    # clinostat parameters
    R = float(config[4]) # outer radius of circular face of clinostat [m]
    r = float(config[5]) # inner radius of circular face of clinostat [m]
    H = float(config[6]) # length of clinostat down the z axis [m]
    
    # system constants
    density = float(config[7]) # density of medium in clinostat [kg/m^3]
    g = float(config[8]) # acceleration due to gravity [m/s^2]
    
    # coefficients
    rotational_diffusion_coefficient = float(config[9]) # inversely proportional to time it takes bacterium to forget direction its travelling in [1/s]
    viscosity_coefficient = float(config[10]) # viscosity coefficient of clinostat medium at room temp [Pa/s] (water during testing)
    diffusion_coefficient = float(config[11]) # diffusion coefficient for medium within clinostat at room temp [m^2/s]

    tumbling_rate = int(config[12]) # how often in a second a bacterium should tumble
    
    centrifugal_force_status = config[13] # does the bacterial dynamics include centrifugal force?
    
    # read in information about intial state of bacteium (position etc...) from initial conditions file
    
    # Determining length of inputfile = no. bacteria in system   
    with open(inputfile_name, 'r') as fp:
        lines = len(fp.readlines()) # no. bacteria present 
        
    fp.close() # close input file to ensure no mess or overwriting
    
    infile = open(inputfile_name, "r")# opening input file to use
    print(infile)

    bacteria = [] # will contain instances of Bacteria3D, i.e initial conditions of each bacterium
    
    # initialising bacteria instances
    for i in range(lines):
        bacteria.append(bac.new_b3d(infile, density))
        
    
    # 2. GENERATING INITIAL CONDITIONS ########################################
    
    initial_positions = np.zeros([lines, 3]) # array that records initial 3D position of bacteria in cartesian coords (XYZ)
    
    # reading inital position of bacteria from initial conditions file
    for i in range(len(bacteria)):        
        initial_positions[i] = bacteria[i].pos

    # Storage for data
    pos_array = np.zeros([lines, numstep, 3]) # xyz position of every bacteria every timestep
    time_array = np.zeros(numstep) 
    tumble_array = np.zeros(numstep) # track if bacterium tumbles or not
    swim_direction = np.zeros([lines, numstep, 3]) # track direction of bacteria swimming
    
    # generate
    for i in range(lines):
        bacteria[i].terminal_vel(viscosity_coefficient, density, g) # terminal velocity - remains a constant in the time integrator
        bacteria[i].rotational_vel(omega) # rotational velocity
        
        # initial planar position (needed for rotational components of the velocity)
        x = initial_positions[i][0]
        y = initial_positions[i][1]
        planar_pos = np.array([x, y, 0])
        
        bacteria[i].centripetal_force(viscosity_coefficient, density, omega, planar_pos, dt, centrifugal_force_status) # centrifugal velocity
        
        bacteria[i].update_vel(dt, diffusion_coefficient) # get intial velocity vector
        

    # 3. beginning of time integration  #######################################
    
    for i in range(numstep):  # Anything that happens per each timestep 
      
        # text progress bar
        if i%100 == 0:
            print(f'Progress: {i} out of {numstep}')
        
        time += dt # add another timestep to current time
        time_array[i] = time # storing current time in simulation
        
        for j in range(lines): # Anything for each bacterium
  
            # radius of bacterium
            a = bacteria[j].rad
            
            planar_pos = np.array([bacteria[j].pos[0], bacteria[j].pos[1], 0])
            #print('inner boundry, new planar magnitude is: ' + str(np.linalg.norm(planar_pos)))
            bacteria[j].centripetal_force(viscosity_coefficient, density, omega, planar_pos, dt, centrifugal_force_status)
            
            bacteria[j].rotational_vel(omega)
            
            # updating velocity of each bacterium [m/s]
           # bacteria[j].update_vel(dt, diffusion_coefficient)
            
            # # Update position of each bacterium, using updated velocity
            bacteria[j].update_pos(dt)
            
            # updating velocity of each bacterium [m/s]
            bacteria[j].update_vel(dt, diffusion_coefficient)
            
            #
            #bacteria[j].rotational_vel(omega)
  
            # updating velocity of each bacterium [m/s]
            #bacteria[j].update_vel(dt, diffusion_coefficient)
            
            # # Update position of each bacterium, using updated velocity
            #bacteria[j].update_pos(dt)
            
            # Update position of each bacterium, using updated velocity
            #bacteria[j].update_pos(dt) 
            #pos_array[j, i] = bacteria[j].pos # store each position
            
            # boundary conditions (x & y)
            
            #position on a 2D circle (set z = 0)
            planar_position = bacteria[j].pos - [0, 0, bacteria[j].pos[2]]
            planar_magnitude = np.linalg.norm(planar_position)
            
            # if planar_magnitude > 0.05:
            #     print(planar_magnitude)
            #print(planar_magnitude)
            
            # zone where boundry conditions are applied
            outer_bc_zone = R - a # outer radius minus bacterium radii
            inner_bc_zone = r + a # inner radius plus bacterium radii
            
            boundry_conditions = True
            
            if boundry_conditions == True:
                
                # if bacterium position, r, within 1 bacterial radii from the edge, apply boundry conditions
                if planar_magnitude >= outer_bc_zone:
                    #print(time)
                    #print('outer boundry, current planar magnitude is :' +str(planar_magnitude))
                    
                    # last updated velocity of bacterium
                    current_vel = bacteria[j].vel # maybe need to use np.copy?
    
                    # radial magnitude and direction of the velocity
                    planar_rad_mag, planar_rad_dir = f.radial_velocity(planar_position, bacteria[j].vel)
                    
                    # removing radial component of velocity, i.e setting velocity to its tangential component
                    bacteria[j].vel = current_vel - (planar_rad_mag*planar_rad_dir)
                                 
                    # saving swimming direction for output later
                    swim_direction[j, i] = bacteria[j].swim_direction # saving swimming direction
                    
                    # radial direction of position
                    rad_dir = np.copy(bacteria[j].pos)
                    rad_dir[2] = 0
                    rad_dir /= np.linalg.norm(rad_dir)
                    
                    # moving bacterium to outside of boundry condition zone
                    frac = 0.1   #Fraction of body size to set inside outer_bc_zone (make parameter later)
                    pos_z = bacteria[j].pos[2] # storing z coord of bacterium
                    bacteria[j].pos = (R - (1.0 + frac)*a)*rad_dir # setting position to some fraction outside the boundry zone but inside the clinostat, this only sets xy parameters
                    bacteria[j].pos[2] = pos_z # setting z parameter back to original
                    #pos_array[j, i] = bacteria[j].pos # storing new position
                    #print(pos_array[j, i])
                    
                    #xs, ys = bacteria[j].pos[0], bacteria[j].pos[1]
                    
                    #print(np.linalg.norm([xs, ys]))
                    
                    # saving variables for output
                    
                    # getting the swimming direciton vector for saving to output.
                    tumble = bac.tumble_probability(dt, tumbling_rate) # does bacterium tumble? 1 = yes, 0 = no
                    tumble_array[i] = tumble
                    bacteria[j].update_swimming_vel(omega, rotational_diffusion_coefficient, dt, tumble) # updating swimming velocity
                    swim_direction[j, i] = bacteria[j].swim_direction # saving swimming direction
                    
                    planar_pos = np.array([bacteria[j].pos[0], bacteria[j].pos[1], 0])
                    #print('outer boundry, new planar magnitude is: ' + str(np.linalg.norm(planar_pos)))
                    
                    #bacteria[j].centripetal_force(viscosity_coefficient, density, omega, planar_pos, dt)
                                   
                    # # apply end boundry conditions if applicable, along side wall conditions
                    
                    # just before one wall
                    if bacteria[j].pos[2] >= (H - a): 
                        
                        # updating z position and re-recording the position
                        bacteria[j].pos[2] = H - (2*a)
                        #pos_array[j, i] = bacteria[j].pos
                        bacteria[j].vel[2] *= 0
                    
                    # close to the other wall
                    if bacteria[j].pos[2] <= (0 + a):
                        
                        # updating z position and re-recording the position
                        bacteria[j].pos[2] = 0 + (2*a)
                        #pos_array[j, i] = bacteria[j].pos
                        bacteria[j].vel[2] *= 0
                
                # apply inner wall boundry conditions, if applicable
                elif planar_magnitude <= inner_bc_zone:
                    #print('inner boundry')
                    
                    #print(planar_magnitude, inner_bc_zone)
                    
                    # last updated velocity of bacterium
                    current_vel = bacteria[j].vel
    
                    # radial magnitude and direction of the velocity
                    planar_rad_mag, planar_rad_dir = f.radial_velocity(planar_position, bacteria[j].vel)
                    
                    # removing radial component of velocity, i.e setting velocity to its tangential component
                    bacteria[j].vel = current_vel - (planar_rad_mag*planar_rad_dir)
                    
                    # saving swimming direction for output later
                    swim_direction[j, i] = bacteria[j].swim_direction # saving swimming direction
                    
                    # radial direction of position
                    rad_dir = np.copy(bacteria[j].pos)
                    rad_dir[2] = 0
                    rad_dir /= np.linalg.norm(rad_dir) # planar radial unit vector of current position
                    
                    # moving bacterium to outside of boundry condition zone
                    frac = 0.1   #Fraction of body size to set inside inner_bc_zone (make parameter later)
                    pos_z = bacteria[j].pos[2] # storing z coord of bacterium
                    bacteria[j].pos = (r + (1.0 + frac)*a)*rad_dir # setting position to some fraction outside the boundry zone but inside the clinostat, this only sets xy parameters
                    bacteria[j].pos[2] = pos_z # setting z parameter back to original
                    #pos_array[j, i] = bacteria[j].pos # storing new position
                    
                    # # saving variables for output
                    
                    # getting the swimming direciton vector for saving to output.
                    tumble = bac.tumble_probability(dt, tumbling_rate) # does bacterium tumble? 1 = yes, 0 = no
                    tumble_array[i] = tumble
                    bacteria[j].update_swimming_vel(omega, rotational_diffusion_coefficient, dt, tumble) # updating swimming velocity
                    swim_direction[j, i] = bacteria[j].swim_direction # saving swimming direction
                    
                    
                    #planar_pos = np.array([bacteria[j].pos[0], bacteria[j].pos[1], 0])
                    #print('inner boundry, new planar magnitude is: ' + str(np.linalg.norm(planar_pos)))
                    #bacteria[j].centripetal_force(viscosity_coefficient, density, omega, planar_pos, dt)
                    
                    # apply end boundry conditions if applicable, along side wall conditions
                    
                    # just before one wall
                    if bacteria[j].pos[2] >= (H - a): 
                        
                        # updating z position and re-recording the position
                        bacteria[j].pos[2] = H - (2*a)
                        #pos_array[j, i] = bacteria[j].pos
                        
                        bacteria[j].vel[2] *= 0 #set z velocity to 0
                    
                    # close to the other wall
                    if bacteria[j].pos[2] <= (0 + a):
                        
                        # updating z position and re-recording the position
                        bacteria[j].pos[2] = 0 + (2*a)
                        #pos_array[j, i] = bacteria[j].pos
                        
                        # setting z velocity to be 0, therefore velocity only in xy plane
                        bacteria[j].vel[2] *= 0
                
                # upper z boundry condition
                elif bacteria[j].pos[2] >= (H - a): 
                    
                   # updating z position and re-recording the position
                   bacteria[j].pos[2] = H - (2*a)
                   #pos_array[j, i] = bacteria[j].pos
                   
                   bacteria[j].vel[2] *= 0
                 
                # lower z boundry condition   
                elif bacteria[j].pos[2] <= (0 + a):
                    
                    # updating z position and re-recording the position
                    bacteria[j].pos[2] = 0 + (2*a)
                    #pos_array[j, i] = bacteria[j].pos
                    
                    bacteria[j].vel[2] *= 0
                    
                # dont apply wall boundry conditions
                else:
                    #print('no boundry, current planar magnitude is ' + str(planar_magnitude))
                    
                    # updating velocity of each bacterium [m/s]
                    #bacteria[j].update_vel(dt, diffusion_coefficient) - this makes the bacterium spiral outwards
                    
                    # # Update position of each bacterium, using updated velocity
                    #bacteria[j].update_pos(dt) 
                    #pos_array[j, i] = bacteria[j].pos # store each position
                    
                    # updating rotational velocity as it depends on position
                    #bacteria[j].rotational_vel(omega)
                    
                    # get planar positions
                    x = bacteria[j].pos[0]
                    y = bacteria[j].pos[1]
                    
                    planar_pos = np.array([x, y, 0])
                    
                    #bacteria[j].centripetal_force(viscosity_coefficient, density, omega, planar_pos, dt)
                    
                    # updating the swimming velocity and saving variables
                    tumble = bac.tumble_probability(dt, tumbling_rate) # does bacterium tumble? 1 = yes, 0 = no
                    tumble_array[i] = tumble
                    
                    bacteria[j].update_swimming_vel(omega, rotational_diffusion_coefficient, dt, tumble) # updating swimming velocity
                    swim_direction[j, i] = bacteria[j].swim_direction # saving swimming direction
           
            
           
            #bacteria[j].rotational_vel(omega)
           
            # record position after all relevant conditions applied
            pos_array[j, i] = bacteria[j].pos
            
            # get planar positions
            x = pos_array[j, i][0]
            y = pos_array[j, i][1]
            
            planar_pos = np.array([x, y, 0])
            
            if np.linalg.norm(planar_pos) > 0.05:
                print(np.linalg.norm(planar_pos))
            #print(pos_array[j, i])
          
           # if statement for inner boundry conditions, including the same thing about the z conditions 
    
    #pos_array[0] = pos_array[0][::100]
    # saving parameters to output file
    np.savetxt(output_file, pos_array[0], delimiter=",")
    #np.savetxt(tumble_file, tumble_array, delimiter=",")
    np.savetxt(time_file, time_array, delimiter=",") # outputting time for ease of analysis
    #np.savetxt(swimming_file, swim_direction[0], delimiter = ",")
    # 4. PLOTTING ################################################################
            
    # plotting path of bacteria
    nopoints = 50
    
    x_coords = pos_array[0,:,0]
    y_coords = pos_array[0,:,1]
    z_coords = pos_array[0,:,2]
    
    
    # fig, ax = plt.subplots(figsize=(10, 10))
    # # cir = plt.Circle((0, 0), R, facecolor='#c7c7c7', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
    # # ax.add_patch(cir)
    # # cir2 = plt.Circle((0, 0), r, facecolor='white', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
    # # ax.add_patch(cir2)
    # ax.scatter(x_coords[::nopoints], y_coords[::nopoints], s = 15, c = 'navy')
    # ax.set_xlabel('x (m)')#, fontsize = 30)
    # ax.set_ylabel('y (m)')#, fontsize = 30)
    # ax.tick_params(axis='x')#, labelsize=20)
    # ax.tick_params(axis='y')#, labelsize=20)
    # fig.suptitle('Diffusion', fontsize = 40)
    # plt.show()
    # ax = fig.add_subplot(projection='3d')
    # ax.set_xlabel('x (m)')
    # ax.set_ylabel('y (m)')
    
    #ax.scatter(x_coords[::nopoints], z_coords[::nopoints], y_coords[::nopoints])
    plt.show()
    
    fig,ax = plt.subplots(1,4, figsize = (35,12))
    ax[0].scatter(x_coords[::nopoints], z_coords[::nopoints], s=10, label = 'Bacteria 1')
    ax[0].set_xlabel('x')#, fontsize = 30)
    ax[0].set_ylabel('z')#, fontsize = 30)
    ax[0].tick_params(axis='x')#, labelsize = 20)
    ax[0].tick_params(axis='y')#, labelsize=20)
    #ax[0].legend()
    
    # view of circular face of clinostat
    cir = plt.Circle((0, 0), R, facecolor='#c7c7c7', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
    ax[1].add_patch(cir)
    cir2 = plt.Circle((0, 0), r, facecolor='white', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
    ax[1].add_patch(cir2)
    ax[1].plot(x_coords[::nopoints], y_coords[::nopoints], 'o',  zorder = 1, label = 'Bacteria 1')
    ax[1].set_xlabel('x')#, fontsize = 30)
    ax[1].set_ylabel('y')#, fontsize = 30)
    ax[1].tick_params(axis='x')#, labelsize=20)
    ax[1].tick_params(axis='y')#, labelsize=20)

    # plotting y position against time
    # ax[2].scatter(time_array[::10]/1e2, x_coords[::10], s = 5, label = 'Bacteria 1')
    ax[2].plot(time_array[::nopoints]/1e2, y_coords[::nopoints], 'o', label = 'Bacteria 1')
    ax[2].set_xlabel('t ($10^2$s)')#, fontsize = 30)
    ax[2].set_ylabel('y')#, fontsize = 30)
    ax[2].tick_params(axis='x')#, labelsize=20)
    ax[2].tick_params(axis='y')#, labelsize=20)
    
    ax[3].scatter(time_array[::nopoints]/1e2, np.sqrt(y_coords[::nopoints]**2 + x_coords[::nopoints]**2))
    
    fig.suptitle('Sim length = ' + str(total_time) + 's' + ', $\Delta$t = ' + str(dt) + 's' + ', RPM = ' + str(clino_rotation_rate), fontsize=30)
    #plt.savefig(figure_output_file, dpi = 300)
    plt.show()
    
    # PUT DPI 
    
# Execute main method, but only when directly invoked
if __name__ == "__main__":
    
    config = 'Swimming'
    
    output_file = f"D:\MPhys\ExpData\AdditionalData\Varying_configurations\Data\{config}\positions.csv"
    time_file = f"D:/MPhys/ExpData/AdditionalData/Varying_configurations/Data/{config}/time.csv"
    figure_output = f"D:\MPhys\ExpData\AdditionalData\Varying_configurations\Data\{config}\trajectory.png"
    
    main('test_config.txt', 'bacpos.csv', 'tumbles.csv', 'time.csv', 'swimming_direction.csv', 'traj.png')