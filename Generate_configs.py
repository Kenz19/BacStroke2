'''
This script generates all of the configuration & initalcondition files
needed for experiment 1 and stores them within this directory

!!!  need to add code to copy defualt files into directory instead of just
# manual copy paste
'''

# imports #####################################################################

# modules
import numpy as np
import os
import shutil

# external files
import Functions as f

# Generate folders ############################################################

# create storage directories for configuration and initialconditions files
directory_path_config = 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Gravity/Configs'
directory_path_conditions = 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Gravity/Initial_conditions'

default_config_path = 'C:/Users/kenzi/Documents/Masters/Summer/test_config.txt'
default_ic_path = 'C:/Users/kenzi/Documents/Masters/Summer/Object_initial_conditions/Ecoli.txt'

if not os.path.exists(directory_path_config): 
      
    # if the directory is not present  
    # then create it. 
    os.makedirs(directory_path_config)
    
    # file path of default config file
    default_path = os.path.abspath(default_config_path)
    
if not os.path.exists(directory_path_conditions): 
      
    # if the directory is not present then create it. 
    os.makedirs(directory_path_conditions)
    
# define g & vs values ########################################################

g = np.array([0, 0.01, 0.03, 0.1, 0.3, 1, 3])*9.8 # [m/s^2] gravities
vs = np.array([0, 0.1, 0.3, 1, 3, 10, 30])*25E-6 # [m/s] swimming speeds

# create configuration files ##################################################

# !!! need to add code to copy defualt files into directory instead of just
# manual copy paste

for i in range(len(g)):
    
    new_config_path = f'{directory_path_config}/g=' + str(round(g[i], 3))
    
    # create new config file by copying default config file 
    shutil.copyfile(default_config_path, new_config_path)
    
    # open new file and change g
    with open(new_config_path, 'r') as file:
        data = file.readlines()       
    
    # setting diffusion coefficient and rotation to 0 in config file
    data[25] = str(g[i]) + '\n' # \n to ensure new line added, gravity

    # and write everything back to config file
    with open(new_config_path, 'w') as file:
        file.writelines(data)
        
# create initial conditions files #############################################

for j in range(len(vs)):
    
    # name if new initial conditions file
    new_ic_path = f'{directory_path_conditions}/vs=' + str(round(vs[j], 8))
    
    # create copy of default ic file under new name
    shutil.copyfile(default_ic_path, new_ic_path)
    
    # change swimming speed in new ic file
    f.change_swimming_speed(new_ic_path, [vs[j]])
