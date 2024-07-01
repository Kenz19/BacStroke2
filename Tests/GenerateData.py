# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 14:38:20 2024

Generate data

"""

import os
import BacStroke as bac
import Functions as f

# configs
config_path = 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Tests/Configs'
config_list = os.listdir(config_path)
config_count = len(config_list)

# initial condition files
conditions_path = 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Tests/Initial_conditions'
conditions_list = os.listdir(conditions_path)
conditions_count = len(conditions_list)

for i in range(config_count):
    
    current_config = f'{config_path}/{config_list[i]}'
    
    # read config file and get clinostat dimensions
    constants = f.read_config(current_config)
    R = float(constants[4])
    r = float(constants[5])
    H = float(constants[6])
    #print(R, r, H)
    
    for j in range(conditions_count):
        
        current_condition = f'{conditions_path}/{conditions_list[j]}'
        print(current_condition)
        
        # change the conditions file in config to correct file
        with open(current_config, 'r') as file:
            data = file.readlines()
            
        # replace name with new initial conditions file
        data[1] = current_condition + '\n' # \n to ensure new line added, gravity

        # and write everything back to config file
        with open(current_config, 'w') as file:
            file.writelines(data)
        
        # create output folder for config and conditions combo
        save_folder = f'C:/Users/kenzi/Documents/Masters/Summer/Studies/Tests/RawData/{config_list[i]},{conditions_list[j]}'
        
        # create save folder if it doesn't exist
        if not os.path.isdir(save_folder): 
            os.makedirs(save_folder)
            
        for k in range(5):
            
            run_output_folder = f'{save_folder}/run{k+1}'
        
            if not os.path.isdir(run_output_folder): 
                os.makedirs(run_output_folder)
            
            output_file = f'{run_output_folder}/output.csv'
            
            # if output.csv hasnt been produced, run bacstroke
            if not os.path.isfile(f'{run_output_folder}/output.csv'): # check if folder exists
                    
                    # change starting coordinates of initial conditions file
                    f.change_starting_coords(current_condition, f.sample_from_hollow_cylinder(r, R, H))
                        
                    # run backstroke
                    bac.main(current_config, output_file)
                    
            

