'''
Calculate the normalized thickness Delta = (R_epi - R_endo)/R_endo for PBPVR
model given chamber volume V_endo (from HeartDeformNets) and thickness.
'''

import numpy as np
import pyvista as pv
import pandas as pd

def calc_Delta(V_endo, thickness):
    '''
    Calculate the normalized thickness Delta = (R_epi - R_endo)/R_endo for PBPVR
    model given chamber volume V_endo (from HeartDeformNets) and thickness.
    Input:
        V_endo: endocardial volume (mL=cm^3)
        thickness: wall thickness (cm)
    '''
    R_endo = (3*V_endo/(4*np.pi))**(1/3) # endocardial radius (cm)
    R_epi = R_endo + thickness           # epicardial radius (cm)
    return (R_epi - R_endo)/R_endo

# Get cardiac chamber volumes
volume_data_file = '/Users/aaronbrown/Library/CloudStorage/GoogleDrive-abrown97@stanford.edu/My Drive/Stanford/Marsden Lab/Papers/Patient-specific coupled BiV Paper/Sims/mesh_and_fibers/volume/whole_heart/volume.txt'
column_names = ['RR_interval', 'V_LA', 'V_LV', 'V_RA', 'V_RV']
data = pd.read_csv(volume_data_file, delim_whitespace=True, nrows=10, comment='#', names=column_names, header=None)

V_LV_70 = data[data['RR_interval'] == 70]['V_LV'].values[0]
V_RV_70 = data[data['RR_interval'] == 70]['V_RV'].values[0]
V_LA_20 = data[data['RR_interval'] == 20]['V_LA'].values[0]
V_RA_20 = data[data['RR_interval'] == 20]['V_RA'].values[0]


# Define cardiac wall thickness
thickness_LV = 1.1 # cm; average thickness of LV free wall measured from CT
thickness_RV = 0.35 # cm; assumed from Strocchi2023
thickness_LA = 0.2 # cm; assumed from Strocchi2023
thickness_RA = 0.2 # cm; assumed from Strocchi2023

# Calculate normalized thickness
D_LV = calc_Delta(V_LV_70, thickness_LV)
D_RV = calc_Delta(V_RV_70, thickness_RV)
D_LA = calc_Delta(V_LA_20, thickness_LA)
D_RA = calc_Delta(V_RA_20, thickness_RA)

print(f'V_LV_70: {V_LV_70:.2f} cm^3,\tt_LV: {thickness_LV:.2f} cm \t-> D_LV: {D_LV:.2f}')
print(f'V_RV_70: {V_RV_70:.2f} cm^3,\tt_RV: {thickness_RV:.2f} cm \t-> D_RV: {D_RV:.2f}')
print(f'V_LA_20: {V_LA_20:.2f} cm^3,\tt_LA: {thickness_LA:.2f} cm \t-> D_LA: {D_LA:.2f}')
print(f'V_RA_20: {V_RA_20:.2f} cm^3,\tt_RA: {thickness_RA:.2f} cm \t-> D_RA: {D_RA:.2f}')