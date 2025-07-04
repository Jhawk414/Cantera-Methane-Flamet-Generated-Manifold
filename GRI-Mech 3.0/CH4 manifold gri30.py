"Course: AE 774 Combustion - Spring 2025"
"Author: Ryan Medlin"
print("This code develops a Methane-Air flamelet-generated manifold (FGM) using the GRI-Mech 3.0 Chemical Reaction Mechanism. Please ensure you have read the included user-guide.\n")

#==== Import packages =====
import cantera as ct
import numpy as np
from scipy import interpolate
from scipy.stats import beta
from scipy import integrate
from scipy.interpolate import RegularGridInterpolator

#==== Define Simple Parameters =====
T_f  = 300  # fuel stream temp
T_ox = 300 # ox stream temp
T_gas = 300 # ambient temp

manifold_range = [0.4, 0.5, 1, 3, 5, 7, 9, 10] # Range of Pressures to simulate each flamelet to build manifold.

# String w/ all 53 species (GRI 3.0 only)
spec_list = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "C", "CH", "CH2", "CH2(S)", "CH3", "CH4", "CO", "CO2", "HCO", "CH2O", "CH2OH", "CH3O", "CH3OH", "C2H", "C2H2", "C2H3", "C2H4", "C2H5", "C2H6", "HCCO", "CH2CO", "HCCOH", "N", "NH", "NH2", "NH3", "NNH", "NO", "NO2", "N2O", "HNO", "CN", "HCN", "H2CN", "HCNN", "HCNO", "HOCN", "HNCO", "NCO", "N2", "AR", "C3H7", "C3H8", "CH2CHO", "CH3CHO"]

for i in manifold_range:
    Pressure = i*ct.one_atm # Pressure to be used in the counterflow combustor
    print(f"A flamelet at P = {i} atmospheres is now being simulated. Please wait.")
    
    #==== Domain Details
    domain_width = 0.02
    gas = ct.Solution('gri30.yaml') # Instantiates a "gas" object
    gas.TP = T_gas, Pressure
    
    #==== Fuel/Ox Stream Details
    m_dot_f = 0.15 # fuel mass flow per area [kg/s*m^2]
    m_dot_ox = 0.75 # oxidizer mass flow per area [kg/s*m^2]
    
    fuel_comp = 'CH4:1' # Fuel stream composition on a molar basis
    ox_comp = 'O2:1, N2:3.76' # Oxidizer stream composition on a molar basis
    
    #==== Set flamelet Object's Parameters
    flamelet = ct.CounterflowDiffusionFlame(gas, width=domain_width)
    flamelet.transport_model = 'mixture-averaged'
    flamelet.set_refine_criteria(ratio=2.25, slope=0.125, curve=0.2, prune=0.04) 
    
    flamelet.fuel_inlet.mdot = m_dot_f      # Sets fuel stream mass flux
    flamelet.fuel_inlet.X = fuel_comp       # Sets fuel stream composition
    flamelet.fuel_inlet.T = T_f             # Sets fuel stream temp
    
    flamelet.oxidizer_inlet.mdot = m_dot_ox # Sets ox inlet mass flux
    flamelet.oxidizer_inlet.X = ox_comp     # Sets ox inlet composition
    flamelet.oxidizer_inlet.T = T_ox        # Sets ox inlet temp
        
    # Radiation
    flamelet.radiation_enabled = False;  flamelet.boundary_emissivities = 0.0, 0.0
    
    #==== Define other solving parameters====
    loglevel = 0
    
    #==== Solve ! ====
    flamelet.solve(loglevel,auto=True, refine_grid=True)
    
    #==== Extract f, PV, norm_PV, Temperature, and Y_k ====
    Mix_Fraction = flamelet.mixture_fraction('C') # Mixture fraction
    PV = flamelet.Y[gas.species_index('CO2')] + flamelet.Y[gas.species_index('H2O')] # Progress Variable = Y_CO2 + Y_H2O
    Temp = flamelet.T 
    density = flamelet.density
    
    #===== Interpolate Output on evenly-spaced Z =====
    new_Z = np.linspace(0.01, 0.999, 75) # Define evenly-spaced mixture fraction f (Z)
        
    # Approximate the properties as a function of the raw mixture fraction (returns a f'n to be called after)
    z_approx_PV = interpolate.interp1d(Mix_Fraction, PV)
    z_approx_T = interpolate.interp1d(Mix_Fraction, Temp)
    z_approx_rho = interpolate.interp1d(Mix_Fraction, density)

    # Interpolate (map) on equidistant Z (returns actual values, not a f'n)
    interp_PV = z_approx_PV(new_Z)
    interp_T = z_approx_T(new_Z)
    interp_rho = z_approx_rho(new_Z)

    # Grab all the PVs into a matrix for normalization
    #if i==min(manifold_range):
    #    group_all_PV = interp_PV
    #else:
    #    group_all_PV = np.column_stack((group_all_PV, interp_PV))

    # Group φ and corresponding PV for each flamelet to built large solution table
    T_C_group   = np.column_stack((interp_T, interp_PV)) # 100x2 array of T and C for the i^th solution
    rho_C_group = np.column_stack((interp_rho, interp_PV))  

    # >> Repeat above by looping through the 53 species' mass fractions.
    Yholder = {}
    z_Yapprox_holder = {}
    interp_Y_hold = {}
    Y_C_holder = {} # Instantiates dictionary for grouping the adjusted Y and PVs @ each flamelet
    
    for k in spec_list:
        Yholder['Y_' + str(k)] = flamelet.Y[gas.species_index(str(k))] # extracts all 53 mass fractions into dictionary
        z_Yapprox_holder['z_approx_' + str(k)] = interpolate.interp1d(Mix_Fraction, Yholder['Y_' + str(k)]) # interp1d finds the relation f'n
        interp_Y_hold['interp_Y_' + str(k)] = z_Yapprox_holder['z_approx_' + str(k)](new_Z) # Calculate Y_k @ new_Z
        Y_C_holder['Y_'+str(k)+'_C_group'] = np.column_stack((interp_Y_hold['interp_Y_'+str(k)] , interp_PV)) # groups arrays of Y and C for the i^th simulation
    
   
   # Generate solution table of [Z | φ C | φ C | φ C | φ C], where φ and C are the mixture property and PV, respectively.
    if i==min(manifold_range):
        ZCT_matrix = np.column_stack((new_Z, T_C_group)) # Adds the len(new_Z)x2 T-C array to the 100x1 Z array on the 1st flamelet.
        ZCrho_matrix = np.column_stack((new_Z, rho_C_group))
        # >> Repeat the above by looping through the mass fractions.
        ZCY_holder = {} # Instantiates dictionary
        for k in spec_list:
            ZCY_holder['ZCY_'+str(k)+'_soln'] = np.column_stack((new_Z, Y_C_holder['Y_'+str(k)+'_C_group']))
        
    else: # Adds the newest [ φ C ] solution array to the library
        ZCT_matrix   = np.column_stack((ZCT_matrix, T_C_group))
        ZCrho_matrix = np.column_stack((ZCrho_matrix,rho_C_group))
        # >> Repeat the above by looping through the mass fractions.
        for k in spec_list:
            ZCY_holder['ZCY_'+str(k)+'_soln'] = np.column_stack((ZCY_holder['ZCY_'+str(k)+'_soln'], Y_C_holder['Y_'+str(k)+'_C_group']))

print('\nAll flamelets have now been simulated.\n')

#====== Interpolate Output on evenly-spaced PV =====
# NOTE: At this point in the code, a table of [Z | φ C | φ C | φ C | φ C] should be generated. We now need to reshape the solution matrix at constant Z-values to use interp1d to find the φ - normalized C relation, then interpolate on a rectilinear C (note that Z is constant for each interpolation — i.e.,  solution interpolation will now occur on "slices" of Z).

new_C = np.linspace(0.01, 0.99, 25)

# Grab the 1D array solution for a "slice" of Z and reshape into 2D array capable of being used by interp1d. # of rows is # of flamelets simulated.
for zee in range(len(new_Z)):
    Soln_2D_T = ZCT_matrix[zee, 1:].reshape(len(manifold_range),2)
    Soln_2D_rho = ZCrho_matrix[zee, 1:].reshape(len(manifold_range),2)
    
# Normalize the PV (1st column) within each constant-Z "slice": PV_norm = (PV - PV_min) / (PV_max - PV_min).
    PV_T_norm = (Soln_2D_T[:,1] - min(Soln_2D_T[:,1]) )/(max(Soln_2D_T[:,1]) - min(Soln_2D_T[:,1]))
    PV_rho_norm = (Soln_2D_rho[:,1] - min(Soln_2D_rho[:,1]) )/(max(Soln_2D_rho[:,1]) - min(Soln_2D_rho[:,1]))

# Then, find φ's relation to the rectilinear grid of the normalized PV.
    c_approx_T   = interpolate.interp1d(PV_T_norm, Soln_2D_T[:,0], fill_value="extrapolate") # Note: PV is in 1st column and φ is in 0th column.
    c_approx_rho = interpolate.interp1d(PV_rho_norm, Soln_2D_rho[:,0], fill_value="extrapolate")

# Evaluate interpolated φ on rectilinear normalized PV.
    T_rectilinear = c_approx_T(new_C)
    rho_rectilinear = c_approx_rho(new_C)
    
    # >> Repeat φ - normalized PV interpolation by looping through the 53 species' mass fractions:
    Soln_2D_Y_hold = {}; PV_Y_norm_hold = {}; c_approx_Y_hold = {}; Y_rect_hold = {} # Instantiates dictionaries
    for k in spec_list:
        Soln_2D_Y_hold['2D_Y_'+str(k)] = ZCY_holder['ZCY_'+str(k)+'_soln'][zee,1:].reshape(len(manifold_range),2)
        PV_Y_norm_hold['PV_Y_'+str(k)] = (Soln_2D_Y_hold['2D_Y_'+str(k)][:,1] - min(Soln_2D_Y_hold['2D_Y_'+str(k)][:,1])) / (max(Soln_2D_Y_hold['2D_Y_'+str(k)][:,1]) - min(Soln_2D_Y_hold['2D_Y_'+str(k)][:,1])) # normalized PV within each Z "slice"
        c_approx_Y_hold['Y_'+str(k)+'_appx'] = interpolate.interp1d(PV_Y_norm_hold['PV_Y_'+str(k)], Soln_2D_Y_hold['2D_Y_'+str(k)][:,0], fill_value="extrapolate")
        Y_rect_hold['Y_'+str(k)+'_rect'] = c_approx_Y_hold['Y_'+str(k)+'_appx'](new_C)
        
    # ==== Form lookup Table ====
    # NOTE: Now that output values are interpolated w.r.t. both Z and a normalized C, we can now form THE lookup table that'll be used for the unfiltered 2D interpolation. Note that only the unfiltered φ will reside in the array, not the Z or C "header" since RegularGridInterpolator doesn't require those.
    
    if zee==0:
         # On first Z "slice", place new_C rectilinear output into the 2D lookup data array.
         T_Lookup   = T_rectilinear 
         rho_Lookup = rho_rectilinear
         
         # >> Repeat above by looping through the 53 species' mass fractions.
         Y_lookup_hold = {} #instantiates dictionary
         for k in spec_list:
             Y_lookup_hold['Y_'+str(k)+'_lookup'] = Y_rect_hold['Y_'+str(k)+'_rect']
    else:
        # Adds newest Interpolated φ to the 2D lookup data array.
        T_Lookup   = np.column_stack((T_Lookup, T_rectilinear))
        rho_Lookup = np.column_stack((rho_Lookup, rho_rectilinear))
        # >> Repeat above by looping through the 53 species' mass fractions.
        for k in spec_list:
            Y_lookup_hold['Y_'+str(k)+'_lookup'] = np.column_stack((Y_lookup_hold['Y_'+str(k)+'_lookup'], Y_rect_hold['Y_'+str(k)+'_rect']))

#==== Post-process the lookup tables ====
# Transpose to have constant Z "slices" as rows
T_Lookup   = np.transpose(T_Lookup.astype(int)) # Converts T to integers to read cleaner
rho_Lookup = np.transpose(rho_Lookup)
# >> Repeat above by looping through the 53 species' mass fractions.
for k in spec_list:
    Y_lookup_hold['Y_'+str(k)+'_lookup'] = np.transpose(Y_lookup_hold['Y_'+str(k)+'_lookup'])

T_reg_interp = RegularGridInterpolator((new_Z, new_C), T_Lookup) # 2D rectilinear interpolating with (new_Z, new_C) as the "headers" and T_Lookup as the 2D data
rho_reg_interp = RegularGridInterpolator((new_Z, new_C), rho_Lookup)

#==== Prompt user for input ====
print("Please provide the following inputs, noting the bounds listed:\n\n• Mixture fraction Z (0.001 ≤ Z ≤ 0.999)\n• Normalized progress variable C (0.001 ≤ C ≤ 0.999)\n• Mean mixture fraction variance Z\" (0.001 ≤ Z\" ≤ 0.1)\n\nNote: Do not separate inputs with a comma, rather only with a single space.\n\n")
Z_usr_input, C_usr_input , Zvar_usr_input = input('Z  C  Z": ').split() # prompts user to input Z & C and splits them into respective variable

# Evaluate unfiltered φ
T_unfiltered = T_reg_interp((float(Z_usr_input), float(C_usr_input)))
rho_unfiltered = rho_reg_interp((float(Z_usr_input), float(C_usr_input)))

# >> Repeat above by looping through the 53 species' mass fractions
Y_reg_interp_hold = {}; Y_unfiltered_hold = {} #instantiates dictionary
for k in spec_list:
    Y_reg_interp_hold['Y_'+str(k)+'_reg_interp'] = RegularGridInterpolator(((new_Z, new_C)), Y_lookup_hold['Y_'+str(k)+'_lookup'])
    Y_reg_grab = Y_reg_interp_hold['Y_'+str(k)+'_reg_interp'] # Grabs the scipy object out of the dictionary since RegularGridInterpolator can't index strings
    Y_unfiltered_hold['Y_'+str(k)+'_unfiltered'] = Y_reg_grab((float(Z_usr_input), float(C_usr_input)))


#===== Apply β-PDF to unfiltered φ =====
Zvar_usr_input = float(Zvar_usr_input) # convert user input to float before use

gma = ( (new_Z*(1-new_Z)) / Zvar_usr_input**2) - 1 # Gamma Function

# Shape Parameters for the β-PDF function
a = new_Z*gma
b = (1-new_Z)*gma

bPDF = beta.pdf(new_Z, a, b) # β-PDF function
bPDF_no_NaN = bPDF[~np.isnan(bPDF)] # Removes NaN entries to enable integration

# Integrate using copmosite Simpson's rule
P_int = integrate.simpson(bPDF_no_NaN)

# Normalize the β-PDF function by the integral and sum
bPDF_norm = sum(bPDF_no_NaN / P_int)

# Apply to unfiltered φ
T_filtered = T_unfiltered*bPDF_norm
print("\nThe filtered temperature is approximately %.1f K." % T_filtered)

rho_filtered = rho_unfiltered*bPDF_norm
print("The filtered density is approximately %.3f kg/m³." % rho_filtered)


# >> Loop through the mass fractions to filter and print()
Y_filtered_hold = {} #instantiates a dictionary

print("\nThe following is a table of the filtered mass fractions.\n\nSpecies  |  Filtered Yₖ")

for k in spec_list:
    Y_filtered_hold['Y_'+str(k)+'_filt'] = Y_unfiltered_hold['Y_'+str(k)+'_unfiltered']*bPDF_norm
    print(f"{k}     {Y_filtered_hold['Y_'+str(k)+'_filt']:.4E}")
