# 1 Geochemistry with kinetically-controlled mineral changesin the aquifer can run in standalone fashion. 
#   couples with porous_flow.i simulation using MultiApps.
# 2 This file receives pf_rate_H pf_rate_Na pf_rate_K... and temperature as AuxVariables from porous_flow.i
# 3 The pf_rate [kg/s] from porousflow.i at each node, but the geochemistry needs rates-of-changes[moles/s]. 
#   Geochemistry.i considers just 1 litre of aqueous solution at every node, the nodal_void_volume is used to convert pf_rate_* into rate_*_per_1l [ mol/s/1_litre_of_aqueous_solution].
# 4  This file sends massfrac_H massfrac_Na ... to porous_flow.i.  
#    These are computed from the corresponding transported_* quantities.

# define kinetic rates

[UserObjects]

  [rate_Diopside]
    type = GeochemistryKineticRate
    kinetic_species_name = Diopside
    intrinsic_rate_constant = 1E-17
    multiply_by_mass = true
    area_quantity = 34
    activation_energy = 40.57E3
    one_over_T0 = 0.003354 # this is room tem T=298.15K, 1/T is constant
  []

  [rate_Hedenbergite]
    type = GeochemistryKineticRate
    kinetic_species_name = Hedenbergite
    intrinsic_rate_constant = 1E-17
    multiply_by_mass = true
    area_quantity = 34
    activation_energy = 40.57E3
    one_over_T0 = 0.003354 
  []
  [rate_Albite]
    type = GeochemistryKineticRate
    kinetic_species_name = Albite
    intrinsic_rate_constant = 1E-17
    multiply_by_mass = true
    area_quantity = 225
    activation_energy = 66.49E3
    one_over_T0 = 0.003354
  []
  [rate_Anorthite]
    type = GeochemistryKineticRate
    kinetic_species_name = Anorthite
    intrinsic_rate_constant = 1.0E-13
    multiply_by_mass = true
    area_quantity = 225
    activation_energy = 66.49E3
    one_over_T0 = 0.003354
  []
  [rate_Antigorite]
    type = GeochemistryKineticRate
    kinetic_species_name = Antigorite
    intrinsic_rate_constant = 1.0E-17
    multiply_by_mass = true
    area_quantity = 126000
    activation_energy = 56.58E3
    one_over_T0 = 0.003354
  []
  [rate_Greenalite]
    type = GeochemistryKineticRate
    kinetic_species_name = Greenalite
    intrinsic_rate_constant = 1.0E-17
    multiply_by_mass = true
    area_quantity = 126000
    activation_energy = 56.58E3
    one_over_T0 = 0.003354
  []
  [rate_Forsterite]
    type = GeochemistryKineticRate
    kinetic_species_name = Forsterite
    intrinsic_rate_constant = 1.0E-17
    multiply_by_mass = true
    area_quantity = 900
    activation_energy = 78.96E3
    one_over_T0 = 0.003354
  []
  [rate_Fayalite]
    type = GeochemistryKineticRate
    kinetic_species_name = Fayalite
    intrinsic_rate_constant = 1.0E-17
    multiply_by_mass = true
    area_quantity = 360
    activation_energy = 94.35E3
    one_over_T0 = 0.003354
  []
  [rate_K-feldspar]
    type = GeochemistryKineticRate
    kinetic_species_name = K-feldspar
    intrinsic_rate_constant = 1E-17
    multiply_by_mass = true
    area_quantity = 190
    activation_energy = 37.97E3
    one_over_T0 = 0.003354
  []
  [rate_Calcite]
    type = GeochemistryKineticRate
    kinetic_species_name = Calcite
    intrinsic_rate_constant = 1.0E-10
    multiply_by_mass = true
    area_quantity = 370
    activation_energy = 23.5E3
    one_over_T0 = 0.003354
  []
  [rate_Magnesite]
    type = GeochemistryKineticRate
    kinetic_species_name = Magnesite
    intrinsic_rate_constant = 1.0E-10
    multiply_by_mass = true
    area_quantity = 662
    activation_energy = 23.5E3
    one_over_T0 = 0.003354
  []
  [rate_Siderite]
    type = GeochemistryKineticRate
    kinetic_species_name = Siderite
    intrinsic_rate_constant = 1.0E-10
    multiply_by_mass = true
    area_quantity = 1050
    activation_energy = 62.73E3
    one_over_T0 = 0.003354
  []
  [rate_Gibbsite]
    type = GeochemistryKineticRate
    kinetic_species_name = Gibbsite
    intrinsic_rate_constant = 1.0E-10
    multiply_by_mass = true
    area_quantity = 190000
    activation_energy = 61.18E3
    one_over_T0 = 0.003354
  []
  [rate_Quartz]
    type = GeochemistryKineticRate
    kinetic_species_name = Quartz
    intrinsic_rate_constant = 1E-18
    multiply_by_mass = true
    area_quantity = 225
    activation_energy = 78.96E3
    one_over_T0 = 0.003354
  []
  [definition]
    type = GeochemicalModelDefinition
    database_file = '../../../../geochemistry/database/moose_geochemdb.json'
    basis_species = 'H2O H+ Na+ K+ Ca++ Mg++ Al+++ Fe++ SiO2(aq) Cl- HCO3-'
    remove_all_extrapolated_secondary_species = true    
    kinetic_minerals = 'Diopside Hedenbergite Albite Anorthite Antigorite Greenalite Forsterite Fayalite K-feldspar Calcite Magnesite Siderite Gibbsite Quartz'
    kinetic_rate_descriptions='rate_Diopside rate_Hedenbergite rate_Albite rate_Anorthite rate_Antigorite rate_Greenalite rate_Forsterite rate_Fayalite rate_K-feldspar rate_Calcite rate_Magnesite rate_Siderite rate_Gibbsite rate_Quartz'
  []

  [nodal_void_volume_uo]
    type = NodalVoidVolume
    porosity = porosity
    execute_on = 'initial timestep_end' # "initial" means this is evaluated properly for the first timestep
  []
[]

[SpatialReactionSolver]
  model_definition = definition
  geochemistry_reactor_name = reactor
  charge_balance_species = 'Cl-'
  constraint_species ='H2O             H+       Na+         K+         Ca++        Mg++       Al+++     Fe++     SiO2(aq)     Cl-         HCO3-'
# Following numbers are from aquifer_eqm_out.csv_ use molar* values
  #constraint_value = '1.0 1.78E-09 0.009956037 1.6600E-05 0.006715892  5.8526E-06 1.97E-19 1.61E-08 0.000168617  0.009717499 1.59E-05'
  constraint_value = '0.40489896507677 1.78E-09 0.009956037 1.6600E-05 0.006715892  5.8526E-06 1.97E-19 1.61E-08 0.000168617  0.009717499 1.59E-05'

  constraint_meaning = 'kg_solvent_water free_concentration       free_concentration    free_concentration      free_concentration     free_concentration       free_concentration      free_concentration      free_concentration  bulk_composition        free_concentration'
  constraint_unit = '   kg               molal               molal            molal              molal             molal               molal              molal            molal                  moles              molal'
  initial_temperature = 40
  temperature = temperature # this from porous flow.i
  kinetic_species_name =          'Diopside    Hedenbergite Albite      Anorthite   Antigorite  Greenalite  Forsterite  Fayalite    K-feldspar  Calcite     Magnesite   Siderite    Gibbsite    Quartz'
  kinetic_species_initial_value = '8.399152671	4.644216691	5.150894374	7.474451831	0.002858613	0.020124774	7.243662937	4.273550334	12.01846239	0.135376618	0.178456706	0.174641984	0.220380818	0.15646514'
  kinetic_species_unit = '         moles        moles        moles       moles       moles       moles       moles      moles       moles        moles      moles        moles       moles       moles'
  evaluate_kinetic_rates_always = true # otherwise will easily "run out" of dissolving species
  source_species_names = 'H2O             H+            Na+            K+             Ca++          Mg++           Al+++          Fe++          SiO2(aq)          Cl-             HCO3-'
  source_species_rates = 'rate_H2O_per_1l rate_H_per_1l rate_Na_per_1l rate_K_per_1l rate_Ca_per_1l rate_Mg_per_1l rate_Al_per_1l rate_Fe_per_1l rate_SiO2_per_1l  rate_Cl_per_1l  rate_HCO3_per_1l'
  # above rate come from porousflow.i after unit conversion.  
  ramp_max_ionic_strength_initial = 0 # max_ionic_strength in such a simple problem does not need ramping
  execute_console_output_on = ''
  add_aux_molal = false # save some memory and reduce variables in output exodus
  add_aux_mg_per_kg = false # save some memory and reduce variables in output exodus
  add_aux_free_mg = false # save some memory and reduce variables in output exodus
  add_aux_activity = false # save some memory and reduce variables in output exodus
  add_aux_bulk_moles = false # save some memory and reduce variables in output exodus
  adaptive_timestepping = true
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 15
    ny = 10
    xmin = -100
    xmax = 200
    ymin = -100
    ymax = 100
  []
  [injection_node]
    input = gen
    type = ExtraNodesetGenerator
    new_boundary = injection_node
    coord = '0 0 0'
  []
[]

[Executioner]
  type = Transient
  [TimeStepper]
    type = FunctionDT
    function = 'max(1E3, 0.3 * t)'
  []
  end_time = 4E8
[]

[AuxVariables]
  [temperature]
    initial_condition = 40.0
  []
  [porosity]
    initial_condition = 0.2
  []    
  [nodal_void_volume]
  []
  [free_cm3_Kfeldspar] # necessary because of the minus sign in K-feldspar which does not parse correctly in the porosity AuxKernel
  [] 
  [pf_rate_H] # change in H mass (kg/s) at each node provided by the porous-flow simulatio; come from porous_flow.i
  []
  [pf_rate_Na]# variables received from Porous_flow.i
  []
  [pf_rate_K]
  []
  [pf_rate_Ca]
  []
  [pf_rate_Mg]
  []
  [pf_rate_Al]
  []
  [pf_rate_Fe]
  []
  [pf_rate_SiO2]
  []
  [pf_rate_Cl]
  []
  [pf_rate_HCO3]
  []
  [pf_rate_H2O] # change in H2O mass (kg/s) at each node provided by the porous-flow simulation
  []

  [rate_H_per_1l] # variables received from Porous_flow.i and convert it into (mole/s/L) for gechemistry reaction use.
  []
  [rate_Na_per_1l]
  []
  [rate_K_per_1l]
  []
  [rate_Ca_per_1l]
  []
  [rate_Mg_per_1l]
  []
  [rate_Al_per_1l]
  [] 
  [rate_Fe_per_1l]
  [] 
  [rate_SiO2_per_1l]
  [] 
  [rate_Cl_per_1l]
  []
   [rate_HCO3_per_1l]
  []
  [rate_H2O_per_1l] # variables received from Porous_flow.i and convert it into (mole/s/L) for gechemistry reaction use
  []

  [transported_H] # variables (mole/s/L) after chemical reaction and will feedback to porous_flow 
  []
  [transported_Na]
  []
  [transported_K]
  []
  [transported_Ca]
  []
  [transported_Mg]
  []
  [transported_Al]
  []  
  [transported_Fe]
  []   
  [transported_SiO2]
  []
  [transported_Cl]
  []
  [transported_HCO3]
  []
  [transported_H2O]
  []# variables (mole/s/L) after chemical reaction and will feedback to porous_flow 


  [transported_mass] # variables (mole/s/L) change unit to kg/s for porous_flow use
  []
  [massfrac_H]
  []
  [massfrac_Na]
  []
  [massfrac_K]
  []
  [massfrac_Ca]
  []
  [massfrac_Mg]
  []
  [massfrac_Al]
  [] 
  [massfrac_Fe]
  []
  [massfrac_SiO2]
  []
  [massfrac_Cl]
  []
  [massfrac_HCO3]
  []
  [massfrac_H2O]
  []
[]

[AuxKernels]
  [free_cm3_Kfeldspar]
    type = GeochemistryQuantityAux
    variable = free_cm3_Kfeldspar
    species = 'K-feldspar'
    quantity = free_cm3
    execute_on = 'timestep_begin timestep_end'
  []
  [porosity_auxk]
    type = ParsedAux
    coupled_variables =             'free_cm3_Diopside free_cm3_Hedenbergite free_cm3_Albite free_cm3_Anorthite free_cm3_Antigorite free_cm3_Greenalite free_cm3_Forsterite free_cm3_Fayalite free_cm3_Kfeldspar free_cm3_Calcite free_cm3_Magnesite free_cm3_Siderite free_cm3_Gibbsite free_cm3_Quartz'
    expression = '1000.0 / (1000.0 + free_cm3_Diopside + free_cm3_Hedenbergite + free_cm3_Albite + free_cm3_Anorthite + free_cm3_Antigorite + free_cm3_Greenalite + free_cm3_Forsterite + free_cm3_Fayalite + free_cm3_Kfeldspar + free_cm3_Calcite + free_cm3_Magnesite + free_cm3_Siderite + free_cm3_Gibbsite + free_cm3_Quartz )'
    variable = porosity
    execute_on = 'timestep_end'
  []
  [nodal_void_volume_auxk]
    type = NodalVoidVolumeAux
    variable = nodal_void_volume
    nodal_void_volume_uo = nodal_void_volume_uo
    execute_on = 'initial timestep_end' # "initial" to ensure it is properly evaluated for the first timestep
  []
  [rate_H_per_1l_auxk] # unit is mole/s/L
    type = ParsedAux
    coupled_variables = 'pf_rate_H nodal_void_volume'
    variable = rate_H_per_1l
    expression = 'pf_rate_H / 1.0079 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Na_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Na nodal_void_volume'
    variable = rate_Na_per_1l
    expression = 'pf_rate_Na / 22.9898 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_K_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_K nodal_void_volume'
    variable = rate_K_per_1l
    expression = 'pf_rate_K / 39.0983 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Ca_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Ca nodal_void_volume'
    variable = rate_Ca_per_1l
    expression = 'pf_rate_Ca / 40.08 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Mg_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Mg nodal_void_volume'
    variable = rate_Mg_per_1l
    expression = 'pf_rate_Mg / 24.305 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Al_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Al nodal_void_volume'
    variable = rate_Al_per_1l
    expression = 'pf_rate_Al / 26.9815 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Fe_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Fe nodal_void_volume'
    variable = rate_Al_per_1l
    expression = 'pf_rate_Fe / 55.847 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_SiO2_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_SiO2 nodal_void_volume'
    variable = rate_SiO2_per_1l
    expression = 'pf_rate_SiO2 / 60.0843 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Cl_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Cl nodal_void_volume'
    variable = rate_Cl_per_1l
    expression = 'pf_rate_Cl / 35.453 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_HCO3_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_HCO3 nodal_void_volume'
    variable = rate_HCO3_per_1l
    expression = 'pf_rate_HCO3 / 61.0171 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_H2O_per_1l_auxk] # unit is mole/s/L
    type = ParsedAux
    coupled_variables = 'pf_rate_H2O nodal_void_volume'
    variable = rate_H2O_per_1l
    expression = 'pf_rate_H2O / 18.01801802 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [transported_H_auxk] # unit is mole/s/L
    type = GeochemistryQuantityAux
    variable = transported_H
    species = 'H+'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Na_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Na
    species = 'Na+'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_K_auxk]
    type = GeochemistryQuantityAux
    variable = transported_K
    species = 'K+'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Ca_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Ca
    species = 'Ca++'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Mg_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Mg
    species = 'Mg++'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Al_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Al
    species = 'Al+++'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Fe_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Fe
    species = 'Fe++'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []

  [transported_SiO2_auxk]
    type = GeochemistryQuantityAux
    variable = transported_SiO2
    species = 'SiO2(aq)'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Cl_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Cl
    species = 'Cl-'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_HCO3_auxk]
    type = GeochemistryQuantityAux
    variable = transported_HCO3
    species = 'HCO3-'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_H2O_auxk] # unit is mole/s/L
    type = GeochemistryQuantityAux
    variable = transported_H2O
    species = 'H2O'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_mass_auxk] # unit is g/s/L
    type = ParsedAux
    coupled_variables = ' transported_H transported_Na transported_K transported_Ca transported_Mg transported_Al transported_Fe transported_SiO2  transported_Cl transported_HCO3 transported_H2O'
    variable = transported_mass
    expression = 'transported_H * 1.0079 + transported_Cl * 35.453 + transported_Fe * 55.847 + transported_HCO3 * 61.0171 + transported_SiO2 * 60.0843 + transported_Al * 26.9815 + transported_Ca * 40.08 + transported_Mg * 24.305 + transported_K * 39.0983 + transported_Na * 22.9898 + transported_H2O * 18.01801802'
    execute_on = 'timestep_end' 
  []
  [massfrac_H_auxk] #fraction dimentionless
    type = ParsedAux
    coupled_variables = 'transported_H transported_mass'
    variable = massfrac_H
    expression = 'transported_H * 1.0079 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Na_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Na transported_mass'
    variable = massfrac_Na
    expression = 'transported_Na * 22.9898 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_K_auxk]
    type = ParsedAux
    coupled_variables = 'transported_K transported_mass'
    variable = massfrac_K
    expression = 'transported_K * 39.0983 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Ca_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Ca transported_mass'
    variable = massfrac_Ca
    expression = 'transported_Ca * 40.08 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Mg_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Mg transported_mass'
    variable = massfrac_Mg
    expression = 'transported_Mg * 24.305 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Al_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Al transported_mass'
    variable = massfrac_Al
    expression = 'transported_Al * 26.9815 / transported_mass'
    execute_on = 'timestep_end'
  []

  [massfrac_Fe_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Fe transported_mass'
    variable = massfrac_Fe
    expression = 'transported_Fe * 55.847 / transported_mass'
    execute_on = 'timestep_end'
  []

  [massfrac_SiO2_auxk]
    type = ParsedAux
    coupled_variables = 'transported_SiO2 transported_mass'
    variable = massfrac_SiO2
    expression = 'transported_SiO2 * 60.0843 / transported_mass'
    execute_on = 'timestep_end'
  []

  [massfrac_Cl_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Cl transported_mass'
    variable = massfrac_Cl
    expression = 'transported_Cl * 35.453 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_HCO3_auxk]
    type = ParsedAux
    coupled_variables = 'transported_HCO3 transported_mass'
    variable = massfrac_HCO3
    expression = 'transported_HCO3 * 61.0171 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_H2O_auxk]
    type = ParsedAux
    coupled_variables = 'transported_H2O transported_mass'
    variable = massfrac_H2O
    expression = 'transported_H2O * 18.01801802 / transported_mass'
    execute_on = 'timestep_end'
  []
[]

[GlobalParams]
  point = '0 0 0'
  reactor = reactor
[]

[Postprocessors]
  [temperature]
    type = PointValue
    variable = 'solution_temperature'
  []
  [porosity]
    type = PointValue
    variable = porosity
  []
  [solution_temperature]
    type = PointValue
    variable = solution_temperature
  []
  [massfrac_H]
    type = PointValue
    variable = massfrac_H
  []
  [massfrac_Na]
    type = PointValue
    variable = massfrac_Na
  []
  [massfrac_K]
    type = PointValue
    variable = massfrac_K
  []
  [massfrac_Ca]
    type = PointValue
    variable = massfrac_Ca
  []
  [massfrac_Mg]
    type = PointValue
    variable = massfrac_Mg
  []
  [massfrac_Al]
    type = PointValue
    variable = massfrac_Al
  []
  [massfrac_Fe]
    type = PointValue
    variable = massfrac_Fe
  []
  
  [massfrac_SiO2]
    type = PointValue
    variable = massfrac_SiO2
  []
  [massfrac_Cl]
    type = PointValue
    variable = massfrac_Cl
  []
  [massfrac_HCO3]
    type = PointValue
    variable = massfrac_HCO3
  []
  [massfrac_H2O]
    type = PointValue
    variable = massfrac_H2O
  []

  [kg_solvent_water]
    type = PointValue
    variable = kg_solvent_H2O
  []
  [free_cm3_Diopside]
    type = PointValue
    variable = free_cm3_Diopside  
  []
  [free_cm3_Hedenbergite]
    type = PointValue
    variable = free_cm3_Hedenbergite  
  []
  [free_cm3_Albite]
    type = PointValue
    variable = free_cm3_Albite  
  []
  [free_cm3_Anorthite]
    type = PointValue
    variable = free_cm3_Anorthite  
  []
  [free_cm3_Antigorite]
    type = PointValue
    variable = free_cm3_Antigorite  
  []
  [free_cm3_Greenalite]
    type = PointValue
    variable = free_cm3_Greenalite  
  []
  [free_cm3_Forsterite]
    type = PointValue
    variable = free_cm3_Forsterite  
  []         
  [free_cm3_Fayalite]
    type = PointValue
    variable = free_cm3_Fayalite  
  []      
  [free_cm3_Kfeldspar]
    type = PointValue
    variable = 'free_cm3_K-feldspar'
  [] 
  
  [free_cm3_Calcite]
    type = PointValue
    variable = free_cm3_Calcite  
  []       
  [free_cm3_Magnesite]
    type = PointValue
    variable = free_cm3_Magnesite  
  [] 
  [free_cm3_Siderite]
    type = PointValue
    variable = free_cm3_Siderite  
  [] 
  [free_cm3_Gibbsite ]
    type = PointValue
    variable = free_cm3_Gibbsite   
  []        
  [free_cm3_Quartz ]
    type = PointValue
    variable = free_cm3_Quartz   
  [] 
  [cm3_mineral]
    type = LinearCombinationPostprocessor
    pp_names =  'free_cm3_Diopside free_cm3_Hedenbergite free_cm3_Albite free_cm3_Anorthite free_cm3_Antigorite free_cm3_Greenalite free_cm3_Forsterite free_cm3_Fayalite free_cm3_Kfeldspar free_cm3_Calcite free_cm3_Magnesite free_cm3_Siderite free_cm3_Gibbsite free_cm3_Quartz'
    pp_coefs = '1 1 1 1 1 1 1 1 1 1 1 1 1 1'
  []

  [pH]
    type = PointValue
    variable = 'pH'
  []
[]

[Outputs]
  # [exo]
  #   type = Exodus
  #   execute_on = final
  # []
  exodus=true
  csv = true
[]
