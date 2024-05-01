# Input file modified from RobPodgorney version
# - 2D instead of 3D with different resolution.
# - Celsius instead of Kelvin
# - a DirichletBC is used instead
# - Use of PorousFlowFullySaturated instead of PorousFlowUnsaturated, and the save_component_rate_in feature to record the change in kg of each species at each node for passing to the Geochem simulation
# - MultiApps and Transfers to transfer information between this simulation and the aquifer_geochemistry.i simulation
# injection_rate = -0.463 # kg/s; since 40t/day=0.463 kg/s,negative because injection as a source
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

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [f_H]
    initial_condition = -1.61189E-07
  []
  [f_Na]
    initial_condition = 0.000228573
  []
  [f_K]
    initial_condition = 6.48158E-07
  []
  [f_Ca]
    initial_condition = 0.000282584
  []
  [f_Mg]
    initial_condition = 1.44718E-07
  []
  [f_Al]
    initial_condition = 6.86074E-07
  []  
  [f_Fe]
    initial_condition = 1.00339E-09
  []  
  [f_SiO2]
    initial_condition = 1.19029E-05
  []

  [f_Cl]
    initial_condition = 0.000849695
  []
  [f_HCO3]
    initial_condition = 1.30872E-06
  []
  [porepressure]
    initial_condition = 7.9E6
  []
  [temperature]
    initial_condition = 40 # degC
    scaling = 1E-6 # fluid enthalpy is roughly 1E6
  []
[]

[BCs]
  [source_temperature]
    type = DirichletBC
    boundary = injection_node
    variable = temperature
    value = 30 # degC
  []
[]

[DiracKernels]

  [inject_H]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux = 9.16942E-07
    variable = f_H
  []
  [inject_Na]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux = 9.25785E-06
    variable = f_Na
  []
  [inject_K]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux = 4.63317E-07
    variable = f_K
  []
  [inject_Ca]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux = 4.63819E-07
    variable = f_Ca
  []
  [inject_Mg]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux = 4.62400E-07
    variable = f_Mg
  []
  [inject_Al]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux =4.63362E-07
    variable = f_Al
  []
  [inject_Fe]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux =4.62736E-07
    variable = f_Fe
  []
  [inject_SiO2]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux = 4.61689E-07
    variable = f_SiO2
  []
  [inject_Cl]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux = 2.46361E-05
    variable = f_Cl
  []
  [inject_HCO3]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux = 4.62893E-05
    variable = f_HCO3
  []
  [inject_H2O]
    type = PorousFlowPointSourceFromPostprocessor
    point = ' 0 0 0'
    mass_flux = 0.462879085
    variable = porepressure
  []
[]

[Postprocessors]
   [dt]
    type = TimestepSize
    execute_on = 'timestep_begin'
  []
[]

[FluidProperties]
  [the_simple_fluid]
    type = SimpleFluidProperties
    thermal_expansion = 2E-4
    bulk_modulus = 2E9
    viscosity = 1E-3
    density0 = 980
    cv = 4000.0
    cp = 4000.0
    porepressure_coefficient = 0
  []
[]

[PorousFlowFullySaturated]
  coupling_type = ThermoHydro
  porepressure = porepressure
  temperature = temperature
  mass_fraction_vars = 'f_H f_Na f_K f_Ca f_Mg f_SiO2 f_Al f_Cl f_Fe f_HCO3'
  save_component_rate_in = 'rate_H rate_Na rate_K rate_Ca rate_Mg rate_SiO2 rate_Al rate_Cl rate_Fe rate_HCO3 rate_H2O' # change in kg at every node / dt
  fp = the_simple_fluid
  temperature_unit = Celsius
[]

[AuxVariables]
  [rate_H]
  []
  [rate_Na]
  []
  [rate_K]
  []
  [rate_Ca]
  []
  [rate_Mg]
  []
  [rate_SiO2]
  []
  [rate_Al]
  []
  [rate_Cl]
  []
  [rate_Fe]
  []
  [rate_HCO3]
  []
  [rate_H2O]
  []
[]

[Materials]
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.2
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1.2E-13 0 0   0 1.2E-13 0   0 0 1.2E-13'
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.3 0 0  0 2.3 0  0 0 2.3'
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    density = 2900.0
    specific_heat_capacity = 950.0
  []
[]

[Preconditioning]
  active = typically_efficient
  [typically_efficient]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = ' hypre    boomeramg'
  []
  [strong]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      ilu           NONZERO                   2'
  []
  [probably_too_strong]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  #end_time = 7776000 # 90 days
  end_time = 120000 # 90 days
  # [TimeStepper]
  #   type = SolutionTimeAdaptiveDT
  #   dt = 500
  # []
  [TimeStepper]
    type = FunctionDT
    function = 'min(3E3, max(1E4, 0.2 * t))'
  []


[]

[Outputs]
  exodus = true
  csv = true
[]

[MultiApps]
  [react]
    type = TransientMultiApp
    input_files = aquifer_geochemistry.i
    clone_parent_mesh = true
    execute_on = 'timestep_end'
  []
[]
[Transfers]
  [changes_due_to_flow]
    type = MultiAppCopyTransfer
    source_variable = 'rate_H rate_Na rate_K rate_Ca rate_Mg rate_SiO2 rate_Al rate_Cl rate_Fe rate_HCO3 rate_H2O temperature'
    variable = 'pf_rate_H pf_rate_Na pf_rate_K pf_rate_Ca pf_rate_Mg pf_rate_SiO2 pf_rate_Al pf_rate_Cl pf_rate_Fe pf_rate_HCO3 pf_rate_H2O temperature'
    to_multi_app = react
  []
  [massfrac_from_geochem]
    type = MultiAppCopyTransfer
    source_variable = 'massfrac_H massfrac_Na massfrac_K massfrac_Ca massfrac_Mg massfrac_SiO2 massfrac_Al massfrac_Cl massfrac_Fe massfrac_HCO3'
    variable = 'f_H f_Na f_K f_Ca f_Mg f_SiO2 f_Al f_Cl f_Fe f_HCO3'
    from_multi_app = react
  []
[]
