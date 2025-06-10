[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 20
    ny = 10
    nz = 1
    xmin = 0
    ymin = -1
    zmin = -1
    xmax = 10
    ymax = 1
    zmax = 1
  []
  second_order = true
[]

[AuxVariables]
  [vf_x]
    family = LAGRANGE
    order = SECOND
  []
  [vf_y]
    family = LAGRANGE
    order = FIRST
  []
  [vf_z]
    family = LAGRANGE
    order = FIRST
  []
  [rho_f]
    family = LAGRANGE
    order = FIRST 
  []
[]

[AuxKernels]
  [evaluate_vf_x]
    type = ParsedAux
    variable = vf_x
    expression = '1-y^2'
    use_xyzt = true
  []
  [evaluate_vf_y]
    type = ParsedAux
    variable = vf_y
    expression = '0'
    use_xyzt = true
  []
  [evaluate_vf_z]
    type = ParsedAux
    variable = vf_z
    expression = '0'
    use_xyzt = true
  []
  [evaluate_rho_f]
    type = ParsedAux
    variable = rho_f
    expression = '1400'
    use_xyzt = true
  []
[]


[UserObjects]
  [particle_tracking_study]
    type = RepeatableRayStudy
    names = 'particle_1'
    start_points = '0.25 0.5000001 0'
    directions = '1 0 0'
  []
[]

[RayKernels]
  # [tracking]
  #   type = ParticleAdvectionKernel
  #   fluid_velocity = 'vf_x vf_y vf_z'
  #   dt = 0.1
  #   ray_refraction = false
  # []
  [particle_tracking]
    type = ParticleTrackingKernel
    particle_density = 2000
    particle_diameter = 1e-6
    initial_position = '0.25 0.5000001 0'
    initial_velocity = '0 0 0'
    fluid_velocity = 'vf_x vf_y vf_z'
    fluid_density = rho_f
    fluid_viscosity = '14'
    gravity = '0 -10 0'
    dt = 0.0001
    ray_refraction = false
  []
[]

[RayBCs]
  [kill_particle]
    type = KillRayBC
    boundary = 'bottom right top left front back'
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]






