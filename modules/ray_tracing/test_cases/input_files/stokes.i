[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 1
    nz = 10
    xmin = 0
    ymin = -1
    zmin = 0
    xmax = 1
    ymax = 1
    zmax = 1
  []
  second_order = true
[]

[AuxVariables]
  [vf_x]
    family = LAGRANGE
    order = FIRST
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
    expression = '0'
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
    expression = '920'
    use_xyzt = true
  []
[]


[UserObjects]
  [particle_tracking_study]
    type = RepeatableRayStudy
    names = 'particle_1'
    start_points = '0.05 0 0.95'
    directions = '0 0 -1'
  []
[]

[RayKernels]
  [particle_tracking]
    type = ParticleTrackingKernel
    particle_density = 7800
    particle_diameter = 1e-3
    initial_position = '0.05 0 0.95'
    initial_velocity = '0 0 0'
    fluid_velocity = 'vf_x vf_y vf_z'
    fluid_density = rho_f
    fluid_viscosity = '0.081'
    gravity = '0 0 -10'
    dt = 0.002
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






