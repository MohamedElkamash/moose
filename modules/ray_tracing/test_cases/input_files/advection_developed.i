[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 20
    nz = 1
    xmin = 0
    xmax = 0.1
    ymin = -0.01
    ymax = 0.01
    zmin = -0.05
    zmax = 0.05
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
[]

[AuxKernels]
  [evaluate_vf_x]
    type = ParsedAux
    variable = vf_x
    expression = '0.1*(1-y^2)'
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
[]


[UserObjects]
  [particle_tracking_study]
    type = RepeatableRayStudy
    names = 'particle_1'
    start_points = '0.0005 0.0005 0'
    directions = '1 0 0'
  []
[]

[RayKernels]
  [tracking]
    type = ParticleAdvectionKernel
    fluid_velocity = 'vf_x vf_y vf_z'
    dt = 0.1
    ray_refraction = false
  []
[]

[RayBCs]
  [kill_particle]
    type = KillRayBC
    boundary = 'bottom right top left'
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






