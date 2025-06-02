[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 1
    xmin = -1
    ymin = -1
    zmin = -1
    xmax = 1
    ymax = 1
    zmax = 1
  []
[]

[AuxVariables]
  [fluid_velocity]
    family = LAGRANGE_VEC
    order = FIRST
  []
[]

[AuxKernels]
  [parsed]
    type = ParsedVectorAux
    variable = fluid_velocity
    expression_x = '-y'
    expression_y = 'x'
    expression_z = '0'
    use_xyzt = true
  []
[]

[UserObjects]
  [particle_tracking_study]
    type = RepeatableRayStudy
    names = 'particle_1'
    start_points = '0.5 0 0'
    directions = '0 0.5 0'
  []
[]

[RayKernels]
  [tracking]
    type = TrackingKernel
    fluid_velocity = fluid_velocity
    dt = 0.1
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






