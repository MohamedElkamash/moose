[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 3
    nz = 1
    xmax = 4
    ymax = 3
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
    expression_x = '2'
    expression_y = '1'
    expression_z = '0'
    use_xyzt = true
  []
[]

[UserObjects]
  [particle_tracking_study]
    type = RepeatableRayStudy
    names = 'particle_1'
    start_points = '0.5 0.5 0'
    directions = '2 1 0'
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






