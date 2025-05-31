[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 4
    ny = 3
    xmax = 4
    ymax = 3
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
    directions = '1 0 0'
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






