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
[]

[AuxKernels]
  [evaluate_vf_x]
    type = ParsedAux
    variable = vf_x
    expression = '-y'
    use_xyzt = true
  []
  [evaluate_vf_y]
    type = ParsedAux
    variable = vf_y
    expression = 'x'
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
    start_points = '0.5 0 0'
    directions = '0 0.5 0'
  []
[]

[RayKernels]
  [tracking]
    type = TrackingKernel
    fluid_velocity = 'vf_x vf_y vf_z'
    dt = 0.1
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