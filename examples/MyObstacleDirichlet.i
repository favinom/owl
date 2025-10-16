[Mesh]
  [./gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 20
    #parallel_type = replicated
  []

  [./nodeset]
    type = BoundingBoxNodeSetGenerator
    input = gmg
    new_boundary = middle_node
    top_right = '0.6 1.01 0'
    bottom_left = '0.4 0.99 0'
  []
[]


[Variables]
[disp_x] order = FIRST []
[disp_y] order = FIRST []
[]



[Materials]
[./material] type = NeoHookean mu = 1.0 lambda = 2.0 disp_x = disp_x disp_y = disp_y []
#[./material] type = LinearElasticityMaterial mu = 1.0 lambda = 20.0 disp_x = disp_x disp_y = disp_y []
[]

[Kernels]
[./linearelasticityx] type = Elasticity variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y [../]
[./linearelasticityy] type = Elasticity variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y [../]
[]

#[NodalKernels]
#[./lagrangemultiplierx] type = ObstacleContactLagrangeMultiplier variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y lambda = lambda   []
#[./lagrangemultipliery] type = ObstacleContactLagrangeMultiplier variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y lambda = lambda   []
#[./obstacle] type = EnforceObstacleConstraint variable = lambda  disp_x = disp_x disp_y = disp_y []
#[]


[BCs]
[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = NeumannBC variable = disp_y value = 0.0 boundary = top  [../]

[./mny] type = FunctionDirichletBC variable = disp_y function = -0.1*t boundary = middle_node  [../]

[./bottomx] type = DirichletBC variable = disp_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = DirichletBC variable = disp_y value = 0.0 boundary = bottom  [../]
#[./bottomy] type = FunctionDirichletBC variable = disp_y function = 0.1*t boundary = bottom  [../]


#[./bottomy] type = FunctionNeumannBC variable = disp_y function =  boundary = bottom  [../]

#[./left] type = NeumannBC   variable = pressure value = 1.0 boundary = left  [../]
#[./left] type = DirichletBC variable = pippo_x value = 0.0 boundary = left [../]
#[./rite] type = DirichletBC variable = pressure value = 0.0 boundary = right [../]
#[./top] type = DirichletBC variable = pressure value = 0.0 boundary = top [../]
[]

#[Preconditioning]
#[./myprec] type = FDP full = true []
#[]

[Executioner]
 type = Transient
 solve_type = NEWTON
 # solve_type = LINEAR
 start_time = 0.0
 end_time = 2
 dt = 0.1
 dtmin = 0.1
 #nl_rel_tol = 1e-1
 # petsc_options = '-ksp_view_pmat'   # Stampa a schermo la matrice Jacobiana
 line_search = none
 # petsc_options_iname = '-pc_type -pc_hypre_type '
 # petsc_options_value = ' hypre    boomeramg     '
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '


#[Quadrature] type = GRID order = ELEVENTH []

[]

[Outputs]
 file_base = prova
 exodus = true
[]


