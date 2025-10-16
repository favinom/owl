#[GlobalParams]
#  my_common_parameter = some_value
#[]




[Mesh]

#[gmg]
 type = GeneratedMesh  # DistributedRectilinearMeshGenerator # GeneratedMeshGenerator
 dim = 2
 nx = 20
 ny = 20
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
 elem_type = QUAD9
#[]
[]

[Variables]
[disp_x] order = SECOND []
[disp_y] order = SECOND []
[lambda] order = SECOND []
[]



[Materials]
[./material] type = NeoHookean mu = 1.0 lambda = 2.0 disp_x = disp_x disp_y = disp_y []
#[./materialf] type = ObstacleFunctions set_g = (4*(x-0.5)*(x-0.5)-0.1*t+1-y) set_gx = 8*(x-0.5) set_gy = -1 set_gxx = 8 set_gxy = 0 set_gyy = 0 disp_x = disp_x disp_y = disp_y []
[./materialf] type = ObstacleFunctions set_g = -(x*x+y*y-x-2*y*(-3+1)-2*3+1)-0.1*t set_gx = 2*x-1 set_gy = 2*y-2*(-3+1) set_gxx = 2 set_gxy = 0 set_gyy = 2 disp_x = disp_x disp_y = disp_y []
#[./material] type = LinearElasticityMaterial mu = 1.0 lambda = 2.0 disp_x = disp_x disp_y = disp_y []
[]

[Kernels]
[./linearelasticityx] type = Elasticity variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y [../]
[./linearelasticityy] type = Elasticity variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y [../]
[]

[NodalKernels]
[./obstacle] type = IdentityNoContact variable = lambda  disp_x = disp_x disp_y = disp_y []
[]



[BCs]
[./lagrangemultiplierx] type = WeakObstacleContactLagMult variable = disp_x component = 0 activate_Nanson = true scoperta_del_secolo = false disp_x = disp_x disp_y = disp_y lambda = lambda boundary = top []
[./lagrangemultipliery] type = WeakObstacleContactLagMult variable = disp_y component = 1 activate_Nanson = true scoperta_del_secolo = false disp_x = disp_x disp_y = disp_y lambda = lambda boundary = top []
[./lagrangemultiplierlambda] type = WeakEnforceObstacleConstraint variable = lambda disp_x = disp_x disp_y = disp_y boundary = top []

[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = NeumannBC variable = disp_y value = 0.0 boundary = top  [../]


#[./leftl] type = DirichletBC variable = lambda value = 0.0 boundary = left  [../]
#[./ritel] type = DirichletBC variable = lambda value = 0.0 boundary = right  [../]
#[./bottoml] type = DirichletBC variable = lambda value = 0.0 boundary = bottom  [../]

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
 #solve_type = NEWTON
 solve_type = LINEAR
 start_time = 0.0
 end_time = 8
 dt = 0.1
 dtmin = 0.1
 nl_rel_tol = 1e-4
 #petsc_options = '-ksp_view_pmat'   # Stampa a schermo la matrice Jacobiana
 line_search = none
 # petsc_options_iname = '-pc_type -pc_hypre_type '
 # petsc_options_value = ' hypre    boomeramg     '
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '

[Quadrature]  type = SIMPSON [] # type = GRID order = ELEVENTH


[]


[Outputs]
 file_base = prova
 exodus = true
[]


