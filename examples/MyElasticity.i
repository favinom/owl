[Mesh]

#[gmg]
 type = GeneratedMesh  # DistributedRectilinearMeshGenerator # GeneratedMeshGenerator
 dim = 2
 nx = 20
 ny = 8
 xmin = 0
 xmax = 10
 ymin = 0
 ymax = 4
 elem_type = QUAD9
#[]
[]

[Variables]
[disp_x] order = SECOND []
[disp_y] order = SECOND []
[]



[Materials]
[./material] type = DeSaintVenant mu = 1.0 lambda = 50.0 disp_x = disp_x disp_y = disp_y []
#[./material] type = LinearElasticityMaterial mu = 1.0 lambda = 50.0 disp_x = disp_x disp_y = disp_y []
[]

[Kernels]
[./linearelasticityx] type = Elasticity variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y [../]
[./linearelasticityy] type = Elasticity variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y [../]
[]

[BCs]
[./leftx] type = DirichletBC variable = disp_x value = 0.0 boundary = left [../]
[./lefty] type = DirichletBC variable = disp_y value = 0.0 boundary = left [../]

[./ritex] type = FunctionDirichletBC variable = disp_x function = t boundary = right [../]
[./ritey] type = DirichletBC variable = disp_y value = 0.0 boundary = right [../]

[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = NeumannBC variable = disp_y value = 0.0 boundary = top  [../]

[./bottomx] type = NeumannBC variable = disp_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = NeumannBC variable = disp_y value = 0.0 boundary = bottom  [../]


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
 start_time = 0.0
 end_time = 2.0
 dt = 0.1
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


