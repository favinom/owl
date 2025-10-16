[Mesh]

#[gmg]
 type = GeneratedMesh  # DistributedRectilinearMeshGenerator # GeneratedMeshGenerator
 dim = 3 
 nx = 5
 ny = 5
 nz = 5
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
 zmin = 0
 zmax = 1
 elem_type = HEX27
#[]
[]

[Variables]
[disp_x] order = SECOND []
[disp_y] order = SECOND []
[disp_z] order = SECOND []
[]


[Materials]
#[./material] type = HolmesMow mu = 0.37 lambda = 2.46 disp_x = disp_x disp_y = disp_y []
[./material] type = HolmesMowPlasticityGradeZero mu = 1 lambda = 10.0 disp_x = disp_x disp_y = disp_y disp_z = disp_z []
[./materialms] type = GradeZeroExplicit lambda_p = 0.1 []
[]

[Kernels]
[./elasticityx] type = ElastoPlasticity variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y disp_z = disp_z[../]
[./elasticityy] type = ElastoPlasticity variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y disp_z = disp_z[../]
[./elasticityz] type = ElastoPlasticity variable = disp_z component = 2 disp_x = disp_x disp_y = disp_y disp_z = disp_z[../]
[]

[BCs]
[./leftx] type = NeumannBC variable = disp_x value = 0.0 boundary = left [../]
[./lefty] type = NeumannBC variable = disp_y value = 0.0 boundary = left [../]
[./leftz] type = NeumannBC variable = disp_z value = 0.0 boundary = left [../]

[./ritex] type = NeumannBC variable = disp_x value = 0.0 boundary = right [../]
[./ritey] type = NeumannBC variable = disp_y value = 0.0 boundary = right [../]
[./ritez] type = NeumannBC variable = disp_z value = 0.0 boundary = right [../]

[./frontx] type = NeumannBC variable = disp_x value = 0.0 boundary = front [../]
[./fronty] type = NeumannBC variable = disp_y value = 0.0 boundary = front [../]
[./frontz] type = NeumannBC variable = disp_z value = 0.0 boundary = front [../]

[./backx] type = NeumannBC variable = disp_x value = 0.0 boundary = back [../]
[./backy] type = NeumannBC variable = disp_y value = 0.0 boundary = back [../]
[./backz] type = NeumannBC variable = disp_z value = 0.0 boundary = back [../]

[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = FunctionNeumannBC variable = disp_y function = -0.1*t*(t<5)+(0.5*(t-4.95)-0.5)*(t<4.95)*(t>5.95) boundary = top  [../]
[./topz] type = NeumannBC variable = disp_z value = 0.0 boundary = top  [../]


[./bottomx] type = DirichletBC variable = disp_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = DirichletBC variable = disp_y value = 0.0 boundary = bottom  [../]
[./bottomz] type = DirichletBC variable = disp_z value = 0.0 boundary = bottom  [../]


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
 end_time = 7
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


