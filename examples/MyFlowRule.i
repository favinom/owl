[Mesh]

#[gmg]
 type = GeneratedMesh  # DistributedRectilinearMeshGenerator # GeneratedMeshGenerator
 dim = 2 #3 
 nx = 10
 ny = 10
 #nz = 10
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
 #zmin = 0
 #zmax = 1
 elem_type = QUAD9
#[]
[]

[Variables]
[disp_x] order = SECOND []
[disp_y] order = SECOND []
[]



[Materials]
#[./material] type = HolmesMow mu = 0.37 lambda = 2.46 disp_x = disp_x disp_y = disp_y []
[./material] type = NeoHookeanPlasticityGradeZero mu = 1 lambda = 50.0 disp_x = disp_x disp_y = disp_y []
[./materialms] type = GradeZeroExplicit lambda_p = 0.1 []
[]

[Kernels]
[./linearelasticityx] type = ElastoPlasticity variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y [../]
[./linearelasticityy] type = ElastoPlasticity variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y [../]
[]

[BCs]
[./leftx] type = NeumannBC variable = disp_x value = 0.0 boundary = left [../]
[./lefty] type = NeumannBC variable = disp_y value = 0.0 boundary = left [../]

[./ritex] type = NeumannBC variable = disp_x value = 0.0 boundary = right [../]
[./ritey] type = NeumannBC variable = disp_y value = 0.0 boundary = right [../]

[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = FunctionNeumannBC variable = disp_y function = -0.1*t*(t<10) boundary = top  [../]

[./bottomx] type = DirichletBC variable = disp_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = DirichletBC variable = disp_y value = 0.0 boundary = bottom  [../]


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
 end_time = 15
 dt = 0.05
 #line_search = none
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


