[Mesh]

[gmg]
 type = DistributedRectilinearMeshGenerator # GeneratedMeshGenerator #
 dim = 2
 nx = 20
 ny = 2
 xmin = 0
 xmax = 10
 ymin = 0
 ymax = 1
# elem_type = QUAD9
[]

[]

[Variables]
[pippo_x] []
[pippo_y] []
[]

[Kernels]
[./linearelasticityx] type = LinearElasticity variable = pippo_x mu = 20000.0 lambda = 10000.0 component = 0 disp_x = pippo_x disp_y = pippo_y [../]
[./linearelasticityy] type = LinearElasticity variable = pippo_y mu = 20000.0 lambda = 10000.0 component = 1 disp_x = pippo_x disp_y = pippo_y [../]
[]

[BCs]
[./leftx] type = DirichletBC variable = pippo_x value = 0.0 boundary = left [../]
[./lefty] type = DirichletBC variable = pippo_y value = 0.0 boundary = left [../]

[./ritex] type = DirichletBC variable = pippo_x value = 0.0 boundary = right [../]
[./ritey] type = DirichletBC variable = pippo_y value = 0.0 boundary = right [../]

[./topx] type = NeumannBC variable = pippo_x value = 0.0 boundary = top  [../]
[./topy] type = NeumannBC variable = pippo_y value = 0.0 boundary = top  [../]

[./bottomx] type = NeumannBC variable = pippo_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = FunctionNeumannBC variable = pippo_y function = 10*x*(x-10) boundary = bottom  [../]

#[./left] type = NeumannBC   variable = pressure value = 1.0 boundary = left  [../]
#[./left] type = DirichletBC variable = pippo_x value = 0.0 boundary = left [../]
#[./rite] type = DirichletBC variable = pressure value = 0.0 boundary = right [../]
#[./top] type = DirichletBC variable = pressure value = 0.0 boundary = top [../]
[]

[Preconditioning]
[./myprec] type = FDP full = true []
[]

[Executioner]
 type = Steady
 solve_type = NEWTON
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


