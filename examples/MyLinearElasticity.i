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
[disp_x] []
[disp_y] []
[]

[Kernels]
[./linearelasticityx] type = LinearElasticity variable = disp_x mu = 20000.0 lambda = 10000.0 component = 0 disp_x = disp_x disp_y = disp_y [../]
[./linearelasticityy] type = LinearElasticity variable = disp_y mu = 20000.0 lambda = 10000.0 component = 1 disp_x = disp_x disp_y = disp_y [../]
[]

[BCs]
[./leftx] type = DirichletBC variable = disp_x value = 0.0 boundary = left [../]
[./lefty] type = DirichletBC variable = disp_y value = 0.0 boundary = left [../]

[./ritex] type = DirichletBC variable = disp_x value = 0.0 boundary = right [../]
[./ritey] type = DirichletBC variable = disp_y value = 0.0 boundary = right [../]

[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = NeumannBC variable = disp_y value = 0.0 boundary = top  [../]

[./bottomx] type = NeumannBC variable = disp_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = FunctionNeumannBC variable = disp_y function = 10*x*(x-10) boundary = bottom  [../]

#[./left] type = NeumannBC   variable = pressure value = 1.0 boundary = left  [../]
#[./left] type = DirichletBC variable = pippo_x value = 0.0 boundary = left [../]
#[./rite] type = DirichletBC variable = pressure value = 0.0 boundary = right [../]
#[./top] type = DirichletBC variable = pressure value = 0.0 boundary = top [../]
[]

#[Preconditioning]
#[./myprec] type = FDP full = true []
#[]

[Executioner]
 type = Steady
 solve_type = NEWTON
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


