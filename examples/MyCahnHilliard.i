[Mesh]

[gmg]
 type = DistributedRectilinearMeshGenerator # GeneratedMeshGenerator #
 dim = 2
 nx = 10
 ny = 10
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
# elem_type = QUAD9
[]

[]

[Variables]
[J_k] []
[mu_k] []
[]

[Kernels]
[./standardtimederivativeJk] type = StandardTimeDerivative variable = J_k     J_k = J_k     mu_k = mu_k [../]
[./diffusivetimederivativeJk] type = DiffusiveTimeDerivative variable = J_k     m = 0.01      kv = 0.0     J_k = J_k      mu_k = mu_k [../]
[./chempotdiff] type = ChemicalPotentialDiffusion variable = J_k     m = 0.01      J_k = J_k       mu_k = mu_k [../]
[./chempotdefinizione] type = ChemicalPotentialReaction variable = mu_k     d = 0.001      J_k = J_k      mu_k = mu_k [../]
#[./potenziale] type = DoubleWellCH variable = mu_k coeff = 10.0   min1 = 0.7   max = 1.2   min2 = 1.3   J_k = J_k   mu_k = mu_k [../]
[]



[ICs]
[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1-0.5*cos(2*pi*x)*cos(2*pi*y)  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 0.2+10*exp(-10*(x^2+y^2))  [../]
#[./mu_k_ic] type = FunctionIC  variable = 'mu_k'  function = pi*pi*cos(pi*x)*cos(pi*y)  [../]
[]



[BCs]
[./leftx] type = NeumannBC variable = J_k value = 0.0 boundary = left [../]
[./lefty] type = NeumannBC variable = mu_k value = 0.0 boundary = left [../]

[./ritex] type = NeumannBC variable = J_k value = 0.0 boundary = right [../]
[./ritey] type = NeumannBC variable = mu_k value = 0.0 boundary = right [../]

[./topx] type = NeumannBC variable = J_k value = 0.0 boundary = top  [../]
[./topy] type = NeumannBC variable = mu_k value = 0.0 boundary = top  [../]

[./bottomx] type = NeumannBC variable = J_k value = 0.0 boundary = bottom  [../]
[./bottomy] type = NeumannBC variable = mu_k value = 0.0 boundary = bottom  [../]

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
 end_time = 100.0
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


