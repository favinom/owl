[Mesh]

[./gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 20
     xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
    #parallel_type = replicated
  []

  [./nodeset]
    type = BoundingBoxNodeSetGenerator
    input = gmg
    new_boundary = middle_node
    top_right = '0.52 0.01 0'
    bottom_left = '0.48 -0.01 0'
  []
[]

[Variables]
[disp_x] order = FIRST []
[disp_y] order = FIRST []
[lambda] order = FIRST []
[J_k] []
[mu_k] []
[]



[Materials]
#[./material] type = HolmesMow mu = 0.37*1000 lambda = 2.46*1000 disp_x = disp_x disp_y = disp_y []
[./material] type = NeoHookeanPlasticityGradeZero mu = 1e3 lambda = 2.0e3 disp_x = disp_x disp_y = disp_y []
[./materialms] type = GradeOneGrowth lambda_p = 0 J_k = J_k []
[]

[Kernels]
[./linearelasticityx] type = ElastoPlasticity variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y [../]
[./linearelasticityy] type = ElastoPlasticity variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y [../]

[./standardtimederivativeJk] type = StandardTimeDerivative variable = J_k     J_k = J_k     mu_k = mu_k [../]
[./diffusivetimederivativeJk] type = DiffusiveTimeDerivative variable = J_k     m = 1e-6      kv = 100     J_k = J_k      mu_k = mu_k [../]
[./chempotdiff] type = ChemicalPotentialDiffusion variable = J_k     m = 1e-6      J_k = J_k       mu_k = mu_k [../]

[./chempotdefinizione] type = ChemicalPotentialReaction variable = mu_k     d = 100      J_k = J_k      mu_k = mu_k [../]
[./couplingmu] type = CouplingGrowthMechanics      variable = mu_k J_k = J_k      mu_k = mu_k  disp_x = disp_x disp_y = disp_y  [../]
[./potenziale] type = DoubleWellCH variable = mu_k coeff = 1000.0   min1 = 1   max = 1.02   min2 = 1.05   J_k = J_k   mu_k = mu_k [../]
[]


[NodalKernels]
[./lagrangemultiplierx] type = ObstacleContactLagrangeMultiplier variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y lambda = lambda   []
[./lagrangemultipliery] type = ObstacleContactLagrangeMultiplier variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y lambda = lambda   []
[./obstacle] type = EnforceObstacleConstraint variable = lambda  disp_x = disp_x disp_y = disp_y []
[]


[ICs]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1-0.5*cos(pi*x)*cos(pi*y)  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 0.2+10*exp(-10*(x^2+y^2))  [../]
[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1  [../]
#[./mu_k_ic] type = FunctionIC  variable = 'mu_k'  function = pi*pi*cos(pi*x)*cos(pi*y)  [../]
[]



[BCs]
[./leftx] type = NeumannBC variable = disp_x value = 0.0 boundary = left [../]
[./lefty] type = NeumannBC variable = disp_y value = 0.0 boundary = left [../]

[./ritex] type = NeumannBC variable = disp_x value = 0.0 boundary = right [../]
[./ritey] type = NeumannBC variable = disp_y value = 0.0 boundary = right [../]

[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = NeumannBC variable = disp_y value = 0.0 boundary = top  [../]

[./bottomx] type = NeumannBC variable = disp_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = DirichletBC variable = disp_y value = 0.0 boundary = bottom  [../]
[./mny] type = DirichletBC variable = disp_y value = 0.0 boundary = middle_node  [../]
[./mnx] type = DirichletBC variable = disp_x value = 0.0 boundary = middle_node  [../]



[./mukall] type = NeumannBC variable = mu_k value = 0.0 boundary = 'bottom left right top'  [../]
[./Jkall] type = NeumannBC variable = J_k value = 0.0 boundary = 'bottom left right top' [../]

#[./leftJk] type = DirichletBC variable = J_k value = 0.0 boundary = left [../]
#[./leftmu] type = NeumannBC variable = mu_k value = 0.0 boundary = left [../]

#[./riteJk] type = NeumannBC variable = J_k value = 0.0 boundary = right [../]
#[./ritemu] type = NeumannBC variable = mu_k value = 0.0 boundary = right [../]

#[./topJk] type = NeumannBC variable = J_k value = 0.0 boundary = top  [../]
#[./topmu] type = NeumannBC variable = mu_k value = 0.0 boundary = top  [../]

#[./bottomJk] type = NeumannBC variable = J_k value = 0.0 boundary = bottom  [../]
#[./bottommu] type = NeumannBC variable = mu_k value = 0.0 boundary = bottom  [../]

[]

[Preconditioning]
#[./myprec] type = FDP full = true []
[]

[Executioner]
 type = Transient
 solve_type = NEWTON
 #solve_type = LINEAR
 start_time = 0.0
 end_time = 60
 dt = 0.1
 line_search = none
 # petsc_options_iname = '-pc_type -pc_hypre_type '
 # petsc_options_value = ' hypre    boomeramg     '
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '


#[Quadrature] type = GRID order = ELEVENTH []

[]


[Postprocessors]
 [./average]    type = ElementAverageValue    variable = mu_k  [../]
[]



[Outputs]
 file_base = prova
 exodus = true
 csv = true
[]


