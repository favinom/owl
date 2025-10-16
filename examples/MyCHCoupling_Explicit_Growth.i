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
[J_k] []
[mu_k] []
[]



[Materials]
#[./material] type = HolmesMow mu = 0.37*1000 lambda = 2.46*1000 disp_x = disp_x disp_y = disp_y []
[./material] type = NeoHookeanPlasticityGradeZero mu = 1e1 lambda = 2.0e1 disp_x = disp_x disp_y = disp_y []
[./materialms] type = GradeOneGrowth lambda_p = 0 J_k = J_k []
[]

[Kernels]
[./linearelasticityx] type = ElastoPlasticity variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y [../]
[./linearelasticityy] type = ElastoPlasticity variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y [../]

[./standardtimederivativeJk] type = StandardTimeDerivative variable = J_k     J_k = J_k     mu_k = mu_k [../]
[./diffusivetimederivativeJk] type = DiffusiveTimeDerivative variable = J_k     m = 1e-4      kv = 1.0e1     J_k = J_k      mu_k = mu_k [../]
[./chempotdiff] type = ChemicalPotentialDiffusion variable = J_k     m = 1e-4      J_k = J_k       mu_k = mu_k [../]

[./growth] type = GrowthRate variable = J_k  Phi_natural = 0.8  alpha = 1  beta = 0.2  J_k = J_k  mu_k = mu_k  disp_x = disp_x  disp_y = disp_y  [../]

[./chempotdefinizione] type = ChemicalPotentialReaction variable = mu_k     d = 0.1      J_k = J_k      mu_k = mu_k [../]
#[./couplingmu] type = CouplingGrowthMechanics      variable = mu_k  J_k = J_k      mu_k = mu_k  disp_x = disp_x disp_y = disp_y  [../]


#[./source] type = BodyForce variable = J_k function = 0.1*(((x-0.5)^2+(y-0.5)^2)<0.01)*t*(t<5)[../]
[]



[ICs]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1-0.5*cos(pi*x)*cos(pi*y)  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 0.2+10*exp(-10*(x^2+y^2))  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1+0.3*(x>0.899)*(y>0.899)  [../]
[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1+0.3*(y<0.199)  [../]
#[./mu_k_ic] type = FunctionIC  variable = 'mu_k'  function = pi*pi*cos(pi*x)*cos(pi*y)  [../]
[]



[BCs]
[./leftx] type = NeumannBC variable = disp_x value = 0.0 boundary = left [../]
[./lefty] type = NeumannBC variable = disp_y value = 0.0 boundary = left [../]

[./ritex] type = NeumannBC variable = disp_x value = 0.0 boundary = right [../]
[./ritey] type = NeumannBC variable = disp_y value = 0.0 boundary = right [../]

[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
#[./topy] type = FunctionNeumannBC variable = disp_y function = -800*0.05*t*(t<10) boundary = top  [../]

[./bottomx] type = DirichletBC variable = disp_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = DirichletBC variable = disp_y value = 0.0 boundary = bottom  [../]

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
 start_time = 0.0
 end_time = 10
 dt = 0.1
 line_search = none
 # petsc_options_iname = '-pc_type -pc_hypre_type '
 # petsc_options_value = ' hypre    boomeramg     '
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '


#[Quadrature] type = GRID order = ELEVENTH []

[]


[Postprocessors]
  [./average]
    type = ElementAverageValue
    variable = J_k
  [../]
[]



[Outputs]
 file_base = prova
 exodus = true
 csv = true
[]


