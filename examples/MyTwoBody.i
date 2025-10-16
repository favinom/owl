[Mesh]

[gmg]
 type = GeneratedMeshGenerator # GeneratedMesh  # DistributedRectilinearMeshGenerator
 dim = 2
 nx = 25
 ny = 10
 xmin = 0
 xmax = 2.5
 ymin = 0
 ymax = 1
elem_type = QUAD9
[]



  [SubdomainBoundingBox1]
    type = SubdomainBoundingBoxGenerator
    input = gmg
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '1 1 0'
  []
  
  
  [SubdomainBoundingBox2]
    type = SubdomainBoundingBoxGenerator
    input = SubdomainBoundingBox1
    block_id = 2
    bottom_left = '1.5 0 0'
    top_right = '2.5 1 0'
  []
  
  
  [SubdomainBoundingBox3]
    type = SubdomainBoundingBoxGenerator
    input = SubdomainBoundingBox2
    block_id = 3
    bottom_left = '1 0 0'
    top_right = '1.5 1 0'
  []

  [interior_sideset1]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = 1
    paired_block = 3
    input = SubdomainBoundingBox3
    new_boundary = right1
  []
  
  
   [interior_sideset2]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = 2
    paired_block = 3
    input = interior_sideset1
    new_boundary = left2
  []
  


  [ed0]
    type = BlockDeletionGenerator
    input = interior_sideset2
    block = 3
  []
  
  
  
    [ed1]
    type = BlockDeletionGenerator
    input = ed0
    block = 2
  []
  
  
   [./createNewSidesetOne]
    type = SideSetsFromPointsGenerator
    input = ed1
    new_boundary= 'top1 left1 bottom1'

    points = '0.5 1 0
    		0 0.5 0
    		0.5 0 0'
  []
  
  
[]

[Variables]
[disp_x] order = SECOND []
[disp_y] order = SECOND []
[]



[Materials]
#[./material] type = DeSaintVenant mu = 1.0 lambda = 50.0 disp_x = disp_x disp_y = disp_y []
[./material] type = NeoHookean mu = 1.0 lambda = 100.0 disp_x = disp_x disp_y = disp_y []
[]

[Kernels]
[./linearelasticityx] type = Elasticity variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y [../]
[./linearelasticityy] type = Elasticity variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y [../]
[]

[BCs]
[./left1x] type = DirichletBC variable = disp_x value = 0.0 boundary = left1 [../]
[./left1y] type = DirichletBC variable = disp_y value = 0.0 boundary = left1 [../]

[./ritex] type = FunctionDirichletBC variable = disp_x function = 0.1*t boundary = right1 [../]
[./ritey] type = DirichletBC variable = disp_y value = 0.0 boundary = right1 [../]

#[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
#[./topy] type = NeumannBC variable = disp_y value = 0.0 boundary = top  [../]

#[./bottomx] type = NeumannBC variable = disp_x value = 0.0 boundary = bottom  [../]
#[./bottomy] type = NeumannBC variable = disp_y value = 0.0 boundary = bottom  [../]


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
 dt = 0.02
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


