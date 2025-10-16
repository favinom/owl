[Mesh]

[gmg]
 type = DistributedRectilinearMeshGenerator # GeneratedMeshGenerator #
 dim = 2
 nx = 50
 ny = 50
# elem_type = QUAD9
[]

[]

[Variables]
[density] []
[]

[Kernels]
[./diffusion] type = SAC_Diffusion variable = density alpha = 0.01 beta = 100[]
[./timederivative] type = TimeDerivative variable = density []
[]

[ICs]
[./density_ic] type = FunctionIC variable = 'density' function = 0.25+0.55*(x<0.2)*(y<0.2)[]
#[./density_ic] type = FunctionIC variable = 'density' function = x*y[]
[]


[BCs]
[./left] type = NeumannBC variable = density value = 0.0 boundary = left  [../]
[./rite] type = NeumannBC variable = density value = 0.0 boundary = right [../]

#[./left] type = DirichletBC variable = pressure value = 0.0 boundary = left [../]
#[./top] type = DirichletBC variable = pressure value = 0.0 boundary = top [../]
[]


#[Preconditioning]
#[./myprec] type = FDP full = true []
#[]


[Executioner]
 type = Transient
 solve_type = NEWTON
 start_time = 0.0
 end_time = 10.0
 dt = 0.1
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


