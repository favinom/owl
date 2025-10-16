[Mesh]

[gmg]
 type = DistributedRectilinearMeshGenerator # GeneratedMeshGenerator #
 dim = 2
 nx = 10
 ny = 10
# elem_type = QUAD9
[]

[]

[Variables]
[pressure] []

[]

[Kernels]
[./diffusion] type = Diffusion variable = pressure [../]
[./bodyforce] type = BodyForce variable = pressure value = 5.0 []
[./reaction] type = Reaction variable = pressure rate = 2.0 []
[./timederivative] type = MyTimeDerivative variable = pressure []
[]

[BCs]
#[./left] type = NeumannBC   variable = pressure value = 1.0 boundary = left  [../]
[./left] type = DirichletBC variable = pressure value = 0.0 boundary = left [../]
[./rite] type = DirichletBC variable = pressure value = 0.0 boundary = right [../]
[./top] type = DirichletBC variable = pressure value = 0.0 boundary = top [../]
[]

[Executioner]
 type = Transient
 solve_type = LINEAR
 start_time = 0.0
 end_time = 2.0
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


