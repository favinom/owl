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
[]

[BCs]
[./left] type = NeumannBC   variable = pressure value = 1.0 boundary = left  [../]
[./rite] type = DirichletBC variable = pressure value = 1.0 boundary = right [../]
[]

[Executioner]
 type = Steady
 solve_type = LINEAR
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


