[Mesh]
  [cyl2d_iga]
    type = FileMeshGenerator
    file = octSphere_TRIAL.xda
  []
  #uniform_refine = 2
[]


[Variables]
[pressure] u = FIRST []

[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[BCs]
[./leftp] type = DirichletBC variable = u value = 0.0 boundary = left [../]
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


