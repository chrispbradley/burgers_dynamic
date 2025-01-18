PROGRAM BurgersDynamic

  USE OpenCMISS
  
  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Test program parameters
  REAL(OC_RP), PARAMETER :: NU_PARAM=1.0_OC_RP !Viscous coefficient
  REAL(OC_RP), PARAMETER :: LENGTH=1.0_OC_RP !Lenth of domain
  REAL(OC_RP), PARAMETER :: U_L=2.0_OC_RP !Dependent value at left hand of domain
  REAL(OC_RP), PARAMETER :: U_R=0.0_OC_RP !Dependent value at right hand of domain
  INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_GLOBAL_X_ELEMENTS=4 !Number of elements

  !Set time parameter
  REAL(OC_RP), PARAMETER :: DYNAMIC_SOLVER_START_TIME=0.0_OC_RP !Start time
  REAL(OC_RP), PARAMETER :: DYNAMIC_SOLVER_STOP_TIME=5.0_OC_RP !Stop time
  REAL(OC_RP), PARAMETER :: DYNAMIC_SOLVER_TIME_INCREMENT=1.0_OC_RP !Time increment

  !Set numerical parameters
  INTEGER(OC_Intg), PARAMETER :: DYNAMIC_SOLVER_OUTPUT_FREQUENCY=1 !Output every time step
  REAL(OC_RP), PARAMETER :: DYNAMIC_SOLVER_THETA=1.0_OC_RP !Time integrator theta value (1.0 = implicit)
  REAL(OC_RP), PARAMETER :: NONLINEAR_SOLVER_ABSOLUTE_TOLERANCE=1.0E-6_OC_RP 
  REAL(OC_RP), PARAMETER :: NONLINEAR_SOLVER_RELATIVE_TOLERANCE=1.0E-6_OC_RP 
  REAL(OC_RP), PARAMETER :: LINEAR_SOLVER_ABSOLUTE_TOLERANCE=1.0E-6_OC_RP 
  REAL(OC_RP), PARAMETER :: LINEAR_SOLVER_DIVERGENCE_TOLERANCE=1.0E5_OC_RP 
  REAL(OC_RP), PARAMETER :: LINEAR_SOLVER_RELATIVE_TOLERANCE=1.0E-6_OC_RP
  INTEGER(OC_Intg), PARAMETER :: LINEAR_SOLVER_MAXIMUM_ITERATIONS=10000
  INTEGER(OC_Intg), PARAMETER :: LINEAR_SOLVER_RESTART_VALUE=3000
  LOGICAL, PARAMETER :: LINEAR_SOLVER_DIRECT_FLAG=.FALSE. !Use a direct linear solver

  !User numbers
  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=6
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=7
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=8
  INTEGER(OC_Intg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=9
  INTEGER(OC_Intg), PARAMETER :: EQUATION_SET_FIELD_USER_NUMBER=10
  INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=11
  INTEGER(OC_Intg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=12
  INTEGER(OC_Intg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=13
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=14
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=15
  INTEGER(OC_Intg), PARAMETER :: DYNAMIC_SOLVER_NUMBER=1

  !Program variables
  INTEGER(OC_Intg) :: err,firstNodeDomain,decompositionIndex,equationsSetIndex,firstNodeNumber,lastNodeDomain,lastNodeNumber, &
    & myComputationalNodeNumber,numberOfComputationalNodes
  LOGICAL :: directoryExists,exportField
  TYPE(OC_BasisType) :: basis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) :: geometricField,equationsSetField,dependentField,materialsField,analyticField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_NodesType) :: nodes
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_ControlLoopType) :: controlLoop
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: dynamicSolver,nonlinearSolver,linearSolver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  !-----------------------------------------------------------------------------------------------------------
  ! INITIALISE
  !-----------------------------------------------------------------------------------------------------------

  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  CALL OC_OutputSetOn("BurgersDynamic",err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! CONTEXT
  !-----------------------------------------------------------------------------------------------------------

  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  !Get the computational nodes information
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)

  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,myComputationalNodeNumber,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a new RC coordinate system
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  !Set the coordinate system to be 1D
  CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,1,err)
  !Finish the creation of the coordinate system
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"BurgersRegion",err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(region,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a basis
  CALL OC_Basis_Initialise(basis,err)
  CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  CALL OC_Basis_TypeSet(basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CALL OC_Basis_NumberOfXiSet(basis,1,err)
  !Set the basis xi interpolation and number of Gauss points
  CALL OC_Basis_InterpolationXiSet(basis,[OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
  CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[2],err)
  !Finish the creation of the basis
  CALL OC_Basis_CreateFinish(basis,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x mesh
  CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL OC_GeneratedMesh_BasisSet(generatedMesh,basis,err)
  !Define the mesh on the region
  CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[LENGTH],err)
  CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[NUMBER_OF_GLOBAL_X_ELEMENTS],err)
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  !Finish the decomposition
  CALL OC_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !Start to create a default (geometric) field on the region
  CALL OC_Field_Initialise(geometricField,err)
  CALL OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the scaling to use
  CALL OC_Field_ScalingTypeSet(geometricField,OC_FIELD_NO_SCALING,err)
  !Set the domain to be used by the field components.
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
  !Finish creating the field
  CALL OC_Field_CreateFinish(geometricField,err)
  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations_set for a dynamic nonlinear burgers equation
  CALL OC_EquationsSet_Initialise(equationsSet,err)
  CALL OC_Field_Initialise(equationsSetField,err)
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[OC_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & OC_EQUATIONS_SET_BURGERS_EQUATION_TYPE,OC_EQUATIONS_SET_BURGERS_SUBTYPE],EQUATION_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSet,err)
  !Finish creating the equations set
  CALL OC_EquationsSet_CreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL OC_Field_Initialise(dependentField,err)
  CALL OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  !Set the mesh component to be used by the field components.
  CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
  CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,1,1,err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)
  !Initialise dependent field
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_OC_RP,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIALS FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set material field variables
  CALL OC_Field_Initialise(materialsField,err)
  CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
  !Finish the equations set material field variables
  CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)
  
  !Initialise materials field
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,-NU_PARAM,err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set analytic field variables
  CALL OC_Field_Initialise(analyticField,err)
  CALL OC_EquationsSet_AnalyticCreateStart(equationsSet,OC_EQUATIONS_SET_BURGERS_EQUATION_ONE_DIM_1, &
   & ANALYTIC_FIELD_USER_NUMBER,analyticField,err)
  !Finish the equations set analytic field variables
  CALL OC_EquationsSet_AnalyticCreateFinish(equationsSet,err)

  !Set u_L to 2.0
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,U_L,err)
  !Set u_R to 0.0
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,U_R,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL OC_Equations_Initialise(equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type (Sparse/Full)
  !CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_FULL_MATRICES,err)
  !Set the equations set output
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
  CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_FLUID_MECHANICS_CLASS, &
    & OC_PROBLEM_BURGERS_EQUATION_TYPE,OC_PROBLEM_DYNAMIC_BURGERS_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! CONTROL LOOP
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem control
  CALL OC_ControlLoop_Initialise(controlLoop,err)
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  !Get the control loop
  CALL OC_Problem_ControlLoopGet(problem,OC_CONTROL_LOOP_NODE,controlLoop,err)
  !Set the times
  CALL OC_ControlLoop_TimesSet(controlLoop,DYNAMIC_SOLVER_START_TIME,DYNAMIC_SOLVER_STOP_TIME,DYNAMIC_SOLVER_TIME_INCREMENT,err)
  !Set the output timing
  CALL OC_ControlLoop_TimeOutputSet(controlLoop,DYNAMIC_SOLVER_OUTPUT_FREQUENCY,err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solvers
  CALL OC_Solver_Initialise(dynamicSolver,err)
  CALL OC_Solver_Initialise(nonlinearSolver,err)
  CALL OC_Solver_Initialise(linearSolver,err)
  CALL OC_Problem_SolversCreateStart(problem,err)

  !Get the dymamic solver
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,DYNAMIC_SOLVER_NUMBER,dynamicSolver,err)
  !Set the output type
  !CALL OC_Solver_OutputTypeSet(dynamicSolver,OC_SOLVER_NO_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(dynamicSolver,OC_SOLVER_MATRIX_OUTPUT,err)
  !Set theta
  CALL OC_Solver_DynamicThetaSet(dynamicSolver,DYNAMIC_SOLVER_THETA,err)
  !Get the dynamic nonlinear solver
  CALL OC_Solver_DynamicNonlinearSolverGet(dynamicSolver,nonlinearSolver,err)
  !Set the nonlinear Jacobian type
  CALL OC_Solver_NewtonJacobianCalculationTypeSet(nonlinearSolver,OC_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
  !CALL OC_Solver_NewtonJacobianCalculationTypeSet(nonlinearSolver,OC_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  !Set the output type
  !CALL OC_Solver_OutputTypeSet(nonlinearSolver,OC_SOLVER_NO_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(nonlinearSolver,OC_SOLVER_MATRIX_OUTPUT,err)
  !Set the solver settings
  CALL OC_Solver_NewtonAbsoluteToleranceSet(nonlinearSolver,NONLINEAR_SOLVER_ABSOLUTE_TOLERANCE,err)
  CALL OC_Solver_NewtonRelativeToleranceSet(nonlinearSolver,NONLINEAR_SOLVER_RELATIVE_TOLERANCE,err)
  !Get the dynamic nonlinear linear solver
  CALL OC_Solver_NewtonLinearSolverGet(nonlinearSolver,linearSolver,err)
  !Set the output type
  CALL OC_Solver_OutputTypeSet(linearSolver,OC_SOLVER_NO_OUTPUT,err)
  !Set the solver settings
  IF(LINEAR_SOLVER_DIRECT_FLAG) THEN
    CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
    CALL OC_Solver_LibraryTypeSet(linearSolver,OC_SOLVER_MUMPS_LIBRARY,err)
  ELSE
    CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,err)
    CALL OC_Solver_LinearIterativeMaximumIterationsSet(linearSolver,LINEAR_SOLVER_MAXIMUM_ITERATIONS,err)
    CALL OC_Solver_LinearIterativeDivergenceToleranceSet(linearSolver,LINEAR_SOLVER_DIVERGENCE_TOLERANCE,err)
    CALL OC_Solver_LinearIterativeRelativeToleranceSet(linearSolver,LINEAR_SOLVER_RELATIVE_TOLERANCE,err)
    CALL OC_Solver_LinearIterativeAbsoluteToleranceSet(linearSolver,LINEAR_SOLVER_ABSOLUTE_TOLERANCE,err)
    CALL OC_Solver_LinearIterativeGMRESRestartSet(linearSolver,LINEAR_SOLVER_RESTART_VALUE,err)
  ENDIF
  !Finish the creation of the problem solver
  CALL OC_Problem_SolversCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem solver equations
  CALL OC_SolverEquations_Initialise(solverEquations,err)
  CALL OC_Problem_SolverEquationsCreateStart(problem,err)
  !Get the dynamic solver equations
  CALL OC_Solver_Initialise(dynamicSolver,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,dynamicSolver,err)
  CALL OC_Solver_SolverEquationsGet(dynamicSolver,solverEquations,err)
  !Set the solver equations sparsity (Sparse/Full)
  !CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES,err)
  CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_FULL_MATRICES,err)
  !Add in the equations set
  CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Set the initial value to u_L+u_R/2.0
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1, &
    & (U_L+U_R)/2.0_OC_RP,err)
  
  !Set up the boundary conditions
  CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  !Set the fixed boundary conditions at the first node and last nodes
  firstNodeNumber=1
  CALL OC_Nodes_Initialise(nodes,err)
  CALL OC_Region_NodesGet(region,nodes,err)
  CALL OC_Nodes_NumberOfNodesGet(nodes,LastNodeNumber,err)
  CALL OC_Decomposition_NodeDomainGet(decomposition,1,firstNodeNumber,firstNodeDomain,err)
  CALL OC_Decomposition_NodeDomainGet(decomposition,1,lastNodeNumber,lastNodeDomain,err)
  IF(firstNodeDomain==myComputationalNodeNumber) THEN
    CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1, &
      & OC_NO_GLOBAL_DERIV,firstNodeNumber,1,OC_BOUNDARY_CONDITION_FIXED,1.0_OC_RP,err)
  ENDIF
  IF(lastNodeDomain==myComputationalNodeNumber) THEN
    CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1, &
      & OC_NO_GLOBAL_DERIV,lastNodeNumber,1,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
  ENDIF
  CALL OC_SolverEquations_BoundaryConditionsAnalytic(solverEquations,err)
  !Finish the creation of the equations set boundary conditions
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)


  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------
  
  INQUIRE(file="./output", EXIST=directoryExists)
  IF (.NOT.directoryExists) CALL EXECUTE_COMMAND_LINE("mkdir ./output")
  
  !Solve the problem
  WRITE(*,'("Solving problem...")')
  CALL OC_Problem_Solve(problem,err)
  WRITE(*,'("Problem solved!")')

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  !Output Analytic analysis
  Call OC_AnalyticAnalysis_Output(dependentField,"BurgersDynamicAnalytics",err)

  !export fields
  exportField=.FALSE.
  IF(exportField) THEN
    CALL OC_Fields_Initialise(fields,err)
    CALL OC_Fields_Create(region,fields,err)
    CALL OC_Fields_NodesExport(fields,"BurgersDynamic","FORTRAN",err)
    CALL OC_Fields_ElementsExport(fields,"BurgersDynamic","FORTRAN",err)
    CALL OC_Fields_Finalise(fields,err)
  ENDIF
  
  !-----------------------------------------------------------------------------------------------------------
  ! FINALISE
  !-----------------------------------------------------------------------------------------------------------

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)
  
  WRITE(*,'("Program successfully completed.")')

END PROGRAM BurgersDynamic
