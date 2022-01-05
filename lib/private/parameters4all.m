function localVariables_=parameters4all(localVariables_)
% Declare input parameters common to the 4 tenscalc functions:
%   cmex2optimizeCS.m
%   class2optimizeCS.m
%   cmex2equilibriumLatentCS.m
%   class2equilibriumLatentCS.m
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    declareParameter(...
        'VariableName','parameters',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the '
            'parameters (must be given, not optimized).'
                      });

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Symbolic computation parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    declareParameter(...
        'VariableName','packOptimizationVariables',...
        'AdmissibleValues',{false,true},...
        'DefaultValue',true,...
        'Description',{
            'When true all optimization variables are packed into a'
            'single column vector before performing symbolic'
            'differentiation. This speeds up symbolic differentiation,'
            'but may result in slower matlab solvers.'
            ' '
            'This parameter should not affect the speed of C solvers.'
                      });

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Interior-point solver stopping criteria parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    declareParameter(...
        'VariableName','gradTolerance',...
        'DefaultValue',1e-4,...
        'Description',{
            'Maximum infinity norm for the gradient below which the first order optimality'
            'conditions assumed to by met.'
                      });

    declareParameter(...
        'VariableName','equalTolerance',...
        'DefaultValue',1e-4,...
        'Description',{
            'Maximum infinite norm for the vector of equality constraints below which the'
            'equalities are assumed to hold.'
                      });

    declareParameter(...
        'VariableName','desiredDualityGap',...
        'DefaultValue',1e-5,...
        'Description',{
            'Value for the duality gap that triggers the end of'
            'the constrained optimization. The overall optimization '
            'terminates at the end of the first Newton step for which'
            'the duality gap becomes smaller than this value.'
                      });

    declareParameter(...
        'VariableName','maxIter',...
        'DefaultValue',200,...
        'Description',{
            'Maximum number of Newton iterations.'
                      });

    declareParameter(...
        'VariableName','LDLthreshold',...
        'DefaultValue',1e-5,...
        'Description',{
            'Pivot threshold for MATLAB''s LDL factorization.'
            ' '
            'From help LDL: LDLthreshold must be a double scalar lying in'
            'the interval [0, 0.5]. The default value for THRESH is 0.01.'
            'Using smaller values of LDLthreshold may give faster factorization'
            'times and fewer entries, but may also result in a less stable factorization.'
            ' '
            'Currently this only affects the generation of matlab code.'
                      });
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Interior-point Scaling parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    declareParameter(...
        'VariableName','scaleInequalities',...
        'AdmissibleValues',{false,true},...
        'DefaultValue',true,...
        'Description',{
            'When true, all inequalities are multiplied by appropriate'
            'constants so that at the first iteration they start as ''1>=0''.'
                      });

    declareParameter(...
        'VariableName','scaleCost',...
        'DefaultValue',0,...
        'Description',{
            'When nonzero, the cost is multiplied by an appropriate constant'
            'so that at the first iteration it starts at ''scaleCost''.';
            'A run-time error will occur if the cost starts at zero.'
                      });

    evalin('caller','scaleEqualities=false;');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Interior-point solver parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % declareParameter(...
    %     'VariableName','method',...
    %     'DefaultValue','primalDual',...
    %     'AdmissibleValues',{'primalDual'},...
    %     'Description',{
    %         'Variable that specifies which method should be used:'
    %         '* |primalDual| - interior point primal-dual method'
    %         '* |barrier|    - interior point barrier method'
    %         '                   (not yet implemented)'
    %                   });

    declareParameter(...
        'VariableName','muFactorAggressive',...
        ...%'DefaultValue',1/3,...
        'DefaultValue',.2,...
        'Description',{
            'Multiplicative factor used to update the barrier parameter |mu|'
            '(must be smaller than 1).'
            'This value is used when there is good progress along the'
            'Newton direction, the equality constraints are satisfied,'
            'and the gradient is small.'
            'Nice convex problems can take as low as 1/100, but '
            'poorly conditioned problems may require as high as 1/3.'
            'This parameter is only used when |skipAffine=true|.'
                      });
    declareParameter(...
        'VariableName','muFactorConservative',...
        'DefaultValue',.95,...
        ...%'DefaultValue',1,...
        'Description',{
            'Multiplicative factor used to update the barrier parameter |mu|'
            '(must be smaller than 1).'
            'This value is used when there is good progress along the'
            'Newton direction, the equality constraints are satisfied,'
            'but the gradient is still large.'
            'A value not much smaller than one (or 1) is preferable.'
                      });

    declareParameter(...
        'VariableName','skipAffine',...
        'DefaultValue',true,...
        'AdmissibleValues',{false,true},...
        'Description',{;
            'When |false| the barrier parameter |mu| is updated based on how much';
            'progress can be achieved in the search direction obtained with |mu=0|.';
            'This is known as the ''affine search direction step''.';
            'This step (obtained by setting to |skipAffine|=|false|) can significantly'
            'speed up convergence by rapidly decreasing the barrier parameter.'
            'However, it can be fragile for tough non-convex problems.'
                      });

    declareParameter(...
        'VariableName','delta',...
        'DefaultValue',3,...
        'AdmissibleValues',{2,3},...
        'Description',{
            'Delta parameter used to determine mu based on the affine search direction step:'
            'Set |delta=3| for well behaved problems (for an aggressive convergence)'
            'and |delta=2| in poorly conditioned problems (for a more robust behavior).'
            'This parameter is only used when |skipAffine=false|.'
                      });

    declareParameter(...
        'VariableName','alphaMin',...
        'DefaultValue',1e-7,...
        'Description',{
            'Minimum value for the scalar gain in Newton''s method line search,'
            'below which a search direction is declared to have failed.'
                      });

    declareParameter(...
        'VariableName','alphaMax',...
        'DefaultValue',1,...
        'Description',{
            'Maximum value for the scalar gain in Newton''s method line search.'
            'Should only be set lower than 1 for very poorly scaled problems.'
                      });

    declareParameter(...
        'VariableName','coupledAlphas',...
        'DefaultValue',true,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true| the same scalar gain is used for the primal and dual variables'
            'in Newton''s method line search.';
            ' '
            'Setting to |false| speeds up convergence in some problems, but appears';
            'to be somewhat fragile, so it should be set to |false| with caution.'
                      });

    declareParameter(...
        'VariableName','smallerNewtonMatrix',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true| the matrix that needs to be inverted to compute a Newton step'
            'is reduced by first eliminating the dual variables associated with inequality'
            'constraints.'
            'However, often the smaller matrix is not as sparse so the computation'
            'may actually increase.'
                      });

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Debugging parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    declareParameter(...
        'VariableName','solverVerboseLevel',...
        'DefaultValue',2,...
        'Description',{
            'Level of verbose for the solver outputs:';
            '* 0  - the solver does not produce any output and does not report timing information';
            '* 1  - the solver does not produce any output but reports timing information';
            '* 2  - the solver only report a summary of the solution status'
            '       when the optimization terminates';
            '* 3  - the solver reports a summary of the solution status '
            '       at each iteration';
            '* >3 - the solver produces several (somewhat unreadable) outputs'
            '       at each iteration step'
                      });

    declareParameter(...
        'VariableName','debugConvergence',...
        'DefaultValue',false,...
        'AdmissibleValues',{true,false},...
        'Description',{
            'Includes additional printed output to help debug failed convergence.';
            ' ';
            'ATTENTION: Not implementd in cmex2optimizeCS().'
            '           debugConvergence should be used in class2optimizeCS().'
                      });

    declareParameter(...
        'VariableName','debugConvergenceThreshold',...
        'DefaultValue',1e5,...
        'Description',{
            'Threshold above which solver warns about large values.';
            'Only used when |debugConvergence=true|';
                      });

    declareParameter(...
        'VariableName','allowSave',...
        'DefaultValue',false,...
        'AdmissibleValues',{true,false},...
        'Description',{
            'Generates code that permit saving the "hessian" matrix, for subsequent'
            'optimization of pivoting, row/column permutations, and scaling'
            'for the hessian''s LU factorization.'
            'The hessian in saved in two files named'
            '%   |{classname}_WW.subscripts| and |{classname}_WW.values|'
            'that store the sparsity structure and the actual values, respectively,'
            'at some desired iteration (see output parameter |solver|).'
            ' ';
            'This parameter is only used for C-code solvers.';
                      });

    declareParameter(...
        'VariableName','profiling',...
        'DefaultValue',false,...
        'AdmissibleValues',{true,false},...
        'Description',{
            'When |true|, adds profiling to the C code.';
            ' '
            'Accumulated profiling information is diplayed on the screen when the'
            'dynamic library is unloaded.'
            ' ';
            'This parameter is only used for C-code solvers.';
                      });

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Code generation parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    declareParameter(...
        'VariableName','classname',...
        'DefaultValue',getFromPedigree(),...
        'Description', {
            'Name of the class to be created.'
            'A matlab class will be created with this name plus a |.m| extension.';
            'The class will have the following methods:';
            '   * |obj=classname()| - creates class and loads the dynamic library';
            '                         containing the C code';
            '   * |delete(obj)|     - deletes the class and unloads the dynamic library';
            '   * |setP_{parameter}(obj,value)|'
            '                   - sets the value of one of the parameters'
            '   * |setV_{variable}(obj,value)|'
            '                   - sets the value of one of the optimization variables';
            '   * |[y1,y2,...]=getOutputs(obj)|';
            '                   - gets the values of the |outputExpressions|';
            '   * |[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter),addEye2Hessian)|'
            'where'
            '   * |mu0|      - initial value for the barrier variable (default=1)'
            '   * |maxIter|  - maximum number of Newton iterations (default=200)'
            '   * |saveIter| - iteration # when to save the "hessian" matrix (default=-1)'
            '                  (for subsequent pivoting/permutations/scaling optimization)'
            '                  only saves when |allowSave| is true.';
            ' '
            '                  When |saveIter=0|, the hessian matrix is saved'
            '                                     at the last iteration; and'
            '                  when |saveIter<0|, the hessian matrix is not saved.';
            ' '
            '                  The "hessian" matrix will be saved regardless of the';
            '                  value of |saveIter|, when the solver exists with |status=4|';
            '   * |addEye2Hessian| - see help for |addEye2Hessian| and |adjustAddEye2Hessian|'
            '   * |status|   - solver exist status';
            '                 *  0  = success';
            '                 *  >0 = solver terminated unexpectedly';
            '                 a nonzero status indicates the reason for termination';
            '                 in a binary format:'
            '                   bit  0 = 1 - (primal) variables violate constraints'
            '                   bit  1 = 1 - dual variables are negative'
            '                   bit  2 = 1 - failed to invert hessian'
            '                   bit  3 = 1 - maximum # of iterations reached';
            '                 when the solver exists because the maximum # of iterations'
            '                 was reached (bit 3 = 1), the remaining bits provide'
            '                 information about the solution returned';
            '                   bit  4 = 1 - gradient larger then |gradTolerance|'
            '                   bit  5 = 1 - equality constraints violate |equalTolerance|';
            '                   bit  6 = 1 - duality gap larger than |desiredDualityGap|';
            '                   bit  7 = 1 - barrier variable larger than minimum value';
            '                   bit  8 = 1 - scalar gain |alpha| in Newton direction';
            '                                smaller than alphaMin'
            '                   bit  9 = 1 - scalar gain |alpha| in Newton direction';
            '                                smaller than .1'
            '                   bit 10 = 1 - scalar gain |alpha| in Newton direction';
            '                                smaller than .5'
            '   * |iter|    - number of iterations'
            '   * |time|    - solver''s compute time (in secs).'
            ' '
            'One can look "inside" this class to find the name of the cmex functions.'
                       });

    declareParameter(...
        'VariableName','folder',...
        'DefaultValue','.',...
        'Description', {
            'Path to the folder where the class and cmex files will be created.';
            'The folder will be created if it does not exist and it will be added';
            'to the begining of the path if not there already.'
                       });

    declareParameter(...
        'VariableName','simulinkLibrary',...
        'DefaultValue','',...
        'Description', {
            'Name of a simulink library to be created with Simulink blocks that can be used'
            'to call the different solver functions.'
            'The blocks are created with direct feedthrough.'
            'No library and no simulink blocks are created if |simulinkLibrary| is an empty string.'
                       });


    if 0
        declareParameter(...
            'VariableName','codeType',...
            'DefaultValue','C',...
            'AdmissibleValues',{'C','C+asmSB','C+asmLB'},...
            'Description',{
                'Type of code produced:'
                '* |C| - all computations done in pure C code. Ideal for final code.'
                '          '
                '          Impact on non-optimized compilation:'
                '            * medium compilation times'
                '            * largest code size'
                '            * slowest run times';
                '          '
                '          Impact on optimized code:'
                '            Gives the most freedom to the compiler for optimization'
                '            * slowest compile optimization times'
                '            * fastest run times'
                '            * smallest code sizes';
                '          '
                '* |C+asmLB| - little C code, with most of the computations done'
                '          by large blocks of inlined assembly code. Ideal for testing.'
                '          '
                '          Impact on non-optimized compilation:'
                '            * fastest compilation times'
                '            * smallest code size'
                '            * fastest run times (for non-optimized code)'
                '          '
                '          Impact on optimized compilation:'
                '            Most of the compiler optimization is restricted to re-ordering'
                '            and/or inlining the large blocks of asm code'
                '            * fastest compile optimization times'
                '            * slowest run times'
                '            * largest optimized code sizes (due to inlining large blocks)';
                '          '
                '          NOT FULLY IMPLEMENTED.';
                '          '
                '* |C+asmSB| - little C code, with most of the computations done'
                '          by small blocks of inlined assembly code'
                '          '
                '          Impact on non-optimized compilation:'
                '            * medium compilation times'
                '            * medium code size'
                '            * medium run times'
                '          '
                '          Impact on optimized code:'
                '            Most of the compiler optimization is restricted to re-ordering'
                '            and/or inlining the small blocks of asm code'
                '            * medium compile optimization times,'
                '            * medium run times'
                '            * medium code sizes';
                '          '
                '          NOT FULLY IMPLEMENTED NOR TESTED.'
                ' ';
                'This parameter is only used for C-code solvers.';
                          });
    end

    declareParameter(...
        'VariableName','scratchbookType',...
        'DefaultValue','double',...
        'AdmissibleValues',{'double','float'},...
        'Description',{
            'C variable type used for the scratchbook. See csparse documentation.';
            ' ';
            'This parameter is only used for C-code solvers.';
                      });
    declareParameter(...
        'VariableName','fastRedundancyCheck',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'when true, very intensive operations (like the lu factorization) '
            'do not check if it is possible to reuse computations.'
            'See csparse documentation.'
            ' ';
            'This parameter is only used for C-code solvers.';
                      });

    declareParameter(...
        'VariableName','minInstructions4loop',...
        'DefaultValue',50,...
        'Description',{
            'Minimum number of similar instruction that will be execute as part of a'
            'for(;;) loop, rather than being executed as independent C commands.'
            'When equal to ''inf'', instructions will never be grouped into foor loops.'
            ' ';
            'This parameter is only used for C-code solvers.';
                      });

    declareParameter(...
        'VariableName','maxInstructionsPerFunction',...
        'DefaultValue',100,...
        'Description',{
            'Maximum number of instructions to be included into a single function.'
            'When equal to ''inf'', there is no limit on the size of a single function.';
            ' ';
            'Large values of |maxInstructionsPerFunction| and therefore large functions';
            'give more opportunities for compiler optimization, but can resul in very';
            'slow compilation (especially with compiler optimization turned on).';
            ' ';
            'This parameter is only used for C-code solvers.';
                      });

    declareParameter(...
        'VariableName','useUmfpack',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true| the system of linear equations needed for the Newton step';
            'is solved using the umfpack package.';
            'Using this package tyically increases the solve time, but results in';
            'code that is more robust.'
            ' ';
            'This parameter is only used for C-code solvers.';
                      });

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compilation/linking parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    declareParameter(...
        'VariableName','compilerOptimization',...
        'DefaultValue','-O0',...
        'AdmissibleValues',{'-O0','-O1','-O2','-O3','-Os','-Ofast'},...
        'Description',{
            'Optimization parameters passed to the C compiler.'
            '* |-Ofast| often generates the fastest code'
            '* |-Os| often generates the smallest code'
            '* |-O1| includes the most important optimizations, while still compiling fast'
            '* |-O0| compiles the fastest'
            'Only used when either |compileGateways|, |compileLibrary|, or '
            '|compileStandalones| are set to |true|.'
            ' ';
            'This parameter is only used for C-code solvers.';
        });

    declareParameter(...
        'VariableName','compileGateways',...
        'DefaultValue',false,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'When |true| the gateway functions are compiled using |cmex|.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });
    declareParameter(...
        'VariableName','compileLibrary',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'When |true| the dynamicLibrary is compiled using |gcc|.'
            'This parameter is used only when |callType=''dynamicLibrary''|.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });
    declareParameter(...
        'VariableName','compileStandalones',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'When |true| the standalone/server executable is compiled using |gcc|.'
            'This parameter is used only when |callType| has one of the values:';
            '   |''standalone''| or |''client-server''|'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','callType',...
        'DefaultValue','dynamicLibrary',...
        'AdmissibleValues',{'dynamicLibrary','client-server'},...
        'Description', {
            'Method used to interact with the solver:';
            '* |dynamicLibrary| - the solver is linked with matlab using a dynamic library';
            '* |client-server|   - the solver runs as a server in independent process,';
            '                     and a socket is used to exchange data.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','targetComputer',...
        'DefaultValue',lower(computer),...
        'AdmissibleValues',{'maci64','glnxa64','pcwin64'},...
        'Description', {
            'OS where the C files will be compiled.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','serverComputer',...
        'DefaultValue',lower(computer),...
        'AdmissibleValues',{'maci64','glnxa64','pcwin64'},...
        'Description', {
            'OS where the server will be compiled.'
            'This parameter is used only when |callType=''client-server''|.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','absolutePath',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'When ''true'' the the cmex functions use an absolute path to open';
            'the dynamic library, which means that the dynamic library cannot';
            'be moved away from the folder where it was created.';
            ' '
            'When ''false'' no path information about the dynamic library is'
            'included in the cmex function, which must then rely on the OS-specific';
            'method used to find dynamic libraries. See documentation of ''dlopen'''
            'for linux and OSX or ''LoadLibrary'' for Microsoft Windows.'
            ' '
            'This parameter is used only when ''callType''=''dynamicLibrary''.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','serverProgramName',...
        'DefaultValue','',...
        'Description', {
            'Name of the executable file for the server executable.'
            'This parameter is used only when |callType=''client-server''|.'
                       });

    declareParameter(...
        'VariableName','serverAddress',...
        'DefaultValue','localhost',...
        'Description', {
            'IP address (or name) of the server.'
            'This parameter is used only when |callType=''client-server''|.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','port',...
        'DefaultValue',1968,...
        'Description', {
            'Port number for the socket that connects client and server.'
            'This parameter is used only when |callType=''client-server''|.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Output parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    declareOutput(...
        'VariableName','classname',...
        'Description', {
            'Name of the class created.'
                       });

    declareOutput(...
        'VariableName','code',...
        'Description', {
            'csparse object with the code generated.'
                       });

    if 0
        declareOutput(...
            'VariableName','debugInfo',...
            'Description',{
                'Data to be passed to debugConvergenceAnalysis.'
                          });
    end

end
