function varargout=class2optimizeCS(varargin)	
% To get help, type class2optimizeCS('help')
%
% Copyright 2012-2017 Joao Hespanha

% This file is part of Tencalc.
%
% TensCalc is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% TensCalc is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.    

    %% Function global help
    declareParameter(...
        'Help', {
            'creates a matlab class for solving optimization problems'
            'of the form:'
            '%  objective(optimizationVariables^*,parameters) ='
            '%       = minimum     objective(optimizationVariables,parameters)'
            '%         w.r.t.      optimizationVariables'
            '%         subject to  constraints(optimizationVariables,parameters)'
            'and returns'
            '%  outputExpressions(optimizationVariables^*,parameters)'
            ' '
            'The solver is accessed through a matlab class.'
            'See |ipm.pdf| for details of the optimization engine.';
                });

    %% Declare input parameters
    declareParameter(...
        'VariableName','classname',...
        'DefaultValue',getFromPedigree(),...
        'Description',{
            'Name of the class to be created.'
            'A matlab class will be created with this name plus a |.m| extension.';
            'The class will have the following methods:';
            '   * |obj=classname()| - creates class';
            '   * |delete(obj)|     - deletes the class';
            '   * |setP_{parameter}(obj,value)|'
            '                   - sets the value of one of the parameters'
            '   * |setV_{variable}(obj,value)|'
            '                   - sets the value of one of the optimization variables'
            '   * |[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter))|'
            'where'
            '   * |mu0|      - initial value for the barrier variable'
            '   * |maxIter|  - maximum number of Newton iterations'
            '   * |saveIter| - Parameter not used (included for compatibility'
            '                  with |cmex2optimizeCS|).'
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
            'Path to the folder where the files will be created.';
            'Needs to be in the Matlab path.'
                       });
    
    declareParameter(...
        'VariableName','objective',...
        'Description',{
            'Scalar Tcalculus symbolic object to be optimized.'
                      });
    declareParameter(...
        'VariableName','optimizationVariables',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'variables to be optimized.'
                      });
    declareParameter(...
        'VariableName','constraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints. Both equality and inequality constraints'
            'are allowed.'
                      });
    declareParameter(...
        'VariableName','parameters',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the '
            'parameters (must be given, not optimized).'
                      });

    declareParameter(...
        'VariableName','outputExpressions',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the '
            'variables to be returned.';
            ' ';
            'The following Tcalculus symbolic variables are assigned special values';
            'and can be using in outputExpressions';
            '* |lambda1_|,|_lambda2_|,... - Lagrangian multipliers associated with';
            '                         the inequalities constraints';
            '                         (in the order that they appear and with';
            '                          the same size as the corresponding constraints)';
            '* |nu1_|,|nu2_|,...    - Lagrangian multipliers associated with';
            '                         the equality constraints';
            '                         (in the order that they appear and with';
            '                          the same size as the corresponding constraints)';
            '* |Hess_| - Hessian matrix used by the (last) newton step to update';
            '            the primal variables (not including |addEye2Hessian|).'
            '* |dHess_| - D factor in the LDL factorization of the Hessian matrix'
            '             used by the (last) newton step to update the primal variables'
            '             (including |addEye2Hessian|, unlike Hess_).'
            ' ';
            'ATTENTION: To be able to include these variables as input parameters,';
            '           they will have to be created outside this function'
            '           *with the appropriate sizes*.'
            '           Eventually, their values will be overridden by the solver'
            '           to reflect the values above.'
                      });

    declareParameter(...
        'VariableName','method',...
        'DefaultValue','primalDual',...
        'AdmissibleValues',{'primalDual'},...
        'Description',{
            'Variable that specifies which method should be used:'
            '* |primalDual| - interior point primal-dual method'
            '* |barrier|    - interior point barrier method'
            '                   (not yet implemented)'
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
            'Should only be set lower to 1 for very poorly scaled problems.'
                      });

    declareParameter(...
        'VariableName','coupledAlphas',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true| the same scalar gain is used for the primal and dual variables'
            'in Newton''s method line search.'
                      });

    declareParameter(...
        'VariableName','skipAffine',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true| the affine search direction step is omitted.'
                      });

    declareParameter(...
        'VariableName','useLDL',...
        'DefaultValue',true,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true| the search directions are computed using an'
            'LDL instead of an LU factorization.'
                      });
    
    declareParameter(...
        'VariableName','umfpack',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                      });
    
    declareParameter(...
        'VariableName','smallerNewtonMatrix',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true| the matrix that needs to be inverted to compute a Newton step'
            'is reduced by first eliminating the dual variables associated with inequality'
            'constrainst.'
            'However, often the smaller matrix is not as sparse so the computation'
            'may actually increase.'
                      });
    
    declareParameter(...
        'VariableName','delta',...
        'DefaultValue',3,...
        'AdmissibleValues',{2,3},...
        'Description',{
            'Delta parameter used to determine mu based on the affine direction.'
            'Set |Delta=3| for well behaved problems (for an aggressive convergence)'
            'and |Delta=2| in poorly conditioned problems (for a more robust behavior).'
            'This parameter is only used when |skipAffine=false|.'
                      });
    declareParameter(...
        'VariableName','muFactorAggressive',...
        'DefaultValue',1/3,...
        'Description',{
            'Multiplicative factor used to update the barrier parameter'
            '(must be smaller than 1).'
            'This value is used when there is good progress along the'
            'Newton direction.'
            'Nice convex problems can take as low as 1/100, but '
            'poorly conditioned problems may require as high as 1/3.'
            'This parameter is only used when |skipAffine=true|.'
                      });
    declareParameter(...
        'VariableName','muFactorConservative',...
        'DefaultValue',.75,...
        'Description',{
            'Multiplicative factor used to update the barrier parameter'
            '(must be smaller than 1).'
            'This value is used when there is poor or no progress along the'
            'Newton direction. A value not much smaller than one is preferable.'
                      });

    declareParameter(...
        'VariableName','gradTolerance',...
        'DefaultValue',1e-4,...
        'Description',{
            'Maximum norm for the gradient below which the first order optimality'
            'conditions assumed to by met.'
                      });

    declareParameter(...
        'VariableName','equalTolerance',...
        'DefaultValue',1e-4,...
        'Description',{
            'Maximum norm for the vector of equality constraints below which the'
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
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                      });

    declareParameter(...
        'VariableName','addEye2Hessian',...
        'DefaultValue',1e-9,...
        'Description',{
            'Add to the Hessian matrix appropriate identity matrices scaled by this constant.'
                      });
    
    declareParameter(...
        'VariableName','scratchbookType',...
        'DefaultValue','double',...
        'AdmissibleValues',{'double','float'},...
        'Description',{
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                      });
    
    declareParameter(...
        'VariableName','fastRedundancyCheck',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                      });
    
    declareParameter(...
        'VariableName','codeType',...
        'DefaultValue','C',...
        'AdmissibleValues',{'C','C+asmSB','C+asmLB'},...
        'Description',{
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                      });

    declareParameter(...
        'VariableName','compilerOptimization',...
        'DefaultValue','-O1',...
        'AdmissibleValues',{'-O0','-O1','-O2','-O3','-Ofast'},...
        'Description',{
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                      });
    
    declareParameter(...
        'VariableName','callType',...
        'DefaultValue','dynamicLibrary',...
        'AdmissibleValues',{'dynamicLibrary','client-server'},...
        'Description', {
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                       });

    declareParameter(...
    'VariableName','serverProgramName',...
    'DefaultValue','',...
    'Description', {
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
        });

    declareParameter(...
        'VariableName','serverAddress',...
        'DefaultValue','localhost',...
        'Description', {
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                       });
    
    declareParameter(...
        'VariableName','port',...
        'DefaultValue',1968,...
        'Description', {
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                       });

    declareParameter(...
        'VariableName','compileGateways',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                       });
    
    declareParameter(...
        'VariableName','compileLibrary',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                       });
    
    declareParameter(...
        'VariableName','compileStandalones',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                       });
    
    declareParameter(...
        'VariableName','targetComputer',...
        'DefaultValue',lower(computer),...
        'AdmissibleValues',{'maci64','glnxa64','pcwin64'},...
        'Description', {
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                       });
    
    declareParameter(...
        'VariableName','serverComputer',...
        'DefaultValue',lower(computer),...
        'AdmissibleValues',{'maci64','glnxa64','pcwin64'},...
        'Description', {
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                       });

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
            'Includes additional printed output to help debug failed convergence.'
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
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                      });

    declareParameter(...
        'VariableName','profiling',...
        'DefaultValue',false,...
        'AdmissibleValues',{true,false},...
        'Description',{
            'When nonzero, adds profiling to the C code.'
                      });
    
    %% Output parameters
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Retrieve parameters and inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [stopNow,params]=setParameters(nargout,varargin);
    if stopNow
        return 
    end

    %% transfer any folder in classname into folder
    [folder,classname]=fileparts(fsfullfile(folder,classname));

    %% create folder if it does not exist
    if ~strcmp(folder,'.') && ~exist(folder,'dir')
        fprintf('class2optimizeCS: outputs folder ''%s'' does not exist, creating it\n',folder);
        if mkdir(folder)==0
            error('Unable to create folder ''%s''\n',folder)
        end
    end
    
    rmpath(folder);
    addpath(folder);
    
    %% Fix class when gotten from pedigree
    classname=regexprep(classname,'+TS=','_TS_');
    classname=regexprep(classname,'-','_');
    classname=regexprep(classname,'+','_');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check input parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~iscell(parameters)
        parameters
        error('parameters must be a cell array of Tcalculus variables');
    end

    for i=1:length(parameters)
        if ~isequal(class(parameters{i}),'Tcalculus')
            parameters{i}
            error('%dth parameter must be of class ''Tcalculus'' (not [%s])\n',...
                  i,class(parameters{i}));
        end
        if ~isequal(type(parameters{i}),'variable')
            parameters{i}
            error('%dth parameter must be of the type ''variable'' (not [%s])\n',...
                  i,type(parameters{i}));
        end
    end

    if ~iscell(optimizationVariables)
        optimizationVariables
        error('optimizationVariables must be a cell array of Tcalculus variables');
    end

    for i=1:length(optimizationVariables)
        if ~isequal(type(optimizationVariables{i}),'variable')
            optimizationVariables{i}
            error('%dth optimizationVariables must be of the type ''variable'' (not [%s])\n',...
                  i,type(optimizationVariables{i}));
        end
    end

    if ~isempty(size(objective))
        error('Minimization criterion must be scalar (not [%s])',index2str(size(objective)));
    end
    
    if ~isempty(constraints) && ~iscell(constraints)
        error('Constraints parameter must be a cell array\n');
    end

    fprintf('class2optimizeCS:... ');
    t_classCS=clock();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Declare the problem-specific variables 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t_csparse=clock();
    debug=false;
    tprod2matlab=true;
    scratchbookType='double';
    code=csparse(scratchbookType,debug,tprod2matlab); % using instructionsTable.c
    classhelp={'% Create object';
               sprintf('obj=%s();',classname)};

    %% Declare 'sets' for initializing parameters
    if length(parameters)>0
        classhelp{end+1}='% Set parameters';
    end
    for i=1:length(parameters)
        declareSet(code,parameters{i},sprintf('setP_%s',name(parameters{i})));
        msize=size(parameters{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setP_%s(obj,{[%s] matrix});',...
                                name(parameters{i}),index2str(msize));
    end

    %% Declare 'sets' for initializing primal variables
    classhelp{end+1}='% Initialize primal variables';
    for i=1:length(optimizationVariables)
        declareSet(code,optimizationVariables{i},...
                   sprintf('setV_%s',name(optimizationVariables{i})));
        msize=size(optimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                name(optimizationVariables{i}),index2str(msize));
    end

    %% Define constraints, dual variables, and declare 'sets' for initializing dual variables

    if verboseLevel>1
        fprintf('Defining constraints and dual variables... ');
    end

    [G,F,nus,lambdas,outputExpressions]=...
        parseConstraints(code,classname,constraints,outputExpressions);
    
    %% Pack primal variables
    [u,packU,unpackU,objective,outputExpressions,F,G]...
        =packVariables(optimizationVariables,'x_',objective,outputExpressions,F,G);
    u0=packExpressions(optimizationVariables);

    src={u0};
    dst={u};
    declareCopy(code,dst,src,'initPrimal__');

    %% Pack dual variables
    if size(G,1)>0
        [nu,~,~,outputExpressions]=packVariables(nus,'nu_',outputExpressions);
        nu0=packExpressions(nus);
        src{end+1}=nu0;
        dst{end+1}=nu;
    else
        nu=Tzeros(0);
    end
    if size(F,1)>0
        [lambda,~,~,outputExpressions]=packVariables(lambdas,'lambda_',outputExpressions);
        lambda0=packExpressions(lambdas);
        src{end+1}=lambda0;
        dst{end+1}=lambda;
    else
        lambda=Tzeros(0);
    end
    declareCopy(code,dst,src,'initPrimalDual__');

    %% Generate the code for the functions that do the raw computation
    t_ipmPD=clock();
    [Hess__,dHess__]=ipmPD_CS(code,objective,u,lambda,nu,F,G,...
                              smallerNewtonMatrix,addEye2Hessian,skipAffine,...
                              useLDL,false,...
                              classname,allowSave,debugConvergence);
    code.statistics.time.ipmPD=etime(clock,t_ipmPD);
    outputExpressions=substitute(outputExpressions,...
                                 Tvariable('Hess_',size(Hess__)),Hess__);
    outputExpressions=substitute(outputExpressions,...
                                 Tvariable('dHess_',size(dHess__)),dHess__);

    %% Declare ipm solver 
    classhelp{end+1}='% Solve optimization';
    classhelp{end+1}=...
        sprintf('[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));');

    defines.nU=size(u,1);
    defines.nG=size(G,1);
    defines.nF=size(F,1);
    defines.gradTolerance=gradTolerance;
    defines.equalTolerance=equalTolerance;
    defines.desiredDualityGap=desiredDualityGap;
    defines.alphaMin=alphaMin;
    defines.alphaMax=alphaMax;
    defines.coupledAlphas=double(coupledAlphas);
    defines.muFactorAggressive=muFactorAggressive;
    defines.muFactorConservative=muFactorConservative;
    defines.delta=delta;
    defines.skipAffine=double(skipAffine);
    defines.allowSave=double(allowSave);
    defines.debugConvergence=double(debugConvergence);
    defines.debugConvergenceThreshold=debugConvergenceThreshold;
    defines.profiling=double(profiling);
    defines.verboseLevel=solverVerboseLevel;
    
    pth=fileparts(which('class2optimizeCS.m'));
    declareFunction(code,fsfullfile(pth,'ipmPD_CSsolver.m'),'solve',defines);
    
    %% Declare 'gets' for output expressions
    classhelp{end+1}='% Get outputs';
    classhelp{end+1}='';
    for i=1:length(outputExpressions)
        classhelp{end}=sprintf('%sy%d,',classhelp{end},i);
    end
    classhelp{end}=sprintf('[%s]=getOutputs(obj);',classhelp{end}(1:end-1));
    declareGet(code,outputExpressions,'getOutputs');

    code.statistics.time.csparse=etime(clock,t_csparse);
    code.statistics.defines=defines;
    
    fprintf('  done creating csparse object (%.3f sec)\n',etime(clock,t_csparse));

    %% Compile code
    fprintf(' creating matlab code... ');
    t_compile2matlab=clock();
    compile2matlab(code,...
                   fsfullfile(folder,sprintf('%s.m',classname)),...
                   fsfullfile(folder,sprintf('%s.log',classname)),...
                   classhelp,...
                   profiling);
    code.statistics.time.compile2matlab=etime(clock,t_compile2matlab);

    fprintf(' done creating matlab code (%.3f sec)\n',etime(clock,t_compile2matlab));

    %% debug info to be passed to debugConvergenceAnalysis
    
    debugInfo.optimizationVariables=optimizationVariables;
    debugInfo.constraints=constraints;
    
    code.statistics.time.classCS=etime(clock,t_classCS);
    fprintf('done class2optimizeCS (%.3f sec)\n',etime(clock,t_classCS));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    varargout=setOutputs(nargout,params);

end