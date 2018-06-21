function varargout=class2equilibriumLatentCS(varargin)	
% To get help, type class2equilibriumLatentCS('help')
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
            'Creates a matlab class for computing a Nash equilibrium'
            'of the form'
            '%  P1objective(P1variables^*,P2variables^*,latentVariables^*,parameters) ='
            '%        = minimize    P1objective(P1variables,P2variables^*,latentVariables,parameters)'
            '%          w.r.t.      P1variables,latentVariables'
            '%          subject to  P1constraints(P1variables,P2variables^*,latentVariables,parameters)'
            '%                      latentConstraints(P1variables,P2variables^*,latentVariables,parameters)'
            '%  P2objective(P1variables^*,P2variables^*,latentVariables^*,parameters) ='
            '%        = minimize    P2objective(P1variables^*,P2variables,latentVariables,parameters)'
            '%          w.r.t.      P2variables,latentVariables'
            '%          subject to  P2constraints(P1variables^*,P2variables,latentVariables,parameters)'
            '%                      latentConstraints(P1variables^*,P2variables,latentVariables,parameters)'
            'and returns'
            '  outputExpressions(P1variables^*,P2variables^*,latentVariables^*,parameters)'
            'See ipm.pdf for details of the optimization engine.'
            ' '
            'The solver is accessed through a matlab class.'
                });

    %% Declare input parameters
    declareParameter(...
        'VariableName','classname',...
        'DefaultValue',getFromPedigree(),...
        'Description',{
            'Name of the class to be created.'
            'A matlab class will be created with this name plus a |.m| extension.';
            'The class will have the following methods:';
            '   * |obj=classname()| - creates class and loads the dynamic library';
            '                         containing the C code';
            '   * |delete(obj)|     - deletes the class and unload the dynamic library';
            '   * |setP_{parameter}(obj,value)|'
            '                   - sets the value of one of the parameters'
            '   * |setV_{variable}(obj,value)|'
            '                   - sets the value of one of the optimization variables'
            '   * |[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter))|'
            'where'
            '   * |mu0|      - initial value for the barrier variable'
            '   * |maxIter|  - maximum number of Newton iterations'
            '   * |saveIter| - Parameter not used (included for compatibility'
            '                  with cmex2equilibriumLatentCS).'
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
            'The solver may override "mu0" if it deems it too low, so by giving'
            '"mu0" equal to 0 one gets the smallest value acceptable to the solver.'
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
        'VariableName','P1objective',...
        'Description',{
            'Scalar Tcalculus symbolic object to be optimized for player 1.'
                      });
    declareParameter(...
        'VariableName','P1optimizationVariables',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'variables to be optimized by player 1.'
                      });
    declareParameter(...
        'VariableName','P1constraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints for player 1. Both equality and inequality'
            'constraints are allowed.'
                      });

    declareParameter(...
        'VariableName','P2objective',...
        'Description',{
            'Scalar Tcalculus symbolic object to be optimized for player 2.'
                      });
    declareParameter(...
        'VariableName','P2optimizationVariables',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'variables to be optimized by player 2.'
                      });
    declareParameter(...
        'VariableName','P2constraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints for player 2. Both equality and inequality'
            'constraints are acceptable.'
                      });

    declareParameter(...
        'VariableName','latentVariables',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'common latent variables.'
                      });
    declareParameter(...
        'VariableName','latentConstraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints that define the latent variables.'
            'Only equality constraints are acceptable.'
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
            'variables to be returned.'
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
            'Maximum number of Newton iterations.'
                      });

    declareParameter(...
        'VariableName','addEye2Hessian',...
        'DefaultValue',0*1e-6,...
        'Description',{
            'Add to the Hessian matrix appropriate identity matrices scaled by this constant.'
            '[Should be avoided since does not seem to have a good justification for Nash'
            'equilibria.]'
                      });
    
    declareParameter(...
        'VariableName','umfpack',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                      });
    
    declareParameter(...
        'VariableName','scratchbookType',...
        'DefaultValue','double',...
        'AdmissibleValues',{'double','float'},...
        'Description',{
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                      });
    
    declareParameter(...
        'VariableName','fastRedundancyCheck',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                      });
    
    declareParameter(...
        'VariableName','codeType',...
        'DefaultValue','C',...
        'AdmissibleValues',{'C','C+asmSB','C+asmLB'},...
        'Description',{
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                      });

    declareParameter(...
        'VariableName','minInstructions4loop',...
        'DefaultValue',100,...
        'Description',{
            'Parameter not used (included for compatibility with |cmex2optimizeCS|).'
                      });
        
    declareParameter(...
        'VariableName','compilerOptimization',...
        'DefaultValue','-O1',...
        'AdmissibleValues',{'-O0','-O1','-O2','-O3','-Ofast'},...
        'Description',{
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                      });
    
    declareParameter(...
        'VariableName','callType',...
        'DefaultValue','dynamicLibrary',...
        'AdmissibleValues',{'dynamicLibrary','client-server'},...
        'Description', {
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                       });

    declareParameter(...
    'VariableName','serverProgramName',...
    'DefaultValue','',...
    'Description', {
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
        });

    declareParameter(...
        'VariableName','serverAddress',...
        'DefaultValue','localhost',...
        'Description', {
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                       });
    
    declareParameter(...
        'VariableName','port',...
        'DefaultValue',1968,...
        'Description', {
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                       });

    declareParameter(...
        'VariableName','compileGateways',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                       });
    
    declareParameter(...
        'VariableName','compileLibrary',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                       });
    
    declareParameter(...
        'VariableName','compileStandalones',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                       });
    
    declareParameter(...
        'VariableName','targetComputer',...
        'DefaultValue',lower(computer),...
        'AdmissibleValues',{'maci64','glnxa64','pcwin64'},...
        'Description', {
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
                       });
    
    declareParameter(...
        'VariableName','serverComputer',...
        'DefaultValue',lower(computer),...
        'AdmissibleValues',{'maci64','glnxa64','pcwin64'},...
        'Description', {
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
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
            'Includes additional output to help debug failed convergence.'
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
            'Parameter not used (included for compatibility with cmex2equilibriumLatentCS).'
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
        fprintf('class2equilibriumLatentCS: outputs folder ''%s'' does not exist, creating it\n',folder);
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
            error('all parameters must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,class(parameters{i}));
        end
        if ~isequal(type(parameters{i}),'variable')
            parameters{i}
            error('all parameters must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(parameters{i}));
        end
    end

    if ~iscell(P1optimizationVariables)
        P1optimizationVariables
        error('P1optimizationVariables must be a cell array of Tcalculus variables');
    end

    for i=1:length(P1optimizationVariables)
        if ~isequal(class(P1optimizationVariables{i}),'Tcalculus')
            P1optimizationVariables{i}
            error('all P1optimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,class(P1optimizationVariables{i}));
        end
        if ~isequal(type(P1optimizationVariables{i}),'variable')
            P1optimizationVariables{i}
            error('all P1optimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(P1optimizationVariables{i}));
        end
    end

    if ~iscell(P2optimizationVariables)
        P2optimizationVariables
        error('P2optimizationVariables must be a cell array of Tcalculus variables');
    end

    for i=1:length(P2optimizationVariables)
        if ~isequal(class(P2optimizationVariables{i}),'Tcalculus')
            P1optimizationVariables{i}
            error('all P1optimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,class(P2optimizationVariables{i}));
        end
        if ~isequal(type(P2optimizationVariables{i}),'variable')
            P2optimizationVariables{i}
            error('all P2optimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(P2optimizationVariables{i}));
        end
    end

    if ~isempty(size(P1objective))
        error('P1''s minimization criterion must be scalar (not [%s])',...
              index2str(size(P1objective)));
    end

    if ~isempty(size(P2objective))
        error('P2''s minimization criterion must be scalar (not [%s])',...
              index2str(size(P2objective)));
    end

    if ~isempty(P1constraints) && ~iscell(P1constraints)
        error('P1''s constraints parameter must be a cell array\n');
    end

    if ~isempty(P2constraints) && ~iscell(P2constraints)
        error('P2''s constraints parameter must be a cell array\n');
    end

    if ~iscell(latentVariables)
        latentVariables
        error('latentOptimizationVariables must be a cell array of Tcalculus variables');
    end

    if ~isempty(latentConstraints) && ~iscell(latentConstraints)
        error('latent constraints parameter must be a cell array\n');
    end

    fprintf('class2equilibriumLatentCS: ...');
    t_classCS=clock();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Declare the problem-specific variables 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf(' declaring sets for parameters and primal variables... ');

    t_csparse=clock();
    debug=false;
    tprod2matlab=true;
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
    % Player 1
    for i=1:length(P1optimizationVariables)
        declareSet(code,P1optimizationVariables{i},...
                   sprintf('setV_%s',name(P1optimizationVariables{i})));
        msize=size(P1optimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                 name(P1optimizationVariables{i}),index2str(msize));
    end
    % Player 2
    for i=1:length(P2optimizationVariables)
        declareSet(code,P2optimizationVariables{i},...
                   sprintf('setV_%s',name(P2optimizationVariables{i})));
        msize=size(P2optimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                 name(P2optimizationVariables{i}),index2str(msize));
    end
    % Latent variables
    for i=1:length(latentVariables)
        declareSet(code,latentVariables{i},sprintf('setV_%s',name(latentVariables{i})));
        msize=size(latentVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                 name(latentVariables{i}),index2str(msize));
    end
    
    %% Define constraints, dual variables, and declare 'sets' for initializing dual variables

    if verboseLevel>1
        fprintf('Defining constraints and dual variables... ');
    end

    % Player 1
    [Gu,Fu,P1nus,P1lambdas,outputExpressions]=...
        parseConstraints(code,classname,P1constraints,outputExpressions,'P1');
    
    % Player 2
    [Gd,Fd,P2nus,P2lambdas,outputExpressions]=...
        parseConstraints(code,classname,P2constraints,outputExpressions,'P2');
    
    % Latent constraints (only one H, but two sets of dual variables)
    [H,~,P1xnus,err,outputExpressions]=...
        parseConstraints(code,classname,latentConstraints,outputExpressions,'P1x');
    [H,~,P2xnus,err,outputExpressions]=...
        parseConstraints(code,classname,latentConstraints,outputExpressions,'P2x');

    if ~isempty(err)
        latentConstraints
        error('latent constraints cannot be inequalities\n');
    end
            
    % HCells={};
    % P1xnus={};
    % P2xnus={};
    % for k=1:length(latentConstraints)
    %     switch (type(latentConstraints{k}))
    %       case {'iszero'}
    %         % remove iszero
    %         op1=Tcalculus(operands(latentConstraints{k}));
    %         HCells{end+1}=op1; 
    %         % create appropriate nu
    %         P1xnus{end+1}=Tvariable(sprintf('P1xnu%d_',length(P1xnus)+1),size(op1));
    %         P2xnus{end+1}=Tvariable(sprintf('P2xnu%d_',length(P2xnus)+1),size(op1));
    %         declareSet(code,P1xnus{end},sprintf('setD_%s',classname,name(P1xnus{end})));
    %         declareSet(code,P2xnus{end},sprintf('setD_%s',classname,name(P2xnus{end})));
    %       otherwise
    %         error('latent constraint of type ''%s'' not implemented\n',...
    %               type(latentConstraints{k}));
    %     end
    % end
    
    %% Pack constraints

    fprintf('packing expressions and variables... ');

    %% Pack primal variables
    % Player 1
    [u,packU,unpackU,P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H]...
        =packVariables(P1optimizationVariables,'u_',...
                       P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H);
    u0=packExpressions(P1optimizationVariables);
    % Player 2
    [d,packD,unpackD,P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H]...
        =packVariables(P2optimizationVariables,'d_',...
                       P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H);
    d0=packExpressions(P2optimizationVariables);
    src={u0,d0};
    dst={u,d};
    % Latent
    if ~isempty(latentVariables)
        [x,packX,unpackX,P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H]...
            =packVariables(latentVariables,'x_',...
                           P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H);
        x0=packExpressions(latentVariables);
        src{end+1}=x0;
        dst{end+1}=x;
    else
        x=Tzeros(0);
    end
    
    if length(H)~=length(x)
        error('number of latent variables (=%d) must match number of latent (equality) constraints (=%d) so that latent variables are uniquely defined',length(x),length(H));
    end
    
    declareCopy(code,dst,src,'initPrimal__');

    %% Pack dual variables
    % Player 1
    if size(Gu,1)>0
        P1nu=packVariables(P1nus,'P1nu_');
        P1nu0=packExpressions(P1nus);
        src{end+1}=P1nu0;
        dst{end+1}=P1nu;
    else
        P1nu=Tzeros(0);
    end
    if size(Fu,1)>0
        P1lambda=packVariables(P1lambdas,'P1lambda_');
        P1lambda0=packExpressions(P1lambdas);
        src{end+1}=P1lambda0;
        dst{end+1}=P1lambda;
    else
        P1lambda=Tzeros(0);
    end
    % Player 2
    if size(Gd,1)>0
        P2nu=packVariables(P2nus,'P2nu_');
        P2nu0=packExpressions(P2nus);
        src{end+1}=P2nu0;
        dst{end+1}=P2nu;
    else
        P2nu=Tzeros(0);
    end
    if size(Fu,1)>0
        P2lambda=packVariables(P2lambdas,'P2lambda_');
        P2lambda0=packExpressions(P2lambdas);
        src{end+1}=P2lambda0;
        dst{end+1}=P2lambda;
    else
        P2lambda=Tzeros(0);
    end
    if size(H,1)>0
        P1xnu=packVariables(P1xnus,'P1xnu_');
        P2xnu=packVariables(P2xnus,'P2xnu_');
        P1xnu0=packExpressions(P1xnus);
        P2xnu0=packExpressions(P2xnus);
        src{end+1}=P1xnu0;
        dst{end+1}=P1xnu;
        src{end+1}=P2xnu0;
        dst{end+1}=P2xnu;
    else
        P1xnu=Tzeros(0);
        P2xnu=Tzeros(0);
    end
    declareCopy(code,dst,src,'initPrimalDual__');

    %% Generate the code for the functions that do the raw computation
    t_ipmPD=clock();
    ipmPDeqlat_CS(code,P1objective,P2objective,u,d,x,P1lambda,P1nu,P1xnu,P2lambda,P2nu,P2xnu,...
                  Fu,Gu,Fd,Gd,H,...
                  smallerNewtonMatrix,addEye2Hessian,skipAffine,...
                  false,...
                  classname,allowSave,debugConvergence,profiling);
    code.statistics.time.ipmPD=etime(clock,t_ipmPD);
    
    %% Declare ipm solver 
    classhelp{end+1}='% Solve optimization';
    classhelp{end+1}='[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));';
    defines.nZ=size(u,1)+size(d,1)+size(x,1);
    defines.nU=size(u,1);
    defines.nD=size(d,1);
    defines.nX=size(x,1);
    defines.nG=size(Gu,1)+size(Gd,1)+size(H,1);
    defines.nF=size(Fu,1)+size(Fd,1);
    defines.nNu=size(Gu,1)+size(Gd,1)+2*size(H,1);
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
    
    pth=fileparts(which('class2equilibriumLatentCS.m'));
    declareFunction(code,fsfullfile(pth,'ipmPDeq_CSsolver.m'),'solve',defines);
    
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

    fprintf(' done creating matlab code (%.3f sec)\n',etime(clock,t_csparse));

    %% debug info to be passed to debugConvergenceAnalysis
    
    debugInfo.P1optimizationVariables=P1optimizationVariables;
    debugInfo.P2optimizationVariables=P2optimizationVariables;
    debugInfo.P1constraints=P1constraints;
    debugInfo.P2constraints=P2constraints;
    
    code.statistics.time.cmexCS=etime(clock,t_classCS);
    fprintf('done class2equilibriumLatentCS (%.3f sec)\n',etime(clock,t_classCS));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    varargout=setOutputs(nargout,params);

end