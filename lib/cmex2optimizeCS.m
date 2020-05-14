function varargout=cmex2optimizeCS(varargin)	
% To get help, type cmex2optimizeCS('help')
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
            'Creates a set of cmex functions for solving optimization problems'
            'of the form:'
            '%  objective(optimizationVariables^*,parameters) ='
            '%       = minimum     objective(optimizationVariables,parameters)'
            '%         w.r.t.      optimizationVariables'
            '%         subject to  constraints(optimizationVariables,parameters)'
            'and returns'
            '%  outputExpressions(optimizationVariables^*,parameters)'
            ' '
            'The solver is accessed through several cmex functions that can be';
            'accessed directly or through a matlab class.'
            'See |ipm.pdf| for details of the optimization engine.';
                });

    %% Declare input parameters
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
        'VariableName','sensitivityVariables',...
        'DefaultValue',{},...        
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'optimization variables with respect to which we want to compute cost';
            'sensitivity.'
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
        'VariableName','outputExpressions',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the '
            'variables to be returned.';
            ' ';
            'The following Tcalculus symbolic variables are assigned special values';
            'and can be using in outputExpressions';
            '* |lambda1_|,|lambda2_|,... - Lagrangian multipliers associated with';
            '                         the inequalities constraints';
            '                         (in the order that they appear and with';
            '                          the same size as the corresponding constraints)';
            '* |nu1_|,|nu2_|,...    - Lagrangian multipliers associated with';
            '                         the equality constraints';
            '                         (in the order that they appear and with';
            '                          the same size as the corresponding constraints)';
            '* |Hess_| - Hessian matrix used by the (last) Newton step to update';
            '            the primal variables (including |addEye2Hessian|).'
            '* |dHess_| - D factor in the LDL factorization of the Hessian matrix'
            '             used by the (last) Newton step to update the primal variables'
            '             (including |addEye2Hessian|).'
            '* |Grad_| - gradient of Lagrangian at the (last) Newton step.'
            '* |mu_|   - barrier parameter at the (last) Newton step.'
            '* |u_|    - vector stacked with all primal variables at the (last) Newton step.'
            '* |F_|    - vector stacked with all equalities at the (last) Newton step.'
            '* |G_|    - vector stacked with all inequalities at the (last) Newton step.'
            '* |nu_|   - vector stacked with all dual equality variables at the (last) Newton step.'
            '* |lambda_|   - vector stacked with all dual inequality variables at the (last) Newton step.'
            ' ';
            'ATTENTION: To be able to include these variables as input parameters,';
            '           they have to be previously created outside *with the appropriate sizes*.'
            '           Eventually, their values will be overridden by the solver'
            '           to reflect the values listted above.'
                      });

    localVariables_=parameters4all(localVariables_);

    declareParameter(...
        'VariableName','addEye2Hessian',...
        'DefaultValue',true,...
        'AdmissibleValues',{false,true},...        
        'Description',{
            'When |true|, adds to the Hessian matrix identity matrices scaled by a small constant,'
            'whose value is set when the solver is called.'
            ' '
            'A larger value for the scaling constant has two main effects:'
            '1) Improves the numerical conditioning of the system of equations that'
            '   finds the Newton search direction.';
            '2) Moves the search direction towards the gradient descent of';
            '   the Lagragian (and away from the Newton direction).';
            'Both effects improve the robustness of the solver, but this is typically';
            'achieved at the expense of slower convergence.'
            'For convex problems, one typically chooses |addEye2Hessian| equal to the';
            'square root of the machine precision.'
            'For non-convex problems, one can try to increase this parameter when';
            'the Newton direction actually causes an increase of the Lagragian.'
                      });
    
    declareParameter(...
        'VariableName','useLDL',...
        'DefaultValue',true,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true| the search directions are computed using an'
            'LDL instead of an LU factorization.';
            ' '
            'In general, the LDL factorization leads to faster code.';
            'However, the current implementation is restricted to a pure';
            'diagonal matrix (no 2x2 blocks in the D factor) so it may';
            'fail with the message ''ldl needs pivoting''. If this happens';
            'either set |useLDL=false| or use a nonzero value for |addEye2Hessian|.';
                      });
    
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
        fprintf('cmex2optimizeCS: outputs folder ''%s'' does not exist, creating it\n',folder);
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
    
    debugConvergence=false; % not implemented for cmex

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check input parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    parameters=checkParameters(parameters);

    if isstruct(optimizationVariables)
        optimizationVariables=struct2cell(optimizationVariables);
    end

    if ~iscell(optimizationVariables)
        optimizationVariables
        error('optimizationVariables must be a cell array of Tcalculus variables');
    end

    for i=1:length(optimizationVariables)
        if ~isequal(class(optimizationVariables{i}),'Tcalculus')
            optimizationVariables{i}
            error('all optimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,class(optimizationVariables{i}));
        end
        if ~isequal(type(optimizationVariables{i}),'variable')
            optimizationVariables{i}
            error('all optimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(optimizationVariables{i}));
        end
        for j=1:length(parameters)
            if isequal(name(optimizationVariables{i}),name(parameters{j}))
                optimizationVariables{i}
                error('optimization variable ''%s'' cannot also be a parameter\n',name(optimizationVariables{i}));
            end
        end
    end

    if ~isempty(size(objective))
        error('Minimization criterion must be scalar (not [%s])',index2str(size(objective)));
    end

    if ~isempty(constraints) && ~iscell(constraints)
        error('Constraints parameter must be a cell array\n');
    end

    [outputExpressions,outputNames]=checkOutputExpressions(outputExpressions);

    fprintf('cmex2optimizeCS:... ');
    t_cmexCS=clock();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Declare the problem-specific variables 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t_csparse=clock();
    debug=false;
    tprod2matlab=false;
    %code=csparse0(scratchbookType,debug);   % using fastTable.m
    %code=csparse1(scratchbookType,debug);  % using fastTable.m, string I_ instruction types
    %code=csparse2(scratchbookType,debug);  % using fastTable.m, integer instruction types
    code=csparse(scratchbookType,debug,tprod2matlab,fastRedundancyCheck); % using instructionsTable.c
    classhelp={'Create object';
               sprintf('obj=%s();',classname)};
    
    % template for createGateway
    template=cmextoolsTemplate();
    %% Declare 'sets' for initializing parameters
    if length(parameters)>0
        classhelp{end+1}='Set parameters';
    end
    for i=1:length(parameters)
        template(end+1,1).MEXfunction=sprintf('%s_set_%s',classname,name(parameters{i}));
        if ~isempty(simulinkLibrary)
            template(end).Sfunction=sprintf('%sS_set_%s',classname,name(parameters{i}));
        end
        template(end).Cfunction=sprintf('%s_set_%s',classname,name(parameters{i}));
        template(end).method=sprintf('setP_%s',name(parameters{i}));
        template(end).inputs =struct('type','double',... % will be overwriten after compile2C
                                     'name',name(parameters{i}),...
                                     'sizes',size(parameters{i})); % will be overwriten after compile2C
        declareSet(code,parameters{i},template(end).Cfunction);
        msize=size(parameters{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setP_%s(obj,{[%s] matrix});',...
                                name(parameters{i}),index2str(msize));
    end

    %% Declare 'sets' for initializing primal variables
    classhelp{end+1}='Initialize primal variables';
    for i=1:length(optimizationVariables)
        template(end+1,1).MEXfunction=sprintf('%s_set_%s',...
                                            classname,name(optimizationVariables{i}));
        if ~isempty(simulinkLibrary)
            template(end).Sfunction=sprintf('%sS_set_%s',classname,name(optimizationVariables{i}));
        end
        template(end).Cfunction=sprintf('%s_set_%s',classname,name(optimizationVariables{i}));
        template(end).method=sprintf('setV_%s',name(optimizationVariables{i}));
        template(end).inputs =struct('type','double',...
                                     'name',name(optimizationVariables{i}),...
                                     'sizes',size(optimizationVariables{i}));
        declareSet(code,optimizationVariables{i},template(end).Cfunction);
        msize=size(optimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                name(optimizationVariables{i}),index2str(msize));
    end

    %% Define constraints, dual variables, and declare 'sets' for
    %% initializing primal and dual variables

    if verboseLevel>0
        fprintf('  Defining primal variables, constraints, and dual variables... ');
    end

    [G,F,nus,lambdas,outputExpressions,tpl]=...
        parseConstraints(code,classname,constraints,outputExpressions);
    template=[template;tpl];    
    
    %% Pack primal variables
    [u,whereVariables,packU,unpackU,objective,outputExpressions,F,G]...
        =packVariables(optimizationVariables,'x_',objective,outputExpressions,F,G);
    u0=packExpressions(optimizationVariables);

    src={u0};
    dst={u};
    declareCopy(code,dst,src,'initPrimal__');

    %% Pack dual variables
    if size(G,1)>0
        [nu,~,~,~,outputExpressions]=packVariables(nus,'nu_',outputExpressions);
        nu0=packExpressions(nus);
        src{end+1}=nu0;
        dst{end+1}=nu;
    else
        nu=Tzeros(0);
    end
    if size(F,1)>0
        [lambda,~,~,~,outputExpressions]=packVariables(lambdas,'lambda_',outputExpressions);
        lambda0=packExpressions(lambdas);
        src{end+1}=lambda0;
        dst{end+1}=lambda;
    else
        lambda=Tzeros(0);
    end
    declareCopy(code,dst,src,'initPrimalDual__');

    %% Find indices of sensitivityVariables in x
    isSensitivity=false(length(u),1);
    for i=1:length(sensitivityVariables)
        if ~strcmp(type(sensitivityVariables{i}),'variable')
            sensitivityVariables{i}
            error('all sensitivityVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(sensitivityVariables{i}));
        end
        found=false;
        for j=1:length(optimizationVariables)
            if isequal(sensitivityVariables{i},optimizationVariables{j})
                fprintf('   sensitivityVariable %s: %d values\n',sensitivityVariables{i}.name,length(whereVariables{j}));
                isSensitivity(whereVariables{j})=true;
                found=true;
                break;
            end
        end
        if ~found
            error('sensitivityVariable %s is not an optimizationVariable\n',name(sensitivityVariables{i}));
        end
    end
    fprintf('      %d sensitivity variables\n',sum(isSensitivity));
    
    %% Generate the code for the functions that do the raw computation
    t_ipmPD=clock();
    Tout=ipmPD_CS(code,objective,u,lambda,nu,F,G,isSensitivity,...
                 smallerNewtonMatrix,addEye2Hessian,skipAffine,...
                 scaleInequalities,scaleCost,scaleEqualities,...
                 useLDL,umfpack,...
                 classname,allowSave,debugConvergence);
    code.statistics.time.ipmPD=etime(clock,t_ipmPD);
    
    % Replace solver variables into output expression
    fn=fields(Tout);
    for i=1:length(fn)
        varname=sprintf('%s_',fn{i});
        outputExpressions=substitute(outputExpressions,...
                                     Tvariable(varname,size(Tout.(fn{i})),true),Tout.(fn{i}));
    end    
    
    %% Declare ipm solver 
    classhelp{end+1}='Solve optimization';
    classhelp{end+1}='[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter),addEye2Hessian);';
    template(end+1,1).MEXfunction=sprintf('%s_solve',classname);
    if ~isempty(simulinkLibrary)
        template(end).Sfunction=sprintf('%sS_solve',classname);
    end
    template(end).Cfunction='ipmPD_CSsolver';
    template(end).method='solve';
    template(end).inputs(1) =struct('type','double','name','mu0','sizes',1,'default',1);
    template(end).inputs(2) =struct('type','int32','name','maxIter','sizes',1,'default',200);
    template(end).inputs(3) =struct('type','int32','name','saveIter','sizes',1,'default',-1);
    template(end).inputs(4) =struct('type','double','name','addEye2Hessian','sizes',1,'default',1e-9);
    template(end).outputs(1)=struct('type','int32','name','status','sizes',1);
    template(end).outputs(2)=struct('type','int32','name','iter','sizes',1);
    template(end).outputs(3)=struct('type','double','name','time','sizes',1);

    %folder='.';
    classFolder=fsfullfile(folder,sprintf('@%s',classname));
    if ~exist(classFolder,'dir')
        fprintf('class classFolder @%s does not exist, creating it... ',classname);
        mkdir(classFolder);
    end
    
    defines.saveNamePrefix=['"',fsfullfile(classFolder,classname),'"'];
    defines.nU=size(u,1);
    defines.nG=size(G,1);
    defines.nF=size(F,1);
    defines.gradTolerance=sprintf('%e',gradTolerance); % to make double
    defines.equalTolerance=sprintf('%e',equalTolerance); % to make double
    defines.desiredDualityGap=sprintf('%e',desiredDualityGap); % to make double
    defines.scaleCost=double(scaleCost~=0);
    defines.scaleInequalities=double(scaleInequalities);
    defines.scaleEqualities=double(scaleEqualities);
    defines.alphaMin=sprintf('%e',alphaMin); % to make double
    defines.alphaMax=sprintf('%e',alphaMax); % to make double
    defines.coupledAlphas=double(coupledAlphas);
    defines.muFactorAggressive=sprintf('%e',muFactorAggressive); % to make double
    defines.muFactorConservative=sprintf('%e',muFactorConservative); % to make double
    defines.delta=delta;
    defines.skipAffine=double(skipAffine);
    defines.smallerNewtonMatrix=double(smallerNewtonMatrix);
    defines.useLDL=double(useLDL);
    defines.umfpack=double(umfpack);
    defines.allowSave=double(allowSave);
    defines.debugConvergence=double(debugConvergence);
    defines.debugConvergenceThreshold=debugConvergenceThreshold;
    defines.profiling=double(profiling);
    defines.verboseLevel=solverVerboseLevel;
    
    pth=fileparts(which('cmex2optimizeCS.m'));
    declareFunction(code,fsfullfile(pth,'ipmPD_CSsolver.c'),'ipmPD_CSsolver',...
                    defines,template(end).inputs,template(end).outputs,template(end).method);

    %% Declare 'gets' for output expressions
    classhelp{end+1}='Get outputs';
    classhelp{end+1}='';
    template(end+1,1).MEXfunction=sprintf('%s_getOutputs',classname);
    if ~isempty(simulinkLibrary)
        template(end).Sfunction=sprintf('%sS_getOutputs',classname);
    end
    template(end).Cfunction=sprintf('%s_getOutputs',classname);
    template(end).method='getOutputs';
    template(end).outputs=struct('type',{},'name',{},'sizes',{});
    for i=1:length(outputExpressions)
        template(end).outputs(i).type='double';% will be overwriten after compile2C
        template(end).outputs(i).name=outputNames{i};
        template(end).outputs(i).sizes=size(outputExpressions{i});% will be overwriten after compile2C
        classhelp{end}=[classhelp{end},outputNames{i},','];
    end
    classhelp{end}=sprintf('[%s]=getOutputs(obj);',classhelp{end}(1:end-1));
    classhelp{end+1}=sprintf('[y (struct)]=getOutputs(obj);');
    declareGet(code,outputExpressions,template(end).Cfunction);

    code.statistics.time.csparse=etime(clock,t_csparse);
    code.statistics.defines=defines;
    
    fprintf('  done creating csparse object (%.3f sec)\n',etime(clock,t_csparse));

    %% Compile code
    fprintf(' Creating C code... ');
    t_compile2C=clock();
    compile2C(code,minInstructions4loop,maxInstructionsPerFunction,...
              sprintf('%s.c',classname),...
              sprintf('%s.h',classname),...
              sprintf('%s.log',classname),...
              classFolder,...
              profiling);
    % update templates in case compile2C made changes (e.g., turning matrix into sparse
    for i=1:length(template)
        for j=1:length(code.template)
            if strcmp(template(i).MEXfunction,code.template(j).MEXfunction)
                %template(i),code.template(j),;
                for l=1:length(template(i).inputs)
                    template(i).inputs(l).type=code.template(j).inputs(l).type;
                    template(i).inputs(l).sizes=code.template(j).inputs(l).sizes;
                end
                for l=1:length(template(i).outputs)
                    template(i).outputs(l).type=code.template(j).outputs(l).type;
                    template(i).outputs(l).sizes=code.template(j).outputs(l).sizes;
                end
                break
            end
        end
    end
    for j=1:length(code.template)
        if strcmp('profilingView',code.template(j).method)
            template(end+1)=code.template(j);
            break
        end
    end
    code.statistics.time.compile2C=etime(clock,t_compile2C);
    
    fprintf('  done creating C code (%.3f sec)\n',etime(clock,t_compile2C));

    classhelp{end+1}='Delete object';
    classhelp{end+1}='clear obj';

    t_createGateway=clock();
    %% Create gateway & compile library
    switch (callType) 
      case 'dynamicLibrary'
        statistics=createGateway('template',template,...
                                 'CfunctionsSource',fsfullfile(classFolder,sprintf('%s.c',classname)),...
                                 'simulinkLibrary',simulinkLibrary,'dummySimulinkIOs',true,...
                                 'callType','dynamicLibrary',...
                                 'dynamicLibrary',classname,'absolutePath',absolutePath,...
                                 'folder',folder,...
                                 'className',classname,'classHelp',classhelp,...
                                 'targetComputer',targetComputer,...
                                 'compilerOptimization',compilerOptimization,...
                                 'compileGateways',compileGateways,...
                                 'compileLibrary',compileLibrary,...
                                 'compileStandalones',compileStandalones,...
                                 'verboseLevel',verboseLevel);
      case 'client-server'
        statistics=createGateway('template',template,...
                                 'CfunctionsSource',fsfullfile(classFolder,sprintf('%s.c',classname)),...
                                 'simulinkLibrary',simulinkLibrary,'dummySimulinkIOs',true,...
                                 'folder',folder,...
                                 'className',classname,'classHelp',classhelp,...
                                 'callType','client-server','serverProgramName',serverProgramName,...
                                 'serverComputer',serverComputer,...
                                 'serverAddress',serverAddress,'port',port,...
                                 'targetComputer',targetComputer,...
                                 'compilerOptimization',compilerOptimization,...
                                 'compileGateways',compileGateways,...
                                 'compileLibrary',compileLibrary,...
                                 'compileStandalones',compileStandalones,...
                                 'verboseLevel',verboseLevel);
    end
    code.statistics.time.createGateway=etime(clock,t_createGateway);
    code.statistics.createGateway=statistics;
    
    fprintf('  done creating & compiling gateways & library (%.2f sec)\n',...
            etime(clock,t_createGateway));

    if verboseLevel>3
        for i=1:length(template)
            fprintf('    mexFunction(%d): %s\n',i,template(i).Cfunction);
            if length(template(i).inputs)>0
                for j=1:length(template(i).inputs);
                    fprintf('          input(%d):\n',j);
                    fprintf('                    type(%d): %s\n',...
                            j,template(i).inputs(j).type);
                    fprintf('                    name(%d): %s\n',...
                            j,template(i).inputs(j).name);
                    fprintf('                    size(%d): %s\n',...
                            j,index2str(template(i).inputs(j).sizes));
                end
            end
            if length(template(i).outputs)>0
                for j=1:length(template(i).outputs);
                    fprintf('          output(%d):\n',j);
                    fprintf('                    type(%d): %s\n',...
                            j,template(i).outputs(j).type);
                    fprintf('                    name(%d): %s\n',...
                            j,template(i).outputs(j).name);
                    fprintf('                    size(%d): %s\n',...
                            j,index2str(template(i).outputs(j).sizes));
                end
            end
        end
    end

    %% debug info to be passed to debugConvergenceAnalysis
    
    debugInfo.optimizationVariables=optimizationVariables;
    debugInfo.constraints=constraints;
    
    code.statistics.time.cmexCS=etime(clock,t_cmexCS);
    fprintf('done cmex2optimizeCS (%.3f sec)\n',etime(clock,t_cmexCS));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    varargout=setOutputs(nargout,params);

end
