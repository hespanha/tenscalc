function varargout=class2optimizeCS(varargin)
% To get help, type class2optimizeCS('help')
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    %% Function global help
    declareParameter(...
        'Help', {
            'Creates a matlab class for solving optimization problems'
            'of the form:'
            '  objective(optimizationVariables^*,parameters) ='
            '       = minimum     objective(optimizationVariables,parameters)'
            '         w.r.t.      optimizationVariables'
            '         subject to  constraints(optimizationVariables,parameters)'
            'and returns'
            '  outputExpressions(optimizationVariables^*,parameters)'
            ' '
            'The solver is accessed through a matlab class.'
            'See |ipm.pdf| for details of the optimization engine.';
                });

    localVariables_=parameters4optimize(localVariables_);

    localVariables_=parameters4all(localVariables_);

    declareParameter(...
        'VariableName','useLDL',...
        'DefaultValue',true,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true| the search directions are computed using an'
            'LDL instead of an LU factorization.';
            ' '
            'In general, the LDL factorization leads to faster code.';
                      });

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Retrieve parameters and inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [stopNow,params]=setParameters(nargout,varargin);
    if stopNow
        return
    end

    if useUmfpack
        warning('umfpack does not include LDL factorization, using LU factorization instead');
        useLDL=false;
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

    fprintf('class2optimizeCS:... ');
    t_classCS=clock();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Declare the problem-specific variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t_csparse=clock();
    debug=false;
    scratchbookType='double';
    tprod2matlab=true;
    code=csparse(scratchbookType,debug,tprod2matlab); % using instructionsTable.c
    code.LDLthreshold=LDLthreshold;
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

    %% Define constraints, dual variables, and declare 'sets' for
    %% initializing primal and dual variables

    if verboseLevel>0
        fprintf('  Defining primal variables, constraints, and dual variables... ');
    end

    [G,F,nus,lambdas,outputExpressions]=...
        parseConstraints(code,classname,constraints,outputExpressions);

    %% Pack primal variables
    Joriginal=objective;
    Foriginal=F;
    Goriginal=G;
    [u,whereVariables,~,~,objective,outputExpressions,F,G]...
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

    %% Get indices of sensitivity variables
    isSensitivity=variableIndices(u,optimizationVariables,whereVariables,sensitivityVariables);
    %% Generate the code for the functions that do the raw computation
    t_ipmPD=clock();
    Tout=ipmPD_CS(struct(...
        'code',code,...
        'u',u,...         % single column vector
        'f',objective,...    % as a function of u
        'F',F,...            % as a function of u
        'G',G,...            % as a function of u
        'uList',{optimizationVariables},... % list of variables
        'f_uList',Joriginal,...              % as a function of uList
        'F_uList',Foriginal,...              % as a function of uList
        'G_uList',Goriginal,...              % as a function of uList
        'lambda',lambda,...
        'nu',nu,...
        'packOptimizationVariables',packOptimizationVariables,...
        'isSensitivity',isSensitivity,...
        'smallerNewtonMatrix',smallerNewtonMatrix,...
        'addEye2Hessian',addEye2Hessian,...
        'skipAffine',skipAffine,...
        'scaleInequalities',scaleInequalities,...
        'scaleCost',scaleCost,...
        'scaleEqualities',scaleEqualities,...
        'useLDL',useLDL,...
        'atomicFactorization',useUmfpack,...
        'cmexfunction',classname,...
        'allowSave',allowSave,...
        'debugConvergence',debugConvergence));
    code.statistics.time.ipmPD=etime(clock,t_ipmPD);

    % Replace solver variables into output expression
    fn=fields(Tout);
    for i=1:length(fn)
        varname=sprintf('%s_',fn{i});
        outputExpressions=substitute(outputExpressions,...
                                     Tvariable(varname,size(Tout.(fn{i})),true),Tout.(fn{i}));
    end

    %% Declare ipm solver
    classhelp{end+1}='% Solve optimization';
    classhelp{end+1}=...
        sprintf('[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter),addEye2Hessian);');

    defines.nU=size(u,1);
    defines.nG=size(G,1);
    defines.nF=size(F,1);
    defines.gradTolerance=gradTolerance;
    defines.addEye2HessianUtolerance=addEye2HessianUtolerance;
    defines.equalTolerance=equalTolerance;
    defines.desiredDualityGap=desiredDualityGap;
    defines.scaleCost=double(scaleCost);
    defines.scaleInequalities=double(scaleInequalities);
    defines.scaleEqualities=double(scaleEqualities);
    defines.alphaMin=alphaMin;
    defines.alphaMax=alphaMax;
    defines.coupledAlphas=double(coupledAlphas);
    defines.setAddEye2Hessian=double(addEye2Hessian~=0);
    defines.adjustAddEye2Hessian=double(adjustAddEye2Hessian~=0);
    defines.useInertia=double(useInertia~=0);
    defines.muFactorAggressive=muFactorAggressive;
    defines.muFactorConservative=muFactorConservative;
    defines.delta=delta;
    defines.skipAffine=double(skipAffine);
    defines.smallerNewtonMatrix=double(smallerNewtonMatrix);
    defines.useLDL=double(useLDL);
    defines.useUmfpack=double(useUmfpack);
    defines.allowSave=double(allowSave);
    defines.debugConvergence=double(debugConvergence);
    defines.debugConvergenceThreshold=debugConvergenceThreshold;
    defines.profiling=double(profiling);
    defines.verboseLevel=solverVerboseLevel;

    pth=fileparts(which('class2optimizeCS.m'));
    declareFunction(code,fsfullfile(pth,'ipmPD_CSsolver.m'),'solve',defines,[],[],'solve');

    %% Declare 'gets' for output expressions
    classhelp{end+1}='% Get outputs';
    classhelp{end+1}='';
    for i=1:length(outputExpressions)
        classhelp{end}=[classhelp{end},outputNames{i},','];
    end
    classhelp{end}=sprintf('[%s]=getOutputs(obj);',classhelp{end}(1:end-1));
    classhelp{end+1}=sprintf('[y (struct)]=getOutputs(obj);',classhelp{end}(1:end-1));

    declareGet(code,cell2struct(outputExpressions,outputNames),'getOutputs');

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