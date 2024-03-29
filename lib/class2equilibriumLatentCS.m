function varargout=class2equilibriumLatentCS(varargin)
% To get help, type class2equilibriumLatentCS('help')
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    %% Function global help
    declareParameter(...
        'Help', {
            'Creates a matlab class for computing a Nash equilibrium'
            'of the form'
            '  P1objective(P1variables^*,P2variables^*,latentVariables^*,parameters) ='
            '        = minimize    P1objective(P1variables,P2variables^*,latentVariables,parameters)'
            '          w.r.t.      P1variables,latentVariables'
            '          subject to  P1constraints(P1variables,P2variables^*,latentVariables,parameters)'
            '                      latentConstraints(P1variables,P2variables^*,latentVariables,parameters)'
            '  P2objective(P1variables^*,P2variables^*,latentVariables^*,parameters) ='
            '        = minimize    P2objective(P1variables^*,P2variables,latentVariables,parameters)'
            '          w.r.t.      P2variables,latentVariables'
            '          subject to  P2constraints(P1variables^*,P2variables,latentVariables,parameters)'
            '                      latentConstraints(P1variables^*,P2variables,latentVariables,parameters)'
            'and returns'
            '  outputExpressions(P1variables^*,P2variables^*,latentVariables^*,parameters)'
            'See ipm.pdf for details of the optimization engine.'
            ' '
            'The solver is accessed through a matlab class.'
                });

    localVariables_=parameters4equilibrium(localVariables_);

    localVariables_=parameters4all(localVariables_);

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

    parameters=checkParameters(parameters);

    if isstruct(P1optimizationVariables)
        P1optimizationVariables=struct2cell(P1optimizationVariables);
    end

    if isstruct(P2optimizationVariables)
        P2optimizationVariables=struct2cell(P2optimizationVariables);
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
        for j=1:length(parameters)
            if isequal(name(P1optimizationVariables{i}),name(parameters{j}))
                P1optimizationVariables{i}
                error('optimization variable ''%s'' cannot also be a parameter\n',name(P1optimizationVariables{i}));
            end
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
        for j=1:length(parameters)
            if isequal(name(P2optimizationVariables{i}),name(parameters{j}))
                P2optimizationVariables{i}
                error('optimization variable ''%s'' cannot also be a parameter\n',name(P2optimizationVariables{i}));
            end
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

    [outputExpressions,outputNames]=checkOutputExpressions(outputExpressions);

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
    code.LDLthreshold=LDLthreshold;
    classhelp={'Create object';
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

    if verboseLevel>0
        fprintf('  Defining primal variables, constraints, and dual variables... ');
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

    %% Pack constraints

    if verboseLevel>1
        fprintf('Packing expressions and variables... ');
    end

    %% Pack primal variables
    % Player 1
    [u,~,~,~,P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H]...
        =packVariables(P1optimizationVariables,'u_',...
                       P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H);
    u0=packExpressions(P1optimizationVariables);
    % Player 2
    [d,~,~,~,P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H]...
        =packVariables(P2optimizationVariables,'d_',...
                       P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H);
    d0=packExpressions(P2optimizationVariables);
    src={u0,d0};
    dst={u,d};
    % Latent
    if ~isempty(latentVariables)
        [x,~,~,~,P1objective,P2objective,outputExpressions,Fu,Gu,Fd,Gd,H]...
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
    if size(Fd,1)>0
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
    Tout=ipmPDeqlat_CS(struct(...
        'code',code,...
        'P1objective',P1objective,...
        'P2objective',P2objective,...
        'u',u,...
        'd',d,...
        'x',x,...
        'P1lambda',P1lambda,...
        'P1nu',P1nu,...
        'P1xnu',P1xnu,...
        'P2lambda',P2lambda,...
        'P2nu',P2nu,...
        'P2xnu',P2xnu,...
        'Fu',Fu,...
        'Gu',Gu,...
        'Fd',Fd,...
        'Gd',Gd,...
        'H',H,...
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
        'debugConvergence',debugConvergence,...
        'profiling',profiling));
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
    classhelp{end+1}='[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter),addEye2Hessian);';
    defines.nZ=size(u,1)+size(d,1)+size(x,1);
    defines.nU=size(u,1);
    defines.nD=size(d,1);
    defines.nX=size(x,1);
    defines.nG=size(Gu,1)+size(Gd,1)+size(H,1);
    defines.nF=size(Fu,1)+size(Fd,1);
    defines.nNu=size(Gu,1)+size(Gd,1)+2*size(H,1);
    defines.gradTolerance=gradTolerance;
    defines.addEye2Hessian1tolerance=addEye2Hessian1tolerance;
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

    pth=fileparts(which('class2equilibriumLatentCS.m'));
    declareFunction(code,fsfullfile(pth,'ipmPDeq_CSsolver.m'),'solve',defines,[],[],'solve');

    %% Declare 'gets' for output expressions
    classhelp{end+1}='% Get outputs';
    classhelp{end+1}='';
    for i=1:length(outputExpressions)
        classhelp{end}=[classhelp{end},outputNames{i},','];
    end
    classhelp{end}=sprintf('[%s]=getOutputs(obj);',classhelp{end}(1:end-1));
    classhelp{end+1}=sprintf('[y (struct)]=getOutputs(obj);');
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
