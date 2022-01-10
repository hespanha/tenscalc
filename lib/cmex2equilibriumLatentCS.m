function varargout=cmex2equilibriumLatentCS(varargin)
% To get help, type cmex2equilibriumLatentCS('help')
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    %% Function global help
    declareParameter(...
        'Help', {
            'Creates a set of cmex functions for computing a Nash equilibrium'
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
            'The solver is accessed through several cmex functions that can be'
            'accessed directly or through a matlab class.'
            'See |ipm.pdf| for details of the optimization engine.';
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
        fprintf('cmex2equilibriumLatentCS: outputs folder ''%s'' does not exist, creating it\n',folder);
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

    fprintf('cmex2equilibriumLatentCS: ...');
    t_cmexCS=clock();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Declare the problem-specific variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf(' declaring sets for parameters and primal variables... ');

    t_csparse=clock();
    debug=false;
    tprod2matlab=false;
    %code=csparse0(scratchbookType,debug);   % using fastTable.m
    %code=csparse1(scratchbookType,debug);  % using fastTable.m, string I_ instruction types
    %code=csparse2(scratchbookType,debug);  % using fastTable.m, integer instruction types
    code=csparse(scratchbookType,debug,tprod2matlab,fastRedundancyCheck); % using instructionsTable.c
    code.LDLthreshold=LDLthreshold;
    classhelp={'% Create object';
               sprintf('obj=%s();',classname)};

    % template for createGateway
    template=cmextoolsTemplate();
    %% Declare 'sets' for initializing parameters
    if length(parameters)>0
        classhelp{end+1}='Set parameters';
    end
    for i=1:length(parameters)
        template(end+1,1).MEXfunction=sprintf('%s_set_%s',classname,name(parameters{i}));
        template(end).Cfunction=sprintf('%s_set_%s',classname,name(parameters{i}));
        template(end).method=sprintf('setP_%s',name(parameters{i}));
        template(end).inputs =struct('type','double',...
                                     'name',name(parameters{i}),...
                                     'sizes',size(parameters{i}));
        declareSet(code,parameters{i},template(end).MEXfunction);
        msize=size(parameters{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setP_%s(obj,{[%s] matrix});',...
                                name(parameters{i}),index2str(msize));
    end

    %% Declare 'sets' for initializing primal variables
    classhelp{end+1}='Initialize primal variables';
    % Player 1
    for i=1:length(P1optimizationVariables)
        template(end+1,1).MEXfunction=sprintf('%s_set_%s',...
                                            classname,name(P1optimizationVariables{i}));
        template(end).Cfunction=sprintf('%s_set_%s',classname,name(P1optimizationVariables{i}));
        template(end).method=sprintf('setV_%s',name(P1optimizationVariables{i}));
        template(end).inputs =struct('type','double',...
                                     'name',name(P1optimizationVariables{i}),...
                                     'sizes',size(P1optimizationVariables{i}));
        declareSet(code,P1optimizationVariables{i},template(end).MEXfunction);
        msize=size(P1optimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                 name(P1optimizationVariables{i}),index2str(msize));
    end
    % Player 2
    for i=1:length(P2optimizationVariables)
        template(end+1,1).MEXfunction=sprintf('%s_set_%s',...
                                            classname,name(P2optimizationVariables{i}));
        template(end).Cfunction=sprintf('%s_set_%s',classname,name(P2optimizationVariables{i}));
        template(end).method=sprintf('setV_%s',name(P2optimizationVariables{i}));
        template(end).inputs =struct('type','double',...
                                     'name',name(P2optimizationVariables{i}),...
                                     'sizes',size(P2optimizationVariables{i}));
        declareSet(code,P2optimizationVariables{i},template(end).MEXfunction);
        msize=size(P2optimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                 name(P2optimizationVariables{i}),index2str(msize));
    end
    % Latent variables
    for i=1:length(latentVariables)
        template(end+1,1).MEXfunction=sprintf('%s_set_%s',classname,name(latentVariables{i}));
        template(end).Cfunction=sprintf('%s_set_%s',classname,name(latentVariables{i}));
        template(end).method=sprintf('setV_%s',name(latentVariables{i}));
        template(end).inputs =struct('type','double',...
                                     'name',name(latentVariables{i}),...
                                     'sizes',size(latentVariables{i}));
        declareSet(code,latentVariables{i},template(end).MEXfunction);
        msize=size(latentVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                 name(latentVariables{i}),index2str(msize));
    end

    %% Define constraints, dual variables, and declare 'sets' for initializing dual variables

    if verboseLevel>0
        fprintf('  Defining primal variables, constraints, and dual variables... ');
    end

    % Player 1
    [Gu,Fu,P1nus,P1lambdas,outputExpressions,tpl]=...
        parseConstraints(code,classname,P1constraints,outputExpressions,'P1');
    template=[template;tpl];


    % Player 2
    [Gd,Fd,P2nus,P2lambdas,outputExpressions,tpl]=...
        parseConstraints(code,classname,P2constraints,outputExpressions,'P2');
    template=[template;tpl];

    % Latent constraints (only one H, but two sets of dual variables)
    [H,~,P1xnus,err,outputExpressions,tpl]=...
        parseConstraints(code,classname,latentConstraints,outputExpressions,'P1x');
    template=[template;tpl];
    [H,~,P2xnus,err,outputExpressions,tpl]=...
        parseConstraints(code,classname,latentConstraints,outputExpressions,'P2x');
    template=[template;tpl];

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
    % Hess__=ipmPDeqlat_CS(code,P1objective,P2objective,u,d,x,P1lambda,P1nu,P1xnu,P2lambda,P2nu,P2xnu,...
    %                      Fu,Gu,Fd,Gd,H,...
    %                      smallerNewtonMatrix,addEye2Hessian,skipAffine,...
    %                      scaleInequalities,scaleCost,scaleEqualities,...
    %                      useUmfpack,...
    %                      classname,allowSave,debugConvergence,profiling);
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
    nZ=size(u,1)+size(d,1)+size(x,1);
    nG=size(Gu,1)+size(Gd,1)+size(H,1);
    nF=size(Fu,1)+size(Fd,1);
    nNu=size(Gu,1)+size(Gd,1)+2*size(H,1);

    classhelp{end+1}='Solve optimization';
    classhelp{end+1}='[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter),addEye2Hessian);';
    template(end+1,1).MEXfunction=sprintf('%s_solve',classname);
    template(end).Cfunction='ipmPDeq_CSsolver';
    template(end).method='solve';
    template(end).inputs(1) =struct('type','double','name','mu0','sizes',1,'default',1);
    template(end).inputs(2) =struct('type','int32','name','maxIter','sizes',1,'default',200);
    template(end).inputs(3) =struct('type','int32','name','saveIter','sizes',1,'default',-1);
    template(end).inputs(4) =struct('type','double','name','addEye2Hessian','sizes',2,'default',[1e-9;1e-9]);
    template(end).outputs(1)=struct('type','int32','name','status','sizes',1);
    template(end).outputs(2)=struct('type','int32','name','iter','sizes',1);
    template(end).outputs(3)=struct('type','double','name','time','sizes',1);
    % template(end).outputs(4)=struct('type','double','name','z','sizes',[nZ,maxIter+1]);
    % template(end).outputs(5)=struct('type','double','name','nu','sizes',[nNu,maxIter+1]);
    % template(end).outputs(6)=struct('type','double','name','lambda','sizes',[nF,maxIter+1]);
    % template(end).outputs(7)=struct('type','double','name','dZ_s','sizes',[nZ,maxIter]);
    % template(end).outputs(8)=struct('type','double','name','dNu_s','sizes',[nNu,maxIter]);
    % template(end).outputs(9)=struct('type','double','name','dLambda_s','sizes',[nF,maxIter]);
    % template(end).outputs(10)=struct('type','double','name','G','sizes',[nG,maxIter+1]);
    % template(end).outputs(11)=struct('type','double','name','F','sizes',[nF,maxIter+1]);
    % template(end).outputs(12)=struct('type','double','name','primaAlpha_s','sizes',[maxIter+1]);
    % template(end).outputs(13)=struct('type','double','name','dualAlpha_s','sizes',[maxIter+1]);
    % template(end).outputs(14)=struct('type','double','name','finalAlpha','sizes',[maxIter+1]);

    %folder='.';
    classFolder=fsfullfile(folder,sprintf('@%s',classname));
    if ~exist(classFolder,'dir')
        fprintf('class classFolder @%s does not exist, creating it... ',classname);
        mkdir(classFolder);
    end

    defines.saveNamePrefix=['"',fsfullfile(classFolder,classname),'"'];
    defines.nZ=size(u,1)+size(d,1)+size(x,1);
    defines.nU=size(u,1);
    defines.nD=size(d,1);
    defines.nX=size(x,1);
    defines.nG=size(Gu,1)+size(Gd,1)+size(H,1);
    defines.nF=size(Fu,1)+size(Fd,1);
    defines.nNu=size(Gu,1)+size(Gd,1)+2*size(H,1);
    defines.gradTolerance=sprintf('%e',gradTolerance); % to make double
    defines.addEye2Hessian1tolerance=sprintf('%e',addEye2Hessian1tolerance); % to make double
    defines.equalTolerance=sprintf('%e',equalTolerance); % to make double
    defines.desiredDualityGap=sprintf('%e',desiredDualityGap); % to make double
    defines.scaleCost=double(scaleCost~=0);
    defines.scaleInequalities=double(scaleInequalities);
    defines.scaleEqualities=double(scaleEqualities);
    defines.alphaMin=sprintf('%e',alphaMin); % to make double
    defines.alphaMax=sprintf('%e',alphaMax); % to make double
    defines.coupledAlphas=double(coupledAlphas);
    defines.setAddEye2Hessian=double(addEye2Hessian~=0);
    defines.adjustAddEye2Hessian=double(adjustAddEye2Hessian~=0);
    defines.muFactorAggressive=sprintf('%e',muFactorAggressive); % to make double
    defines.muFactorConservative=sprintf('%e',muFactorConservative); % to make double
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

    pth=fileparts(which('cmex2equilibriumLatentCS.m'));
    declareFunction(code,fsfullfile(pth,'ipmPDeq_CSsolver.c'),'ipmPDeq_CSsolver',...
                    defines,template(end).inputs,template(end).outputs,template(end).method);

    %% Declare 'gets' for output expressions
    classhelp{end+1}='Get outputs';
    classhelp{end+1}='';
    template(end+1,1).MEXfunction=sprintf('%s_getOutputs',classname);
    template(end).Cfunction=sprintf('%s_getOutputs',classname);
    template(end).method='getOutputs';
    template(end).outputs=struct('type',{},'name',{},'sizes',{});
    for i=1:length(outputExpressions)
        template(end).outputs(i).type='double';
        template(end).outputs(i).name=outputNames{i};
        template(end).outputs(i).sizes=size(outputExpressions{i});% will be overwriten after compile2C
        classhelp{end}=[classhelp{end},outputNames{i},','];
    end
    classhelp{end}=sprintf('[%s]=getOutputs(obj);',classhelp{end}(1:end-1));
    classhelp{end+1}=sprintf('[y (struct)]=getOutputs(obj);');
    declareGet(code,outputExpressions,template(end).MEXfunction);

    code.statistics.time.csparse=etime(clock,t_csparse);
    code.statistics.defines=defines;

    fprintf('  done creating csparse object (%.3f sec)\n',etime(clock,t_csparse));

    %% Compile code
    fprintf(' creating C code... ');
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
            fprintf('    mexFunction(%d): %s\n',i,template(i).MEXfunction);
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

    debugInfo.P1optimizationVariables=P1optimizationVariables;
    debugInfo.P2optimizationVariables=P2optimizationVariables;
    debugInfo.P1constraints=P1constraints;
    debugInfo.P2constraints=P2constraints;

    code.statistics.time.cmexCS=etime(clock,t_cmexCS);
    fprintf('done cmex2equilibriumLatentCS (%.3f sec)\n',etime(clock,t_cmexCS));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    varargout=setOutputs(nargout,params);

end
