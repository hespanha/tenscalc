function varargout=cmex2minmaxCS(varargin)
% To get help, type cmex2minmaxCS('help')
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha, Raphael Chinchilla).  All rights reserved.

    %% Function global help
    declareParameter(...
        'Help', {
            'Creates a matlab class for solving a min-max optimization'
            '(Stackelberg equilibrium) of the form'
            '  objective(minVariables^*,maxVariables^*,parameters) ='
            '        = minimize                                             maximize    objective(minVariables,maxVariables,parameters)'
            '          w.r.t.      minVariables                             w.r.t.       maxVariables'
            '          subject to  minConstraints(minVariables,parameters)  siubject to  maxConstraints(minVariables,maxVariables,parameters)'
            'and returns'
            '  outputExpressions(minVariables^*,maxVariables^*,parameters)'
            ' '
            'See TO_BE_FIXED for details of the optimization engine.'
            ' '
            'The solver is accessed through several cmex functions that can be'
            'accessed directly or through a matlab class.'
                });

    localVariables_=parameters4minmax(localVariables_);

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
        fprintf('cmex2minmaxCS: outputs folder ''%s'' does not exist, creating it\n',folder);
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

    if isstruct(minOptimizationVariables)
        minOptimizationVariables=struct2cell(minOptimizationVariables);
    end

    if isstruct(maxOptimizationVariables)
        maxOptimizationVariables=struct2cell(maxOptimizationVariables);
    end

    if ~iscell(minOptimizationVariables)
        minOptimizationVariables
        error('minOptimizationVariables must be a cell array of Tcalculus variables');
    end

    for i=1:length(minOptimizationVariables)
        if ~isequal(class(minOptimizationVariables{i}),'Tcalculus')
            minOptimizationVariables{i}
            error('all minOptimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,class(minOptimizationVariables{i}));
        end
        if ~isequal(type(minOptimizationVariables{i}),'variable')
            minOptimizationVariables{i}
            error('all minOptimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(minOptimizationVariables{i}));
        end
        for j=1:length(parameters)
            if isequal(name(minOptimizationVariables{i}),name(parameters{j}))
                minOptimizationVariables{i}
                error('optimization variable ''%s'' cannot also be a parameter\n',name(minOptimizationVariables{i}));
            end
        end
    end

    if ~iscell(maxOptimizationVariables)
        maxOptimizationVariables
        error('maxOptimizationVariables must be a cell array of Tcalculus variables');
    end

    for i=1:length(maxOptimizationVariables)
        if ~isequal(class(maxOptimizationVariables{i}),'Tcalculus')
            minOptimizationVariables{i}
            error('all minOptimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,class(maxOptimizationVariables{i}));
        end
        if ~isequal(type(maxOptimizationVariables{i}),'variable')
            maxOptimizationVariables{i}
            error('all maxOptimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(maxOptimizationVariables{i}));
        end
        for j=1:length(parameters)
            if isequal(name(maxOptimizationVariables{i}),name(parameters{j}))
                maxOptimizationVariables{i}
                error('optimization variable ''%s'' cannot also be a parameter\n',name(maxOptimizationVariables{i}));
            end
        end
    end

    if ~isempty(size(objective))
        error('optimization criterion must be scalar (not [%s])',...
              index2str(size(objective)));
    end

    if ~isempty(minConstraints) && ~iscell(minConstraints)
        error('minimizer''s constraints parameter must be a cell array\n');
    end

    if ~isempty(maxConstraints) && ~iscell(maxConstraints)
        error('maximizer''s constraints parameter must be a cell array\n');
    end

    [outputExpressions,outputNames]=checkOutputExpressions(outputExpressions);

    fprintf('cmex2minmaxCS: ...');
    t_cmexCS=clock();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Declare the problem-specific variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf(' declaring sets for parameters and primal variables... ');

    t_csparse=clock();
    debug=false;
    tprod2matlab=false;
    code=csparse(scratchbookType,debug,tprod2matlab,fastRedundancyCheck); % using instructionsTable.c
    code.LDLthreshold=LDLthreshold;
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
    % Minimizer
    for i=1:length(minOptimizationVariables)
        template(end+1,1).MEXfunction=sprintf('%s_set_%s',...
                                            classname,name(minOptimizationVariables{i}));
        template(end).Cfunction=sprintf('%s_set_%s',classname,name(minOptimizationVariables{i}));
        template(end).method=sprintf('setV_%s',name(minOptimizationVariables{i}));
        template(end).inputs =struct('type','double',...
                                     'name',name(minOptimizationVariables{i}),...
                                     'sizes',size(minOptimizationVariables{i}));
        declareSet(code,minOptimizationVariables{i},template(end).MEXfunction);
        msize=size(minOptimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                 name(minOptimizationVariables{i}),index2str(msize));
    end
    % Maximizer
    for i=1:length(maxOptimizationVariables)
        template(end+1,1).MEXfunction=sprintf('%s_set_%s',...
                                            classname,name(maxOptimizationVariables{i}));
        template(end).Cfunction=sprintf('%s_set_%s',classname,name(maxOptimizationVariables{i}));
        template(end).method=sprintf('setV_%s',name(maxOptimizationVariables{i}));
        template(end).inputs =struct('type','double',...
                                     'name',name(maxOptimizationVariables{i}),...
                                     'sizes',size(maxOptimizationVariables{i}));
        declareSet(code,maxOptimizationVariables{i},template(end).MEXfunction);
        msize=size(maxOptimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                 name(maxOptimizationVariables{i}),index2str(msize));
    end

    %% Define constraints, dual variables, and declare 'sets' for initializing dual variables

    if verboseLevel>0
        fprintf('  Defining primal variables, constraints, and dual variables... ');
    end

    % Minimizer
    [Gu,Fu,minNus,minLambdas,outputExpressions,tpl]=...
        parseConstraints(code,classname,minConstraints,outputExpressions,'minimizer');
    template=[template;tpl];


    % Maximizer
    [Gd,Fd,maxNus,maxLambdas,outputExpressions,tpl]=...
        parseConstraints(code,classname,maxConstraints,outputExpressions,'maximizer');
    template=[template;tpl];

    %% Pack constraints

    if verboseLevel>1
        fprintf('Packing expressions and variables... ');
    end

    %% Pack primal variables
    % Minimizer
    [u,~,~,~,objective,outputExpressions,Fu,Gu,Fd,Gd]...
        =packVariables(minOptimizationVariables,'u_',...
                       objective,outputExpressions,Fu,Gu,Fd,Gd);
    u0=packExpressions(minOptimizationVariables);
    % Maximizer
    [d,~,~,~,objective,outputExpressions,Fu,Gu,Fd,Gd]...
        =packVariables(maxOptimizationVariables,'d_',...
                       objective,outputExpressions,Fu,Gu,Fd,Gd);
    d0=packExpressions(maxOptimizationVariables);
    src={u0,d0};
    dst={u,d};

    declareCopy(code,dst,src,'initPrimal__');

    %% Pack dual variables
    % Minimizer
    if size(Gu,1)>0
        minNu=packVariables(minNus,'minNu_');
        minNu0=packExpressions(minNus);
        src{end+1}=minNu0;
        dst{end+1}=minNu;
    else
        minNu=Tzeros(0);
    end
    if size(Fu,1)>0
        minLambda=packVariables(minLambdas,'minLambda_');
        minLambda0=packExpressions(minLambdas);
        src{end+1}=minLambda0;
        dst{end+1}=minLambda;
    else
        minLambda=Tzeros(0);
    end
    % Maximizer
    if size(Gd,1)>0
        maxNu=packVariables(maxNus,'maxNu_');
        maxNu0=packExpressions(maxNus);
        src{end+1}=maxNu0;
        dst{end+1}=maxNu;
    else
        maxNu=Tzeros(0);
    end
    if size(Fd,1)>0
        maxLambda=packVariables(maxLambdas,'maxLambda_');
        maxLambda0=packExpressions(maxLambdas);
        src{end+1}=maxLambda0;
        dst{end+1}=maxLambda;
    else
        maxLambda=Tzeros(0);
    end
    declareCopy(code,dst,src,'initPrimalDual__');

    %% Generate the code for the functions that do the raw computation
    t_ipmPD=clock();
    Tout=ipmPDminmax_CS(struct(...
        'code',code,...
        'objective',objective,...
        'u',u,...
        'd',d,...
        'minLambda',minLambda,...
        'minNu',minNu,...
        'maxLambda',maxLambda,...
        'maxNu',maxNu,...
        'Fu',Fu,...
        'Gu',Gu,...
        'Fd',Fd,...
        'Gd',Gd,...
        'addEye2Hessian',addEye2Hessian,...
        'scaleInequalities',scaleInequalities,...
        'scaleCost',scaleCost,...
        'scaleEqualities',scaleEqualities,...
        'atomicFactorization',useUmfpack,...
        'cmexfunction',classname,...
        'allowSave',allowSave,...
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
    classhelp{end+1}='Solve optimization';
    classhelp{end+1}='[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter),addEye2Hessian);';
    template(end+1,1).MEXfunction=sprintf('%s_solve',classname);
    template(end).Cfunction='ipmPDminmax_CSsolver';
    template(end).method='solve';
    template(end).inputs(1) =struct('type','double','name','mu0','sizes',1,'default',1);
    template(end).inputs(2) =struct('type','int32','name','maxIter','sizes',1,'default',200);
    template(end).inputs(3) =struct('type','int32','name','saveIter','sizes',1,'default',-1);
    template(end).inputs(4) =struct('type','double','name','addEye2Hessian','sizes',2,'default',[1e-9;1e-9]);
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
    defines.nZ=size(u,1)+size(d,1);
    defines.nU=size(u,1);
    defines.nD=size(d,1);
    defines.nG=size(Gu,1)+size(Gd,1);
    defines.nF=size(Fu,1)+size(Fd,1);
    defines.nGu=size(Gu,1);
    defines.nFu=size(Fu,1);
    defines.nGd=size(Gd,1);
    defines.nFd=size(Fd,1);
    defines.nNu=size(Gu,1)+size(Gd,1);
    defines.gradTolerance=sprintf('%e',gradTolerance); % to make double
    defines.addEye2HessianUtolerance=sprintf('%e',addEye2HessianUtolerance); % to make double
    defines.addEye2HessianDtolerance=sprintf('%e',addEye2HessianDtolerance); % to make double
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
    defines.profiling=double(profiling);
    defines.verboseLevel=solverVerboseLevel;

    pth=fileparts(which('cmex2minmaxCS.m'));
    declareFunction(code,fsfullfile(pth,'ipmPDminmax_CSsolver.c'),'ipmPDminmax_CSsolver',...
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

    code.statistics.time.cmexCS=etime(clock,t_cmexCS);
    fprintf('done cmex2minmaxCS (%.3f sec)\n',etime(clock,t_cmexCS));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    varargout=setOutputs(nargout,params);

end
