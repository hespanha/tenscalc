function varargout=class2minmaxCS(varargin)
% To get help, type class2minmaxCS('help')
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
            'The solver is accessed through a matlab class.'
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
        fprintf('class2minmaxCS: outputs folder ''%s'' does not exist, creating it\n',folder);
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

    fprintf('class2minmaxCS: ...');
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
        classhelp{end+1}='Set parameters';
    end
    for i=1:length(parameters)
        declareSet(code,parameters{i},sprintf('setP_%s',name(parameters{i})));
        msize=size(parameters{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setP_%s(obj,{[%s] matrix});',...
                                name(parameters{i}),index2str(msize));
    end

    %% Declare 'sets' for initializing primal variables
    classhelp{end+1}='Initialize primal variables';
    % Minimizer
    for i=1:length(minOptimizationVariables)
        declareSet(code,minOptimizationVariables{i},...
                   sprintf('setV_%s',name(minOptimizationVariables{i})));
        msize=size(minOptimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                 name(minOptimizationVariables{i}),index2str(msize));
    end
    % Maximizer
    for i=1:length(maxOptimizationVariables)
        declareSet(code,maxOptimizationVariables{i},...
                   sprintf('setV_%s',name(maxOptimizationVariables{i})));
        msize=size(maxOptimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                 name(maxOptimizationVariables{i}),index2str(msize));
    end

    %% Define constraints, dual variables, and declare 'sets' for initializing dual variables

    if verboseLevel>0
        fprintf('  Defining primal variables, constraints, and dual variables... ');
    end

    % Minimizer
    [Gu,Fu,nusU,lambdasU,outputExpressions]=...
        parseConstraints(code,classname,minConstraints,outputExpressions,'minimizer');

    % Maximizer
    [Gd,Fd,nusD,lambdasD,outputExpressions]=...
        parseConstraints(code,classname,maxConstraints,outputExpressions,'maximizer');

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
        nuU=packVariables(nusU,'nuU_');
        nuU0=packExpressions(nusU);
        src{end+1}=nuU0;
        dst{end+1}=nuU;
    else
        nuU=Tzeros(0);
    end
    if size(Fu,1)>0
        lambdaU=packVariables(lambdasU,'lambdaU_');
        lambdaU0=packExpressions(lambdasU);
        src{end+1}=lambdaU0;
        dst{end+1}=lambdaU;
    else
        lambdaU=Tzeros(0);
    end
    % Maximizer
    if size(Gd,1)>0
        nuD=packVariables(nusD,'nuD_');
        nuD0=packExpressions(nusD);
        src{end+1}=nuD0;
        dst{end+1}=nuD;
    else
        nuD=Tzeros(0);
    end
    if size(Fd,1)>0
        lambdaD=packVariables(lambdasD,'lambdaD_');
        lambdaD0=packExpressions(lambdasD);
        src{end+1}=lambdaD0;
        dst{end+1}=lambdaD;
    else
        lambdaD=Tzeros(0);
    end
    declareCopy(code,dst,src,'initPrimalDual__');

    %% Generate the code for the functions that do the raw computation
    t_ipmPD=clock();
    Tout=ipmPDminmax_CS(struct(...
        'code',code,...
        'objective',objective,...
        'u',u,...
        'd',d,...
        'Fu',Fu,...
        'Gu',Gu,...
        'Fd',Fd,...
        'Gd',Gd,...
        'lambdaU',lambdaU,...
        'nuU',nuU,...
        'lambdaD',lambdaD,...
        'nuD',nuD,...
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
    defines.nZ=size(u,1)+size(d,1);
    defines.nU=size(u,1);
    defines.nD=size(d,1);
    defines.nG=size(Gu,1)+size(Gd,1);
    defines.nF=size(Fu,1)+size(Fd,1);
    defines.nNu=size(Gu,1)+size(Gd,1);
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
    defines.useUmfpack=double(useUmfpack);
    defines.allowSave=double(allowSave);
    defines.profiling=double(profiling);
    defines.verboseLevel=solverVerboseLevel;

    pth=fileparts(which('class2minmaxCS.m'));
    declareFunction(code,fsfullfile(pth,'ipmPDminmax_CSsolver.m'),'solve',defines,[],[],'solve');

    %% Declare 'gets' for output expressions
    classhelp{end+1}='Get outputs';
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

    code.statistics.time.cmexCS=etime(clock,t_classCS);
    fprintf('done class2minmaxCS (%.3f sec)\n',etime(clock,t_classCS));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    varargout=setOutputs(nargout,params);

end
