function varargout=TClasso(varargin);
% To get help, type TClasso('help')
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

declareParameter(...
    'Help', {
    'This script generates a tenscal optimization engine that does'
    'lasso fit of data to a affine/linear function of'
    'the form:'
    '    f(x) = c+x*w'
    '  where';
    '    c is a scalar (optional)'
    '    x is an n (row) vector';
    '    w is an n (column) vector';
    'The vector w is found using an optimization of the form'
    '  minimize sum_i (f(x_i)-y_i)^2 + l1weight sum_i |w_i|'
    });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    '   * |setP_y(obj,value)| - sets the value of y'
    '   * |setP_l1weight(obj,value)| - sets the value of l1weight'
    '   * |setV_c(obj,value)| - sets the value of c (optional)'
    '   * |setV_W(obj,value)| - sets the value of W'
    '   * |setV_absW(obj,value)| - sets upper values for |W|'
    '   * |[W,c,J]=solve(obj,mu0,int32(maxIter),int32(saveIter))|'
    'where'
    '   * |mu0|      - initial value for the barrier variable'
    '   * |maxIter|  - maximum number of Newton iterations'
    '   * |saveIter| - iteration # when to save the "hessian" matrix'
    '                  (for subsequent pivoting/permutations/scaling optimization)'
    '                  only saves when |allowSave| is true.';
    ' '
    '                  When |saveIter=0|, the hessian matrix is saved'
    '                                     at the last iteration; and'
    '                  when |saveIter<0|, the hessian matrix is not saved.';
    ' '
    '                  The "hessian" matrix will be saved regardless of the';
    '                  value of |saveIter|, when the solver exists with |status=-2|'
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
    'VariableName','dimension',...
    'DefaultValue',[],...
    'Description', {
    'Dimension of feature vectors'
    });

declareParameter(...
    'VariableName','nTrainingPoints',...
    'DefaultValue',[],...
    'Description', {
    'Number of points used for training. Needed when ''x'' is an optimization parameter.'
    });

declareParameter(...
    'VariableName','x',...
    'DefaultValue',[],...
    'Description', {
    'Vector ''x'' for the optimization.'
    ' When empty, the vector ''x'' is passed as an optimization parameter.'
    });

declareParameter(...
    'VariableName','addConstant',...
    'DefaultValue',true,...
    'AdmissibleValues',{true,false},...
    'Description', {
    'Include the constan term c in the regression function';
    '    f(x) = c+x*w'
    'For well conditioned problems it is faster to exclude the sqrt(),';
    'but for worst conditioned problems the sqrt() helps.'
    });

declareParameter(...
    'VariableName','preMultiplied',...
    'DefaultValue',false,...
    'AdmissibleValues',{true,false},...
    'Description', {;
    'When false, the inputs to the solver are'
    '  y (nTrainingPoints column vector) with dependent variables (one entry per sample)'
    '  X (nTrainingPoints x dimension matrix) with indepdent variables (one row per sample)'
    'When true, the inputs to the solver are the matrix'
    '   X  = X''*X  (dimension x dimension matrix)'
    '   yX = 2*X''*y (dimension vector)';
    'Can only be true when the vector ''x'' is passed as an optimization parameter.'
    });

declareParameter(...
    'VariableName','useSqrt',...
    'DefaultValue',false,...
    'AdmissibleValues',{true,false},...
    'Description', {
    'When ''true''a sqrt() is included in the l2-term of the cost:'
    '  minimize sqrt(sum_i (f(x_i)-y_i)^2) + l1weight sum_i |w_i|';
    ' '
    'For well conditioned problems it is faster to exclude the sqrt(),';
    'but for worst conditioned problems the may sqrt() help.'
    });

declareParameter(...
    'VariableName','smoothSqrt',...
    'DefaultValue',false,...
    'AdmissibleValues',{true,false},...
    'Description', {
    'Use a slack variable and a smooth constraint to emulate the sqrt():'
    '  minimize sqrte + l1weight sum_i |w_i|'
    '  subject to    sqrte>=0, sqrte^2>sum_i (f(x_i)-y_i)^2'
    ' ';
    'One used when useSqrt=true';
    ' ';
    'Does not seem to be advantageous'
    });

declareParameter(...
    'VariableName','solverType',...
    'AdmissibleValues',{'matlab','C','Copt'},...
    'DefaultValue','C',...
    'Description', {
    'Type of engine generated';
    '* ''matlab'' - Matlab class';
    '* ''C''      - C code (not optimized)'
    '* ''Copt''   - optimized C code'
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
    'VariableName','coupledAlphas',...
    'DefaultValue',false,...
    'AdmissibleValues',{false,true},...
    'Description',{
    'When |true| the same scalar gain is used for the primal and dual variables'
    'in Newton''s method line search.'
    });

declareParameter(...
    'VariableName','solverVerboseLevel',...
    'DefaultValue',1,...
    'Description',{;
    'See cmex2optimizeCS help'
    });

declareParameter(...
    'VariableName','addEye2Hessian',...
    'DefaultValue',true,...
    'AdmissibleValues',{false,true},...
    'Description',{
    'See cmex2optimizeCS help'
    });

declareParameter(...
    'VariableName','adjustAddEye2Hessian',...
    'DefaultValue',true,...
    'AdmissibleValues',{false,true},...
    'Description',{
    'See cmex2optimizeCS help'
    });

declareParameter(...
    'VariableName','gradTolerance',...
    'DefaultValue',1e-5,...
    'Description',{
    'See cmex2optimizeCS help'
    });

declareParameter(...
    'VariableName','equalTolerance',...
    'DefaultValue',1e-5,...
    'Description',{
    'See cmex2optimizeCS help'
    });

declareParameter(...
    'VariableName','desiredDualityGap',...
    'DefaultValue',1e-5,...
    'Description',{
    'See cmex2optimizeCS help'
    });

declareOutput(...
    'VariableName','classname',...
    'Description', {
    'Name of the class used to access the optimizatio engine.'
    'Use help(classname) to see how to use it.'
    });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,parameters__]=setParameters(nargout,varargin);
if stopNow
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('TClasso: starting...');
t0=clock;

if ~isempty(x)
    if ~isempty(nTrainingPoints) && nTrainingPoints~=size(x,1)
        error('TClasso: ''nTrainingPoints'' does not match 1st dimension of x [%dx%d]\n',size(x));
    end
    if ~isempty(dimension) && dimension~=size(x,1)
        error('TClasso: ''dimension'' does not match 1st dimension of x [%dx%d]\n',size(x));
    end
    [nTrainingPoints,dimension]=size(x);
end

%% Setup optimization

% L1-weight
Tvariable l1weight [];

% Optimization variables
Tvariable W [dimension];
optimizationVariables={W};
outputVariables={W};

constraints={};


% Training data
if preMultiplied
    Tvariable y2 [];
    Tvariable XX   [dimension,dimension];
    Tvariable yX   [dimension];

    parameters={y2;XX;yX;l1weight};
    e2=tprod(W,[-1],XX,[-1,-2],W,[-2])+y2-2*yX*W;
    if addConstant
        Tvariable c [];
        optimizationVariables{end+1}=c;
        outputVariables{end+1}=c;
        Tvariable sumY [];
        Tvariable sumX [dimension];
        parameters(end+1:end+2)={sumY;sumX};
        e2=e2-2*sumY*c+2*c*sumX*W+nTrainingPoints*norm2(c);
    end
else
    if isempty(x)
        Tvariable y [nTrainingPoints];
        Tvariable X [nTrainingPoints,dimension];
        parameters={y,X,l1weight};
    else
        Tvariable y [nTrainingPoints];
        X=x;
        parameters={y,l1weight};
    end

    e=X*W-y;
    if addConstant
        Tvariable c [];
        optimizationVariables{end+1}=c;
        outputVariables{end+1}=c;
        e=e+c;
    end
    e2=norm2(e);
end


if useSqrt
    if smoothSqrt
        Tvariable sqrte [];
        optimizationVariables{end+1}=sqrte;
        J=sqrte;
        constraints(end+1:end+2)={sqrte>=0;norm2(sqrte)>=e2};

        skipAffine=true;
        muFactorAggressive=.1;
        muFactorConservative=.99;
    else
        J=sqrt(e2);

        skipAffine=true;
        muFactorAggressive=.05;
        muFactorConservative=.99;
    end
else
    J=e2;
    skipAffine=false;
    muFactorAggressive=.1;
    muFactorConservative=.9;
end

Tvariable absW [dimension];
optimizationVariables{end+1}=absW;

constraints(end+1:end+2)={
    W<=absW;
    W>=-absW
    };

J=J+l1weight*sum(absW,1);

outputVariables{end+1}=J;

allowSave=false;
debugConvergence=false;

switch solverType
    case 'matlab'
        optimizeCS=@class2optimizeCS;
        %classname='tmpM';% matlab;
        compilerOptimization='-O0';
        allowSave=true;
        debugConvergence=false;
    case 'C'
        optimizeCS=@cmex2optimizeCS;
        %classname='tmpC';  % c
        compilerOptimization='-O0';
        allowSave=false;
        debugConvergence=false;
    case 'Cfast'
        optimizeCS=@cmex2optimizeCS;
        compilerOptimization='-Ofast';
        %classname='tmpC';  % c;
        debugConvergence=false;
    otherwise
        error('Unkown ''solverType'' = ''%s''\n',solverType);
end
%classname=sprintf('%s_%x_%x',classname,nTrainingPoints,dimension); %toolong

[classname,code]=optimizeCS(...
    'classname',classname,...
    ...%'pedigreeClass',classname,...
    ...%'folder',folder,...
    'executeScript','asneeded',...
    'objective',J,...
    'optimizationVariables',optimizationVariables,...
    'constraints',constraints,...
    'outputExpressions',outputVariables,...
    'parameters',parameters,...
    'muFactorConservative',muFactorConservative,...
    'muFactorAggressive',muFactorAggressive,...
    'gradTolerance',gradTolerance,...
    'equalTolerance',equalTolerance,...
    'desiredDualityGap',desiredDualityGap,...
    'skipAffine',skipAffine,...
    'useLDL',useLDL,...
    'smallerNewtonMatrix',smallerNewtonMatrix,...
    'coupledAlphas',coupledAlphas,...
    'debugConvergence',debugConvergence,...
    'debugConvergenceThreshold',1e4,...
    'profiling',false,...
    'addEye2Hessian',addEye2Hessian,...
    'scratchbookType','double',...
    'compilerOptimization',compilerOptimization,...
    'allowSave',allowSave,...
    'solverVerboseLevel',solverVerboseLevel);
%profile viewer

if isa(classname,'outputWithPedigree')
    classname=getValue(classname);
end

fprintf('done TClasso (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,parameters__);

end
