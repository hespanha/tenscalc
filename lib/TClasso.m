function varargout=TClasso(varargin);
% To get help, type TClasso('help')
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
            '   * |setP_{parameter}(obj,value)|'
            '                   - sets the value of one of the parameters'
            '   * |setV_{variable}(obj,value)|'
            '                   - sets the value of one of the optimization variables'
            '   * |[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter))|'
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
    'Description', {
        'Dimension of feature vectors'
                   });
    
declareParameter(...
    'VariableName','nTrainingPoints',...
    'Description', {
        'Number of points used for training.'
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
        '  X (nTrainingPoints x dimensions matrix) with indepdent variables (one row per sample)'
        'When true, the inputs to the solver are the matrix'
        '   X  = X''*X  (dimensions x dimensions matrix)'
        '   yX = 2*X''*y (dimensions vector)'
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
    'AdmissibleValues',{'matlab','C','Cfast'},...
    'DefaultValue','C',...
    'Description', {
        'Type of engine generated';
        '* ''matlab'' - Matlab class';
        '* ''C''      - C code (not optimized)'
        '* ''Cfast''  - optimized C code'
                   });
    
declareParameter(...
    'VariableName','solverVerboseLevel',...
    'DefaultValue',1,...
    'Description',{;
                   'See cmex2optimizeCS help'
                  });

declareParameter(...
    'VariableName','addEye2Hessian',...
    'DefaultValue',1e-6,...
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
    Tvariable y [nTrainingPoints];
    Tvariable X [nTrainingPoints,dimension];

    parameters={y,X,l1weight};
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
    classname='tmpM';% matlab;
    compilerOptimization='-O0';
    allowSave=true;
    debugConvergence=false;
  case 'C'
    optimizeCS=@cmex2optimizeCS;
    classname='tmpC';  % c
    compilerOptimization='-O0';
    allowSave=true;
    debugConvergence=false;
  case 'Cfast'
    optimizeCS=@cmex2optimizeCS;
    compilerOptimization='-Ofast';
    classname='tmpC';  % c;
    debugConvergence=false;
  otherwise
    error('Unkown ''solverType'' = ''%s''\n',solverType);
end
classname=sprintf('%s_%d_%d',classname,nTrainingPoints,dimension);

[classname,code]=optimizeCS(...
    'pedigreeClass',classname,...
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
    'debugConvergence',debugConvergence,...
    'debugConvergenceThreshold',1e4,...
    'profiling',false,...
    'addEye2Hessian',addEye2Hessian,...
    'scratchbookType','double',...
    'compilerOptimization',compilerOptimization,...
    'allowSave',allowSave,...
    'solverVerboseLevel',solverVerboseLevel);
%profile viewer

classname=getValue(classname);

fprintf('done TClasso (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,parameters__);


