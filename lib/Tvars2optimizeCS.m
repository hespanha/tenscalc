function varargout=Tvars2optimizeCS(varargin)
% To get help, type Tvars2optimizeCS('help')
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

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
    'VariableName','constraints',...
    'DefaultValue',{},...
    'Description',{
    'Cell-array of Tcalculus symbolic objects representing the'
    'constraints. Both equality and inequality constraints'
    'are allowed.'
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
    'VariableName','sensitivityVariables',...
    'DefaultValue',{},...
    'Description',{
    'Cell-array of Tcalculus symbolic objects representing the'
    'optimization variables with respect to which we want to compute cost';
    'sensitivity.'
    });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
declareOutput(...
    'VariableName','Tout',...
    'Description', {
    'Structure with several Tcalculus symbolic objects relevant to the optimization.'
    });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters and inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,params]=setParameters(nargout,varargin);
if stopNow
    return
end

classname='dummy';
scratchbookType='double';
fastRedundancyCheck='false';
outputExpressions={};
skipAffine=true;
useLDL=true;
allowSave=false;
debugConvergence=false;
umfpack=false;
scaleInequalities=false;
scaleCost=false;
scaleEqualities=false;
addEye2Hessian=false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
end

if ~isempty(size(objective))
    error('Minimization criterion must be scalar (not [%s])',index2str(size(objective)));
end

if ~isempty(constraints) && ~iscell(constraints)
    error('Constraints parameter must be a cell array\n');
end

fprintf('Tvars2optimizeCS:... ');
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

%% Define constraints, dual variables

if verboseLevel>0
    fprintf('  Defining primal variables, constraints, and dual variables... ');
end

[G,F,nus,lambdas,outputExpressions,tpl]=...
    parseConstraints(code,classname,constraints,outputExpressions);

%% Pack primal variables
[u,whereVariables,~,~,objective,outputExpressions,F,G]...
    =packVariables(optimizationVariables,'x_',objective,outputExpressions,F,G);
u0=packExpressions(optimizationVariables);

%% Pack dual variables
if size(G,1)>0
    [nu,~,~,~,outputExpressions]=packVariables(nus,'nu_',outputExpressions);
    nu0=packExpressions(nus);
else
    nu=Tzeros(0);
end
if size(F,1)>0
    [lambda,~,~,~,outputExpressions]=packVariables(lambdas,'lambda_',outputExpressions);
    lambda0=packExpressions(lambdas);
else
    lambda=Tzeros(0);
end

%% Get indices of sensitivity variables
isSensitivity=variableIndices(u,optimizationVariables,whereVariables,sensitivityVariables);

%% Generate the code for the functions that do the raw computation
t_ipmPD=clock();
Tout=ipmPD_CS(struct(...
    'code',code,...
    'f',objective,...
    'u',u,...
    'lambda',lambda,...
    'nu',nu,...
    'F',F,...
    'G',G,...
    'packOptimizationVariables',packOptimizationVariables,...
    'isSensitivity',isSensitivity,...
    'smallerNewtonMatrix',smallerNewtonMatrix,...
    'addEye2Hessian',addEye2Hessian,...
    'skipAffine',skipAffine,...
    'scaleInequalities',scaleInequalities,...
    'scaleCost',scaleCost,...
    'scaleEqualities',scaleEqualities,...
    'useLDL',useLDL,...
    'atomicFactorization',umfpack,...
    'cmexfunction',classname,...
    'allowSave',allowSave,...
    'debugConvergence',debugConvergence));
code.statistics.time.ipmPD=etime(clock,t_ipmPD);
code.statistics.time.cmexCS=etime(clock,t_cmexCS);

%% Copy variables with appropriate size to output structure
fn=fields(Tout_);
for i=1:length(fn)
    varname=sprintf('%s_',fn{i});
    Tout.(varname)=Tvariable(varname,size(Tout_.(fn{i})),true);
end

%% constraints and associated dual variables
kl=0;
kn=0;
for k=1:length(constraints)
    switch type(constraints{k})
        case 'ispositive'
            kl=kl+1;
            fld=sprintf('ispositive%d__',kl);
            Tout.(fld)=Tcalculus(operands(constraints{k}));
            fld=name(lambdas{kl});
            Tout.(fld)=lambdas{kl};
        case 'iszero'
            kn=kn+1;
            fld=sprintf('iszero%d__',kn);
            Tout.(fld)=Tcalculus(operands(constraints{k}));
            fld=name(nus{kn});
            Tout.(fld)=nus{kn};
        otherwise
            error('constraint of type ''%s'' not implemented\n',type(constraints{k}));
    end
end

fprintf('done Tvars2optimizeCS (%.3f sec)\n',etime(clock,t_cmexCS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);

end
