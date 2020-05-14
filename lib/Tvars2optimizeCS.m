function varargout=Tvars2optimizeCS(varargin)	
% To get help, type Tvars2optimizeCS('help')
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
        'VariableName','constraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints. Both equality and inequality constraints'
            'are allowed.'
                      });
    
    declareParameter(...
        'VariableName','addEye2Hessian',...
        'DefaultValue',1e-9,...
        'Description',{
            'Add to the Hessian matrix appropriate identity matrices scaled by this constant.';
            'A larger value for |addEye2Hessian| has two main effects:'
            '1) Improves the numerical conditioning of the system of equations that'
            '   finds the Newton search direction.';
            '2) Moves the search direction towards the gradient descent of';
            '   the Lagragian (and away from the Newton direction).';
            'Both effects improve the robustness of the solver, but this is typically';
            'achieved at the expense of slower convergence.'
            'For convex problems, one typically chooses |addEye2Hessian| equal to the';
            'square root of the machine percision.'
            'For non-convex problems, one can try to increase this parameter when';
            'the Newton direction actually causes and increase of the Lagragian.'
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
    [u,whereVariables,packU,unpackU,objective,outputExpressions,F,G]...
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
    Tout_=ipmPD_CS(code,objective,u,lambda,nu,F,G,isSensitivity,...
                  smallerNewtonMatrix,addEye2Hessian,skipAffine,...
                  scaleInequalities,scaleCost,scaleEqualities,...
                  useLDL,umfpack,...
                  classname,allowSave,debugConvergence);
    code.statistics.time.ipmPD=etime(clock,t_ipmPD);
    code.statistics.time.cmexCS=etime(clock,t_cmexCS);
    fprintf('done Tvars2optimizeCS (%.3f sec)\n',etime(clock,t_cmexCS));

    fn=fields(Tout_);
    for i=1:length(fn)
        varname=sprintf('%s_',fn{i});
        Tout.(varname)=Tvariable(varname,size(Tout_.(fn{i})));
    end

    for i=1:length(constraints)
        switch type(constraints{i})
          case 'ispositive'
            fld=sprintf('ispositive%d',i);
            Tout.(fld)=Tcalculus(operands(constraints{i}));
          case 'iszero'
            fld=sprintf('iszero%d',i);
            Tout.(fld)=Tcalculus(operands(constraints{i}));
          otherwise
            Tout.(fld)=constraints{i};
        end
    end
    %Tout.constraints=constraints;
    for i=1:length(lambdas)
        fld=name(lambdas{i});
        Tout.(fld)=lambdas{i};
    end
    for i=1:length(nus)
        fld=name(nus{i});
        Tout.(fld)=nus{i};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    varargout=setOutputs(nargout,params);

end
