function varargout=class2equilibriumLatentCS(varargin)	
% To get help, type class2equilibriumLatentCS('help')
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
            'Creates a matlab class for computing a Nash equilibrium'
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
            'The solver is accessed through a matlab class.'
                });

    %% Declare input parameters

    declareParameter(...
        'VariableName','P1objective',...
        'Description',{
            'Scalar Tcalculus symbolic object to be optimized for player 1.'
                      });
    declareParameter(...
        'VariableName','P1optimizationVariables',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'variables to be optimized by player 1.'
                      });
    declareParameter(...
        'VariableName','P1constraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints for player 1. Both equality and inequality'
            'constraints are allowed.'
                      });

    declareParameter(...
        'VariableName','P2objective',...
        'Description',{
            'Scalar Tcalculus symbolic object to be optimized for player 2.'
                      });
    declareParameter(...
        'VariableName','P2optimizationVariables',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'variables to be optimized by player 2.'
                      });
    declareParameter(...
        'VariableName','P2constraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints for player 2. Both equality and inequality'
            'constraints are acceptable.'
                      });

    declareParameter(...
        'VariableName','latentVariables',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'common latent variables.'
                      });
    declareParameter(...
        'VariableName','latentConstraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints that define the latent variables.'
            'Only equality constraints are acceptable.'
                      });

    declareParameter(...
        'VariableName','outputExpressions',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the '
            'variables to be returned.'
            ' ';
            'The following Tcalculus symbolic variables are assigned special values';
            'and can be using in outputExpressions';
            '* |P1lambda1_|,|P1lambda2_|,... |P2lambda1_|,|P2lambda2_|,...';
            '      - Lagrangian multipliers associated with the inequalities constraints';
            '        for player 1 and 2 (in the order that they appear and with the same size';
            '        as the corresponding constraints).';
            '* |P1nu1_|,|P1nu2_|,... ,|P2nu1_|,|P2nu2_|,...'
            '* |P1xnu1_|,|P1xnu2_|,... ,|P2xnu1_|,|P2xnu2_|,...'    
            '      - Lagrangian multipliers associated with the equality constraints';
            '        for player 1 and 2 (in the order that they appear and with the same size'
            '        as the corresponding constraints). The P1x and P2x variables correspond';
            '        to the latentConstraints.'
            '* |Hess_| - Hessian matrix used by the (last) Newton step to update';
            '            the primal variables (not including |addEye2Hessian|).'
            ' ';
            'ATTENTION: To be able to include these variables as input parameters,';
            '           they have to be previously created outside *with the appropriate sizes*.'
            '           Eventually, their values will be overridden by the solver'
            '           to reflect the values listted above.'
                      });

    localVariables_=parameters4all(localVariables_);

    declareParameter(...
        'VariableName','addEye2Hessian',...
        'DefaultValue',0*1e-6,...
        'Description',{
            'Add to the Hessian matrix appropriate identity matrices scaled by this constant.'
            '[Should be avoided since does not seem to have a good justification for Nash'
            'equilibria.]'
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
            
    % HCells={};
    % P1xnus={};
    % P2xnus={};
    % for k=1:length(latentConstraints)
    %     switch (type(latentConstraints{k}))
    %       case {'iszero'}
    %         % remove iszero
    %         op1=Tcalculus(operands(latentConstraints{k}));
    %         HCells{end+1}=op1; 
    %         % create appropriate nu
    %         P1xnus{end+1}=Tvariable(sprintf('P1xnu%d_',length(P1xnus)+1),size(op1));
    %         P2xnus{end+1}=Tvariable(sprintf('P2xnu%d_',length(P2xnus)+1),size(op1));
    %         declareSet(code,P1xnus{end},sprintf('setD_%s',classname,name(P1xnus{end})));
    %         declareSet(code,P2xnus{end},sprintf('setD_%s',classname,name(P2xnus{end})));
    %       otherwise
    %         error('latent constraint of type ''%s'' not implemented\n',...
    %               type(latentConstraints{k}));
    %     end
    % end
    
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
    Hess__=ipmPDeqlat_CS(code,P1objective,P2objective,u,d,x,P1lambda,P1nu,P1xnu,P2lambda,P2nu,P2xnu,...
                         Fu,Gu,Fd,Gd,H,...
                         smallerNewtonMatrix,addEye2Hessian,skipAffine,...
                         scaleInequalities,scaleCost,scaleEqualities,...
                         false,...
                         classname,allowSave,debugConvergence,profiling);
    code.statistics.time.ipmPD=etime(clock,t_ipmPD);
    outputExpressions=substitute(outputExpressions,...
                                 Tvariable('Hess_',size(Hess__),true),Hess__);
    
    %% Declare ipm solver 
    classhelp{end+1}='Solve optimization';
    classhelp{end+1}='[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));';
    defines.nZ=size(u,1)+size(d,1)+size(x,1);
    defines.nU=size(u,1);
    defines.nD=size(d,1);
    defines.nX=size(x,1);
    defines.nG=size(Gu,1)+size(Gd,1)+size(H,1);
    defines.nF=size(Fu,1)+size(Fd,1);
    defines.nNu=size(Gu,1)+size(Gd,1)+2*size(H,1);
    defines.gradTolerance=gradTolerance;
    defines.equalTolerance=equalTolerance;
    defines.desiredDualityGap=desiredDualityGap;
    defines.scaleCost=double(scaleCost);
    defines.scaleInequalities=double(scaleInequalities);
    defines.scaleEqualities=double(scaleEqualities);
    defines.alphaMin=alphaMin;
    defines.alphaMax=alphaMax;
    defines.coupledAlphas=double(coupledAlphas);
    defines.muFactorAggressive=muFactorAggressive;
    defines.muFactorConservative=muFactorConservative;
    defines.delta=delta;
    defines.skipAffine=double(skipAffine);
    defines.smallerNewtonMatrix=double(smallerNewtonMatrix);
    defines.useUmfpack=double(useUmfpack);
    defines.allowSave=double(allowSave);
    defines.debugConvergence=double(debugConvergence);
    defines.debugConvergenceThreshold=debugConvergenceThreshold;
    defines.profiling=double(profiling);
    defines.verboseLevel=solverVerboseLevel;
    
    pth=fileparts(which('class2equilibriumLatentCS.m'));
    declareFunction(code,fsfullfile(pth,'ipmPDeq_CSsolver.m'),'solve',defines,[],[],'solve');
    
    %% Declare 'gets' for output expressions
    classhelp{end+1}='Get outputs';
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
