function [classname,pars]=createSolver(obj,classPrefix,iterativeSolvers,solverPars)
%   [classname,pars]=createSolver(obj,classPrefix,iterativeSolvers,solverPars)
% creates solver for the model
%
% Parameters
%    iterativeSolvers = true to create an iterative solver (only used when marginalizing variables ?)
%    solverPars = structure with parameters to pass to tenscalc solver generator
%
% Returns
%    classname = class used to access tenscalc's solver
%    pars = parameters used to call cmex2optimizeCS /class2optimizeCS
    if nargin<3
        iterativeSolvers=false;
    end

    %% # unknows
    nDOF=sum(obj.nParameters2optimize)+sum(obj.nStates2optimize)...
         +sum(obj.nStochasticInputsKnownVariance)+sum(obj.nStochasticInputsUnknownVariance)...
         -sum(obj.nStateEquations)-sum(obj.nMeasurements)...
         -length(obj.nStochasticInputsUnknownVariance);

    fprintf('Creating solver:\n');
    fprintf('   %4d parameters to optimize [%s] (including %d variances)\n',...
            sum(obj.nParameters2optimize),index2str(obj.nParameters2optimize),length(obj.nStochasticInputsUnknownVariance));
    fprintf('   %4d states     to optimize [%s]\n',...
            sum(obj.nStates2optimize),index2str(obj.nStates2optimize));
    fprintf('   %4d state eqs. to optimize [%s]\n',...
            sum(obj.nStateEquations),index2str(obj.nStateEquations));
    fprintf('   %4d measurements [%s]\n',...
            sum(obj.nMeasurements),index2str(obj.nMeasurements));
    fprintf('   %4d stochastic inputs with known variance [%s]\n',...
            sum(obj.nStochasticInputsKnownVariance),index2str(obj.nStochasticInputsKnownVariance));
    fprintf('   %4d stochastic inputs with unknown variance [%s]\n',...
            sum(obj.nStochasticInputsUnknownVariance),index2str(obj.nStochasticInputsUnknownVariance));
    fprintf('   %4d DOF for solver (excluding variances)\n',nDOF);

    if any(nDOF>obj.nStochasticInputsUnknownVariance)
        warning('*** ATTENTION: excessive # of DOFs (%d) may allow one or more stochastic inputs with unknown variance (%s) to be zero, leading to unbounded (-infty) likelihood\n',nDOF,index2str(obj.nStochasticInputsUnknownVariance));
    end

    %% Split variables between optimizations & determine marginalization variables
    optimizationVariables1=struct();
    optimizationVariables2=struct();
    marginalizationVariables={};
    marginalizationVariablesScaling={};
    nOptimizationVariables1=0;
    nOptimizationVariables2=0;
    solverParameters1=struct();
    solverParameters2=struct();
    fn=fields(obj.optimizationVariables);
    for f=1:length(fn)
        switch(obj.optimizationVariablesType.(fn{f}))
          case {'state1','state01'}
            % state with equality constraints
            % (currently cannot be marginalized since derivative does not take into account equality constraints)
            % . optimization variable in optimization 1
            optimizationVariables1.(fn{f})=obj.optimizationVariables.(fn{f});
            % . parameter in optimization 2
            solverParameters2.(fn{f})=obj.optimizationVariables.(fn{f});
          case {'state01_toMarginalize','state1_toMarginalize'}
            % state as random variable to be marginalized
            % . optimization variable in optimization 1, to be marginalized
            optimizationVariables1.(fn{f})=obj.optimizationVariables.(fn{f});
            marginalizationVariables{end+1,1}=obj.optimizationVariables.(fn{f});
            [variable,scalingVariable]=getScalingVariables(obj,obj.states.(fn{f}));
            marginalizationVariablesScaling{end+1,1}=scalingVariable;
            % . parameter in optimization 2
            solverParameters2.(fn{f})=obj.optimizationVariables.(fn{f});
          case 'parameter'
            if obj.marginalizeParameters
                % . optimization variable in optimization 1, to marginalize
                optimizationVariables1.(fn{f})=obj.optimizationVariables.(fn{f});
                marginalizationVariables{end+1,1}=obj.optimizationVariables.(fn{f});
                [variable,scalingVariable]=getScalingVariables(obj,obj.parameters.(fn{f}));
                marginalizationVariablesScaling{end+1,1}=scalingVariable;
                % . parameter in optimization 2
                solverParameters2.(fn{f})=obj.optimizationVariables.(fn{f});
            else
                % . parameter in optimization 1
                solverParameters1.(fn{f})=obj.optimizationVariables.(fn{f});
                % . optimization variable in optimization 2
                optimizationVariables2.(fn{f})=obj.optimizationVariables.(fn{f});
            end
          case {'disturbanceInvVariance','noiseInvVariance'}
            % . parameters in optimization 1
            solverParameters1.(fn{f})=obj.optimizationVariables.(fn{f});
            % . optimization variable in optimization 2
            optimizationVariables2.(fn{f})=obj.optimizationVariables.(fn{f});
          otherwise
            error('unknown type ''%s'' for optimization variable ''%s''\n',obj.optimizationVariablesType.(fn{f}),fn{f});
        end
    end
    % Add output forecasts
    fn=fields(obj.forecasts);
    for f=1:length(fn)
        marginalizationVariables{end+1,1}=obj.forecasts.(fn{f}).variable;
        marginalizationVariablesScaling{end+1}=inf;
    end

    % count optimization variables
    nOptimizationVariables1=sum(cellfun(@numel,struct2cell(optimizationVariables1)));
    nOptimizationVariables2=sum(cellfun(@numel,struct2cell(optimizationVariables2)));

    obj.marginalizationNames=cellfun(@name,marginalizationVariables,'UniformOutput',false);
    obj.marginalizationSizes=cellfun(@msize,marginalizationVariables,'UniformOutput',false);

    %% Compute log-det margialization correction
    logPDFs=struct2cell(obj.logPDF);
    logPDFs=cat(1,logPDFs{:});
    logJoint=sum(logPDFs,'all');
    logPDFsForecasts=struct2cell(obj.logPDFforecasts);
    logPDFsForecasts=cat(1,logPDFsForecasts{:});
    logJointForecasts=sum(logPDFsForecasts,'all');

    if ~isempty(marginalizationVariables)
        fprintf('marginalizing %d variables to reduce DOFs:\n',length(marginalizationVariables));
        for m=1:length(marginalizationVariables)
            fprintf('   %d: %-20s [%s] - %d scalars\n',...
                    m,obj.marginalizationNames{m},...
                    index2str(obj.marginalizationSizes{m}),prod(obj.marginalizationSizes{m}));
        end
        if false && obj.addEye2marginalization>0
            obj.addEye2marginalization=obj.addEye2marginalization*1e4;
            warning('*** ATTENTION: adding %e/2*norm2()/norm2(typical) to logJoint to improve numerical conditioning\n',obj.addEye2marginalization);
            for i=1:length(marginalizationVariables)
                logJoint=logJoint+(obj.addEye2marginalization/2)/norm2(marginalizationVariablesScaling{i})*norm2(marginalizationVariables{i});
            end
        end
        % Compute H matrix
        if 0
            [x2marginalize,whereVariables,packCmd,unpackCmd,logJoint1]= ...
                packVariables(marginalizationVariables,'x2marginalize',logJoint);
        else
            [x2marginalize,whereVariables,packCmd,unpackCmd,logJoint1]= ...
                packVariables(marginalizationVariables,'x2marginalize',logJoint+logJointForecasts);
        end
        H=hessian(logJoint1,x2marginalize);
        % back to original variables
        packMV=packExpressions(marginalizationVariables);
        H=substitute(H,x2marginalize,packMV);
        nH=size(H,1);
        if true && obj.addEye2marginalization>0
            H=H+obj.addEye2marginalization*Teye(size(H));
            warning('*** ATTENTION: adding %e to maginalization matrix to improve numerical conditioning\n',obj.addEye2marginalization);
        end
        logdetHmarginalization=logdet(ldl(H));

        % compute linearizations to estimate confidence intervals
        fn=fields(obj.outputs);
        for f=1:length(fn)
            y=substitute(obj.outputExpressions.(fn{f}),marginalizationVariables,x2marginalize);
            Dy=gradient(y,x2marginalize);
            Dy=substitute(Dy,x2marginalize,packMV);
            obj.outputExpressions.([fn{f},'_jacobian'])=Dy;
        end

        obj.outputExpressions.H=H; %addOutput(obj,'H',H);
        addOutput(obj,'logdetHmarginalization',logdetHmarginalization);
        logPi=nH*log(2*pi);
        addOutput(obj,'logPi',logPi);
        logMarginal=-.5*logPi+.5*logdetHmarginalization+logJoint;
        addOutput(obj,'logMarginal',logMarginal);
    else
        addOutput(obj,'logMarginal',logJoint);
    end
    addOutput(obj,'logJoint',logJoint);
    addOutput(obj,'logPDFs',logPDFs);

    %% replace forecast variables by their optimal values
    fn=fields(obj.forecasts);
    for f=1:length(fn)
        obj.outputExpressions=substitute(obj.outputExpressions,...
                                         obj.forecasts.(fn{f}).variable,obj.forecasts.(fn{f}).optimalValue);
    end

    %% Scale all variables and expressions
    scaleAllVariables(obj);

    %% Split constraints between optimizations -- must appear after scaling!
    constraints1={};
    constraints2={};
    for f=1:length(obj.constraints)
        switch(obj.constraintsType{f})
          case {'state','noise','disturbance','userDefined'}
            % goes to optimization 1
            constraints1{end+1,1}=obj.constraints{f};
          case {'disturbanceInvVariance','noiseInvVariance'}
            % goes to optimization 2
            constraints2{end+1,1}=obj.constraints{f};
          case {'parameter'}
            if obj.marginalizeParameters
                % goes to optimization 1
                constraints1{end+1,1}=obj.constraints{f};
            else
                % goes to optimization 2
                constraints2{end+1,1}=obj.constraints{f};
            end
          otherwise
            obj.constraints{f}
            error('unknown type ''%s'' for %d constraint\n',obj.constraintsType{f},f);
        end
    end

    nConstraints=sum(cellfun(@numel,obj.constraints));
    nConstraints1=sum(cellfun(@numel,constraints1));
    nConstraints2=sum(cellfun(@numel,constraints2));
    equalityConstraints1=cellfun(@(c)strcmp(type(c),'iszero'),constraints1);
    equalityConstraints2=cellfun(@(c)strcmp(type(c),'iszero'),constraints2);
    nInequalityConstraints1=sum(cellfun(@numel,constraints1(~equalityConstraints1)));
    nInequalityConstraints2=sum(cellfun(@numel,constraints2(~equalityConstraints2)));

    if any(equalityConstraints1) && ~isempty(marginalizationVariables)
        error('current implementation cannot handle marginalization with equality constraints');
    end

    %% Common to all optimizations
    if nargin<4
        pars=struct();
    else
        pars=solverPars;
    end
    pars.gradTolerance=5e-7;
    pars.equalTolerance=1e-7;
    pars.desiredDualityGap=1e-8;

    pars.scaleCost=1e1;
    pars.scaleInequalities=true;
    pars.pedigreeClass=classPrefix;
    pars.executeScript='asneeded';
    pars.solverVerboseLevel=4;

    if ~iterativeSolvers || isempty(marginalizationVariables)
        fprintf('Building single full solver\n');
        fprintf('            %4d optimization variables\n',obj.nOptimizationVariables);
        fprintf('            %4d constraints\n',nConstraints);
        %% Full optimization
        pars.objective=obj.outputExpressions.logMarginal;

        pars.optimizationVariables=obj.optimizationVariables;
        pars.parameters=struct2cell(obj.solverParameters);
        pars.constraints=obj.constraints;

        if obj.useTvars
            tpars.objective=pars.objective;
            tpars.optimizationVariables=pars.optimizationVariables;
            tpars.constraints=pars.constraints;
            outVars=Tvars2optimizeCS(tpars);
        end

        pars.outputExpressions=obj.outputExpressions;
        if obj.useTvars
            pars.outputExpressions=mergestruct(pars.outputExpressions,outVars);
        else
            pars.outputExpressions.Hess_=Tvariable('Hess_',[1,1]...
                                                   *(obj.nOptimizationVariables...
                                                     +nConstraints),true,true);
            pars.outputExpressions.u_=Tvariable('u_',obj.nOptimizationVariables,true,true);
            pars.outputExpressions.addEye2Hessian1=Tvariable('addEye2HessianU__',[],true,false);
            pars.outputExpressions.addEye2HessianEq=Tvariable('addEye2HessianEq__',[],true,false);
        end

        pars.outputExpressions.cost=pars.objective;

        %pars.executeScript='yes';
        %pars.debugConvergence=true;
        classname=cmex2optimizeCS(pars);
        obj.solverObjFull=feval(classname);

    else
        fprintf('Building 2 iterative solvers\n');
        %% 1st optimizations
        fprintf('  Solver 1: %4d optimization variables\n',nOptimizationVariables1);
        fprintf('            %4d constraints\n',nConstraints1);
        pars1=pars;

        pars1.muFactorAggressive=.5;

        pars1.objective=obj.outputExpressions.logJoint;

        pars1.optimizationVariables=optimizationVariables1;
        pars1.parameters=struct2cell(mergestruct(obj.solverParameters,solverParameters1));
        pars1.constraints=constraints1;

        if obj.useTvars
            tpars.objective=pars1.objective;
            tpars.optimizationVariables=pars1.optimizationVariables;
            tpars.constraints=pars1.constraints;
            outVars=Tvars2optimizeCS(tpars);
        end

        pars1.outputExpressions=obj.outputExpressions;
        if obj.useTvars
            pars1.outputExpressions=mergestruct(pars1.outputExpressions,outVars);
        elseif 1
            pars1.outputExpressions.Hess_=Tvariable('Hess_',[1,1]...
                                                    *(nOptimizationVariables1...
                                                      +nConstraints1),true,true);
            pars1.outputExpressions.u_=Tvariable('u_',nOptimizationVariables1,true,true);
            pars1.outputExpressions.addEye2HessianU=Tvariable('addEye2HessianU__',[],true,false);
            pars1.outputExpressions.addEye2HessianEq=Tvariable('addEye2HessianEq__',[],true,false);
            pars1.outputExpressions.grad=Tvariable('Lf_u_',nOptimizationVariables1,true,true);
            pars1.outputExpressions.lambda=Tvariable('lambda_',nInequalityConstraints1,true,true);
        end

        pars1.outputExpressions.cost=pars1.objective;
        pars1.addEye2HessianUtolerance=100;
        %pars1.debugConvergence=true;
        pars1.solverVerboseLevel=2;
        %pars1.executeScript='yes';

        %% 2nd opptimizations
        fprintf('  Solver 2: %4d optimization variables\n',nOptimizationVariables2);
        fprintf('            %4d constraints\n',nConstraints2);
        pars2=pars;
        pars2.objective=obj.outputExpressions.logMarginal;

        pars2.optimizationVariables=optimizationVariables2;
        pars2.parameters=struct2cell(mergestruct(obj.solverParameters,solverParameters2));
        pars2.constraints=constraints2;

        if obj.useTvars
            tpars.objective=pars2.objective;
            tpars.optimizationVariables=pars2.optimizationVariables;
            tpars.constraints=pars2.constraints;
            outVars=Tvars2optimizeCS(tpars);
        end

        pars2.outputExpressions=obj.outputExpressions;
        if obj.useTvars
            pars2.outputExpressions=mergestruct(pars2.outputExpressions,outVars);
        elseif 1
            pars2.outputExpressions.Hess_=Tvariable('Hess_',[1,1]...
                                                    *(nOptimizationVariables2 ...
                                                      +nConstraints2),true,true);
            pars2.outputExpressions.u_=Tvariable('u_',nOptimizationVariables2,true,true);
            pars2.outputExpressions.addEye2HessianU=Tvariable('addEye2HessianU__',[],true,false);
            pars2.outputExpressions.addEye2HessianEq=Tvariable('addEye2HessianEq__',[],true,false);
            pars2.outputExpressions.grad=Tvariable('Lf_u_',nOptimizationVariables2,true,true);
            pars2.outputExpressions.lambda=Tvariable('lambda_',nInequalityConstraints2,true,true);
        end

        pars2.outputExpressions.cost=pars2.objective;

        pars2.scaleCost=1e1;
        pars2.addEye2HessianUtolerance=100;
        %pars2.executeScript='yes';
        %pars2.muFactorAggressive=.25
        %pars2.debugConvergence=true;
        %pars2.solverVerboseLevel=5;

        %% create solvers
        classname=cmex2optimizeCS(pars1);
        obj.solverObj1=feval(classname);

        classname=class2optimizeCS(pars2);
        obj.solverObj2=feval(classname);
    end
end
