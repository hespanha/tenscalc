function [outSolver,outExpression,kbest]=callSolver(obj,maxIter,mu0,maxSolversIter,plotIterations)
%    [outSolver,outExpression,kbest]=callSolver(obj,maxIter,mu0,maxSolversIter,plotIterations)
% call a solver
%   maxSolversIter - number of solver calls for iterative solvers
%   plotIterations - function that shows plots as iterative solvers progress
%
    
    if nargin<5
        plotIterations=@internalPlotIterations;
    end

    if isempty(obj.solverObjFull) && isempty(obj.solverObj1)
        error('solver not yet created, call ''createSolver()'' before calling ''solver()''');
    end

    setAllParameters(obj);
    initializeAllStates(obj);

    if nargin<4
        maxSolversIter=100;
    end
    if nargin<3
        mu0=.01;
    end
    if nargin<2
        maxIter=410;
    end
    saveIter=-1;

    checkSavedBounds(obj);

    if ~isempty(obj.solverObjFull)
        % single solver
        fprintf('Calling solver:\n');
        startTime=datetime('now');
        outSolver=solve(obj.solverObjFull,mu0(1),int32(maxIter),int32(saveIter));
        fprintf('Getting outputs:\n');
        outExpression=getOutputs(obj.solverObjFull);
        outSolver.startTime=startTime;
        outSolver.endTime=datetime('now');
        kbest=1;
    end

    if ~isempty(obj.solverObj1)
        tol=1e-3;
        if length(mu0)<2
            mu0(2)=mu0(1);
        end
        if length(maxIter)<2
            maxIter(2)=maxIter(1);
        end
        for i=1:maxSolversIter
            %fprintf('paused\n');pause;
            % call 1st solver
            fprintf('%4d/%4d Calling 1st solver\n',i,maxSolversIter);
            startTime=datetime('now');
            sol=solve(obj.solverObj1,mu0(1),int32(maxIter(1)),int32(saveIter));
            out=getOutputs(obj.solverObj1);
            out=addStddev(obj,out);
            sol.startTime=startTime;
            sol.endTime=datetime('now');
            if i==1
                outSolver=sol;
                outExpression=out;
            else
                outSolver(end+1,1)=sol;
                outExpression(end+1,1)=out;
            end
            fprintf('   %4d/%4d 1st solver returned: log marginal pdf = %.3f\n',...
                    i,maxSolversIter,outExpression(end).logMarginal);
            report(obj,outSolver(end),outExpression(end));
            plotIterations(outSolver,outExpression);
            if outSolver(end).status
                %outExpression(end)
                outSolver(end),
                if length(outExpression)>1
                    % erase last full step if there is a good previous step
                    outSolver(end)=[];
                    outExpression(end)=[];
                end
                break;
            end
            % pass outputs to 2nd solver
            fn=fields(obj.optimizationVariables);
            for f=1:length(fn)
                value=outExpression(end).(fn{f});
                switch(obj.optimizationVariablesType.(fn{f}))
                  case {'state01','state01_toMarginalize'}
                    % never optimized by 2nd solver so no need to check bounds
                    setState(obj,fn{f},value,false);
                  case {'state1','state1_toMarginalize'}
                    % never optimized by 2nd solver so no need to check bounds
                    setState(obj,fn{f},value(2:end),false);
                  case {'parameter'}
                    if ~obj.marginalizeParameters && ...
                            checkBound(obj,obj.parameters.(fn{f}),value,sprintf('parameter ''%s''',fn{f}),tol)
                        % back to default to avoid using a value that hits the boundaries
                        setParameter(obj,fn{f});
                    else
                        setParameter(obj,fn{f},value,false);
                    end
                  case {'disturbanceInvVariance','noiseInvVariance'}
                    if checkBound(obj,obj.parameters.(fn{f}),value,sprintf('parameter ''%s''',fn{f}),tol)
                        % back to default to avoid using a value that hits the boundaries
                        setParameter(obj,fn{f});
                    else
                        setParameter(obj,fn{f},value,false);
                    end
                end
            end
            % call 2nd solver
            fprintf('%4d/%4d Calling 2d solver\n',i,maxSolversIter);
            addEye2Hessian=[1e-6;1e-6];
            startTime=datetime('now');
            sol=solve(obj.solverObj2,mu0(2),int32(maxIter(2)),int32(saveIter),addEye2Hessian);
            out=getOutputs(obj.solverObj2);
            out=addStddev(obj,out);
            sol.startTime=startTime;
            sol.endTime=datetime('now');
            outSolver(end+1,1)=sol;
            outExpression(end+1,1)=out;
            fprintf('   %4d/%4d 2nd solver returned: log marginal pdf = %.3f\n',...
                            i,maxSolversIter,outExpression(end).logMarginal);
            %outExpression(end)
            report(obj,outSolver(end),outExpression(end));
            plotIterations(outSolver,outExpression);
            if outSolver(end).status
                % erase last full step if there is a good previous step
                outSolver(end),
                if  length(outExpression)>2
                    % erase last full step
                    outSolver(end-1:end)=[];
                    outExpression(end-1:end)=[];
                end
                break;
            end
            % pass outputs to 1st solver
            for f=1:length(fn)
                value=outExpression(end).(fn{f});
                switch(obj.optimizationVariablesType.(fn{f}))
                  case {'state01','state01_toMarginalize'}
                    if true %% checkBound(obj,obj.states.(fn{f}),value,sprintf('state ''%s''',fn{f}),tol)
                            % back to default to avoid using a value that hits the boundaries
                        setState(obj,fn{f});
                    else
                        % dangerous if there are user-defined state constraints
                        setState(obj,fn{f},value,false);
                    end
                  case {'state1','state1_toMarginalize'}
                    if true %% checkBound(obj,obj.states.(fn{f}),value(2:end),sprintf('state ''%s''',fn{f}),tol)
                            % back to default to avoid using a value that hits the boundaries
                        setState(obj,fn{f});
                    else
                        % dangerous if there are user-defined state constraints
                        setState(obj,fn{f},value(2:end),false);
                    end
                  case {'parameter'}
                    if obj.marginalizeParameters && ...
                            checkBound(obj,obj.parameters.(fn{f}),value,sprintf('parameter ''%s''',fn{f}),tol)
                        % back to default to avoid using a value that hits the boundaries
                        setParameter(obj,fn{f});
                    else
                        setParameter(obj,fn{f},value,false);
                    end
                  case {'disturbanceInvVariance','noiseInvVariance'}
                    % never optimized by 1st solver so no need to check bounds
                    setParameter(obj,fn{f},value,false);
                end
            end
            outSolver(end,1).startTime=startTime;
            outSolver(end,1).endTime=datetime('now');
        end

        [~,kbest]=min([outExpression.logMarginal]);

        % Compute a-posterior for last and also for best
        % for k=union(kbest,length(outExpression))
        %     outExpression(k)=addStddev(obj,outExpression(k));
        % end
    end
end

function internalPlotIterations(outSolver,outExpression)
    figure(200);clf;
    nIter=cat(1,outSolver.iter);
    logMarginals=cat(1,outExpression.logMarginal);
    addEye2Hessian1=cat(1,outExpression.addEye2Hessian1);
    addEye2Hessian2=cat(1,outExpression.addEye2Hessian2);
    subplot(3,1,1);
    yyaxis left
    plot(1:2:length(logMarginals),nIter(1:2:end),'o-');
    ylabel('# iterations')
    yyaxis right
    plot(2:2:length(logMarginals),nIter(2:2:end),'x-');
    grid on
    subplot(3,1,2);
    yyaxis left
    semilogy(1:2:length(logMarginals),addEye2Hessian1(1:2:end),'o-');
    ylabel('addEye2Hessian1')
    yyaxis right
    semilogy(2:2:length(logMarginals),addEye2Hessian1(2:2:end),'x-');
    grid on
    subplot(3,1,3)
    plot(1:2:length(logMarginals),logMarginals(1:2:end),'o-',...
         2:2:length(logMarginals),logMarginals(2:2:end),'x-');
    grid on
    ylabel('log-marginal');
    title(sprintf('best log marginal pdf = %.3f',min(logMarginals)));
    xlabel('solver iteration');
    drawnow;
end

function out=addStddev(obj,out)
%% Compute a-posteriori covariance matrix for marginalized variables from H
    try
        [out.posterioriCov,posterioriVar,~,out.kernelH]=invSparsePD(out.H,obj.addEye2marginalization);
    catch
        out.posterioriCov=nan(size(out.H));
        posterioriVar=nan(size(out.H,1),1);
        out.kernelH=nan(size(out.H,1),0);  
    end
    out.H=out.H-obj.addEye2marginalization*eye(size(out.H));
    posterioriStd=sqrt(max(posterioriVar,0));
    kNeg=posterioriVar<=-1e-4;
    if any(kNeg)
        disp(posterioriVar(kNeg)')
        posterioriStd(kNeg)=nan;
        warning('*** ATTENTION: %d/%d negative variances in inv(H): %s\n',...
                sum(kNeg),numel(kNeg),summarizeValues(obj,posterioriVar));
    end
    
    %% Compute a-posteriori covariance matrices for outputs
    fn=fields(obj.outputs);
    for f=1:length(fn)
        jname=[fn{f},'_jacobian'];
        if isfield(out,jname)
            Dy=out.(jname);
            nd=length(size(obj.outputExpressions.(fn{f})));
            switch nd
              case 0 % scalar
                out.([fn{f},'_posterioriCov'])=Dy'*out.posterioriCov*Dy;
                var=diag(out.([fn{f},'_posterioriCov']));
                kNeg=var<0;
                if any(kNeg)
                    warning('*** ATTENTION: %d/%d negative variances in %s: %s\n',...
                            sum(kNeg),numel(kNeg),fn{f},summarizeValues(obj,var));
                    var(kNeg)=nan;
                end
                std=sqrt(var);
                out.([fn{f},'_posterioriStd'])=std;
                        out.([fn{f},'_CI2sigma'])=out.(fn{f})+std*[-2,2];
              case 1 % vector
                     % covariance matrix -- only valid if output is a vector (otherwise tprod would be needed)
                out.([fn{f},'_posterioriCov'])=Dy*out.posterioriCov*Dy';
                var=diag(out.([fn{f},'_posterioriCov']));
                kNeg=var<0;
                if any(kNeg)
                    warning('*** ATTENTION: %d/%d negative variances in %s: %s\n',...
                            sum(kNeg),numel(kNeg),fn{f},summarizeValues(obj,var));
                    var(kNeg)=nan;
                end
                std=sqrt(var);
                out.([fn{f},'_posterioriStd'])=std;
                out.([fn{f},'_CI2sigma'])=out.(fn{f})+std*[-2,2];
              otherwise
                out.([fn{f},'_posterioriCov'])=tprod(Dy,[1:nd,-1],Dy,[nd+1:2*nd,-2],...
                                                                  out.posterioriCov,[-1,-2]);
                warning('confidence intervals not computed for tensor outputs')
            end
        end
    end
    
    if 0
        % old code
        nMarg=0;
        for f=1:length(obj.marginalizationNames)
            sz=obj.marginalizationSizes{f};
            nl=prod(sz);
            value=out.(obj.marginalizationNames{f});
            std=reshape(posterioriStd(nMarg+(1:nl)),sz);
            if size(value,1)==size(std,1)+1 % state with given initial conditions
                std=[zeros(1,size(std,2));
                     std];
            end
            name=sprintf('%s_posterioriStd',obj.marginalizationNames{f});
            out.(name)=std;
            if sz(2)==1
                % for 1-D states/parameters, computed 2-sigma confidence intervals
                name=sprintf('%s_CI2sigma',obj.marginalizationNames{f});
                out.(name)=out.(obj.marginalizationNames{f})+std*[-2,2];
            end
            nMarg=nMarg+nl;
        end
        if nMarg~=size(out.H,1)
            nMarg,size(out.H,1)
            error('bug in code: unexpected number of marginalized variable');
        end
    end
    
end
