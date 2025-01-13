classdef TCsysid <handle;

    properties
       T;  % time horizon for state estimation

       %% System model

       % structure of structures of the form
       % struct('name',{},...                  % char array
       %        'toOptimize',{},...            % boolean variable (true if variable will be optimized)
       %        'optimizationVariable',{},...  % Tcalculus optimization variable, empty if not an optimization variable
       %        'TCVariable',{},...            % Tcalculus variable
       %        'whereInX',{},...              % indices of where variable entries appear
                                               %  in the solver's optimization vector
       %        'lowerBound',{},...            % index of parameter that corresponds to bound -inf for no bound
       %        'upperBound',{},...            % index of parameter that corresponds to bound +inf for no bound
       %        'scalingValue',{},...          % index of parameter that contains typical value for scaling
       %                                ...    %   empty means no scaling required
       %        'initializationValue',{});     % initial value for solver;
       parameters=struct();
       nParameters2optimize=[];

       % structure of structures of the form
       % struct('name',{},...                  % char array
       %        'toOptimize',{},...            % boolean variable (true if variable will be optimized)
       %        'optimizationVariable',{},...  % Tcalculus optimization variable
       %        'initializationValue',{},...   % initial value for solver;
       %        'lowerBound',{},...            % index of parameter that corresponds to bound -inf for no bound
       %        'upperBound',{},...            % index of parameter that corresponds to bound +inf for no bound
       %        'scalingValue',{},...          % index of parameter that contains typical value for scaling
       %                                ...    %   empty means no scaling required
       %        'disturbanceInvVariance',{},...% index of parameter that corresponds to 1/disturbance variance
       %                                    ...%   or inf when there is no disturbance
       %        'constraint',{});              % index of equality constraint when variance is zero
       states=struct();
       nStates2optimize=[];
       nStateEquations=[];

       % structure of structures of the form
       % struct('name',{},...                  % char array
       %        'timeInstants',{},...          % time instants at which measurements are available
       %        'noiseInvVariance',{},...      % index of parameter that corresponds to 1/noise variance
       %                                    ...%   or inf when there is no noise
       %        'constraint',{});              % index of equality constraint when variance is zero
       measurements=struct();
       nMeasurements=[];

       % structure of structures of the form
       % struct('name',{},...                  % char array
       %        'timeInstants',{},...          % time instants for which forecasts are wanted
       forecasts=struct();


       nStochasticInputsKnownVariance=[];
       nStochasticInputsUnknownVariance=[];

       % structure of structures of the form
       % struct('name',{});                  % char array
       outputs=struct();

       % Parameters used to build the solver
       constraints={};
       constraintsType={};
       nEqualityConstraints=0;
       nLatentStates=0;
       nLatentEqualityConstraints=0;
       nInequalityConstraints=0;
       latentVariables;
       optimizationVariables;
       optimizationVariablesType;
       nOptimizationVariables=0;
       solverParameters;
       outputExpressions;
       logPDF=struct();

       constraintsStates;
       constraintsParameters;
       logPDFstates;
       logPDFmeasurements;
       logPDFforecasts=struct();

       solverObjFull;
       solverObj1;
       solverObj2;

       % Parameters involved in marginalizing random variables
       marginalizeParameters=true;
       marginalizationNames={};
       marginalizationSizes={};
       addEye2marginalization=1e-4;

       % Parameters for code generation
       scaleVariables=true;
       alreadyScaledVariables=false;
       useTvars=false;
    end


    methods;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create object & utilities
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=TCsysid(T,marginalizeParameters)
            % obj=TCsysid(T,marginalizeParameters)
            %   T = total horizon length for the state variables (including initial condition)
            %   marginalize parameters ??? (see covid19)
            if nargin>1
                obj.marginalizeParameters=marginalizeParameters;
            end
            Tcalculus.clear;
            obj.T=T;
        end

        function str=summarizeValues(obj,value)
            value=full(value);
            if numel(value)<3
                str=sprintf(' %10.3f (%10.2e)',value,value);
            else
                mn=min(value(:));
                mx=max(value(:));
                if mn==mx
                    str=sprintf(' %10.3f (%10.2e)',mn,mn);
                else
                    str=sprintf('[%10.3f,%10.3f] ([%10.2e,%10.2e])',mn,mx,mn,mx);
                end
            end
        end

        function name=findState(obj,TCvariable)
            fn=fields(obj.states);
            for f=1:length(fn)
                if isequal(TCvariable,obj.states.(fn{f}).variable)
                    name=fn{f};
                    return
                end
            end
            disp(TCvariable);
            error('unable to find state variable');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check bounds
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [lower,upper]=hasBounds(obj,entity)
        % [lower,upper]=hasBounds(obj,entity)
        %
        %    returns nan if parameter/state does not have
        %    upper/lower bound, otherwise returns the value of the
        %    bound (as stored in the corresponding parameter).
            if ~isempty(entity.lowerBound) && ~isequal(entity.lowerBound,-inf)
                lower=obj.parameters.(entity.lowerBound).initializationValue;
            else
                lower=-inf;
            end
            if ~isempty(entity.upperBound) && ~isequal(entity.upperBound,inf)
                upper=obj.parameters.(entity.upperBound).initializationValue;
            else
                upper=inf;
            end
        end

        function [msg,kLower,kUpper]=hitBounds(obj,value,lower,upper,tol)
            if nargin<5
                tol=1e-5;
            end
            if isfinite(lower)
                if lower==0
                    kLower=value<tol;
                else
                    kLower=(value<lower+tol*lower);
                end
            else
                kLower=[];
            end
            if isfinite(upper)
                if upper==0
                    kUpper=value>-tol;
                else
                    kUpper=(value>upper-tol*upper);
                end
            else
                kUpper=[];
            end
            if any(kLower)
                if numel(kLower)>1
                    msg=sprintf('hitting lower at %d/%d points',sum(kLower),numel(kLower));
                else
                    msg=sprintf('hitting lower');
                end
            else
                msg='';
            end
            if any(kUpper)
                if numel(kUpper)>1
                    msg=sprintf('%shitting upper at %d/%d points',msg,sum(kUpper),numel(kUpper));
                else
                    msg=sprintf('%shitting upper',msg);
                end
            end
        end

        function err=checkBound(obj,entity,value,name,tol)
            [lower,upper]=hasBounds(obj,entity);
            if nargin<5
                [msg,kLower,kUpper]=hitBounds(obj,value,lower,upper);
            else
                [msg,kLower,kUpper]=hitBounds(obj,value,lower,upper,tol);
            end
            err=false;
            if any(kLower)
                err=true;
                fprintf('**%s %s violates lower bound %s: %s\n',...
                        name,summarizeValues(obj,value),...
                        summarizeValues(obj,lower),msg);
            else
                fprintf('  %s %s satisfies lower bound %s\n',...
                        name,...
                        summarizeValues(obj,value),...
                        summarizeValues(obj,lower));
            end
            if any(kUpper)
                err=true;
                fprintf('**%s %s violates upper bound %s: %s\n',...
                        name,...
                        summarizeValues(obj,value),...
                        summarizeValues(obj,upper),msg);
            else
                fprintf('  %s %s satisfies upper bound %s\n',...
                        name,...
                        summarizeValues(obj,value),...
                        summarizeValues(obj,upper));
            end
        end

        function checkSavedBounds(obj,tol)
            if nargin<2
                tol=1e-9;
            end
            err=false;
            % parameters
            fprintf('Checking parameter bounds:\n');
            fn=fields(obj.parameters);
            for f=1:length(fn)
                if obj.parameters.(fn{f}).toOptimize
                    if checkBound(obj,obj.parameters.(fn{f}),...
                                  obj.parameters.(fn{f}).initializationValue,...
                                  sprintf('parameter ''%s''',fn{f}));
                        err=true;
                    end
                end
            end
            fprintf('Checking state bounds:\n');
            fn=fields(obj.states);
            for f=1:length(fn)
                if obj.states.(fn{f}).toOptimize
                    if checkBound(obj,obj.states.(fn{f}),...
                                  obj.states.(fn{f}).initializationValue,...
                                  sprintf('state ''%s''',fn{f}));
                        err=true;
                    end
                end
            end
            if err
                error('some variables violate bounds');
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Scaling
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [variable,scalingVariable]=getScalingVariables(obj,entity)
        % [variable,scalingVariable]=getScalingVariables(obj,entity)
        %
        %   returns optimization and typical value variables; or
        %   returns variable=[], if entity requires no scaling
            if entity.toOptimize && ~isempty(entity.scalingValue)
                variable=entity.optimizationVariable;
                scalingVariable=obj.parameters.(entity.scalingValue).TCvariable;
            else
                variable=[];
                scalingVariable=[];
            end
        end

        function value=getScalingValue(obj,entity)
        % value=getScalingValue(obj,entity)
        %
        %   returns typical value or 1, if entiry requires no scaling
            if entity.toOptimize && ~isempty(entity.scalingValue)
                value=obj.parameters.(entity.scalingValue).initializationValue;
                if numel(value)==1 && numel(entity.initializationValue)>1
                    value=value*ones(size(entity.initializationValue));
                end
            else
                value=ones(size(entity.initializationValue));
            end
        end

        function varargout=substituteVariable(obj,var,expr,varargin)
            fn=fields(obj.outputExpressions);
            for f=1:length(fn)
                obj.outputExpressions.(fn{f})=substitute(obj.outputExpressions.(fn{f}),var,expr);
            end
            for f=1:length(obj.constraints)
                obj.constraints{f}=substitute(obj.constraints{f},var,expr);
            end
            % fn=fields(obj.logPDF);
            % for f=1:length(fn)
            %     obj.logPDF.(fn{f})=substitute(obj.logPDF.(fn{f}),var,expr);
            % end
            for f=1:length(varargin)
                varargout{f}=substitute(varargin{f},var,expr);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add elemens to sysid problem
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function y=logNormal(obj,x,invVar)
            y=+.5*log(2*pi)*numel(x)-.5*numel(x)*log(invVar)+.5*invVar*norm2(x);
        end

        function [variable,scalingVariable]=addParameter(obj,name,sz,lowerBound,upperBound,typicalValue,latent)
        % addParameter(obj,name,sz,typicalValue)
        %    parameters not optimized
        % addParameter(obj,name,sz,lowerBound,upperBound,typicalValue)
        %     both bounds can be reset at solve time, empty for no bound
            if ~ismember(nargin,[4,5,6,7])
                error('wrong  number of input parameters for addParameter');
            end
            variable=Tvariable(name,sz);
            obj.parameters.(name).TCvariable=variable; % create new entry
            if nargin==4
                % fixed parameter (not to optimize)
                obj.parameters.(name).toOptimize=false;
                obj.parameters.(name).scalingValue=[];
                obj.parameters.(name).initializationValue=lowerBound;
                obj.parameters.(name).lowerBound=[];
                obj.parameters.(name).upperBound=[];
                obj.solverParameters.(name)=variable;
                obj.parameters.(name).optimizationVariable=[];
                obj.outputExpressions.(name)=variable;
            else
                msz=sz;
                while length(msz)<2
                    msz(end+1)=1;
                end
                if nargin<6
                    typicalValue=ones(msz);
                end
                if nargin<7
                    latent=false;
                end
                if ~isequal(msz,size(typicalValue))
                    error('size [%s] does not match size of typicalValue [%s]\n',...
                          index2str(sz),index2str(size(typicalValue)));
                end
                obj.parameters.(name).toOptimize=true;
                obj.parameters.(name).isLatent=latent;
                obj.parameters.(name).optimizationVariable=variable;
                obj.optimizationVariables.(name)=variable;
                obj.optimizationVariablesType.(name)='parameter';
                % Attention: the following assumes that optimization variables appear
                % in the solver in the order they are added to
                % obj.optimizationVariables
                obj.parameters.(name).whereInX=obj.nOptimizationVariables+(1:numel(variable))';
                obj.nOptimizationVariables=obj.parameters.(name).whereInX(end);
                obj.outputExpressions.(name)=variable;
                obj.nParameters2optimize(end+1,1)=numel(variable);

                if all(typicalValue==1) || any(typicalValue==0)
                    obj.parameters.(name).scalingValue=[];
                    scalingVariable=[];
                else
                    scalingName=sprintf('%s_scaling',name);
                    scalingVariable=addParameter(obj,scalingName,sz,typicalValue);
                    obj.parameters.(name).scalingValue=scalingName;
                end
                obj.parameters.(name).initializationValue=typicalValue;
                if all(lowerBound==-inf)
                    obj.parameters.(name).lowerBound=-inf;
                else
                    boundName=sprintf('%s_lowerBound',name);
                    boundVariable=addParameter(obj,boundName,sz,lowerBound);
                    assignin('caller',boundName,boundVariable);
                    obj.parameters.(name).lowerBound=boundName;
                    obj.constraints{end+1,1}=variable>=boundVariable;
                    obj.constraintsType{end+1,1}='parameter';
                    obj.nInequalityConstraints=obj.nInequalityConstraints+numel(variable);
                    obj.constraintsParameters{end+1,1}=variable>=boundVariable;
                    obj.outputExpressions.(sprintf('%s_lowerBound_constraint',name))=variable-boundVariable;
                end
                if all(upperBound==inf)
                    obj.parameters.(name).upperBound=inf;
                else
                    boundName=sprintf('%s_upperBound',name);
                    boundVariable=addParameter(obj,boundName,sz,upperBound);
                    assignin('caller',boundName,boundVariable);
                    obj.parameters.(name).upperBound=boundName;
                    obj.constraints{end+1,1}=variable<=boundVariable;
                    obj.constraintsType{end+1,1}='parameter';
                    obj.nInequalityConstraints=obj.nInequalityConstraints+numel(variable);
                    obj.constraintsParameters{end+1,1}=variable<=boundVariable;
                    obj.outputExpressions.(sprintf('%s_upperBound_constraint',name))=boundVariable-variable;
                end
            end
            assignin('caller',name,variable);
        end

        function [variable,scalingVariable]=addState(obj,name,lowerBound,upperBound,initialValue,typicalValue)
        %  addState(obj,stateName,lowerBound,upperBound,initialValue,typicalValue)
        %     both bounds can be reset at solve time, empty for no bound
        %     initial state value can be set at solve time, nan for unknown initial state
            if nargin<6
                typicalValue=1;
            end
            obj.states.(name).toOptimize=true; % create new entry
            if isnan(initialValue)
                % estimate initial state
                obj.states.(name).variable=Tvariable(name,obj.T);
                obj.states.(name).variable_0=obj.states.(name).variable(1);
                obj.states.(name).variable_1=obj.states.(name).variable(2:end);
                obj.states.(name).optimizationVariable=obj.states.(name).variable;
                obj.states.(name).initializationValue=typicalValue*ones(obj.T,1);
                obj.optimizationVariables.(name)=obj.states.(name).variable;
                obj.optimizationVariablesType.(name)='state01';
                % Attention: the following assumes that optimization variables appear
                % in the solver in the order they are added to
                % obj.optimizationVariables
                obj.states.(name).whereInX=obj.nOptimizationVariables+(1:numel(obj.states.(name).variable))';
                obj.nOptimizationVariables=obj.states.(name).whereInX(end);
                obj.nStates2optimize(end+1,1)=numel(obj.states.(name).variable);
            else
                % given initial state
                obj.states.(name).variable_0=addParameter(obj,sprintf('%s_0',name),[],initialValue);
                obj.states.(name).variable_1=Tvariable(name,obj.T-1);
                obj.states.(name).variable=[obj.states.(name).variable_0;obj.states.(name).variable_1];
                obj.states.(name).optimizationVariable=obj.states.(name).variable_1;
                obj.optimizationVariables.(name)=obj.states.(name).variable_1;
                obj.optimizationVariablesType.(name)='state1';
                % Attention: the following assumes that optimization variables appear
                % in the solver in the order they are added to
                % obj.optimizationVariables
                obj.states.(name).whereInX=obj.nOptimizationVariables+(1:numel(obj.states.(name).variable_1))';
                obj.nOptimizationVariables=obj.states.(name).whereInX(end);
                obj.states.(name).initializationValue=typicalValue*ones(obj.T-1,1);
                obj.nStates2optimize(end+1,1)=numel(obj.states.(name).variable_1);
            end
            obj.outputExpressions.(name)=obj.states.(name).variable;
            if ~isequal(typicalValue,1) && ~isequal(typicalValue,0)
                typicalName=sprintf('%s_typical',name);
                typicalVariable=addParameter(obj,typicalName,[],typicalValue);
                obj.states.(name).scalingValue=typicalName;
            else
                    obj.states.(name).scalingValue=[];
            end
            if lowerBound==-inf
                obj.states.(name).lowerBound=lowerBound;
            else
                boundName=sprintf('%s_lowerBound',name);
                boundVariable=addParameter(obj,boundName,[],lowerBound);
                assignin('caller',boundName,boundVariable);
                obj.states.(name).lowerBound=boundName;
                obj.constraints{end+1,1}=obj.optimizationVariables.(name)>=boundVariable;
                obj.constraintsType{end+1,1}='state';
                obj.nInequalityConstraints=obj.nInequalityConstraints+numel(obj.optimizationVariables.(name));
                obj.constraintsStates{end+1,1}=obj.optimizationVariables.(name)>=boundVariable;
                obj.outputExpressions.(sprintf('%s_lowerBound_constraint',name))=obj.optimizationVariables.(name)-boundVariable;
            end
            if upperBound==inf
                obj.states.(name).upperBound=upperBound;
            else
                boundName=sprintf('%s_upperBound',name);
                boundVariable=addParameter(obj,boundName,[],upperBound);
                assignin('caller',boundName,boundVariable);
                obj.states.(name).upperBound=boundName;
                obj.constraints{end+1,1}=obj.optimizationVariables.(name)<=boundVariable;
                obj.constraintsType{end+1,1}='state';
                obj.nInequalityConstraints=obj.nInequalityConstraints+numel(obj.optimizationVariables.(name));
                obj.constraintsStates{end+1,1}=obj.optimizationVariables.(name)<=boundVariable;
                obj.outputExpressions.(sprintf('%s_upperBound_constraint',name))=boundVariable-obj.optimizationVariables.(name);
            end
            variable=obj.states.(name).variable;
            assignin('caller',name,obj.states.(name).variable);
        end

        function addMeasurement(obj,name,outputEquation,timeInstants,noiseStDev,typicalNoise)
        %  addMeasurement(obj,name,outputEquation,timeInstants,noiseStDev,typicalNoise)
        %     add measurement of the form
        %             output = outputEquation + noise
        %     time instants indexed from 1:horizon length
        %     additive noise std. dev. 0 means no noise, nan means unknown std. dev.
            if nargin<6
                typicalNoise=1;
            end
            obj.measurements.(name).timeInstants=timeInstants; % create new entry
            obj.measurements.(name).outputEquation=outputEquation;

            variable=addParameter(obj,name,numel(timeInstants),nan(numel(timeInstants),1));
            if noiseStDev==0
                % no noise, output is an equality constraint
                obj.constraints{end+1,1}=(variable==outputEquation(timeInstants));
                obj.constraintsType{end+1,1}='noise';
                obj.nEqualityConstraints=obj.nEqualityConstraints+numel(variable);
                obj.constraintsStates{end+1,1}=(variable==outputEquation(timeInstants));
                obj.measurements.(name).constraint=numel(obj.constraints);
                obj.measurements.(name).noiseInvVariance=inf;
                obj.outputExpressions.(sprintf('%s_noiseStdDev',name))=0;
            else
                noise=variable-outputEquation(timeInstants);
                noiseIVname=sprintf('%s_noiseInvVariance',name);
                if isnan(noiseStDev)
                    % unkown std. dev.
                    noiseInvVariance=addParameter(obj,noiseIVname,[],0,nan,1/typicalNoise^2);
                    % upper and lower bounds for addParameter above
                    obj.constraintsType{end-1}='noiseInvVariance';
                    obj.constraintsType{end}='noiseInvVariance';
                    obj.optimizationVariablesType.(noiseIVname)='noiseInvVariance';
                    obj.nStochasticInputsUnknownVariance(end+1,1)=numel(timeInstants);
                else
                    % known std. dev.
                    noiseInvVariance=addParameter(obj,noiseIVname,[],1/noiseStDev^2);
                    obj.nStochasticInputsKnownVariance(end+1,1)=numel(timeInstants);
                end
                assignin('caller',noiseIVname,noiseInvVariance);
                obj.measurements.(name).noiseInvVariance=noiseIVname;
                obj.measurements.(name).noiseInvVarianceVariable=noiseInvVariance;
                obj.outputExpressions.(sprintf('%s_noise',name))=noise;
                obj.logPDF.(sprintf('%s_noise',name))=logNormal(obj,noise,noiseInvVariance);
                obj.logPDFmeasurements.(sprintf('%s_noise',name))=logNormal(obj,noise,noiseInvVariance);
                obj.outputExpressions.(sprintf('%s_noiseStdDev',name))=1/sqrt(noiseInvVariance);
            end
            obj.nMeasurements(end+1,1)=numel(timeInstants);
            assignin('caller',name,variable);
        end

        function addMeasurementForecast(obj,output,timeInstants)
        %  addMeasurementForecast(obj,name,timeInstants)
        %     'output' must be an output created using addMeasurement()
        %     time instants indexed fromm 1:horizon length
            if ischar(output)
                name=output;
            elseif isa(output,'Tcalculus')
                name=output.name;
            end
            obj.forecasts.(name).timeInstants=timeInstants;
            variable=Tvariable([name,'_forecast'],numel(timeInstants));
            obj.forecasts.(name).variable=variable;
            addOutput(obj,variable.name,variable);

            % retrieve outputEquation and noiseInvVarianceVariable from measurement
            noiseInvVariance=obj.measurements.(name).noiseInvVarianceVariable;
            % variable will be replaced by optimal value at the end
            obj.forecasts.(name).optimalValue=obj.measurements.(name).outputEquation(timeInstants);

            noise=variable-obj.forecasts.(name).optimalValue;
            obj.logPDFforecasts.(sprintf('%s_noise',name))=logNormal(obj,noise,noiseInvVariance);
            assignin('caller',name,variable);
        end

        function addDynamics(obj,stateVariable,stateEquation,disturbanceStDev,typicalDisturbance)
        % addDynamics(obj,stateVariable,stateEquation,disturbanceStDev,typicalDisturbance)
        %    add discrete-time dynamics of the form
        %             stateVariable^+ = stateEquation + disturbance
        %    additive disturbance std. dev. 0 means no disturbance, nan means unknown std. dev.
            if nargin<5
                typicalDisturbance=1;
            end

            stateName=findState(obj,stateVariable);
            if disturbanceStDev==0
                % no disturbance, dynamics is an equality constraint
                obj.constraints{end+1,1}=(obj.states.(stateName).variable_1==stateEquation(1:end-1));
                obj.constraintsType{end+1,1}='disturbance';
                obj.nEqualityConstraints=obj.nEqualityConstraints+numel(obj.states.(stateName).variable_1);
                obj.nLatentEqualityConstraints=obj.nLatentEqualityConstraints+numel(obj.states.(stateName).variable_1);
                obj.nLatentStates=obj.nLatentStates+numel(obj.states.(stateName).variable_1);
                obj.constraintsStates{end+1,1}=(obj.states.(stateName).variable_1==stateEquation(1:end-1));
                obj.states.(stateName).constraint=numel(obj.constraints);
                obj.states.(stateName).disturbanceInvVariance=inf;
                obj.outputExpressions.(sprintf('%s_disturbanceStdDev',stateName))=0;
            else
                disturbance=obj.states.(stateName).variable_1-stateEquation(1:end-1);
                disturbanceIVname=sprintf('%s_disturbanceInvVariance',stateName);
                if isnan(disturbanceStDev)
                    % unkown std. dev.
                    disturbanceInvVariance=addParameter(obj,disturbanceIVname,...
                                                        [],0,nan,1/typicalDisturbance^2);
                    % upper and lower bounds
                    obj.constraintsType{end-1}='disturbanceInvVariance';
                    obj.constraintsType{end}='disturbanceInvVariance';
                    obj.optimizationVariablesType.(disturbanceIVname)='disturbanceInvVariance';
                    obj.nStochasticInputsUnknownVariance(end+1,1)=numel(disturbance);
                else
                    % known std. dev.
                    disturbanceInvVariance=addParameter(obj,disturbanceIVname,...
                                                        [],1/disturbanceStDev^2);
                    obj.nStochasticInputsKnownVariance(end+1,1)=numel(disturbance);
                end
                assignin('caller',disturbanceIVname,disturbanceInvVariance);
                obj.states.(stateName).disturbanceInvVariance=disturbanceIVname;
                obj.optimizationVariablesType.(stateName)=[obj.optimizationVariablesType.(stateName),'_toMarginalize'];
                obj.outputExpressions.(sprintf('%s_disturbance',stateName))=disturbance;
                obj.logPDF.(sprintf('%s_disturbance',stateName))=...
                    logNormal(obj,disturbance,disturbanceInvVariance);
                obj.logPDFstates.(sprintf('%s_disturbance',stateName))=...
                    logNormal(obj,disturbance,disturbanceInvVariance);
                obj.outputExpressions.(sprintf('%s_disturbanceStdDev',stateName))=...
                    1/sqrt(disturbanceInvVariance);
            end
            obj.nStateEquations(end+1,1)=numel(obj.states.(stateName).variable_1);
        end

        function outputEquation=addOutput(obj,name,outputEquation)
            obj.outputs.(name)=name;
            obj.outputExpressions.(name)=outputEquation;
            assignin('caller',name,outputEquation);
        end

        function addConstraint(obj,constraint)
            obj.constraints{end+1,1}=constraint;
            obj.constraintsType{end+1,1}='userDefined';
            switch (type(constraint))
              case 'ispositive'
                obj.nInequalityConstraints=obj.nInequalityConstraints+numel(constraint);
              case 'iszero'
                obj.nEqualityConstraints=obj.nEqualityConstraints+numel(constraint);
              otherwise
                error('contrainst if type ''%s'' not yet implements\n',type(constraint));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create Solver
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function scaleAllVariables(obj)
        % Rescales
        %  . optimization parameters
        %  . states
        % by calling substituteVariable to update
        %  . outputExpressions
        %  . constraints
            if obj.scaleVariables && ~obj.alreadyScaledVariables
                fprintf('Scaling optimization variables:\n');

                %% Scale parameters
                fn=fields(obj.parameters);
                for f=1:length(fn)
                    [variable,scalingVariable]=getScalingVariables(obj,obj.parameters.(fn{f}));
                    if ~isempty(variable)
                        fprintf('  scaling %s by *= %s\n',name(variable),name(scalingVariable));
                        substituteVariable(obj,variable,scalingVariable.*variable);
                    end
                end
                %% Scale state variables
                fn=fields(obj.states);
                for f=1:length(fn)
                    % scale state
                    [variable,scalingVariable]=getScalingVariables(obj,obj.states.(fn{f}));
                    if ~isempty(variable)
                        fprintf('  scaling %s by *= %s\n',name(variable),name(scalingVariable));
                        substituteVariable(obj,variable,scalingVariable*variable);
                        if isinf(obj.states.(fn{f}).disturbanceInvVariance)
                            % scale equalities for dynamics
                            constr=obj.constraints{obj.states.(fn{f}).constraint};
                            if ~strcmp(type(constr),'iszero')
                                disp(constr)
                                error('unexpected constraint type');
                            else
                                operand=Tcalculus(operands(constr));
                                fprintf('  scaling %s equation by /= %s\n',name(variable),name(scalingVariable));
                                obj.constraints{obj.states.(fn{f}).constraint}=(operand/scalingVariable==0);
                            end
                        end
                    end
                end

            end
            obj.alreadyScaledVariables=true; % only do this once
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function setParameter(obj,name,value,saveValue)
        %    setParameter(obj,name,value)
        % or
        %    setParameter(obj,name)
        % set parameter with given value or from typical values
            if obj.parameters.(name).toOptimize
                type='optimization variable';
                fun=sprintf('setV_%s',name);
            else
                type='parameter            ';
                fun=sprintf('setP_%s',name);
            end
            if nargin>2
                if nargin<4
                    saveValue=true;
                end
                oldValue=obj.parameters.(name).initializationValue;
                if ~myisequal(size(oldValue),size(value))
                    error('value size mismatch: expected [%s], given [%s]\n',...
                          index2str(size(oldValue)),index2str(size(value)))
                end
                if saveValue
                    obj.parameters.(name).initializationValue=value;
                end
            else
                value=obj.parameters.(name).initializationValue;
            end
            if any(isnan(value))
                error('  error:  %s ''%s'' has NaN',type,name);
            end
            fprintf('  setting %s ''%s''  (%s)',type,name,index2str(size(value)));
            fprintf(' to %s',summarizeValues(obj,value))
            if obj.scaleVariables
                scaling=getScalingValue(obj,obj.parameters.(name));
                if any(scaling~=1)
                    fprintf('\n               scaling by %s',summarizeValues(obj,scaling));
                    value=value./scaling;
                    if numel(value)==1
                        fprintf(' to %.2e',value)
                    end
                end
            end
            if ~isempty(obj.solverObjFull)
                %fprintf(' calling ''%s''',fun);
                feval(fun,obj.solverObjFull,value);
            end
            if ~isempty(obj.solverObj1)
                if isfield(obj.optimizationVariablesType,name)
                    switch (obj.optimizationVariablesType.(name))
                      case {'disturbanceInvVariance','noiseInvVariance'}
                        fun=sprintf('setP_%s',name);
                        %fprintf('\n        calling ''%s(1)''',fun);
                        feval(fun,obj.solverObj1,value);
                        fun=sprintf('setV_%s',name);
                        %fprintf('\n        calling ''%s(2)''',fun);
                        feval(fun,obj.solverObj2,value);
                      case {'parameter'}
                        if obj.marginalizeParameters
                            fun=sprintf('setV_%s',name);
                        else
                            fun=sprintf('setP_%s',name);
                        end
                        %fprintf('\n        calling ''%s(1)''',fun);
                        feval(fun,obj.solverObj1,value);
                        if obj.marginalizeParameters
                            fun=sprintf('setP_%s',name);
                        else
                            fun=sprintf('setV_%s',name);
                        end
                        %fprintf('\n        calling ''%s(2)''',fun);
                        feval(fun,obj.solverObj2,value);
                      case {'state01_toMarginalize','state1_toMarginalize','state01','state1'}
                        fun=sprintf('setV_%s',name);
                        %fprintf('\n        calling ''%s(1)''',fun);
                        feval(fun,obj.solverObj1,value);
                        fun=sprintf('setP_%s',name);
                        %fprintf('\n        calling ''%s(2)''',fun);
                        feval(fun,obj.solverObj2,value);
                      otherwise
                        error('unknown type ''%s'' for optimization variable ''%s''\n',obj.optimizationVariablesType.(name),name);
                    end
                else
                    %fprintf('\n        calling ''%s(1)''',fun);
                    feval(fun,obj.solverObj1,value);
                    %fprintf('\n        calling ''%s(2)''',fun);
                    feval(fun,obj.solverObj2,value);
                end
            end
            fprintf('\n');
        end

        function setAllParameters(obj)

            fprintf('Setting all parameters/variables:\n');
            fn=fields(obj.parameters);
            for f=1:length(fn)
                setParameter(obj,fn{f});
            end
        end

        function setState(obj,name,value,saveValue)
        % setParameter(obj,name,value)
        % or
        % setParameter(obj,name)
            if nargin>2
                if nargin<4
                    saveValue=true;
                end
                oldValue=obj.states.(name).initializationValue;
                if ~myisequal(size(oldValue),size(value))
                    error('value size mismatch: state ''%s'' expected [%s], given [%s]',...
                          name,index2str(size(oldValue)),index2str(size(value)))
                end
                if saveValue
                    obj.states.(name).initializationValue=value;
                end
            else
                value=obj.states.(name).initializationValue;
            end
            if any(isnan(value))
                error('  error: initialization for state ''%s'' has NaN',name);
            end
            fprintf('  setting state ''%s''  (%s)',name,index2str(size(value)));
            fprintf(' to %s',summarizeValues(obj,value))
            if obj.scaleVariables
                scaling=getScalingValue(obj,obj.states.(name));
                if any(scaling~=1)
                    fprintf('\n               scaling by %s',summarizeValues(obj,scaling));
                    value=value./scaling;
                    if numel(value)==1
                        fprintf(' to %.2e',value);
                    else
                        fprintf(' to %s',summarizeValues(obj,value));
                    end
                end
            end
            fun=sprintf('setV_%s',name);
            if ~isempty(obj.solverObjFull)
                %fprintf(' calling ''%s''',fun);
                feval(fun,obj.solverObjFull,value);
            end
            if ~isempty(obj.solverObj1)
                fun=sprintf('setV_%s',name);
                %fprintf('\n        calling ''%s(1)''',fun);
                feval(fun,obj.solverObj1,value);
                fun=sprintf('setP_%s',name);
                %fprintf('\n        calling ''%s(2)''',fun);
                feval(fun,obj.solverObj2,value);
            end
            fprintf('\n');
        end

        function initializeAllStates(obj)

            fprintf('Initializing all states:\n');
            fn=fields(obj.states);
            for f=1:length(fn)
                setState(obj,fn{f})
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Call Solver
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function figs=plotCost(obj,outSolver,outExpression,fig0)

            function fig=inspectNoise(date,signal,name)

                subplot(1,2,1);
                plot(date,signal,'.-');grid on;
                if isa(date,'datetime')
                    labels2dates();
                end
                yl=ylim;
                ylabel(name);
                subplot(1,2,2);
                h=histogram(signal,20);grid on
                set(h,'orientation','horizontal');
                ylim(yl);
                title(sprintf('\\mu=%.5f, \\sigma=%.5f',...
                              mean(signal),std(signal)),'interpreter','tex');
                %xlabel(name);
end
            %% Solver output
            if outSolver.status==0
                fprintf('Solver succeed at iteration %3d in %7.3f ms, cost=%.3f\n',...
                        outSolver.iter,1e3*outSolver.time,outExpression.cost);
            else
                fprintf('Solver **failed** at iteration %3d in %7.3f ms, status = 0x%x\n',...
                        outSolver.iter,1e3*outSolver.time,outSolver.status);
            end

            %% Cost
            fn=fields(obj.logPDF);
            figs=[];
            for f=1:length(fn)
                sample=outExpression.(fn{f});
                figs(end+1,1)=clearFigure('figureNumber',fig0,'figureName',sprintf('%s',fn{f}));
                inspectNoise(1:length(sample),sample,fn{f});
                subplot(1,2,1)
                title(sprintf('-loglk=%5.0f, \\sigma=%.1e',...
                              outExpression.logPDFs(f),...
                              outExpression.(sprintf('%sStdDev',fn{f}))),'interpreter','tex');
                subplot(1,2,2)
                title(sprintf('sample mse^{1/2}=%.1e, \\mu=%.1e, \\sigma=%.1e',...
                              sqrt(mean(sample.^2)),mean(sample),std(sample)),'interpreter','tex');
                fig0=figs(end)+1;
            end

        end

        function report(obj,outSolver,outExpression,outStdErr)
        % Report results on parameters, states, and outputs
            reportCost(obj,outSolver,outExpression);
            if nargin<4
                reportParameters(obj,outSolver,outExpression);
                reportStates(obj,outSolver,outExpression);
            else
                reportParameters(obj,outSolver,outExpression,outStdErr);
                reportStates(obj,outSolver,outExpression,outStdErr);
            end
            reportOutputs(obj,outSolver,outExpression);

        end


        function reportCost(obj,outSolver,outExpression)

            %% Solver output
            if outSolver.status==0
                fprintf('Solver succeed at iteration %3d in %7.3f ms, cost=%.3f\n',...
                        outSolver.iter,1e3*outSolver.time,outExpression.cost);
            else
                fprintf('Solver **failed** at iteration %3d in %7.3f ms, status = 0x%x\n',...
                        outSolver.iter,1e3*outSolver.time,outSolver.status);
            end

            %% Cost
            fprintf('  Cost = %.3f:\n',outExpression.cost);
            fn=fields(obj.logPDF);
            for f=1:length(fn)
                sample=outExpression.(fn{f});
                fprintf('    %-25s: -likelihood = %8.3f (model std = %8.2e, sample mse^1/2=%8.2e, sample mean=%8.1e, sample std=%8.2e)\n',...
                        fn{f},outExpression.logPDFs(f),...
                        outExpression.(sprintf('%sStdDev',fn{f})),...
                        sqrt(mean(sample.^2)),mean(sample),std(sample));
            end
            if isfield(outExpression,'logPi')
                fprintf('   -#/2 log-det 2pi              :           %8.3f\n',full(-outExpression.logPi/2));
                fprintf('    1/2 log-det Hmarginalization :           %8.3f\n',full(outExpression.logdetHmarginalization/2));
            end

        end

        function reportParameters(obj,outSolver,outExpression,outStdErr)

            fprintf('  Parameter estimates:\n');
            fn=fields(obj.parameters);
            boundTolerance=1e-03;
            for f=1:length(fn)
                if obj.parameters.(fn{f}).toOptimize
                    value=outExpression.(fn{f});
                    scaling=getScalingValue(obj,obj.parameters.(fn{f}));
                    for jj=1:length(value)
                        if jj>2 && jj<length(value)-2
                            if jj==2
                                fprintf('        ...\n');
                            end
                            continue;
                        end
                        if jj==1
                            fprintf('    %-25s: %s',fn{f},summarizeValues(obj,value(jj)));
                        else
                            fprintf('                       (%4d): %s',jj,summarizeValues(obj,value(jj)));
                        end
                        if isfield(outExpression,[fn{f},'_posterioriStd'])
                            svalue=outExpression.([fn{f},'_posterioriStd']);
                            fprintf(' [std = %s]',summarizeValues(obj,svalue(jj)));
                        end
                        [lower,upper]=hasBounds(obj,obj.parameters.(fn{f}));
                        if isfinite(lower)
                            lower=outExpression.(sprintf('%s_lowerBound',fn{f}))(jj);
                        end
                        if isfinite(upper)
                            upper=outExpression.(sprintf('%s_upperBound',fn{f}))(jj);
                        end
                        fprintf(', constrained to [%9.2e,%9.2e] %s\n',lower,upper,...
                                hitBounds(obj,value,lower,upper));
                        if nargin>=4 && isfield(outStdErr,['errStd_',fn{f}])
                            errStd=outStdErr.(['errStd_',fn{f}]);
                            fprintf('                     err std : %s\n',summarizeValues(obj,errStd(jj)));
                        end
                    end
                    % Check whereInX
                    if false && norm2(value-outExpression.u_(obj.parameters.(fn{f}).whereInX).*scaling)>1e-10
                        value,outExpression.u_(obj.parameters.(fn{f}).whereInX).*scaling
                        error('value does not match solver vector');
                    end
                end
            end
        end

        function reportStates(obj,outSolver,outExpression,outStdErr)

            fprintf('  State estimates:\n');
            fn=fields(obj.states);
            boundTolerance=1e-03;
            for f=1:length(fn)
                value=outExpression.(fn{f});
                switch obj.optimizationVariablesType.(fn{f})
                  case {'state01','state01_toMarginalize'}
                  case {'state1','state1_toMarginalize'}
                    value=value(2:end);
                  otherwise
                    error('unexpected type ''%s''\n',obj.optimizationVariablesType.(fn{f}));
                end
                scaling=getScalingValue(obj,obj.states.(fn{f}));
                fprintf('    %-25s: %s',fn{f},summarizeValues(obj,value));
                [lower,upper]=hasBounds(obj,obj.states.(fn{f}));
                if isfinite(lower)
                    lower=outExpression.(sprintf('%s_lowerBound',fn{f}));
                end
                if isfinite(upper)
                    upper=outExpression.(sprintf('%s_upperBound',fn{f}));
                end
                fprintf(', contrained to [%9.2e,%9.2e] %s\n',lower,upper,...
                        hitBounds(obj,value,lower,upper));
                if nargin>=4
                    errStd=outStdErr.(['errStd_',fn{f}]);
                    fprintf('                     err std : %s\n',summarizeValues(obj,errStd));
                end
                % Check whereInX
                if false && norm2(value-outExpression.u_(obj.states.(fn{f}).whereInX).*scaling)>1e-10
                    value,outExpression.u_(obj.states.(fn{f}).whereInX)*scaling
                    error('value does not match solver vector');
                end
            end

        end

        function reportOutputs(obj,outSolver,outExpression)

            fprintf('  Outputs:\n');
            fn=fields(obj.outputs);
            for f=1:length(fn)
                value=outExpression.(fn{f});
                fprintf('    %-25s: %s\n',fn{f},summarizeValues(obj,value));
            end

        end

    end
end
