classdef Tmpc < handle
% Get help with
%   Tmpc help;
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

    properties (SetAccess = 'immutable');

        % properties set by constructor

        sampleTimeName=''  % When sampleTime is a symbolic
                           % variable, this holds the name of the
                           % variable; otherwise empty string.
        stateDerivativeFunction  % anonymous function with the
                                 % nominal state derivative
        controlDelay       % controlDelay (in sampleTime units) before a control can be applied

        solverName         % name of solver class
        solverObject       % solver object
        parameters         % solver symbolic parameters
                           % (including the parameters ones added
                           % by Tmpc);
        objective          % MPC cost
        nOutputExpressions % number of output expressions
                           % (exclusing the ones added by Tmpc)

        nStates            % length of state vector
        nControls          % length of control input vector
        nHorizon           % length of foward horizon

        parameterNames={}  % cell array with the names of all the parameters
        stateName=''       % name of symbolic state variable
        delayedControlName='' % name of symbolic past control variable
        futureControlName=''  % name of symbolic future control variable

        currentStateName   % name of symbolic variable with the current state
        delayedControls    % symbolic variable with the future controls
                           % that cannot be modified (only for controlDelay>0)
        optimizedControls  % symbolic variable with the future controls
                           % that need to be optimized

        times2add=1000;    % number of times to add to history
                           % matrices, for each call to extendHistory
    end

    properties %(GetAccess = {},SetAccess = 'immutable');
        history=struct('time',{[]},...      % vector of time instants
                       'state',{[]},...     % matrix with state history
                                    ...     % (one state per row, one time instant per column)
                       'control',{[]},...   % matrix with control history
                                      ...   % (one control per row, one time instant per column);
                       'objective',{[]},... % vector of MPC costs
                       'status',{[]},...    % vector of solver status
                       'iter',{[]},...      % vector of # of solver iterations
                       'stime',{[]},...     % vector of solver times
                       'currentIndex',{0}); % index of current time (in History arrays)
                         %   time has valid data from 1:currentIndex
                         %   state has valid data from (:,1:currentIndex)
                         %   control has valid data from (:,1:currentIndex+controlDelay-1)
                         %      the control (:,currentIndex-1) is held
                         %      from time(currentIndex-1) to time(currentIndex)
                         %   status/iter/stime have valid data from
                         %                (:,1:currentIndex+controlDelay-1)
                         %   status/iter/stime indices match the index of the computed control

        parameterValues={};
        parametersSet=[];
        controlSet=false;
        stateSet=false;
        sampleTimeValue=NaN % numerical value for sampling time;

    end

    methods

        %%%%%%%%%%%%%%%%%%%%%%%
        %% Class constructor %%
        %%%%%%%%%%%%%%%%%%%%%%%

        function obj=Tmpc(varargin)
        % Creates a matlab class for
        % 1) creating a Tcalculus solver for an MPC controller
        % 2) running the MPC optimization
        % 3) simulating the MPC controller

            declareParameter(...
                'Help', {
                    'Creates a matlab class for';
                    '1) creating a Tcalculus solver for an MPC controller';
                    '2) calling the MPC solver';
                    '3) simulating the MPC controller';
                    ' ';
                    'The MPC solver is designed to be executed after the "current state"'
                    'x(t) has been measured in order to determined the optimal';
                    'sequence of future states:';
                    '   [ x(t+Ts), x(t+2*Ts),  ... , x(t+nHorizon * Ts) ] ';
                    'and the corresponding future controls'
                    '   [ u(t+controlDelay*Ts),  ... , u(t+(nHorizon-1)*Ts) ]';
                    'where ';
                    '  Ts           = sampling time';
                    '  controlDelay = # of time steps until the fist control can be applied';
                    'The case controlDelay=0 corresponds to a scenario where';
                    'the optimization runs sufficiently fast so that u(t) can basically be'
                    'applied "simultaneously" with the measurement of x(t).';
                    ' '
                    'For positive values of controlDelay, the solver must be provided';
                    'with the control input values';
                    '   [ u(t), u(t+Ts),  ... , u(t+(controlDelay-1)*Ts) ]';
                    'that are scheduled to be applied before u(t+controlDelay*Ts).'

                        });

            %% Declare input parameters
            declareParameter(...
                'VariableName','solverClassname',...
                'Description',{
                    'Name of the solver class to be created.';
                              });

            declareParameter(...
                'VariableName','solverType',...
                'DefaultValue','matlab',...
                'AdmissibleValues',{'matlab','C'},...
                'Description',{
                    'Type of solver to be generated.';
                              });
            declareParameter(...
                'VariableName','reuseSolver',...
                'DefaultValue',true,...
                'AdmissibleValues',{true,false},...
                'Description',{
                    'When true, tries to reuse a previously generated solver';
                    'that was called with the same input parameters.';
                    'In this case, a unique suffix is appended to the solverClassname';
                              });

            declareParameter(...
                'VariableName','sampleTime',...
                'Description',{
                    'Sample time for the time discretization of the dynamics.'
                    'Can be either a double or a Tcalculus symbolic variable'
                    });
            declareParameter(...
                'VariableName','stateVariable',...
                'Description',{
                    'Symbolic Tcalculus variable representing the system''s future state.';
                    'Should be a matrix with size';
                    '     [ # states , horizon length ]';
                    'with value at time t equal to';
                    '   [ x(t+Ts), x(t+2*Ts),  ... , x(t+nHorizon * Ts) ] ';
                    'where Ts is the sample time.'
                    ' ';
                    'This variable can/should be used in the expressions';
                    '  . objective';
                    '  . constraints';
                    '  . outputExpressions';
                    ' '
                    'ATTENTION: note the "mismatch" between the times of'
                    '           stateVariable and controlVariable'
                    });
            declareParameter(...
                'VariableName','controlVariable',...
                'Description',{
                    'Symbolic Tcalculus variable representing the current and future'
                    'controls. Should be a matrix with size';
                    '     [ # controls , horizon length ]';
                    'with value at time t equal to';
                    '   [ u(t), u(t+Ts),  ... , u(t+(nHorizon-1)*Ts) ]';
                    'where u(s) is applied from time s to time s+Ts';
                    ' '
                    'This variable can/should be used in the expressions';
                    '  . objective';
                    '  . constraints';
                    '  . outputExpressions';
                    'but care must be placed in not having constraints on';
                    '   [ u(t), u(t+Ts),  ... , u(t+(controlDelay-1)*Ts) ]';
                    'that are scheduled to be applied before u(t+controlDelay*Ts)'
                    'and are therefore not optimization variables.'
                    ' '
                    'ATTENTION: note the "mismatch" between the times of'
                    '           stateVariable and controlVariable'
                    });
            declareParameter(...
                'VariableName','otherOptimizationVariables',...
                'DefaultValue',{},...
                'Description',{
                    'Cell array with any auxiliary optimization variables'
                    });
            declareParameter(...
                'VariableName','controlDelay',...
                'DefaultValue',0,...
                'Description',{
                    'Delay in sampleTime units before a control can be applied,';
                    'i.e., the first control that uses x(t) is ';
                    '          u(t + controlDelay * Ts),';
                    'This means that controlVariable has two distinct components:';
                    '1) one component that must have been determined before time t'
                    '   and cannot be changed by the optimization executed at time t'
                    '   that uses x(t):'
                    '   [ u(t), u(t+Ts),  ... , u(t+(controlDelay-1)*Ts) ]';
                    '2) one component that will be determined by the optimization';
                    '   executed at time t that uses x(t):'
                    '   [ u(t+controlDelay*Ts),  ... , u(t+(nHorizon-1)*Ts) ]';
                    });
            declareParameter(...
                'VariableName','stateDerivativeFunction',...
                'Description',{
                    'Anonymous function';
                    '     f(x,u,parameter1,parameter2,...)';
                    'that, for any integer N, must take as inputs'
                    '  . the state,'
                    '        [ x(t), x(t+Ts),  ... , x(t+(N-1) * Ts)] ';
                    '  . the control,'
                    '        [ u(t), u(t+Ts),  ... , u(t+(N-1) * Ts)] ';
                    '  . any parameters defined in ''parameters''';
                    'and returns the corresponding state derivatives';
                    '        [ dx(t), dx(t+Ts),  ... , dx(t+(N-1) * Ts)] ';
                    'that is used to obtain the stateVariable'
                    '        [ x(t+Ts), x(t+2*Ts),  ... , x(t+nHorizon * Ts) ] ';
                    'through trapesoidal integration with ZOH for u.'
                    ' '
                    'The function must accept all inputs to be either'
                    '  . matrices of doubles of appropriate sizes, or'
                    '  . Tcalculus symbolic matrices also of appropriate size.';
                    ' ';
                    'This derivative is used by MPC and should correspond to the'
                    'nominal dynamics.'
                    });
            declareParameter(...
                'VariableName','objective',...
                'Description',{
                    'Scalar Tcalculus symbolic expression with the MPC criteria'
                    'to be minimized.';
                              });
            declareParameter(...
                'VariableName','constraints',...
                'DefaultValue',{},...
                'Description',{
                    'Cell-array of Tcalculus symbolic objects representing MPC';
                    'constraints. Both equality and inequality constraints';
                    'are allowed.';
                              });

            declareParameter(...
                'VariableName','parameters',...
                'DefaultValue',{},...
                'Description',{
                    'Cell-array of Tcalculus symbolic objects representing the MPC ';
                    'parameters (must be given, not optimized).';
                      });
            declareParameter(...
                'VariableName','outputExpressions',...
                'Description',{
                    'Cell-array of Tcalculus symbolic objects representing the ';
                    'variables to be returned.';
                      });
            declareParameter(...
                'VariableName','solverParameters',...
                'DefaultValue',{},...
                'Description',{
                    'Cell-array of parameter names/values to be passed to the Tcalculus';
                    'solver generator.'
                      });


            %% Retrieve parameters and inputs

            [stopNow,params]=setParameters(nargout,varargin);
            if stopNow
                return
            end

            %% Constructor body

            if ~iscell(parameters)
                parameters={parameters};
            end
            % save parameters
            obj.stateDerivativeFunction=stateDerivativeFunction;
            obj.objective=objective;
            obj.nStates=size(stateVariable,1);
            obj.nControls=size(controlVariable,1);
            obj.nHorizon=size(controlVariable,2);
            obj.controlDelay=controlDelay;
            % check sizes
            if size(stateVariable,2)~=size(controlVariable,2)
                size(stateVariable),size(controlVariable),
                error('inconsistent time horizons in state and control variables');
            end
            if controlDelay<0
                error('controlDelay %d must be positive or zero',controlDelay);
            end
            if controlDelay>=obj.nHorizon
                error('controlDelay %d must be smaller than horizon length %s',...
                      controlDelay,obj.nHorizon);
            end
            % check sampleTime
            switch class(sampleTime)
              case 'double'
                if ~isequal(size(sampleTime),[1,1])
                    sampleTime,
                    error('sampleTime must be a scalar');
                end
                obj.sampleTimeValue=sampleTime;
              case 'Tcalculus'
                if ~isempty(size(sampleTime))
                    sampleTime,
                    error('sampleTime must be a scalar');
                end
                if ~strcmp(type(sampleTime),'variable')
                    sampleTime,
                    error('symbolic sampleTime must be a variable');
                end
                obj.sampleTimeName=name(sampleTime);
                found=false;
                for i=1:length(parameters)
                    if ~strcmp(class(parameters{i}),'Tcalculus')
                        disp(parameters{i})
                        error('Tmpc: parameter %d not of class Tcalculus',i);
                    end
                    if ~strcmp(type(parameters{i}),'variable')
                        disp(parameters{i})
                        error('Tmpc: parameter %d not a Tcalculus variable',i);
                    end
                    if strcmp(name(parameters{i}),obj.sampleTimeName)
                        found=true;
                        break;
                    end
                end
                if ~found
                    sampleTime,parameters
                    error('symbolic sampleTime must be one of the parameters');
                end
              otherwise
                sampleTime,
                error('sampleTime must be a double or a Tcalculus object');
            end
            % check sample time

            extendHistory(obj)

            %% Save names of parameters, control, and state
            % no need to remember values of parameters added by Tmpcmhe
            % (since not included in anonymous function)
            obj.parameterValues=cell(size(parameters));
            obj.parametersSet=false(size(parameters));
            obj.parameterNames=cell(size(parameters));
            for i=1:length(parameters)
                if ~strcmp(class(parameters{i}),'Tcalculus')
                    disp(parameters{i})
                    error('Tmpc: parameter %d not of class Tcalculus',i);
                end
                if ~strcmp(type(parameters{i}),'variable')
                    disp(parameters{i})
                    error('Tmpc: parameter %d not a Tcalculus variable',i);
                end
                obj.parameterNames{i}=name(parameters{i});
            end
            obj.stateName=name(stateVariable);
            obj.parameters=parameters;

            %% create dynamics
            % currentState      = x(t)
            % stateVariable     = [x(t+Ts) x(t+Ts) ... x(t+T*Ts) ]
            % thisState         = [x(t) x(t+Ts) ... x(t+(T-1)*Ts)]
            % delayedControls   = [u(t) ... u(t+(controlDelay-1)*Ts)
            % optimizedControls = [u(t+controlDelay*Ts ... u(t+(T-1)*Ts)]
            % thisControl       = [u(t) u(t+Ts) ... u(t+(T-1)*Ts)]
            %                   = [delayedControls optimizedControls]
            currentState=Tvariable([obj.stateName,'_initial'],[obj.nStates,1]);
            thisState=[currentState,stateVariable(:,1:end-1)];
            if obj.controlDelay>0
                % first controlDelay controls are parameters, whereas remaining are
                % optimization variables
                obj.delayedControls=Tvariable([name(controlVariable),'_delayed'],...
                                              [obj.nControls,obj.controlDelay]);
                obj.delayedControlName=name(obj.delayedControls);
                obj.optimizedControls=Tvariable([name(controlVariable),'_optimized'],...
                                                [obj.nControls,obj.nHorizon-obj.controlDelay]);
                obj.futureControlName=name(obj.optimizedControls);
                thisControl=[obj.delayedControls,obj.optimizedControls];
                obj.objective=substitute(obj.objective,controlVariable,thisControl);
                constraints=substitute(constraints,controlVariable,thisControl);
                outputExpressions=substitute(outputExpressions,controlVariable,thisControl);
                obj.parameters{end+1}=obj.delayedControls;
                obj.parameterNames{end+1}=name(obj.delayedControls);
                obj.parametersSet(end+1)=false;
            else
                thisControl=controlVariable;
                obj.optimizedControls=controlVariable;
                obj.futureControlName=name(controlVariable);
            end
            obj.parameters{end+1}=currentState;
            obj.currentStateName=name(currentState);
            obj.parameterNames{end+1}=name(currentState);
            obj.parametersSet(end+1)=false;

            try
                if 1
                    %% Trapesoidal integration for u with ZOH
                    %    (x_{k+1} - x_k)/Ts == (f(x_{k+1},u_k)+f(x_k,u_k))/2
                    dynamics=(stateVariable-thisState)==.5*sampleTime*(...
                        stateDerivativeFunction(stateVariable,...
                                                thisControl,...
                                                parameters{:})...
                        +stateDerivativeFunction(thisState,...
                                                 thisControl,...
                                                 parameters{:}));
                else
                    %% Forward Euler integration
                    dynamics=(stateVariable==thisState+...
                              sampleTime*stateDerivativeFunction(thisState,...
                                                                 thisControl,...
                                                                 parameters{:}));
                end

                % add output expressions for debug
                % outputExpressions{end+1}=thisState;
                % outputExpressions{end+1}=stateDerivativeFunction(stateVariable,...
                %                                          controlVariable,...
                %                                          parameters{:});
                % outputExpressions{end+1}=sampleTime;
            catch me
                fprintf('error in appyling ''stateDerivativeFunction'' to symbolic variables\n');
                rethrow(me);
            end


            % add dynamics to constraints
            constraints(end+1)={dynamics};

            % add state & control to outputs
            obj.nOutputExpressions=length(outputExpressions);
            outputExpressions=[outputExpressions,{obj.optimizedControls,stateVariable,obj.objective}];

            %% Call solver
            switch solverType
              case 'matlab'
                solver=@class2optimizeCS;
              case 'C'
                solver=@cmex2optimizeCS;
              otherwise
                error('unknown solver type "%s"',solverType);
            end

            if reuseSolver
                executeScript='asneeded';
            else
                executeScript='yes';
            end

            if ~isempty(solverParameters)
                if ~iscell(solverParameters)
                    solverParameters,
                    error('solverParameters must be a cell array');
                end
                if size(solverParameters,1)~=1 && size(solverParameters,2)~=1
                    solverParameters,
                    error('solverParameters must be a vector of cells');
                end
                for i=1:2:length(solverParameters)
                    if ~ischar(solverParameters{i})
                        solverParameters,
                        error('solverParameters must be a vector of cells of the form {parameter name, parameter value, ...}');
                    end
                end
            end

            [classname,code]=solver('pedigreeClass',solverClassname,...
                                    'executeScript',executeScript,...
                                    'objective',obj.objective,...
                                    'optimizationVariables',...
                                    {obj.optimizedControls,stateVariable,otherOptimizationVariables{:}},...
                                    'constraints',constraints,...
                                    'parameters',obj.parameters,...
                                    'outputExpressions',outputExpressions,...
                                    solverParameters{:});

            obj.solverName=getValue(classname);
            obj.solverObject=feval(obj.solverName);

        end

        %%%%%%%%%%%%%%%%%%%
        %% Other methods %%
        %%%%%%%%%%%%%%%%%%%

        function extendHistory(obj)
        % extendHistory(obj)
        %
        % adds time-steps to the History matrices

            obj.history.time=[obj.history.time;nan(obj.times2add,1)];
            obj.history.state=[obj.history.state,nan(obj.nStates,obj.times2add)];
            obj.history.control=[obj.history.control,nan(obj.nControls,obj.times2add)];
            obj.history.objective=[obj.history.objective;nan(obj.times2add,1)];
            obj.history.status=[obj.history.status;nan(obj.times2add,1)];
            obj.history.iter=[obj.history.iter;nan(obj.times2add,1)];
            obj.history.stime=[obj.history.stime;nan(obj.times2add,1)];
        end

        function setParameter(obj,name,value)
        % setParameter(obj,name,value);
        %
        % Sets the value of parameter called 'name' with the given value

            k=find(strcmp(name,obj.parameterNames));
            if isempty(k)
                obj.parameterNames,
                error('unknown parameter "%s"',name);
            elseif length(k)>1
                obj.parameterNames,
                error('multiple parameters with the same name "%s"',name);
            else
                if 0
                    fprintf('Setting %s\n',['setP_',name]);
                    disp(value)
                end
                feval(['setP_',name],obj.solverObject,value);
                if k<=length(obj.parameterValues)
                    % no need to remember values of parameters added by Tmpc
                    % (since not included in anonymous function)
                    obj.parameterValues{k}=value;
                end
                obj.parametersSet(k)=true;
                if ~isempty(obj.sampleTimeName) && strcmp(name,obj.sampleTimeName)
                    if ~isnan(obj.sampleTimeValue)
                        error('sampleTime value can only bet set once');
                    end
                    obj.sampleTimeValue=value;
                end
            end
        end

        function setSolverInputStart(obj,value)
        % setSolverInputStart(obj,value)
        %
        % Sets the value for the input that will be used to initialize
        % the solver.
        %
        % A warm start from the previous MPC optimization often helps.

            feval(['setV_',obj.futureControlName],obj.solverObject,...
                  value(:,obj.controlDelay+1:end));
            obj.controlSet=true;
        end

        function setSolverStateStart(obj,value)
        % setSolverStateStart(obj,value)
        %
        % Sets the value for the state that will be used to
        % initialize the solver.
        %
        % A warm start compatible with the value given to
        % |setInputStart| often helps.


            feval(['setV_',obj.stateName],obj.solverObject,value(:,2:end));
            obj.stateSet=true;
        end


        function setInitialState(obj,tinit,xinit,uinit);
        % setInitialState(obj,tinit,xinit);
        %  or
        % setInitialState(obj,tinit,xinit,uinit);
        %
        % Initilizes the system with
        %   tinit - initial time
        %   xinit - state at time tinit
        %   uinit - control at times tinit..tinit+Ts*(controlDelay-1)
        %           (only needed if controlDelay>0)
            if ~isequal(size(tinit),[1,1]);
                error('time value value must be a scalar');
            end
            if ~isequal(size(xinit),[obj.nStates,1])
                error('state value must be a column %dx1 array',obj.nStates);
            end
            if nargin<4
                uinit=zeros(obj.nControls,obj.controlDelay);
            end
            if ~isequal(size(uinit),[obj.nControls,obj.controlDelay])
                error('control value must be a column %dx%d array',...
                      obj.nControls,obj.controlDelay);
            end
            obj.history.time(1)=tinit;
            obj.history.state(:,1)=xinit;
            obj.history.control(:,1:obj.controlDelay)=uinit;
            obj.history.currentIndex=1;
        end

        function state=setSolverWarmStart(obj,control)
        % state=setSolverWarmStart(obj,control)
        %
        % Given
        %    control = Sequence of future controls (past the delay).
        %            Should be a matrix with size
        %                  # controls x (horizon length - controlDelay)
        %            with values
        %                 [ u(t+controlDelay*Ts),  ... , u(t+(nHorizon-1)*Ts) ]';
        %            where t is the current time and u(s) is applied from time s to time s+Ts
        % Computes
        %    state = Sequence of future states.
        %            Should be a matrix with size
        %                  # states x (horizon length)
        %            with values
        %               [ x(t), x(t+Ts),  ... , x(t+(nHorizon-1)*Ts) ]% ;
        % and initializes the solver with these sequences.
        %
        % The differential equation is integrated using forward
        % Euler.
        %
        % ATTENTION: This function does not enforce state
        % constraints. If state constraints are being used one
        % should move the output |state| away from the constraints
        % and reset it using |setSolverStateStart|

            if ~isequal(size(control),[obj.nControls,obj.nHorizon-obj.controlDelay])
                error('size of control [%d,%d] does not match expected size [%d,%d]',...
                      size(control,1),size(control,2),obj.nControls,obj.nHorizon-obj.controlDelay);
            end
            if ~all(obj.parametersSet(1:length(obj.parameterValues)))
                fprintf('did not set values for:');
                obj.parametersSet
                disp(obj.parameterNames(~obj.parametersSet &...
                                        (1:length(obj.parameterNames))<=length(obj.parameterValues)));
                error('must set all parameter values before calling setSolverWarmStart');
            end
            if obj.history.currentIndex<1
                error('must initialize state');
            end
            state=nan(obj.nStates,obj.nHorizon+1);
            state(:,1)=obj.history.state(:,obj.history.currentIndex);
            if obj.controlDelay>0
                control=[obj.history.control(:,obj.history.currentIndex:obj.history.currentIndex+obj.controlDelay-1),control];
            end
            setParameter(obj,obj.currentStateName,state(:,1));
            for k=1:obj.nHorizon
                try
                    % Forward Euler integration for the dynamics
                    state(:,k+1)=state(:,k)+obj.sampleTimeValue*...
                        obj.stateDerivativeFunction(state(:,k),control(:,k),...
                                                    obj.parameterValues{:});
                catch me
                    fprintf('error in appyling ''stateDerivativeFunction'' to numerical values\n');
                    rethrow(me);
                end
            end
            feval(['setV_',obj.stateName],obj.solverObject,state(:,2:end));
            obj.stateSet=true;
            if obj.controlDelay>0
                setParameter(obj,obj.delayedControlName,...
                                 control(:,1:obj.controlDelay));
            end
            feval(['setV_',obj.futureControlName],obj.solverObject,...
                  control(:,obj.controlDelay+1:end));
            obj.controlSet=true;
        end

        function [solution,varargout]=solve(obj,mu0,maxIter,saveIter)
        % [solution,y1,y2,...]=solve(obj,mu0,maxIter,saveIter)
        %
        % Given
        %   mu0      - initial value for the gap variable
        %   maxIter  - maximum number of iterations
        %   saveIter - iteration # when to save the "hessian" matrix
        %              (see class2optimizeCS/cmex2optimizeCS)
        % Calls the solver and returns:
        %   solution.control  - optimal control signals (past controlDelay)
        %   solution.state  - optimal state
        %   solution.status - solver exit status
        %   solution.iter   - number of solver iterations
        %   solution.time   - time taken by solver
        %   y1,y2,...       - output expressions (passed to the constructor)
        %
        % Should check if all parameters have been set
            if ~all(obj.parametersSet)
                fprintf('did not set values for:');
                disp(obj.parameterNames(~obj.parametersSet));
                error('must set all parameter values before calling solver');
            end
            if ~obj.controlSet
                error('must set initial value for control before calling solver');
            end
            if ~obj.stateSet
                error('must set initial value for state before calling solver');
            end

            [solution.status,solution.iter,solution.time]=...
                solve(obj.solverObject,mu0,int32(maxIter),int32(saveIter));

            varargout=cell(obj.nOutputExpressions,1);
            [varargout{:},solution.control,solution.state,solution.objective]=getOutputs(obj.solverObject);
        end

        function [t,u0_warm,u_applied]=applyControls(obj,solution,ufinal,stateDerivativeFunctionReal);
        % [t,u0_warm,u_applied]=applyControls(obj,solution,ufinal,stateDerivativeFunctionReal);
        %
        % Given a solver solution
        % . applies the first control
        % . computed the state evolution using ODE23
        % and returns
        % . time for the first state value for the next MPC optimization
        %   (usefull to set up reference signals)
        % . warm start for the next MPC optimization,
        %   inserting the control 'ufinal' at the end of the interval.

            if (nargin<4)
                stateDerivativeFunctionReal=obj.stateDerivativeFunction;
            end

            t=obj.history.time(obj.history.currentIndex);
            u_applied=solution.control(:,1);
            t=t+obj.sampleTimeValue;
            obj.history.currentIndex=obj.history.currentIndex+1;
            obj.history.time(obj.history.currentIndex)=t;
            obj.history.control(:,obj.history.currentIndex+obj.controlDelay-1)=...
                solution.control(:,1);
            obj.history.objective(obj.history.currentIndex+obj.controlDelay-1)=...
                solution.objective;
            obj.history.status(obj.history.currentIndex+obj.controlDelay-1)=...
                solution.status;
            obj.history.iter(obj.history.currentIndex+obj.controlDelay-1)=...
                solution.iter;
            obj.history.stime(obj.history.currentIndex+obj.controlDelay-1)=...
                solution.time;
            if size(obj.history.state,2)<obj.history.currentIndex
                extendHistory(obj)
            end
            try
                if 1
                    % Integration for the dynamics using ode23
                    [tout,yout]=ode23(...
                        @(t,x)stateDerivativeFunctionReal(x,...
                                                          obj.history.control(:,obj.history.currentIndex-1),...
                                                          obj.parameterValues{:}),...
                        [obj.history.time(obj.history.currentIndex-1),obj.history.time(obj.history.currentIndex)],...
                        obj.history.state(:,obj.history.currentIndex-1));
                    obj.history.state(:,obj.history.currentIndex)=yout(end,:)';
                else
                    % Forward Euler integration for the dynamics
                    obj.history.state(:,obj.history.currentIndex)=...
                        obj.history.state(:,obj.history.currentIndex-1)+obj.sampleTimeValue*...
                        stateDerivativeFunctionReal(obj.history.state(:,obj.history.currentIndex-1),...
                                                    obj.history.control(:,obj.history.currentIndex-1),...
                                                    obj.parameterValues{:});
                end
            catch me
                fprintf('error in appyling ''stateDerivativeFunctionReal'' to numerical values\n');
                rethrow(me);
            end

            u0_warm=[solution.control(:,2:end),ufinal*ones(1,1)];
            t=t+obj.sampleTimeValue;

            obj.stateSet=false;
            obj.controlSet=false;

        end

        function history=getHistory(obj);
        % [t,x,u,status,iter,stime]=getHistory(obj);
        %
        % Returns the history of
        %    history.t      - sampling times
        %    history.u      - control values
        %    history.x      - state values
        %    history.status - solver status
        %    history.iter   - # of solver iterations
        %    history.stime  - solver times

            history.t=obj.history.time(1:obj.history.currentIndex-1);
            history.u=obj.history.control(:,1:obj.history.currentIndex-1);
            history.x=obj.history.state(:,1:obj.history.currentIndex-1);
            history.objective=obj.history.objective(1:obj.history.currentIndex-1);
            history.status=obj.history.status(1:obj.history.currentIndex-1);
            history.iter=obj.history.iter(1:obj.history.currentIndex-1);
            history.stime=obj.history.stime(1:obj.history.currentIndex-1);
        end
    end
end
