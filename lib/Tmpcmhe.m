classdef Tmpcmhe < handle
% Get help with 
%   Tmpcmhe help;
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
        stateDerivativeFunction % anonymous function with the state derivative
        outputFunction     % anonymous function with the state derivative
        controlDelay       % controlDelay (in sampleTime units) before a control can be applied
        
        solverName         % name of solver class
        solverObject       % solver object
        parameters         % solver symbolic parameters
                           % (including the parameters ones added by Tmpc);
        objective          % MPC cost
        nOutputExpressions % number of output expressions
                           % (exclusing the ones added by Tmpc)

        nStates            % length of state vector
        nControls          % length of control input vector
        nDisturbances      % length of disturbance input vector
        nOutputs           % length of output vector
        nBackwardHorizon   % length of backward horizon (L)
        nForwardHorizon    % length of foward horizon (T)
        
        parameterNames={}    % cell array with the names of all the parameters
        initialStateName=''  % name of symbolic initial state variable
        latentStateName=''   % name of symbolic latent state variable
        pastOutputName=''    % name of symbolic past output variable
        pastControlName=''   % name of symbolic past control variable
        futureControlName='' % name of symbolic future control variable
        disturbanceName=''   % names of symbolic disturbance variable

        currentStateName   % name of symbolic variable with the current state 
        delayedControls    % symbolic variable with the future controls
                           % that cannot be modified (only for controlDelay>0)

        times2add=1000;    % number of times to add to history
                           % matrices, for each call to extendHistory
    end
        
    properties %(GetAccess = {},SetAccess = 'immutable');
        history=struct('time',{[]},...      % vector of time instants
                       'state',{[]},...     % matrix with state history 
                                    ...     % (one state per row, one time instant per column)
                       'control',{[]},...   % matrix with control history
                                      ...   % (one control per row, one time instant per column);
                       'disturbance',{[]},... % matrix with disturbance history
                                      ...   % (one disturbance per row, one time instant per column);
                       'objective',{[]},... % vector of MPC-MHE costs
                       'output',{[]},...    % matrix with output history
                                      ...   % (one output per row, one time instant per column);
                       'status',{[]},...    % vector of solver status
                       'iter',{[]},...      % vector of # of solver iterations
                       'stime',{[]},...     % vector of solver times
                       'currentIndex',{0}); % index of current time (in History arrays)
                         %   time has valid data from 1:currentIndex
                         %   state has valid data from (:,1:currentIndex)
                         %   control has valid data from (:,1:currentIndex+controlDelay-1)
                         %      the control (:,currentIndex-1) is held
                         %      from time(currentIndex-1) to time(currentIndex)
                         %   disturbance has valid data from (:,1:currentIndex+controlDelay-1)
                         %      the disturbance (:,currentIndex-1) is held
                         %      from time(currentIndex-1) to time(currentIndex)
                         %   status/iter/stime have valid data from
                         %                (:,1:currentIndex+controlDelay-1)
                         %   status/iter/stime indices match the index of the computed control
        
        parameterValues={};
        parametersSet=[];
        initialStateSet=false;
        latentStateSet=false;
        futureControlSet=false;
        disturbanceSet=false;
        stateSet=false;
        sampleTimeValue=NaN % numerical value for sampling time;
        
    end
    
    methods
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %% Class constructor %%
        %%%%%%%%%%%%%%%%%%%%%%%
        
        function obj=Tmpcmhe(varargin)	
        % Creates a matlab class for
        % 1) creating a Tcalculus solver for an MPC controller
        % 2) running the MPC optimization
        % 3) simulating the MPC controller
                    
            declareParameter(...
                'Help', {
                    'Creates a matlab class for';
                    '1) creating a Tcalculus solver for an MPC-MHE controller';
                    '2) calling the MPC-MHE solver';
                    '3) simulating the MPC-MHE controller';
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
                    '     [ # states , L + T + 1 ]';
                    'with value at time t equal to';
                    '     [ x(t-L*Ts), ..., x(t), ..., x(t+T*Ts) ]';
                    'where T is the forward horizon,'
                    '      L is the backwards horizon,'
                    '      Ts is the sample time.'
                    ' ';
                    'This variable can/should be used in the expressions';
                    '  . objective';
                    '  . constraints';
                    '  . outputExpressions';
                    });
            declareParameter(...
                'VariableName','pastOutputVariable',...
                'Description',{
                    'Symbolic Tcalculus variable representing the past'
                    'outputs. Should be a matrix with size';
                    '     [ # outputs , L + 1 ]';
                    'with value at time t equal to';
                    '     [ y(t-L*Ts), ..., y(t) ]'
                    ' '
                    'This variable can/should be used in the expressions';
                    '  . objective';
                    '  . constraints';
                    '  . outputExpressions';
                    });
            declareParameter(...
                'VariableName','pastControlVariable',...
                'Description',{
                    'Symbolic Tcalculus variable representing the past'
                    'controls. Should be a matrix with size';
                    '     [ # controls , L + controlDelay ]';
                    'with value at time t equal to';
                    '     [ u(t-L*Ts), ..., u(t+(delay-1)*Ts) ]';
                    'where u(s) is applied from time s to time s+Ts';
                    ' '
                    'This variable can/should be used in the expressions';
                    '  . objective';
                    '  . constraints';
                    '  . outputExpressions';
                    });
            declareParameter(...
                'VariableName','futureControlVariable',...
                'Description',{
                    'Symbolic Tcalculus variable representing the current and future'
                    'controls. Should be a matrix with size';
                    '     [ # controls , T - controlDelay ]';
                    'with value at time t equal to';
                    '     [ u(t+delay*Ts), ..., u(t+(T-1)*Ts) ]'
                    'where u(s) is applied from time s to time s+Ts';
                    ' '
                    'This variable can/should be used in the expressions';
                    '  . objective';
                    '  . constraints';
                    '  . outputExpressions';
                    });
            declareParameter(...
                'VariableName','disturbanceVariable',...
                'Description',{
                    'Symbolic Tcalculus variable representing the past, current, and future'
                    'disturbances. Should be a matrix with size';
                    '     [ # disturbances , L + T ]';
                    'with value at time t equal to';
                    '     [ d(t-L*Ts, ..., x(t), ..., x(t+(T-1)*Ts) ]'
                    'where d(s) is applied from time s to time s+Ts';
                    ' '
                    'This variable can/should be used in the expressions';
                    '  . objective';
                    '  . constraints';
                    '  . outputExpressions';
                    });
            declareParameter(...
                'VariableName','stateDerivativeFunction',...
                'Description',{
                    'Anonymous function';
                    '     f(x,u,d,parameter1,parameter2,...)';
                    'that, for any integer N, must take as inputs'
                    '  . the state,'
                    '        [ x(t), x(t+Ts),  ... , x(t+(N-1) * Ts)] ';
                    '  . the control,'
                    '        [ u(t), u(t+Ts),  ... , u(t+(N-1) * Ts)] ';
                    '  . the disturbance,'
                    '        [ d(t), d(t+Ts),  ... , d(t+(N-1) * Ts)] ';
                    '  . any parameters defined in ''parameters''';
                    'and returns the corresponding state derivatives';
                    '        [ dx(t), dx(t+Ts),  ... , dx(t+(N-1) * Ts)] ';
                    'through trapesoidal integration with ZOH for u and d.'
                    '        [ x(t+Ts), x(t+2*Ts),  ... , x(t+nHorizon * Ts) ] ';
                    ' ';
                    'The function must accept all inputs can be either'
                    '  . matrices of doubles of appropriate sizes, or'
                    '  . Tcalculus symbolic matrices also of appropriate size.';
                    });
            declareParameter(...
                'VariableName','outputFunction',...
                'Description',{
                    'Anonymous function';
                    '     g(x,u,d,parameter1,parameter2,...)';
                    'that, for any integer N, must take as inputs'
                    '  . the state,'
                    '        [ x(t), x(t+Ts),  ... , x(t+(N-1) * Ts)] ';
                    '  . the control,'
                    '        [ u(t), u(t+Ts),  ... , u(t+(N-1) * Ts)] ';
                    '  . any parameters defined in ''parameters''';
                    'and returns the corresponding output';
                    '        [ y(t), y(t+Ts),  ... , y(t+(N-1) * Ts)] ';
                    'The function must accept all inputs to be either'
                    '  . matrices of doubles of appropriate sizes, or'
                    '  . Tcalculus symbolic matrices also of appropriate size.';
                    ' '
                    'ATTENTION: When controlDelay=0, the output at time y(t) cannot'
                    '           depend on the current values of u(t) and d(t).'

                    });
            declareParameter(...
                'VariableName','objective',...
                'Description',{
                    'Scalar Tcalculus symbolic expression with the MPC-MHE criteria.';
                              });
            declareParameter(...
                'VariableName','controlConstraints',...
                'DefaultValue',{},...
                'Description',{
                    'Cell-array of Tcalculus symbolic objects representing the';
                    'control constraints. Both equality and inequality constraints';
                    'are allowed.';
                              });
            declareParameter(...
                'VariableName','disturbanceConstraints',...
                'DefaultValue',{},...
                'Description',{
                    'Cell-array of Tcalculus symbolic objects representing the';
                    'disturbance constraints. Both equality and inequality constraints';
                    'are allowed.';
                              });
            
            declareParameter(...
                'VariableName','parameters',...
                'DefaultValue',{},...
                'Description',{
                    'Cell-array of Tcalculus symbolic objects representing the MPC-MHE';
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
            obj.outputFunction=outputFunction;
            obj.nStates=size(stateVariable,1);
            obj.nControls=size(futureControlVariable,1);
            obj.nDisturbances=size(disturbanceVariable,1);
            obj.nOutputs=size(futureControlVariable,1);
            
            obj.controlDelay=size(pastControlVariable,2)-size(pastOutputVariable,2)+1;
            obj.nBackwardHorizon=size(pastOutputVariable,2)-1;
            obj.nForwardHorizon=size(futureControlVariable,2)+obj.controlDelay;

            % check sizes
            if size(stateVariable,2)~=obj.nForwardHorizon+obj.nBackwardHorizon+1
                obj
                error('inconsistent size of state variable [%s]',...
                      index2str(size(stateVariable)));
            end
            if size(pastControlVariable,2)~=obj.nBackwardHorizon+obj.controlDelay
                obj
                error('inconsistent length of past control variable [%s]',...
                      index2str(size(pastControlVariable)));
            end
            if size(disturbanceVariable,2)~=obj.nForwardHorizon+obj.nBackwardHorizon
                obj
                error('inconsistent length of disturbance variable [%s]',...
                      index2str(size(size(disturbanceVariable))));
            end
            if size(pastControlVariable,1)~=obj.nControls
                error('inconsistent sizes of pastControlVariable [%s] and futureControlVariable [%s]',...
                      size(pastControlVariable),size(futureControlVariable));;
            end

            if obj.controlDelay<0
                error('controlDelay %d must be positive or zero',obj.controlDelay);
            end
            if obj.controlDelay>=obj.nForwardHorizon
                error('controlDelay %d must be smaller than horizon length %s',...
                      obj.controlDelay,obj.nForwardHorizon);
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
                
            extendHistory(obj);

            %% Save names of parameters, control, and state
            obj.pastOutputName=name(pastOutputVariable);
            obj.pastControlName=name(pastControlVariable);
            obj.futureControlName=name(futureControlVariable);
            obj.disturbanceName=name(disturbanceVariable);

            % no need to remember values of parameters added by Tmpcmhe
            % (since not included in anonymous function)
            obj.parameterValues=cell(size(parameters));
            obj.parameters=[parameters(:);
                            {pastOutputVariable;
                             pastControlVariable}];
            obj.parametersSet=false(size(obj.parameters));
            obj.parameterNames=cell(size(obj.parameters));
            for i=1:length(obj.parameters)
                obj.parameterNames{i}=name(obj.parameters{i});
            end
            
            %% create dynamics
            % separate initial state (P2's optimization variable)
            % from the rest of the state (latent variable)
            initialState=Tvariable([name(stateVariable),'_initial'],[obj.nStates,1]);
            obj.initialStateName=name(initialState);
            nextState=Tvariable([name(stateVariable),'_next'],...
                                [obj.nStates,obj.nForwardHorizon+obj.nBackwardHorizon]);
            obj.latentStateName=name(nextState);
            allState=[initialState,nextState];
            
            obj.objective=substitute(obj.objective,stateVariable,allState);
            controlConstraints=substitute(controlConstraints,stateVariable,allState);
            disturbanceConstraints=substitute(disturbanceConstraints,stateVariable,allState);
            outputExpressions=substitute(outputExpressions,stateVariable,allState);

            previousState=[initialState,nextState(:,1:end-1)];
            previousControl=[pastControlVariable,futureControlVariable];
            
            try
                if 1
                    %% Trapesoidal integration for u with ZOH
                    %    (x_{k+1} - x_k)/Ts == (f(x_{k+1},u_k,d_k)+f(x_k,u_k,d_k))/2
                    dynamics=(nextState-previousState)==.5*sampleTime*(...
                        stateDerivativeFunction(previousState,...
                                                previousControl,...
                                                disturbanceVariable,...
                                                parameters{:})...
                        +stateDerivativeFunction(nextState,...
                                                 previousControl,...
                                                 disturbanceVariable,...
                                                 parameters{:}));
                else
                    % Forward Euler integration for the dynamics
                    dynamics=(nextState==previousState+...
                              sampleTime*stateDerivativeFunction(previousState,...
                                                                 previousControl,...
                                                                 disturbanceVariable,...
                                                                 parameters{:}));
                end
            catch me
                fprintf('error in appyling ''stateDerivativeFunction'' to symbolic variables\n');
                rethrow(me);
            end
                
            disturbanceVariable={disturbanceVariable;initialState};
                            
            % add state & control to outputs
            obj.nOutputExpressions=length(outputExpressions);
            outputExpressions=[outputExpressions,...
                               {futureControlVariable,...
                                disturbanceVariable{:},allState,obj.objective}];

            %% Call solver
            switch solverType
              case 'matlab'
                solver=@class2equilibriumLatentCS;
              case 'C'
                solver=@cmex2equilibriumLatentCS;
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
            if 0
                objective,
                futureControlVariable,
                disturbanceVariable{:}
                dynamics
                %objective=objective-norm2(disturbanceVariable{2}-Tconstant([50;50],[2,1]))
            end
            [classname,code]=solver('pedigreeClass',solverClassname,...
                                    'executeScript',executeScript,...
                                    'P1objective',obj.objective,...
                                    'P2objective',-obj.objective,...
                                    'P1optimizationVariables',{futureControlVariable},...
                                    'P2optimizationVariables',disturbanceVariable,...
                                    'P1constraints',controlConstraints,...
                                    'P2constraints',disturbanceConstraints,...
                                    'latentVariables',{nextState},...
                                    'latentConstraints',{dynamics},...
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
            obj.history.disturbance=[obj.history.disturbance,nan(obj.nDisturbances,obj.times2add)];
            obj.history.output=[obj.history.output,nan(obj.nOutputs,obj.times2add)];
            obj.history.objective=[obj.history.objective;nan(obj.times2add,1)];
            obj.history.status=[obj.history.status;nan(obj.times2add,1)];
            obj.history.iter=[obj.history.iter;nan(obj.times2add,1)];
            obj.history.stime=[obj.history.stime;nan(obj.times2add,1)];
        end
    
        function setParameter(obj,name,value)
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
                    % no need to remember values of parameters added by Tmpcmhe
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
    
        function [t,t_y_past,y_past,t_u_past,u_past,t_state,state]=...
                setInitialState(obj,t_init,x_init,u_init,d_init,n_init);
        % [t,t_y_past,y_past,t_u_past,u_past,t_state,state]=...
        %        setInitialState(obj,t_init,x_init,u_init,d_init,n_init)
        %
        % Initilizes the system with 
        %   t_init - initial time t==t_init+L*Ts
        %   x_init - initial state x(t_init)
        %   u_init - past control     [u(t_init),...,u(t_init+(L+delay-1)*Ts)]
        %   d_init - past disturbance [d(t_init),...,d(t_init+(L-1)*Ts)]
        %   n_init - part additive measurement noise [n(t_init),...,n(t_init+L*Ts)]
        % Returns
        %   t      - current time t=(t_init+L*Ts)
        %   u_past = Sequence of past control controls.
        %            Should be a matrix with size
        %                 [ # controls , L + controlDelay ]
        %            with value at (the new) time t equal to
        %                 [ u(t_init), ..., u(t_init+(L+delay-1)*Ts) ]
        %            where u(s) is applied from time s to time s+Ts
        %   t_u_past = Sequence of times for u_past.
        %   y_past = Sequence of past measurements.
        %            Should be a matrix with size
        %                 [ # outputs , L + 1 ]
        %            with value at (the new) time t equal to
        %                 [y(t_init), ..., y(t_init+L*Ts) ]
        %   t_y_past = Sequence of times for y_past.
        %   x_warm = Initial state for a warm start at the next iteration
        %            Should be a matrix with size
        %                 [ # states , 1 ]
        %            with value at (the new) time t equal to
        %                 [ x(t-L*Ts) ]
        %            where T is the forward horizon,
        %                  L is the backwards horizon,
        %                  Ts is the sample time.
        %   u_warm = Sequence of future control for a warm start at the next iteration
        %            Should be a matrix with size
        %                 [ # controls , T - controlDelay ]
        %            with value at (the new) time t equal to
        %                 [ u(t+delay*Ts), ..., u(t+(T-1)*Ts) ]
        %            where u(s) is applied from time s to time s+Ts
        %   d_warm = Sequence of disturbance for a warm start at the next iteration
        %            Should be a matrix with size
        %                 [ # disturbances , L + T ]
        %            with value at (the new) time t equal to
        %                 [ d(t-L*Ts, ..., d(t), ..., d(t+(T-1)*Ts) ]
        %            where d(s) is applied from time s to time s+Ts
        %   state = Sequence of true states
        %             which is a matrix with size
        %                  [ # states , L + 1 ]
        %             with value at time t equal to
        %                  [ x(t-L*Ts), ..., x(t), ..., x(t) ]
        %             where T is the forward horizon,
        %                   L is the backwards horizon,
        %                   Ts is the sample time.
        %   t_state = Sequence of times for state.
            if ~isequal(size(t_init),[1,1]);
                error('time value value must be a scalar');
            end
            if ~isequal(size(x_init),[obj.nStates,1])
                error('state value must be a column %dx1 array (not size=[%s])',...
                      obj.nStates,index2str(size(x_init)));
            end
            if ~isequal(size(u_init),[obj.nControls,obj.nBackwardHorizon+obj.controlDelay])
                error('control value must be a column %dx%d array (not size=[%s])',...
                      obj.nControls,obj.nBackwardHorizon+obj.controlDelay,...
                      index2str(size(u_init)));
            end
            if ~isequal(size(d_init),[obj.nDisturbances,obj.nBackwardHorizon])
                error('control value must be a column %dx%d array (not size=[%s])',...
                      obj.nControls,obj.nBackwardHorizon,...
                      index2str(size(d_init)));
            end
            if ~isequal(size(n_init),[obj.nOutputs,obj.nBackwardHorizon+1])
                error('control value must be a column %dx%d array (not size=[%s])',...
                      obj.nControls,obj.nBackwardHorizon+1,...
                      index2str(size(n_init)));
            end
            obj.history.time(1:obj.nBackwardHorizon+obj.controlDelay+1)=t_init+...
                obj.sampleTimeValue*(0:obj.nBackwardHorizon+obj.controlDelay)';
            obj.history.state(:,1)=x_init;
            obj.history.control(:,1:obj.nBackwardHorizon+obj.controlDelay)=u_init;
            obj.history.disturbance(:,1:obj.nBackwardHorizon)=d_init;
            for k=1:obj.nBackwardHorizon+1
                if k<=obj.nBackwardHorizon
                    try
                        % Forward Euler integration for the dynamics
                        obj.history.state(:,k+1)=obj.history.state(:,k)+obj.sampleTimeValue*...
                            obj.stateDerivativeFunction(obj.history.state(:,k),...
                                                        obj.history.control(:,k),...
                                                        obj.history.disturbance(:,k),...
                                                        obj.parameterValues{:});
                    catch me
                        fprintf('error in appyling ''stateDerivativeFunction'' to numerical values\n');
                        rethrow(me);
                    end
                end
                try
                    % output
                    obj.history.output(:,k)=obj.outputFunction(...
                        obj.history.state(:,k),...
                        obj.history.control(:,k),...
                        obj.history.disturbance(:,k),...
                        obj.parameterValues{:})+n_init(:,k);
                catch me
                    fprintf('error in appyling ''outputFunction'' to numerical values\n');
                    rethrow(me);
                end
            end
            obj.history.currentIndex=obj.nBackwardHorizon+1;

            t=obj.history.time(obj.history.currentIndex);
            
            k_u_past=obj.history.currentIndex-obj.nBackwardHorizon:...
                     obj.history.currentIndex+obj.controlDelay-1;
            u_past=obj.history.control(:,k_u_past);
            t_u_past=obj.history.time(k_u_past);

            k_y_past=obj.history.currentIndex-obj.nBackwardHorizon:...
                     obj.history.currentIndex;
            y_past=obj.history.output(:,k_y_past);
            t_y_past=obj.history.time(k_y_past);
            
            k_state=obj.history.currentIndex-obj.nBackwardHorizon:...
                    obj.history.currentIndex;
            t_state=obj.history.time(k_state);
            state=obj.history.state(:,k_state);
        end
        
        function [t_state,state]=setSolverWarmStart(obj,u_past,y_past,x_warm,u_warm,d_warm)
        % [t_state,state]=setSolverWarmStart(obj,u_past,y_past,x_warm,u_warm,d_warm)
        %
        % Inputs: 
        %    u_past = Sequence of past control controls.
        %             Should be a matrix with size
        %                  [ # controls , L + controlDelay ]
        %             with value at time t equal to
        %                  [ u(t-L*Ts), ..., u(t+(delay-1)*Ts) ]
        %             where u(s) is applied from time s to time s+Ts
        %    y_past = Sequence of past measurements.
        %             Should be a matrix with size
        %                  [ # outputs , L + 1 ]
        %             with value at time t equal to
        %                  [ y(t-L*Ts), ..., y(t)]
        %    x_warm = Initial state
        %             Should be a matrix with size
        %                  [ # states , 1 ]
        %             with value at time t equal to
        %                  [ x(t-L*Ts) ]
        %             where T is the forward horizon,
        %                   L is the backwards horizon,
        %                   Ts is the sample time.
        %    u_warm = Sequence of future control controls.
        %             Should be a matrix with size
        %                  [ # controls , T - controlDelay ]
        %             with value at time t equal to
        %                  [ u(t+delay*Ts), ..., u(t+(T-1)*Ts) ]
        %             where u(s) is applied from time s to time s+Ts
        %    d_warm = Sequence of disturbance controls.
        %             Should be a matrix with size
        %                  [ # disturbances , L + T ]
        %             with value at time t equal to
        %                  [ d(t-L*Ts, ..., x(t), ..., x(t+(T-1)*Ts) ]
        %             where d(s) is applied from time s to time s+Ts
        % Computes 
        %    state = Sequence of states compatible with "warm" inputs
        %             which is a matrix with size
        %                  [ # states , L + T + 1 ]
        %             with value at time t equal to
        %                  [ x(t-L*Ts), ..., x(t), ..., x(t+T*Ts) ]
        %             where T is the forward horizon,
        %                   L is the backwards horizon,
        %                   Ts is the sample time.
        %   t_state = Sequence of times for state.
        %
        % and initializes the solver with these sequences.
        %
        % The differential equation is integrated using forward Euler.

            if ~isequal(size(x_warm),[obj.nStates,1])
                error('size of x_warm [%d,%d] does not match expected size [%d,%d]',...
                      size(x_warm),obj.nStates,1);
            end
            if ~isequal(size(u_warm),[obj.nControls,obj.nForwardHorizon-obj.controlDelay])
                error('size of u_warm [%d,%d] does not match expected size [%d,%d]',...
                      size(u_warm),obj.nControls,obj.nHorizon);
            end
            if ~isequal(size(d_warm),[obj.nDisturbances,obj.nForwardHorizon+obj.nBackwardHorizon])
                error('size of d_warm [%d,%d] does not match expected size [%d,%d]',...
                      size(d_warm),obj.nDisturbances,obj.nForwardHorizon+obj.nBackwardHorizon);
            end
            if obj.history.currentIndex<obj.nBackwardHorizon+1
                error('must initialize state');
            end
            state=nan(obj.nStates,obj.nForwardHorizon+obj.nBackwardHorizon+1);
            state(:,1)=x_warm;
            u=[u_past,u_warm];
            for k=1:obj.nForwardHorizon+obj.nBackwardHorizon
                try
                    % Forward Euler integration for the dynamics
                    state(:,k+1)=state(:,k)+...
                        obj.sampleTimeValue*obj.stateDerivativeFunction(...
                            state(:,k),u(:,k),d_warm(:,k),...
                            obj.parameterValues{:});
                catch me
                    fprintf('error in appyling ''stateDerivativeFunction'' to numerical values\n');
                    rethrow(me);
                end
            end
            
            k_state=obj.history.currentIndex-obj.nBackwardHorizon:...
                    obj.history.currentIndex+obj.nForwardHorizon;
            t_state=obj.history.time(k_state);
            
            setParameter(obj,obj.pastControlName,u_past);
            setParameter(obj,obj.pastOutputName,y_past);
            feval(['setV_',obj.initialStateName],obj.solverObject,state(:,1));
            obj.initialStateSet=true;
            feval(['setV_',obj.latentStateName],obj.solverObject,state(:,2:end));
            obj.latentStateSet=true;
            feval(['setV_',obj.futureControlName],obj.solverObject,u_warm);
            obj.futureControlSet=true;
            feval(['setV_',obj.disturbanceName],obj.solverObject,d_warm);
            obj.disturbanceSet=true;
        end

        function [solution,varargout]=solve(obj,mu0,maxIter,saveIter)
        % [solution,y1,y2,...]=solve(obj,mu0,maxIter,saveIter)
        %
        % Inputs:
        %   mu0      - initial value for the gap variable
        %   maxIter  - maximum number of iterations
        %   saveIter - iteration # when to save the "hessian" matrix
        %              (see class2optimizeCS/cmex2optimizeCS)
        % Calls the solver and returns:
        %   solution.state = Optimal state,
        %            which is a matrix with size
        %                  [ # states , L + T + 1 ]
        %            with value at time t equal to
        %                  [ x(t-L*Ts), ..., x(t), ..., x(t+T*Ts) ]
        %            where T is the forward horizon,
        %                  L is the backwards horizon,
        %                  Ts is the sample time.
        %   solution.futureControl = Sequence of optimal future control,
        %            which is a matrix with size
        %                 [ # controls , T - controlDelay ]
        %            with value at time t equal to
        %                 [ u(t+delay*Ts), ..., u(t+(T-1)*Ts) ]
        %            where u(s) is applied from time s to time s+Ts
        %   solution.disturbance = Sequence of optimal disturbance,
        %            which is a matrix with size
        %                 [ # disturbances , L + T ]
        %            with value at time t equal to
        %                 [ d(t-L*Ts, ..., x(t), ..., x(t+(T-1)*Ts) ]
        %            where d(s) is applied from time s to time s+Ts
        %   solution.status - solver exit status
        %   solution.iter   - number of solver iterations
        %   solution.time   - time taken by solver
        %   y1,y2,...       - output expressions (passed to the constructor)

            if ~all(obj.parametersSet)
                obj.parameterNames,
                obj.parametersSet,
                error('must set all parameter values before calling solver');
            end
            if ~obj.futureControlSet
                error('must set initial value for the future control before calling solver');
            end
            if ~obj.disturbanceSet
                error('must set initial value for the disturbance before calling solver');
            end
            if ~obj.initialStateSet
                error('must set initial value for the initial state before calling solver');
            end
            if ~obj.latentStateSet
                error('must set initial value for the latent state before calling solver');
            end
            
            [solution.status,solution.iter,solution.time]=...
                solve(obj.solverObject,mu0,int32(maxIter),int32(saveIter));
            
            varargout=cell(obj.nOutputExpressions,1);
            [varargout{:},solution.futureControl,...
             solution.disturbance,solution.initialState,...
             solution.state,solution.objective]=getOutputs(obj.solverObject);
            solution.t=obj.history.time(obj.history.currentIndex);
            solution.t_futureControl=solution.t+...
                obj.sampleTimeValue*(obj.controlDelay:obj.nForwardHorizon-1);
            solution.t_disturbance=solution.t+...
                obj.sampleTimeValue*(-obj.nBackwardHorizon:obj.nForwardHorizon-1);
            solution.t_state=solution.t+...
                obj.sampleTimeValue*(-obj.nBackwardHorizon:obj.nForwardHorizon);
        end
        
        function [t,t_y_past,y_past,t_u_past,u_past,x_warm,u_warm,d_warm,t_state,state]=...
                applyControl(obj,u_past,y_past,solution,disturbance,noise,u_final,d_final);
        % [t,t_y_past,y_past,t_u_past,u_past,x_warm,u_warm,d_warm,t_state,state]=...
        %                applyControls(obj,u_past,y_past,solution,disturbance,noise,u_final,d_final)
        %
        % Applied the first step of an MPCMHE control.
        %
        % Inputs:
        %   u_past = Sequence of past control controls.
        %            Should be a matrix with size
        %                 [ # controls , L + controlDelay ]
        %            with value at time t equal to
        %                 [ u(t-L*Ts), ..., u(t+(delay-1)*Ts) ]
        %            where u(s) is applied from time s to time s+Ts
        %   y_past = Sequence of past measurements.
        %            Should be a matrix with size
        %                 [ # outputs , L + 1 ]
        %            with value at time t equal to
        %                 [ y(t-L*Ts), ..., y(t)]
        %   solution = MPCMHE solution returned by solver()
        %   disturbance = Vector of input disturbance at time 
        %                 [ t ]
        %                 if empty, the worst-case disturbance from
        %                 the last solution is used
        %   noise    = Vector of measurement noise that will
        %              appear in the output at time 
        %                 [ t+Ts ]
        %   u_final  = Sequence of control values to be inserted at the end of u_warm
        %              Should be a matrix with size
        %                 [ # controls , 1 ]
        %   d_final  = Sequence of disturbance values to be inserted at the end of d_warm
        %              Should be a matrix with size
        %                 [ # disturbances , 1 ]
        % Computes the resulting states and outputs and returns
        %   t = new current time (essentially the previous time + Ts)
        %   u_past = Sequence of past control controls.
        %            Should be a matrix with size
        %                 [ # controls , L + controlDelay ]
        %            with value at (the new) time t equal to
        %                 [ u(t-L*Ts), ..., u(t+(delay-1)*Ts) ]
        %            where u(s) is applied from time s to time s+Ts
        %   t_u_past = Sequence of times for u_past.
        %   y_past = Sequence of (noisy) past measurements.
        %            Should be a matrix with size
        %                 [ # outputs , L + 1 ]
        %            with value at (the new) time t equal to
        %                 [ y(t-L*Ts), ..., y(t)]
        %   t_y_past = Sequence of times for y_past.
        %   x_warm = Initial state for a warm start at the next iteration
        %            Should be a matrix with size
        %                 [ # states , 1 ]
        %            with value at (the new) time t equal to
        %                 [ x(t-L*Ts) ]
        %            where T is the forward horizon,
        %                  L is the backwards horizon,
        %                  Ts is the sample time.
        %   u_warm = Sequence of future control for a warm start at the next iteration
        %            Should be a matrix with size
        %                 [ # controls , T - controlDelay ]
        %            with value at (the new) time t equal to
        %                 [ u(t+delay*Ts), ..., u(t+(T-1)*Ts) ]
        %            where u(s) is applied from time s to time s+Ts
        %   d_warm = Sequence of disturbance for a warm start at the next iteration
        %            Should be a matrix with size
        %                 [ # disturbances , L + T ]
        %            with value at (the new) time t equal to
        %                 [ d(t-L*Ts, ..., d(t), ..., d(t+(T-1)*Ts) ]
        %            where d(s) is applied from time s to time s+Ts
        %   state = Sequence of true states
        %             which is a matrix with size
        %                  [ # states , L + 1 ]
        %             with value at time t equal to
        %                  [ x(t-L*Ts), ..., x(t), ..., x(t) ]
        %             where T is the forward horizon,
        %                   L is the backwards horizon,
        %                   Ts is the sample time.
        %   t_state = Sequence of times for state.

            t=obj.history.time(obj.history.currentIndex);
            if size(obj.history.state,2)<obj.history.currentIndex+1
                extendHistory(obj)
            end

            t=t+obj.sampleTimeValue;
            obj.history.time(obj.history.currentIndex+1)=t;
            
            % control & disturbance
            obj.history.control(:,obj.history.currentIndex+obj.controlDelay)=...
                solution.futureControl(:,1);
            if isempty(disturbance) 
                obj.history.disturbance(:,obj.history.currentIndex)=...
                    solution.disturbance(:,obj.nBackwardHorizon+1);
            else
                obj.history.disturbance(:,obj.history.currentIndex)=disturbance;
            end
            try
                if 1
                    % Integration for the dynamics using ode23
                    [tout,yout]=ode23(...
                        @(t,x)obj.stateDerivativeFunction(x,...
                                                          obj.history.control(:,obj.history.currentIndex),...
                                                          obj.history.disturbance(:,obj.history.currentIndex),...
                                                          obj.parameterValues{:}),...
                        [obj.history.time(obj.history.currentIndex),obj.history.time(obj.history.currentIndex+1)],...
                        obj.history.state(:,obj.history.currentIndex));
                    obj.history.state(:,obj.history.currentIndex+1)=yout(end,:)';                    
                else
                    % Forward Euler integration for the dynamics
                    obj.history.state(:,obj.history.currentIndex+1)=...
                        obj.history.state(:,obj.history.currentIndex)+...
                        obj.sampleTimeValue*obj.stateDerivativeFunction(...
                            obj.history.state(:,obj.history.currentIndex),...
                            obj.history.control(:,obj.history.currentIndex),...
                            obj.history.disturbance(:,obj.history.currentIndex),...
                            obj.parameterValues{:});
                end
            catch me
                fprintf('error in appyling ''stateDerivativeFunction'' to numerical values\n');
                rethrow(me);
            end
            % save status etc.
            obj.history.objective(obj.history.currentIndex+obj.controlDelay)=...
                solution.objective;
            obj.history.status(obj.history.currentIndex+obj.controlDelay)=...
                solution.status;
            obj.history.iter(obj.history.currentIndex+obj.controlDelay)=...
                solution.iter;
            obj.history.stime(obj.history.currentIndex+obj.controlDelay)=...
                solution.time;
            
            %% Attention: if output depends on the current u & d and
            %%            controlDelay=0 this will results in NaN's
            try
                % output
                obj.history.output(:,obj.history.currentIndex+1)=obj.outputFunction(...
                    obj.history.state(:,obj.history.currentIndex+1),...
                    obj.history.control(:,obj.history.currentIndex+1),...
                    obj.history.disturbance(:,obj.history.currentIndex+1),...
                    obj.parameterValues{:})+noise;
            catch me
                fprintf('error in appyling ''outputFunction'' to numerical values\n');
                rethrow(me);
            end
            
            obj.history.currentIndex=obj.history.currentIndex+1;
            
            u_warm=[solution.futureControl(:,2:end),u_final];
            d_warm=[solution.disturbance(:,2:end),d_final];
            x_warm=solution.state(:,2);

            t=obj.history.time(obj.history.currentIndex);
            
            k_u_past=obj.history.currentIndex-obj.nBackwardHorizon:...
                     obj.history.currentIndex+obj.controlDelay-1;
            u_past=obj.history.control(:,k_u_past);
            t_u_past=obj.history.time(k_u_past);

            k_y_past=obj.history.currentIndex-obj.nBackwardHorizon:...
                     obj.history.currentIndex;
            y_past=obj.history.output(:,k_y_past);
            t_y_past=obj.history.time(k_y_past);
            
            k_state=obj.history.currentIndex-obj.nBackwardHorizon:...
                    obj.history.currentIndex;
            t_state=obj.history.time(k_state);
            state=obj.history.state(:,k_state);

            obj.initialStateSet=false;
            obj.latentStateSet=false;
            obj.futureControlSet=false;
            obj.disturbanceSet=false;            
        end

        function history=getHistory(obj);
        % [t,x,u,d,y,status,iter,stime]=getHistory(obj);
        %        
        % Returns the history of
        %    history.t      - sampling times
        %    history.u      - control values
        %    history.x      - state values
        %    history.objective - MPC/MHE cost
        %    history.status - solver status
        %    history.iter   - # of solver iterations
        %    history.stime  - solver times
            
            history.t=obj.history.time(1:obj.history.currentIndex-1);
            history.u=obj.history.control(:,1:obj.history.currentIndex-1);
            history.d=obj.history.disturbance(:,1:obj.history.currentIndex-1);
            history.x=obj.history.state(:,1:obj.history.currentIndex-1);
            history.y=obj.history.output(:,1:obj.history.currentIndex-1);
            history.objective=obj.history.objective(1:obj.history.currentIndex-1);
            history.status=obj.history.status(1:obj.history.currentIndex-1);
            history.iter=obj.history.iter(1:obj.history.currentIndex-1);
            history.stime=obj.history.stime(1:obj.history.currentIndex-1);
        end
    end
end
