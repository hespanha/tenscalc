classdef nlss < matlab.mixin.Copyable % abstract handle class with a copy method

% Class used to store a non-linear state space dynamical system and
% perform numerical and symbolic (Tcalculus) simulations
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    properties
        % system dynamics:
        % . for a continous time system \dot x = f(x,u,t)
        % . for a discrete time system  x(k+1) = f(x(k),u(k),k)
        f=NaN;

        % output equation: y=g(x,u,t)
        g=NaN;

        % discrete or continuous time?
        discrete=true;

        % Initial condition/last value of the state
        x0=[]; % column vector
        t0=[];

        % State variable name
        stateName='x';

        % sizes
        nStates=NaN;
        nInputs=NaN;
        nOutputs=NaN;

    end % properties

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Object creation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=nlss(f,g,discrete,stateName,varargin)
        % obj=nlss(f,g,discrete,stateName,x0,t0)
        %
        % Creates a non-linear state-space dynamical system with
        %  . dynamics defined by the functions f and g
        %    (output g is optional)
        %  . discrete-time when discrete=true, continous-time otherwise
        %  . initial condition x(t0)=x0 (column vector)
        %  . "stateName" is a string with the name of the state
        %    (must be a valid matlab variable name, which will be used
        %     to represent the state as a symbolic variable in Tcalculus)
        %    When omitted "x" will be used.
        %
        % For a continous-time system, the dynamics and output map are
        %    \dot x(t) = f(x(t),u(t),t)  in  R^n
        %         y(t) = g(x(t),u(t),t)  in  R^k
        % and for a discrete-time system
        %    x(k+1) = f(x(k),u(k),k)     in  R^n
        %    y(k)   = g(x(k),u(k),k)     in  R^m
        %
        % The input u, state x, and output y are should be column vectors,
        % but the function f & g must be vectorzed in the sense that
        % they should be able to take as inputs matrices
        %     x [n x N], u [k x N], t [N x 1]
        % and
        %     f(x,u,t) must return an [n x N] matrix with the ith column
        %              equal to f(x(:,i),u(:,i),t(i))
        %     g(x,u,t) must return an [m x N] matrix with the ith column
        %              equal to g(x(:,i),u(:,i),t(i))
        %
        % For example, for a time-invariant linear system one would set
        %     f=@(x,u,t)A*x+B*u;
        %     g=@(x,u,t)C*x+D*u;

            obj.f=f;
            if nargin>=2
                obj.g=g;
            end
            if nargin>=3
                obj.discrete=discrete;
            end
            if nargin>=4
                obj.stateName=stateName;
            end
            if nargin>=5 && ~isempty(varargin{1})
                setInitialState(obj,varargin{:});
            end

            if false % obj.discrete==false
                error('continuout-time has bugs')
            end

        end % function nlss()

        function sz=size(obj,dim)
        % sz=size(obj,dim)
        %
        % returns a row vectors containing: [# outputs, # inputs, # states]
        %
        % sz=size(obj,dim)
        %
        % returns entry 'dim' of the row evctor above.

            sz=[obj.nOutputs,obj.nInputs,obj.nStates];
            if nargin>1
                sz=sz(dim);
            end
        end

        function setInitialState(obj,x0,t0)
        % setInitialState(obj,x0,t0)
        % Sets the initial state of a non-linear state-space dynamical
        % system to x0 (column vector), at time t0.
        % When t0 is omitted t0=0 is assumed.
            obj.x0=x0;
            if length(size(x0))~=2 || size(x0,2)~=1
                error('setInitialState: initial state for a nlss system must be an [nx1] matrix ([%s] matrix instead)\n',index2str(size(x0)));
            end
            if nargin>=3
                obj.t0=t0;
            else
                obj.t0=0;
            end
            if isnan(obj.nStates)
                obj.nStates=size(x0,1);
            else
                if obj.nStates ~= size(x0,1)
                    obj.nStates
                    x0
                    error('setInitialState: incorrect size for initial state (%d different than expected %d)',...
                          size(x0,1),obj.nStates);
                end
            end
        end % function setInitialState()

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Object display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function s=char(obj)
        % char(obj)
        % Produces a string describing a non-linear state-space
        % dynamical system.
            if obj.discrete
                s=sprintf('  x^+=%s;\n',char(obj.f));
            else
                s=sprintf('  dot x=%s;\n',char(obj.f));
            end
            if strcmp(class(obj.g),'function_handle')
                s=sprintf('  y=%s;\n',char(obj.g));
            end
            if isnumeric(obj.t0)
                s=sprintf('%s  t0=%g;\n',s,obj.t0);
            else
                s=sprintf('%s  t0=Tcalculus(%s);\n',s,str(obj.t0));
            end
            if isnumeric(obj.x0)
                s=sprintf('%s  x0=%s;\n',s,mat2str(obj.x0));
            else
                s=sprintf('%s  x0=Tcalculus(%s);\n',s,str(obj.x0));
            end
        end % function char()

        function disp(obj)
        % disp(obj)
        % Displays a non-linear state-space dynamical system.
            disp(char(obj))
        end % function disp()

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Simulation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [y,x,constraint]=sim(obj,u,t,x0)
        % [y,x,constraint]=sim(obj,u,t,x0)
        % Simulates the time response of a non-linear state-space
        % dynamical system for
        %   t  [N] or [N x 1] = times at which the input was sampled
        %   u  [k x N]   = input vector samples at times t(1), t(2), ..., t(end)
        %   x0 [n x 1]   = initial value of the state at time t(1) (optional)
        % provides
        %   y  [m x N]   = values of the output at times t(1), t(2), ..., t(end)
        %   x  [n x N-1] = for continuous-time:
        %                    values of the state at times t(2), ..., t(end)
        %   x  [n x N]   = for discrete-time
        %                    values of the state at times t(2), ..., t(end),t(end)+1
        % constraint  =   TC object expression a constraint on the state vector
        %                 that enforces the dynamics.
        % upon return the systems initial condition is set to
        %   t0 = t(end), x0 = x(:,end) for continous-time
        % and
        %   t0 = t(end)+1, x0 = f(x(:,end),u(:,end),t(end)) for discrete-time

            if isnumeric(t) && size(t,2)~=1
                error('sim: time must be a [Nx1] numeric column vector ([%s] matrix instead)\n',...
                      index2str(size(t)));
            end
            if ~isnumeric(t) && length(size(t))~=1
                error('sim: time must be a [N] symbolic vector ([%s] matrix instead)\n',...
                      index2str(size(t)));
            end
            nT=size(t,1);
            if size(u,2)~=nT
                error('sim: inconsistent sizes of time [%s] and input [%s]\n',...
                      index2str(size(t)),index2str(size(u)));
            end

            %% Update/check size of input
            if isnan(obj.nInputs)
                obj.nInputs=size(u,1);
            else
                if obj.nInputs ~= size(u,1)
                    obj.nInputs
                    u
                    error('setInitialState: incorrect size for input (%d different than expected %d)',...
                          size(u,1),obj.nInputs);
                end
            end


            %% set initial condition
            if nargin>=4
                setInitialState(obj,x0,t(1));
            end

            y=nan;

            if obj.discrete
                %% Discrete-time Simulation
                if isnumeric(obj.x0) && isnumeric(u)
                    %% numeric simulation
                    xall=[obj.x0,nan(size(obj.x0,1),nT)];
                    for i=2:nT+1
                        xall(:,i)=obj.f(xall(:,i-1),u(:,i-1),t(i-1));
                    end
                    x=xall(:,2:end);
                    if strcmp(class(obj.g),'function_handle')
                        y=obj.g(xall(:,1:end-1),u,t);
                    end
                    % update state
                    obj.x0=xall(:,end);
                    obj.t0=t(end)+1;
                else
                    %% symbolic
                    x=Tvariable(obj.stateName,[size(obj.x0,1),nT]);
                    xpre=[obj.x0,x(:,1:end-1)];
                    constraint=(x==obj.f(xpre,u,t));
                    if strcmp(class(obj.g),'function_handle')
                        y=obj.g(xpre,u,t);
                    end
                end
            else
                %% Continuous-time Simulation
                if isnumeric(obj.x0) && isnumeric(u)
                    %%fprintf('Numeric simulation\n')
                    %% Numerical simulation
                    if 1
                        % ODE23
                        fhold=@(t,yi,ti)yi(:,max(find(ti<=t+sqrt(eps))));
                        % clf,plot(t,u,'b+');hold on
                        % tt=min(t):.1:2*max(t);
                        % for i=1:length(tt)
                        %     plot(tt(i),fhold(tt(i),u,t),'r.')
                        % end
                        fode=@(tt,xx)obj.f(xx,fhold(tt,u,t),t);
                        [tout,x]=ode23(fode,t,obj.x0);
                        if length(t)==2
                            % if T has only 2 elements, ode23
                            % assumes it is just start and final
                            % times and provides intermediate points
                            x=x([1,end],:);
                        end
                        if strcmp(class(obj.g),'function_handle')
                            y=obj.g(x',u,t);
                        end
                        x=x(2:end,:)';
                        % update state
                        obj.x0=x(:,end);
                        obj.t0=t(end);
                    else
                        % Euler
                        xall=[obj.x0,nan(size(obj.x0,1),nT-1)];
                        for i=2:nT
                            xall(:,i)=xall(:,i-1)+(t(i)-t(i-1))*obj.f(xall(:,i-1),u(:,i-1),t(i-1));
                        end
                        % update state
                        obj.x0=xall(:,end);
                        obj.t0=t(end);
                        if strcmp(class(obj.g),'function_handle')
                            y=obj.g(xall,u,t);
                        end
                        x=xall(:,2:end);
                    end
                else
                    %% fprintf('Symbolic simulation\n')
                    %% Symbolic "simulation" - Euler
                    x=Tvariable(obj.stateName,[size(obj.x0,1),nT-1]);
                    xpre=[obj.x0,x(:,1:end-1)];
                    tpre=t(1:end-1);
                    upre=u(:,1:end-1);
                    dx=obj.f(xpre,upre,tpre);
                    constraint=(x==xpre+tprod(t(2:end)-tpre,[2],dx,[1,2]));
                    xall=[obj.x0,x];
                    if strcmp(class(obj.g),'function_handle')
                        y=obj.g(xall,u,t);
                    end
                end
            end

            %% Update/check size of input
            if isnan(obj.nOutputs)
                obj.nOutputs=size(y,1);
            else
                if obj.nOutputs ~= size(y,1)
                    obj.nOutputs
                    u
                    error('setInitialState: incorrect size for initial state (%d different than expected %d)',...
                          size(x0,1),obj.nOutputs);
                end
            end

        end % function sim()

    end % methods

end % classdef

function test()

    clear all
    A=[0,1;0,0];
    B=[0;1];
    C=[1,0];
    D=0;
    f=@(x,u,t)A*x+B*u;
    g=@(x,u,t)C*x+D*u;
    sys=nlss(f,g,0,[],[],'x')

    u=(-2).^(1:5);
    t=(1:4)';

    x0=[.45;.54];
    t0=0;
    [y,t,x]=sim(sys,u,t,x0,t0)

    Tvariable x0 [2];
    t0=0;
    [y,t,x]=sim(sys,u,t,x0,t0)


end