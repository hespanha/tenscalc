% MPC control of a brushed DC motor with the following 2nd order state-space model
%
%   [dot x1]=[0  1][x1]+[0]u
%   [dot x2] [0  p][x2] [k]
%
% where
%   theta = x1 = shaft angle
%   omega = x2 = shaft angular velocity
%   p          = pole (in addition to the pole at the origin)
%   k          = gain
%
% The goal of MPC is to minimize a criterion of the form
%
%         J = \int_0^T ( theta(t)-ref )^2 dt + lambda_u \int_0^T u(t)^2 dt
%
% where
%   ref = desired angle (reference)
%
% subject to the following state and input constraints
%
%        min_theta <= theta(t) <= max_theta, \forall t
%        min_omega <= omega(t) <= max_omega, \forall t
%        min_u     <= u(t)     <= max_u, \forall t
%
% The solver can be generated once and then called multiple
% times. Creating a solver typically takes much longer than solving
% the optimization by creating the solver.
%
% This script gives you the option of skipping creating the solver,
% but you can only skip generating a solver if you have done it once.

% This file is part of Tencalc.
%
% Copyright (C) 2012-22 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all
% remove previous solvers
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

help mpc_dcmotor

createSolvers=input('Do you want to create the solvers (Y/n)? ','s');
createSolvers=~isequal(createSolvers,'n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=30;                        % horizon length in time-steps
nX=2;                        % state size
nU=1;                        % control size

if createSolvers

    %%%%%%%%%%%%%%%%%%%
    %% Generate solver
    %%%%%%%%%%%%%%%%%%%

    Tcalculus.clear();

    Tvariable Ts [];              % sampling time
    Tvariable x [nX,T];           % state:              [x(1) x(2) ... x(T)]
    Tvariable u [nU,T-1];         % control:            [u(1) u(2) ... u(T-1)] (ZOH)

    Tvariable xinit [nX,1];       % initial state:      x(1)

    Tvariable p [];               % pole
    Tvariable k [];               % gain

    A=[0,1;0,p];                  % matrices for the dynamics
    B=[0;k];

    % System dynamics
    dynamics={
        % option 1) forward Euler
        x(:,2:end)==x(:,1:end-1)+Ts*(A*x(:,1:end-1)+B*u) 
        % option 2) trapesoidal, with ZOH for u       
        %x(:,2:end)==x(:,1:end-1)+Ts*(A*(x(:,1:end-1)+x(:,2:end))/2+B*u) 
        % initial state
        x(:,1)==xinit; 
             };

    % Constraints
    Tvariable min_x [nX,1];       % minimum state
    Tvariable max_x [nX,1];       % maximum state
    Tvariable min_u [nU,1];       % minimum inout
    Tvariable max_u [nU,1];       % maximum input

    constraints={
        x(:,2:end)>=repmat(min_x,[1,T-1]);  % constraints only on future states
        x(:,2:end)<=repmat(max_x,[1,T-1]);  % constraints only on future states
        u>=repmat(min_u,[1,T-1]);           % constraints on all controls
        u<=repmat(max_u,[1,T-1]);           % constraints on all controls
                };

    % Criterion
    Tvariable ref [1,T];           % reference for angle
    Tvariable lambda_u [];         % weigth for input

    Jx2=tsIntegral(sum((x(1,:)-ref).^2,1),Ts);     % integral of position error^2
    Ju2=tsIntegral(sum(u.^2,1),Ts);             % integral of control^2

    J = Jx2+lambda_u*Ju2;

    % Warm start for next optimization (shift and  move away from constrainsts)

    xWarm=[x(:,2:end),zeros(nX,1)];
    uWarm=[u(:,2:end),zeros(nU,1)];
    xWarm=max(xWarm,repmat(min_x+.05*(max_x-min_x),[1,T])); % move away from min_x
    xWarm=min(xWarm,repmat(max_x-.05*(max_x-min_x),[1,T])); % move away from max_x
    uWarm=max(uWarm,repmat(min_u+.05*(max_u-min_u),[1,T-1])); % move away from min_u
    uWarm=min(uWarm,repmat(max_u-.05*(max_u-min_u),[1,T-1])); % move away from max_u

    % Solver output
    outputExpressions=struct('J',J,...
                             'Jx2',Jx2,...
                             'Ju2',Ju2,...
                             't',(0:T-1)*Ts,...
                             'u',u,...
                             'x',x,...
                             'ref',ref,...
                             'theta',x(1,:),...
                             'omega',x(2,:),...
                             'xWarm',xWarm,...
                             'uWarm',uWarm);

    classname=cmex2optimizeCS('classname','tmp_dcmotor',...
                              'objective',J,...
                              'optimizationVariables',{u,x},...
                              'constraints',[dynamics;constraints],...
                              'outputExpressions',outputExpressions,...
                              'parameters',{...
                                  Ts;p;k;...
                                  xinit;...
                                  ref;...
                                  min_x;max_x;min_u;max_u;...
                                  lambda_u},...
                              'addEye2Hessian',true,...
                              'adjustAddEye2Hessian',true,...
                              'scaleInequalities',true,...
                              'solverVerboseLevel',2);

end

%%%%%%%%%%%%%%%
%% Use solver
%%%%%%%%%%%%%%%

%% Create object
obj=tmp_dcmotor();
mu0=1e-3;       % low value will speed-up convergence (up to a point where it causes numerical issues)
maxIter=100;
saveIter=-1;

s = RandStream('mt19937ar','Seed',0); % fix seed for reproducibility
RandStream.setGlobalStream(s);

%% Set parameters

ref=@(t)-.35*sign(sin(.5*t)); % reference signals

Ts=.1;
p=2;k=1;
A=[0,1;0,p];                  % matrices for the dynamics
B=[0;k];
min_x=[-.4;-.3];
max_x=[.4;.3];
min_u=-1;
max_u=1;
lambda_u=1/50;

setP_Ts(obj,Ts);
setP_p(obj,p);
setP_k(obj,k);
setP_min_x(obj,min_x);
setP_max_x(obj,max_x);
setP_min_u(obj,min_u);
setP_max_u(obj,max_u);
setP_lambda_u(obj,lambda_u);

xinit=[.2;.2];
xWarm=xinit+.01*rand(nX,T); % random initialization close to initial condition
uWarm=.01*rand(nU,T-1);
t=0;
clear history;   % structure where system's state, control, etc. will be stored
for k=1:150

    %% Set parameters specific to this optimization
    setP_ref(obj,ref(t+(0:T-1)*Ts));
    setP_xinit(obj,xinit);

    %% Solver initialization (cold for k==1 and then warm)
    setV_u(obj,uWarm);
    setV_x(obj,xWarm);

    %% Solve optimization
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    out=getOutputs(obj);
    k
    disp(out)

    %% Save applied sate, control, etc. in `history` structure
    history.t(k)=t;
    history.x(1:nX,k)=xinit;
    history.theta(1:nX,k)=xinit(1);
    history.omega(1:nX,k)=xinit(2);
    history.ref(k)=out.ref(1);
    history.u(1:nU,k)=out.u(:,1);
    history.J(k)=out.J;
    history.Jx2(k)=out.Jx2;
    history.Ju2(k)=out.Ju2;
    history.iter(k)=iter;
    history.time(k)=time;

    %% Apply (constant) control to update initial position and velocity
    [tout,yout]=ode23(@(t,x)A*x+B*out.u(:,1),[0,Ts],xinit);
    xinit=yout(end,:)';
    t=t+Ts;
    xWarm=out.xWarm;
    uWarm=out.uWarm;
    
    if mod(k,10)==1 || status
        %% Plot MPC solution
        fig=clearFigure('figureNumber',10,'figureName','Latest MPC solution');
        plotData(out,Ts)
        %% Plot actual control
        fig=clearFigure('figureNumber',fig+1,'figureName','Actual control');
        plotData(history,Ts)
        fig=clearFigure('figureNumber',fig+1,'figureName','Solver');
        subplot(2,1,1);
        plot(history.t,history.J,'.-',...
             history.t,history.Jx2,'.-',...
             history.t,history.Ju2,'.-')
        legend('J','Jx2','Ju2','location','best');grid on;
        subplot(2,1,2);
        yyaxis left;plot(history.t,history.iter,'.-');ylabel('# solver iter');grid on;
        yyaxis right;plot(history.t,1000*history.time,'.-');ylabel('solver time [ms]');grid on;
        drawnow;
    end

    if status
        error('optimization failed\n');
    end

end

% Delete object
clear obj

function plotData(data,Ts)
subplot(3,1,1);
plot(data.t,data.theta,'.-',...
     data.t,data.ref,'-'); grid on
ylabel('\theta');
subplot(3,1,2);
plot(data.t,data.omega,'.-'); grid on
ylabel('\omega');
subplot(3,1,3);
plot(data.t(1:size(data.u,2)),data.u,'.-'); grid on
ylabel('u');

end
