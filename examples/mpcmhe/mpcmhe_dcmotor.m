% MPC-MHE output-feedback control of a brushed DC motor with the
% following simplified 2nd order continous-time state-space model:
%
%   [dot x1] = [0  1] [x1]+[0] (u+d)
%   [dot x2]   [0  p] [x2] [k]
%   y        = [1  0] [x1] + n
%                     [x2]
% where
%   theta = x1 = shaft angle
%   omega = x2 = shaft angular velocity
%   u          = input voltage
%   d          = (unmeasured) additive input disturbance
%   y     = x1 = measured output (shaft angle)
%   n          = measurement noise
%   p          = pole (in addition to the pole at the origin)
%   k          = gain
%
% The goal of MPC-MHE is to solve the following saddle-point optimization
%
%         min_uFuture max_{x(-L),d,n} \int_0^T ( theta(t)-ref )^2 dt + lambda_u \int_0^T uFuture(t)^2 dt
%                               - lambda_d int_{-L}^T d^2 dt
%                               - lambda_n int_{-L}^0 n^2 dt
%
% where
%   ref     = [r(1),...,r(T)]     = desired angle (reference)
%   uPast   = [u(-L),....,u(-1)]  = past control inputs (selected previously and used with ZHO)
%   uFuture = [u(0),....,u(T-1)]  = future control inputs (to be selected and to be used with ZOH)
%
%   x(-L)                         = unknown initial condition (at the start of backward horizon)
%   d       = [d(-L),....,d(T-1)] = unknown disturance (assumed with ZOH)
%   y       = [y(-L) ... y(-1)]   = past measured output
%   n       = [n(-L) ... n(-1)]   = past measurement noise
%
%
% subject to the following constraints
%
%        |u(t)| <= max_u,   \forall t    [minimizer constraint]
%        |d(t)| <= max_d,   \forall t    [maximizer constraint]
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

clear all;
% remove previous solvers
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

help mpcmhe_dcmotor

createSolvers=input('Do you want to create the solvers (Y/n)? ','s');
createSolvers=~isequal(createSolvers,'n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nX=2;                        % state size
nU=1;                        % control size
nD=1;                        % disturbance size
nY=1;                        % output size

T=60;                        % forward horizon length in time-steps
L=40;                        % backward horizon length in time-steps

if createSolvers

    %%%%%%%%%%%%%%%%%%%
    %% Generate solver
    %%%%%%%%%%%%%%%%%%%

    Tcalculus.clear();

    Tvariable Ts [];              % sampling time
    Tvariable x0      [nX,1];     % initial state:             x(-L)
    Tvariable x1      [nX,L+T];   % (remaining) states:        [x(-L+1) ... x(0) ... x(T)]
    Tvariable uPast   [nU,L];     % past control (ZOH):        [u(-L) ... u(-1)]
    Tvariable uFuture [nU,T];     % future control ZOH):       [u(0)  ... u(T-1)]
    Tvariable d       [nD,L+T];   % past & future disturbance: [d(-L) ... d(0) ... d(T-1)]
    Tvariable yPast   [nY,L];     % past measurement:          [y(-L) ... y(-1)]


    Tvariable p [];               % pole
    Tvariable k [];               % gain

    A=[0,1;0,p];                  % matrices for the continuous-time dynamics
    B=[0;k];
    C=[1,0];

    % System dynamics
    x=[x0,x1];
    u=[uPast,uFuture];
    dynamics={
        % option 1) forward Euler
        x(:,2:end)==x(:,1:end-1)+Ts*(A*x(:,1:end-1)+B*(u+d))
        % option 2) trapesoidal, with ZOH for u
        %x(:,2:end)==x(:,1:end-1)+Ts*(A*(x(:,1:end-1)+x(:,2:end))/2+B*u)
             };

    % Constraints
    Tvariable max_u [nU,1];       % maximum control input
    Tvariable max_d [nD,1];       % maximum disturbance input

    P1constraints={
        uFuture>=repmat(-max_u,[1,T]);
        uFuture<=repmat(max_u,[1,T]);
                };
    P2constraints={
        d>=repmat(-max_d,[1,L+T]);
        d<=repmat(max_d,[1,L+T]);
                };

    % Criterion
    Tvariable ref [1,T];           % reference for angle
    Tvariable lambda_u [];         % weigth for control input
    Tvariable lambda_d [];         % weigth for input disturbance
    Tvariable lambda_n [];         % weigth for measurement noise


    Jx2=tsIntegral(sum((x(1,end-T+1:end)-ref).^2,1),Ts); % integral of position error^2
    Ju2=tsIntegral(sum(uFuture.^2,1),Ts);                % integral of control^2
    Jd2=tsIntegral(sum(d.^2,1),Ts);                      % integral of disturbance^2
    n=C*x(:,1:L)-yPast;                                  % measurement noise
    Jn2=tsIntegral(sum(n.^2,1),Ts);                      % integral of noise^2

    J = Jx2+lambda_u*Ju2-lambda_d*Jd2-lambda_n*Jn2;

    % Warm start for next optimization (shift and  move away from constrainsts)

    x0Warm=[x1(:,1)];
    x1Warm=[x1(:,[2:end,end])];                % repeat last state
    uWarm=[uFuture(:,2:end),zeros(nU,1)];      % append zero control
    dWarm=[d(:,2:end),zeros(nD,1)];            % append zero disturbance
    uWarm=max(uWarm,repmat(-.95*max_u,[1,T])); % move away from -max_u
    uWarm=min(uWarm,repmat(+.95*max_u,[1,T])); % move away from +max_u
    dWarm=max(dWarm,repmat(-.95*max_d,[1,T])); % move away from -max_d
    dWarm=min(dWarm,repmat(+.95*max_d,[1,T])); % move away from +max_d

    % Solver output
    outputExpressions=struct('J',J,...
                             'Jx2',Jx2,...
                             'Ju2',Ju2,...
                             'Jd2',Jd2,...
                             'Jn2',Jn2,...
                             't',(-L:T)*Ts,...
                             'uFuture',uFuture,...
                             'd',d,...
                             'n',n,...
                             'x',x,...
                             'xEst',x1(:,L),...  % estimate of state at time 0
                             'ref',ref,...
                             'test',[d;
                                     repmat(-max_d,[1,L+T]);
                                     d-repmat(-max_d,[1,L+T]);
                                     repmat(max_d,[1,L+T])-d],...
                             'x0Warm',x0Warm,...
                             'x1Warm',x1Warm,...
                             'uWarm',uWarm,...
                             'dWarm',dWarm);

    classname=cmex2equilibriumLatentCS('classname','tmp_mpcmhe_dcmotor',...
                                       'P1objective',J,...
                                       'P2objective',-J,...
                                       'P1optimizationVariables',{uFuture},...
                                       'P1constraints',P1constraints,...
                                       'P2optimizationVariables',{x0,d},...
                                       'P2constraints',P2constraints,...
                                       'latentVariables',{x1},...
                                       'latentConstraints',dynamics,...
                                       'outputExpressions',outputExpressions,...
                                       'parameters',{...
                                           Ts;p;k;...
                                           uPast;yPast;...
                                           ref;...
                                           max_u;max_d;...
                                           lambda_u;lambda_d;lambda_n},...
                                       'scaleCost',0,...
                                       'scaleInequalities',false,...
                                       'solverVerboseLevel',4);  % use 4 to see solver iterations

end

%%%%%%%%%%%%%%%
%% Use solver
%%%%%%%%%%%%%%%

%% Create object
obj=tmp_mpcmhe_dcmotor();
mu0=1e-3;       % low value will speed-up convergence (up to a point where it causes numerical issues)
maxIter=100;
saveIter=-1;

s = RandStream('mt19937ar','Seed',0); % fix seed for reproducibility
RandStream.setGlobalStream(s);

%% Set parameters

ref=@(t)1*sign(sin(.5*t)); % reference signals

Ts=.05;
p=-2;k=1;
A=[0,1;0,p];                  % matrices for the dynamics
B=[0;k];
C=[1,0];
max_u=5;
max_d=10;
lambda_u=1/50;
lambda_d=50;
lambda_n=5;

setP_Ts(obj,Ts);
setP_p(obj,p);
setP_k(obj,k);
setP_max_u(obj,max_u);
setP_max_d(obj,max_d);
setP_lambda_u(obj,lambda_u);
setP_lambda_d(obj,lambda_d);
setP_lambda_n(obj,lambda_n);

% initial condition
xinit=[.2;.2];

% cold-start: random initialization close to initial condition
x0Warm=.01*rand(nX,1);
x1Warm=.01*rand(nX,T+L);
uWarm=.01*rand(nU,T);
dWarm=.01*rand(nD,T+L);

t=0;
clear history;   % structure where system's state, control, etc. will be stored
nSteps=500;
uPast=[];
yPast=[];
for k=1:nSteps

    if size(yPast,2)<L
        % not enough outputs - apply zero input
        out=struct('t',0,...
                   'x',xinit,...
                   'xEst',nan(nX,1),...
                   'uFuture',zeros(nU,1),...
                   'ref',ref(t),...
                   'J',NaN,...
                   'Jx2',NaN,...
                   'Ju2',NaN,...
                   'Jd2',NaN,...
                   'Jn2',NaN);
        iter=NaN;
        time=NaN;
        status=0;
    else

        %% Set parameters specific to this optimization
        setP_ref(obj,ref(t+(0:T-1)*Ts));
        setP_uPast(obj,uPast);
        setP_yPast(obj,yPast);

        %% Initialization primal variables (cold for k==1 and then warm)
        setV_x0(obj,x0Warm);
        setV_x1(obj,x1Warm);
        setV_uFuture(obj,uWarm);
        setV_d(obj,dWarm);

        full(dWarm)

        %% Solve optimization
        [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
        % Get outputs
        out=getOutputs(obj);

        k
        disp(out)
        
        if status
            out.d
            error('optimization failed\n');
        end

        x0Warm=out.x0Warm;
        x1Warm=out.x1Warm;
        uWarm=full(out.uWarm);
        dWarm=full(out.dWarm);
    end

    %% Save applied sate, control, etc. in `history` structure
    history.t(k)=t;
    history.x(1:nX,k)=xinit;
    history.xEst(1:nX,k)=out.xEst;
    history.ref(k)=out.ref(1);
    history.uFuture(1:nU,k)=out.uFuture(:,1);
    history.J(k)=out.J;
    history.Jx2(k)=out.Jx2;
    history.Ju2(k)=out.Ju2;
    history.Jd2(k)=out.Jd2;
    history.Jn2(k)=out.Jn2;
    history.iter(k)=iter;
    history.time(k)=time;

    %% Apply (constant) control to update state and output with random distrubance/noise
    d=.1*randn(nD,1);
    n=.2*randn(nD,1);
    [tout,yout]=ode23(@(t,x)A*x+B*(out.uFuture(:,1)+d),[0,Ts],xinit);

    %% Add data to past inputs and outputs
    yPast=[yPast,C*yout(1,:)'+n];     % output has delay of one sampling time
    history.y(k)=yPast(:,end);
    uPast=[uPast,out.uFuture(:,1)];
    % erase old data
    yPast(1:end-L)=[];     %
    uPast(1:end-L)=[];

    % update current time & state
    xinit=yout(end,:)';
    t=t+Ts;

    if mod(k,10)==1 || k==nSteps || status
        %% Plot MPC solution
        fig=clearFigure('figureNumber',20,'figureName','Latest MPC solution');
        plotData(out,Ts)
        %% Plot actual control
        fig=clearFigure('figureNumber',fig+1,'figureName','Actual control');
        plotData(history,Ts)
        %% Plot solver statistics
        fig=clearFigure('figureNumber',fig+1,'figureName','Solver statistics');
        subplot(2,1,1);
        plot(history.t,history.J,'.-',...
             history.t,history.Jx2,'.-',...
             history.t,history.Ju2,'.-',...
             history.t,history.Jd2,'.-',...
             history.t,history.Jn2,'.-')
        legend('J','Jx2','Ju2','Jd2','Jn2');grid on;
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
    plot(data.t,data.x(1,:),'b.-',...
         data.t(end-size(data.ref,2)+1:end),data.ref,'-'); grid on
    if isfield(data,'y')
        hold on
        plot(data.t,data.xEst(1,:),'g.',...
             data.t,data.y,'r.');
        legend('true','reference','estimated','measured');
    else
        legend('true','measured');
    end
    ylabel('\theta');
    subplot(3,1,2);
    plot(data.t,data.x(2,:),'b.-'); grid on
    if isfield(data,'y')
        hold on
        plot(data.t,data.xEst(2,:),'g.');
        legend('true','estimated');
    end
    ylabel('\omega');
    subplot(3,1,3);
    plot(data.t(1:size(data.uFuture,2)),data.uFuture,'.-'); grid on
    ylabel('u');

end
