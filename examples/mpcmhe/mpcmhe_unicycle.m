% MPC-MHE output-feedback control for a unicycle pursuer trying to
% catch a velocity controlled evader.
%
% Pursuer continuous-time model (unicycle):
%       dot x1 = v cos x3
%       dot x2 = v sin x3
%       dot x3 = u            u in [-uMax,uMax]
%
%   where
%       x1,x2 = pursuer's 2D position
%       x3    = pursuer's orientation [rad]
%       v     = pursuer's (constant) speed
%       u     = pursuer's turning rate -- pursuer's control signal
%
%       the pursuer's position x1,x2 can be measured (with noise), but
%       not the pursuer's orientation x3
%
%  Target continuous-time model:
%       dot x4 = d1
%       dot x5 = d2
%
%    where
%       d1,d2 = evader's 2D velocity -- evader's control signal
%
%       the current evader's position x1,x2 can be measured (with
%       noise), but not the velocity.
%
% The goal of MPC-MHE is to solve the following saddle-point optimization
%
%         min_uFuture max_{x(-L),d,n}
%                \int_0^T ( x1(t)-x4(t) )^2 + ( x2(t)-x5(t) )^2 dt
%                               + lambda_u int_{-L}^T u^2 dt
%                               - lambda_d int_{-L}^T d^2 dt
%                               - lambda_n int_{-L}^0 n^2 dt
%
% where
%   uPast   = [u(-L),....,u(-1)]  = past control inputs (selected previously and used with ZHO)
%   uFuture = [u(0),....,u(T-1)]  = future control inputs (to be selected and to be used with ZOH)
%
%   x(-L)                         = unknown initial condition (at the start of backward horizon),
%                                   which includes the states or pursue and evader.
%   d       = [d(-L),....,d(T-1)] = unknown pursuer's velocity (assumed with ZOH)
%   y       = [y(-L) ... y(-1)]   = past measured output, which includes noisy measurememnts
%                                   of the the positions of pursuer and evader
%   n       = [n(-L) ... n(-1)]   = past measurement noise
%
% subject to the following constraints
%
%        |u(t)| <= max_u,     \forall t    [pursuer's constraint]
%        ||d(t)|| <= max_d,   \forall t    [evader's constraint]
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

help mpcmhe_unicycle

createSolvers=input('Do you want to (re)create the solvers (Y/n)? ','s');
createSolvers=~isequal(createSolvers,'n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nX=5;                        % state size
nU=1;                        % control size
nD=2;                        % disturbance size
nY=4;                        % output size

T=80;                        % forward horizon length in time-steps
L=30;                        % backward horizon length in time-steps

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

    Tvariable v [];               % pursuer's velocity

    % System dynamics
    x=[x0,x1];
    u=[uPast,uFuture];
    dynamics={
        % option 1) forward Euler
        x(:,2:end)==x(:,1:end-1)+Ts*[v*cos(x(3,1:end-1));...
                                     v*sin(x(3,1:end-1));...
                                     u;
                                     d];
        % option 2) trapesoidal, with ZOH for u
        %x(:,2:end)==x(:,1:end-1)+Ts*[.5*v*(cos(x(3,1:end-1))+cos(x(3,2:end)));...
        %                             .5*v*(sin(x(3,1:end-1))+sin(x(3,2:end)));...
        %                             u;
        %                             d];
             };

    % Constraints
    Tvariable max_u [];       % maximum control input
    Tvariable max_d [];           % maximum disturbance input

    P1constraints={
        uFuture.^2<=max_u^2
                };
    P2constraints={
        sum(d.^2,1)<=max_d^2;
                };

    % Criterion
    Tvariable lambda_u [];         % weigth for control input (pursuer turning rate)
    Tvariable lambda_d [];         % weigth for input disturbance (evader velocity)
    Tvariable lambda_n [];         % weigth for measurement noise

    errFuture=x(1:2,end-T+1:end)-x(4:5,end-T+1:end);  % position between pursuer and evader (future)
    Jerr2=tsIntegral(sum(errFuture.^2,1),Ts);         % integral of (future) relative position^2

    Ju2=tsIntegral(sum(uFuture.^2,1),Ts);          % integral of (future) (dot u)^2

    Jd2=tsIntegral(sum(d.^2,1),Ts);                % integral of disturbance^2

    n=x([1,2,4,5],1:L)-yPast;                         % measurement noise
    Jn2=tsIntegral(sum(n.^2,1),Ts);                   % integral of noise^2

    J = Jerr2+lambda_u*Ju2-lambda_d*Jd2-lambda_n*Jn2;

    % Warm start for next optimization (shift and  move away from constrainsts)

    x0Warm=[x1(:,1)];
    x1Warm=[x1(:,[2:end,end])];                % repeat last state
    uWarm=[uFuture(:,2:end),zeros(nU,1)];      % append zero control
    dWarm=[d(:,2:end),zeros(nD,1)];            % append zero disturbance
    uWarm=max(uWarm,-.9*max_u);                % move away from -max_u
    uWarm=min(uWarm,+.9*max_u);                % move away from +max_u
    dWarm=max(dWarm,-.9*max_d/sqrt(2));        % move away from -max_d
    dWarm=min(dWarm,+.9*max_d/sqrt(2));        % move away from +max_d

    % Solver output
    outputExpressions=struct('J',J,...
                             'Jerr2',Jerr2,...
                             'Ju2',Ju2,...
                             'Jd2',Jd2,...
                             'Jn2',Jn2,...
                             't',(-L:T)*Ts,...
                             'uFuture',uFuture,...
                             'd',d,...
                             'n',n,...
                             'x',x,...
                             'x0Warm',x0Warm,...
                             'x1Warm',x1Warm,...
                             'uWarm',uWarm,...
                             'dWarm',dWarm);

    classname=cmex2equilibriumLatentCS('classname','tmp_mpcmhe_unicycle',...
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
                                           Ts;v;...
                                           uPast;yPast;...
                                           max_u;max_d;...
                                           lambda_u;lambda_d;lambda_n},...
                                       ...
                                       'gradTolerance',1e-7,...
                                       'equalTolerance',1e-8,...
                                       'desiredDualityGap',1e-5,...
                                       ...
                                       'scaleCost',0,...
                                       ...%'scaleInequalities',false,...
                                       'muFactorConservative',.99,...
                                       ...
                                       'allowSave',true,...  % needed to save hessian for pivoting
                                       'solverVerboseLevel',4);  % use 4 to see solver iterations

end

%%%%%%%%%%%%%%%
%% Use solver
%%%%%%%%%%%%%%%

%% Create object
obj=tmp_mpcmhe_unicycle();
mu0=1e-1;       % low value will speed-up convergence (up to a point where it causes numerical issues)
maxIter=500;
saveIter=-1; % save on error

s = RandStream('mt19937ar','Seed',0); % fix seed for reproducibility
RandStream.setGlobalStream(s);

%% Set parameters

Ts=.1;

v=1;
max_u=1.5;
max_d=.5;

% needs some penalty on u and d to regularize solution and prevent singularity (constant changes in direction)
lambda_u=1;
lambda_d=1;

lambda_n=1e3;

setP_Ts(obj,Ts);
setP_v(obj,v);
setP_max_u(obj,max_u);
setP_max_d(obj,max_d);
setP_lambda_u(obj,lambda_u);
setP_lambda_d(obj,lambda_d);
setP_lambda_n(obj,lambda_n);

% initial condition
xinit=[0;0;0;     % pursuer at origin, facing right 
       2+L*Ts*v;  % evader in front of pursuer
       2];        % evader above pursuer

t=0;
clear history;   % structure where system's state, control, etc. will be stored
nSteps=200;
uPast=[];
yPast=[];
for k=1:40/Ts

    if size(yPast,2)<L
        % not enough outputs - apply zero input
        out=struct('t',(-L:T)*Ts,...
                   'x',nan(nX,L+T+1),...
                   'uFuture',zeros(nU,T),...
                   'd',zeros(nD,L+T),...
                   'J',NaN,...
                   'Jerr2',NaN,...
                   'Ju2',NaN,...
                   'Jd2',NaN,...
                   'Jn2',NaN,...
                   'x0Warm',[xinit(4:5);0;xinit(4:5)]+.01*rand(nX,1),...
                   'x1Warm',[xinit(4:5);0;xinit(4:5)]+.01*rand(nX,T+L),...
                   'uWarm',max_u/6*randn(nU,T),...
                   'dWarm',max_d/6*randn(nD,T+L));
        iter=NaN;
        time=NaN;
        status=0;
    
        % cold-start: random initialization close to evader
        x0Warm=[xinit(4:5);0;xinit(4:5)]+.01*rand(nX,1);
        x1Warm=[xinit(4:5);0;xinit(4:5)]+.01*rand(nX,T+L);
        uWarm=max_u/6*randn(nU,T);
        dWarm=max_d/6*randn(nD,T+L);

    else

        %% Set parameters specific to this optimization
        setP_uPast(obj,uPast);
        setP_yPast(obj,yPast);

        %% Initialization primal variables (cold for k==1 and then warm)
        setV_x0(obj,x0Warm);
        setV_x1(obj,x1Warm);
        setV_uFuture(obj,uWarm);
        setV_d(obj,dWarm);

        %% Solve optimization
        [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
        k

        if status
            if status==4 && iter==1
                fprintf('Solver needs to be recreated with saved matrix values to improve pivoting!\n');
                return
            else
                warning('optimization failed, re-using results from last optimization\n');
                out.x(:,1)=[];
                out.uFuture(:,1)=[];
                out.d(:,1)=[];
            end
        else
            % Get outputs
            out=getOutputs(obj);
            disp(out);
        end

        x0Warm=out.x0Warm;
        x1Warm=out.x1Warm;
        uWarm=full(out.uWarm);
        dWarm=full(out.dWarm);
    end

    %% Save applied sate, control, etc. in `history` structure
    history.t(k)=t;
    history.x(1:nX,k)=xinit;
    history.uFuture(1:nU,k)=out.uFuture(:,1);
    history.d(1:nD,k)=out.d(:,1);
    history.J(k)=out.J;
    history.Jerr2(k)=out.Jerr2;
    history.Ju2(k)=out.Ju2;
    history.Jd2(k)=out.Jd2;
    history.Jn2(k)=out.Jn2;
    history.iter(k)=iter;
    history.time(k)=time;

    %% Apply (constant) control to update state and output with random distrubance/noise
    u=out.uFuture(:,1);
    if k<L+50*0
        d=[max_d;0];     % start's by moving left
    else
        d=out.d(:,L+1);  % actually evade
    end
    n=.005*randn(nY,1);
    [tout,yout]=ode23(@(t,x)[v*cos(x(3,:));
                             v*sin(x(3,:));
                             u;
                             d],[0,Ts],xinit);

    %% Add data to past inputs and outputs
    yPast=[yPast,yout(1,[1,2,4,5])'+n];     % output has delay of one sampling time
    history.y(1:nY,k)=yPast(:,end);
    uPast=[uPast,out.uFuture(:,1)];
    % erase old data
    yPast(:,1:end-L)=[];
    uPast(:,1:end-L)=[];

    % update current time & state
    xinit=yout(end,:)';
    t=t+Ts;

    if mod(k,5)==0 || k==nSteps || status
        %% Plot MPC solution
        fig=clearFigure('figureNumber',30,'figureName','Latest MPC solution');
        plotData(out,L+1,Ts,max_u,max_d)
        %% Plot actual control
        fig=clearFigure('figureNumber',fig+1,'figureName','Actual control');
        plotData(history,k,Ts,max_u,max_d)
        %% Plot solver statistics
        fig=clearFigure('figureNumber',fig+1,'figureName','Solver statistics');
        subplot(2,1,1);
        plot(history.t,history.J,'.-',...
             history.t,history.Jerr2,'.-',...
             history.t,history.Ju2,'.-',...
             history.t,history.Jd2,'.-',...
             history.t,history.Jn2,'.-')
        legend('J','Jerr2','Ju2','Jd2','Jn2');grid on;
        subplot(2,1,2);
        yyaxis left;plot(history.t,history.iter,'.-');ylabel('# solver iter');grid on;
        yyaxis right;plot(history.t,1000*history.time,'.-');ylabel('solver time [ms]');grid on;
        drawnow;
    end

end

% Delete object
clear obj

function plotData(data,now,Ts,max_u,max_d)
    subplot(3,1,1:2);cla
    plot(data.x(1,:),data.x(2,:),'k.-',...
         data.x(4,:),data.x(5,:),'b.-',...
         data.x(1,now),data.x(2,now),'r*',...
         data.x(4,now),data.x(5,now),'r*');
    axis equal;
    grid on;
    legend('pursuer','evader','location','best');
    subplot(3,1,3);
    yyaxis left
    plot(data.t(end-size(data.uFuture,2)+1:end),data.uFuture,'.-');grid on
    ylim(1.1*[-max_u,max_u])
    ylabel('u');
    yyaxis right
    plot(data.t(1:size(data.d,2)),sqrt(sum(data.d.^2,1)),'.-');grid on
    ylim(1.1*[0,max_d]);
    ylabel('||d||');
end
