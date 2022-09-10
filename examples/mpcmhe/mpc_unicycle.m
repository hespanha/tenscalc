% MPC state-feedback control for a vehicle modeled as a unicycle
% trying to catch a kinematic target.
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
%       u     = pursuer's turning rate
%
%  Target continuous-time model:
%       dot x4 = d1
%       dot x5 = d2
%
%    where
%       d1,d2 = target's "known" constant 2D velocity, which can be viewed as
%               known constant distance
%
% The goal of MPC is to minimize the integral of the distance between
% pursuer and evader, which corresponds to a criterion of the form:
%
%         J = \int_0^T ( x1(t)-x4(t) )^2 + ( x2(t)-x5(t) )^2 dt
%
% subject to turning rate constraints
%
%         |u| <= uMax
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-22 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all;
% remove previous solvers
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

help mpc_unicyle

createSolvers=input('Do you want to create the solvers (Y/n)? ','s');
createSolvers=~isequal(createSolvers,'n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=200;                       % horizon length in time-steps
nX=5;                        % state size
nU=1;                        % control size
nD=2;                        % disturbance size
if createSolvers

    %%%%%%%%%%%%%%%%%%%
    %% Generate solver
    %%%%%%%%%%%%%%%%%%%

    Tcalculus.clear();

    Tvariable Ts [];              % sampling time
    Tvariable x [nX,T];           % state:              [x(1) x(2) ... x(T)]
    Tvariable u [nU,T-1];         % control:            [u(1) u(2) ... u(T-1)] (ZOH)
    Tvariable d [nD,1];           % (constant) disturbance

    Tvariable xinit [nX,1];       % initial state:      x(1)

    Tvariable v [];      % pursuer's velocity

    % System dynamics
    dynamics={
        % option 1) forward Euler
        %x(:,2:end)==x(:,1:end-1)+Ts*[v*cos(x(3,1:end-1));
        %                             v*sin(x(3,1:end-1));
        %                             u;
        %                             repmat(d,[1,T-1])];
        % option 2) trapesoidal, with ZOH for u & d
        x(:,2:end)==x(:,1:end-1)+Ts*[v*(cos(x(3,1:end-1))+cos(x(3,2:end)))/2;
                                     v*(sin(x(3,1:end-1))+sin(x(3,2:end)))/2;
                                     u;
                                     repmat(d,[1,T-1])];
        % initial state
        x(:,1)==xinit;
             };

    % Constraints
    Tvariable max_u [];   % max turning rate

    constraints={ u >= -max_u;
                  u <= +max_u };

    % Criterion
    J = tsIntegral(sum((x(1:2,:)-x(4:5,:)).^2,1),Ts);     % integral of distance^2

    % Warm start for next optimization (shift and  move away from constrainsts)

    uWarm=[u(:,2:end),zeros(nU,1)];   % append zero control
    xWarm=[x(:,[2:end,end])];         % repeat last state
    uWarm=max(uWarm,-.9*max_u);       % move away from -max_u
    uWarm=min(uWarm,.9*max_u);        % move away from max_u

    % Solver output
    outputExpressions=struct('J',J,...
                             't',(0:T-1)*Ts,...
                             'u',u,...
                             'x',x,...
                             'eq',x(:,2:end)-(x(:,1:end-1)+Ts*[v*cos(x(3,1:end-1));
                        v*sin(x(3,1:end-1));
                        u;
                        repmat(d,[1,T-1])]),...
                             'xWarm',xWarm,...
                             'uWarm',uWarm);

    classname=cmex2optimizeCS('classname','tmp_unicycle_pursuit',...
                              'objective',J,...
                              'optimizationVariables',{u,x},...
                              'constraints',[dynamics;constraints],...
                              'outputExpressions',outputExpressions,...
                              'parameters',{...
                                  Ts;
                                  v;d;...
                                  xinit;...
                                  max_u},...
                              'addEye2Hessian',true,...
                              'adjustAddEye2Hessian',true,...
                              'scaleInequalities',true,...
                              'useInertia',true,...
                              ...%'scaleCost',1e4,...
                              ...%'muFactorAggressive',.25,...
                              'solverVerboseLevel',4);  % use 4 to see solver iterations

end

%%%%%%%%%%%%%%%
%% Use solver
%%%%%%%%%%%%%%%

%% Create object
obj=tmp_unicycle_pursuit();
mu0=1e0;       % low value will speed-up convergence (up to a point where it causes numerical issues)
maxIter=300;
saveIter=-1;
addEye2Hessia=[1e-8;1e-8];

s = RandStream('mt19937ar','Seed',0); % fix seed for reproducibility
RandStream.setGlobalStream(s);

%% Set parameters

Ts=.01;
d=[1;0]; % evader's (constant) velocity
v=2;     % pursuer's speed
max_u=2;

setP_Ts(obj,Ts);
setP_d(obj,d);
setP_v(obj,v);
setP_max_u(obj,max_u);

% initial condition
t=0;
xinit=[-1;-1;0;0;0];

% cold-start: random initialization close to initial condition
xWarm=xinit+.01*rand(nX,T);
uWarm=.01*rand(nU,T-1);

t=0;
clear history;   % structure where system's state, control, etc. will be stored
nSteps=1000;
for k=1:nSteps

    %% Set parameters specific to this optimization
    setP_xinit(obj,xinit);

    %% Initialization primal variables (cold for k==1 and then warm)
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
    history.u(1:nU,k)=out.u(:,1);
    history.J(k)=out.J;
    history.iter(k)=iter;
    history.time(k)=time;

    %% Apply (constant) control to update initial position and velocity
    [tout,yout]=ode23(@(t,x)[v*cos(x(3,:));
                             v*sin(x(3,:));
                             out.u(:,1);
                             d],[0,Ts],xinit);
    xinit=yout(end,:)';
    t=t+Ts;
    xWarm=out.xWarm;
    uWarm=out.uWarm;
    
    %xinit,full(xWarm(:,1:2))

    if mod(k,10)==1 || k==nSteps || status
        %% Plot MPC solution
        fig=clearFigure('figureNumber',10,'figureName','Latest MPC solution');
        plotData(out,Ts)
        %% Plot actual control
        fig=clearFigure('figureNumber',fig+1,'figureName','Actual control');
        plotData(history,Ts)
        %% Plot solver statistics
        fig=clearFigure('figureNumber',fig+1,'figureName','Solver statistics');
        subplot(2,1,1);
        plot(history.t,history.J,'.-')
        ylabel('J');grid on;
        subplot(2,1,2);
        yyaxis left;plot(history.t,history.iter,'.-');ylabel('# solver iter');grid on;
        yyaxis right;plot(history.t,1000*history.time,'.-');ylabel('solver time [ms]');grid on;
        drawnow;
    end

    if status
        error('optimization failed\n');
    end

    %fprintf('paused\n');pause;
end

% Delete object
clear obj

function plotData(data,Ts)

    plot(data.x(1,:),data.x(2,:),'k.-',...
         data.x(4,:),data.x(5,:),'b.-',...
         data.x(1,1),data.x(2,1),'k*',...
         data.x(4,1),data.x(5,1),'b*');
    axis equal;
    grid on;
    legend('pursuer','evader','location','best');
end