% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all
% remove previous solvers
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Generate solver

% Create symbolic optimization

nx=2;    % state size;
nu=1;    % control size;
nd=1;    % disturbance size;
ny=1;    % output sizes
T=60;    % forward horizon
L=20;    % backward horizon
delay=1;

Tvariable Ts [];
Tvariable x        [nx,L+T+1];    % [x(t-L*Ts),...,x(t),...,x(t+T*Ts)]
Tvariable y_past   [ny,L+1];      % [y(t-L*Ts),...,y(t)]
Tvariable u_past   [nu,L+delay];  % [u(t-L*Ts),...,u(t+(delay-1)*Ts)]
Tvariable u_future [nu,T-delay];  % [u(t+delay*Ts), ...,u(t+(T-1)*Ts)]
Tvariable d        [nd,L+T];      % [d(t-L*Ts,...,x(t),...,x(t+(T-1)*Ts)]

Tvariable p  [1,1];
Tvariable k  [1,1];
Tvariable r  [1,T-delay];         % [r(t+(delay+1)*Ts), ...,u(t+T*Ts)]
Tvariable uMax [];
Tvariable dMax [];

% DC-motor-like transfer function
% [dot x1]=[0  1][x1]+[0]u
% [dot x2] [0  p][x2] [k]
%        y=x1     

dxFun=@(x,u,d,p,k,Ts,r,cc,uMax,dMax)[0,1;0,p]*x+[0;k]*(u+d);
yFun=@(x,u,d,p,k,Ts,r,cc,uMax,dMax)x(1,:);    % since yFun only uses 1st argument, it can be called with just one argument

u=[u_past,u_future];

JJ=[norm2(yFun(x(:,L+delay+2:L+T+1))-r)/(T-delay);
    norm2(u_future)/(T-delay);
    norm2(d)/(L+T);
    norm2(yFun(x(:,1:L+1))-y_past)/(L+1)];

Tvariable cc size(JJ);    % optimization weights

Jcc=cc.*JJ;
J=cc*JJ;

% create mpc object
mpcmhe=Tmpcmhe('reuseSolver',true,...
               'solverType','C',...
               'solverClassname','tmp1',...
               'sampleTime',Ts,...
               'stateVariable',x,...
               'pastControlVariable',u_past,...
               'pastOutputVariable',y_past,...
               'futureControlVariable',u_future,...
               'disturbanceVariable',d,...
               'stateDerivativeFunction',dxFun,...
               'outputFunction',yFun,...
               'objective',J,...
               'controlConstraints',{u_future<=uMax, u_future>=-uMax},...
               'disturbanceConstraints',{d<=dMax, d>=-dMax},...
               'outputExpressions',{J,JJ,cc,Jcc,x,u_past,u_future,d},...
               'parameters',{p,k,Ts,r,cc,uMax,dMax},...
               'solverParameters',{;
                    'scaleCost',0e4,...
                    'alphaMin',1e-4,...
                    'debugConvergence',false,...
                    'skipAffine',true,...
                    'solverVerboseLevel',3 ...
                   });

%% Simulate system

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% set parameter values
setParameter(mpcmhe,'p',2);
setParameter(mpcmhe,'k',5);
setParameter(mpcmhe,'uMax',50);
setParameter(mpcmhe,'dMax',10);
Ts=.05;
setParameter(mpcmhe,'Ts',Ts);

cc=[10;.05;-5e2;-5e2];
setParameter(mpcmhe,'cc',cc);

% set process initial condition, inputs, disturbances, and noise and
% get the corresponding measurements
t0=0;
x0=[.2;.2];
u0=.1*randn(nu,L+delay);
d0=0*.1*randn(nu,L);
n0=0*randn(nu,L+1);
[t,t_y_past,y_past,t_u_past,u_past,t_state,state]=setInitialState(mpcmhe,t0,x0,u0,d0,n0);

verbose=false;

if verbose
    fprintf('setInitialState()\n');
    fprintf('[t;u_past]\n');
    disp([t_u_past';u_past]);
    fprintf('[t;y_past]\n');
    disp([t_y_past';y_past]);
    fprintf('[t;real state]\n');
    disp([t_state';state]);
end

fig=1;
figure(fig);clf;

mu0=1;
maxIter=100;
saveIter=false;

ref=@(t)1*sign(sin(2*t));

% cold start
x_warm=.1*randn(nx,1);
u_warm=.1*randn(nu,T-delay);
d_warm=0*.1*randn(nd,L+T);

%x_warm=state(:,1); % cheating

for i=1:100
    % move warm-start away from constraints
    u_warm=min(u_warm,.95);
    u_warm=max(u_warm,-.95);
    d_warm=min(d_warm,.95);
    d_warm=max(d_warm,-.95);
    [t_state,state]=setSolverWarmStart(mpcmhe,u_past,y_past,x_warm,u_warm,d_warm);

    % set reference signal
    t_r=t+(delay+1:T)*Ts;
    r=ref(t_r);
    setParameter(mpcmhe,'r',r);
    
    if verbose
        fprintf('setSolverWarmStart()\n');
        fprintf('[t;u_past]\n');
        disp([t_u_past';u_past]);
        fprintf('[t;y_past]\n');
        disp([t_y_past';y_past]);
        fprintf('x0_warm\n');
        disp(x_warm);
        fprintf('x_warm\n');
        disp(state);
        fprintf('u_warm\n');
        disp(u_warm);
        fprintf('d_warm\n');
        disp(d_warm);
        fprintf('[t;warm state]\n');
        disp([t_state';state]);
    end
    
    [solution,J,JJ,cc,Jcc,x,u_past,u_future,d]=solve(mpcmhe,mu0,maxIter,saveIter);
    %u_future,d,
    if verbose
        fprintf('solve()\n');
        fprintf('x_opt\n');
        disp(x);
        fprintf('u_opt\n');
        disp(u_future);
        fprintf('d_opt\n');
        disp(d);
    end
    
    fprintf('J     =%10.3e computed in %g iterations & %g secs\n',J,solution.iter,solution.time);
    fprintf('JJ    =');fprintf('%10.3e ',JJ);fprintf('\n');
    fprintf('cc    =');fprintf('%10.3e ',cc);fprintf('\n');
    fprintf('JJ*cc =');fprintf('%10.3e ',Jcc);fprintf('\n');
    
    % apply optimal controls, get time, measurements, and warm start for the next iteration
    u_final=zeros(nu,1);
    d_final=zeros(nd,1);
    noise=.1*randn(ny,1);
    %disturbance=[];  % worst-case
    disturbance=.1*randn(nu,1);
    [t,t_y_past,y_past,t_u_past,u_past,x_warm,u_warm,d_warm,t_state,state]=...
        applyControl(mpcmhe,u_past,y_past,solution,disturbance,noise,u_final,d_final); 
    
    if verbose
        fprintf('applyControl()\n');
        fprintf('[t;u_past]\n');
        disp([t_u_past';u_past]);
        fprintf('[t;y_past]\n');
        disp([t_y_past';y_past]);
        fprintf('[t;real state]\n');
        disp([t_state';state]);
    
        fprintf('x0_warm\n');
        disp(x_warm);
        fprintf('x_warm\n');
        disp(state);
        fprintf('u_warm\n');
        disp(u_warm);
        fprintf('d_warm\n');
        disp(d_warm);

    end

    if true || mod(i,5)==0
        if 0
            plot(solution.t_state,solution.state,'.-',...
                 t+(-L:delay-1)*Ts,u_past,'.-',...
                 solution.t_futureControl,solution.futureControl,'.-',...
                 solution.t_disturbance,solution.disturbance,'.-',...
                 t_r,r,'.-');grid on;
            legend('x1','x2','u_past','u_future','d','r');
        else
            plot(t_r,r,'k.-',...
                 solution.t_state,solution.state(1,:),'.-');grid on;
            legend('r','x1');
        end
        hold on;
        %pause
    end
    
    if solution.status
        fprintf('solver failed at time %g (paused)\n',t);
        pause
    end

    %fprintf('paused\n');pause
end

history=getHistory(mpcmhe);

fig=fig+1;figure(fig);clf;
set(fig,'Name','Closed-loop (Tmpc)');
subplot(2,1,1)
plot(history.t,ref(history.t),'.-',...
     history.t,history.x(1,:),'.-',...
     history.t,history.y,'.-');grid on;
legend('r','x1','y');
subplot(2,1,2)
plot(history.t,history.x(2,:),'.-',...
     history.t,history.u,'.-',...
     history.t,history.d,'.-');grid on;
legend('x2','u','d');
xlabel('t');

fig=fig+1;figure(fig);clf;
set(fig,'Name','Solver (Tmpc)');
subplot(2,1,1);
plot(history.t,history.objective,'.-');grid on;
ylabel('MPCMHE cost');
subplot(2,1,2);
yyaxis left
plot(history.t,history.iter,'.-');grid on;
ylabel('# iterations');
yyaxis right
plot(history.t(2:end),1000*history.stime(2:end),'.-');grid on;
ylabel('solver time [ms]');
xlabel('t');




