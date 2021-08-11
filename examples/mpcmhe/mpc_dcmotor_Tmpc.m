% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all
% remove previous solvers
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Generate solver

% Create symbolic optimization

nx=2;  % state size
nu=1;  % control size
T=30;  % forward horizon

Tvariable Ts [];
Tvariable x [nx,T];  % [x(t+Ts), x(t+2*Ts), ..., x(t+T*Ts)]
Tvariable u [nu,T];  % [u(t), u(t+Ts), ..., u(t+(T-1)*Ts) ]
Tvariable p [1,1];
Tvariable k [1,1];

% parameters for state and input constraints
Tvariable uMax [];
Tvariable xmax [2,1];

% DC-motor-like transfer function
% [dot x1]=[0  1][x1]+[0]u
% [dot x2] [0  p][x2] [k]
%        y=x1     

dxFun=@(x,u,p,k,Ts,r,uMax,xmax)[0,1;0,p]*x+[0;k]*u;

% Criterium

Tvariable r [1,T];   % reference
J=norm2(x(1,:)-r)+norm2(u)/50;

% create mpc object
mpc=Tmpc('reuseSolver',true,...
         'solverType','C',...
         'solverClassname','tmp1',...
         'sampleTime',Ts,...
         'controlVariable',u,...
         'stateVariable',x,...
         'stateDerivativeFunction',dxFun,....
         'objective',J,...
         'constraints',{u<=uMax, u>=-uMax,...
                        x<=repmat(xmax,[1,T]),x>=-repmat(xmax,[1,T])},...
         'outputExpressions',{J,x,u},...
         'parameters',{p,k,Ts,r,uMax,xmax},...
         'solverParameters',{;
                    'compilerOptimization','-O0',...
                    'solverVerboseLevel',2 ...
                   });

%% Simulate system

% set parameter values
p=2;k=1;
setParameter(mpc,'p',p);
setParameter(mpc,'k',k);
Ts=.1;
setParameter(mpc,'Ts',Ts);

uMax=1;
xmax=[.4;.3];
setParameter(mpc,'uMax',uMax);
setParameter(mpc,'xmax',xmax);

refmax=.35;
refomega=.5;

% set process initial condition
t0=0;
x0=[.2;.2];
setInitialState(mpc,t0,x0);

% reference signals
ref=@(t)-refmax*sign(sin(refomega*t));

% cold start
u_warm=.1*randn(nu,T);               

mu0=1;
maxIter=50;
saveIter=false;

nSteps=150;

fig=10;
figure(fig);clf;
set(fig,'Name','MPC solutions (Tmpc)');

t=t0;
for i=1:nSteps
    % set reference signal
    r=ref(t+(0:T-1)*Ts);
    setParameter(mpc,'r',r);
    
    % move warm start away from constraints
    u_warm=min(u_warm,uMax-.05);
    u_warm=max(u_warm,-uMax+.05);
    x_warm=setSolverWarmStart(mpc,u_warm);
    % move state away from constraints
    x_warm=min(x_warm,xmax-.05);
    x_warm=max(x_warm,-xmax+.05);
    setSolverStateStart(mpc,x_warm);
    
    [solution,J,x,u]=solve(mpc,mu0,maxIter,saveIter);
    
    if solution.status==0
        fprintf('t=%g, J=%g computed in %g iterations & %g ms\n',t,J,solution.iter,1e3*solution.time);
    else
        fprintf('t=%g, J=%g computed in %g iterations & %g ms\n',t,J,solution.iter,1e3*solution.time);
        disp(solution);
        warning('solver failed')
    end
    
    if true %mod(t/Ts,5)==0
        if 0
            plot(t:Ts:t+Ts*(T-1),x,'.-',t:Ts:t+Ts*(T-1),u,'.-',t:Ts:t+Ts*(T-1),r,'.-');grid on;
            legend('x1','x2','u','r');
        else
            plot(t:Ts:t+Ts*(T-1),r,'k.-',t:Ts:t+Ts*(T-1),x(1,:),'.-');grid on;
            legend('x1','r');
        end
        hold on;
    end
    
    ufinal=zeros(nu,1);
    [t,u_warm]=applyControls(mpc,solution,ufinal); 
end

history=getHistory(mpc);

fig=fig+1;figure(fig);clf;
set(fig,'Name','Closed-loop (Tmpc)');
plot(history.t,history.x,'.-',...
     history.t,history.u,'.-',history.t,ref(history.t),'.-');grid on;
legend('x1','x2','u','r');
xlabel('t');

fig=fig+1;figure(fig);clf;
set(fig,'Name','Solver (Tmpc)');
subplot(2,1,1);
plot(history.t,history.objective,'.-');grid on;
ylabel('MPC cost');
subplot(2,1,2);
yyaxis left
plot(history.t,history.iter,'.-');grid on;
ylabel('# iterations');
yyaxis right
plot(history.t(2:end),1000*history.stime(2:end),'.-');grid on;
ylabel('solver time [ms]');
xlabel('t');



