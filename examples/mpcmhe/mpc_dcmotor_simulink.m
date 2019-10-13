% Copyright 2012-2019 Joao Hespanha

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

clear all
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Generate solver

% Create symbolic optimization

nx=2;  % state size;
nu=1;  % control size
T=30;  % forward horizon;
nU=1; % number of controls to apply

Tvariable Ts [];
Tvariable x0 [nx,1];  % [x(t)]
Tvariable x1 [nx,T]; % [x(t+Ts), x(t+2*Ts), ..., x(t+T*Ts)];
Tvariable u  [nu,T];  % [u(t), u(t+Ts), ..., u(t+(T-1)*wTs) ] -- held constant between sample times
Tvariable p  [1,1];
Tvariable k  [1,1];

% parameters for state and input constraints
Tvariable umax [];
Tvariable xmax [2,1];

% DC-motor-like transfer function
% [dot x1]=[0  1][x1]+[0]u
% [dot x2] [0  p][x2] [k]
%        y=x1     

x=[x0,x1];  % whole state;

% ODE: (x(k+1)-x(k)) /Ts == [ f (x(k+1),u(k)) + f (x(k),u(k)) ]/2
ode={x1-x(:,1:end-1)==Ts*( .5*[0,1;0,p]*(x1+x(:,1:end-1))+[0;k]*u)};

% Criterium

Tvariable r [1,T];   % reference
J=norm2(x1(1,:)-r)+norm2(u)/50;

% Warm start

x1warm=[x1(:,nU+1:end),zeros(nx,nU)];
uwarm=[u(:,nU+1:end),zeros(nu,nU)];
x1warm=min(x1warm,.9*repmat(xmax,[1,T]));
x1warm=max(x1warm,-.9*repmat(xmax,[1,T]));
uwarm=min(uwarm,.9*umax);
uwarm=max(uwarm,-.9*umax);

% inequality constraints
constraints={u<=umax, u>=-umax,...
             x1<=repmat(xmax,[1,T]),x1>=-repmat(xmax,[1,T])};

[classname,code]=cmex2optimizeCS(...
    'classname','tmp_mpc_dcmotor',...
    ...
    'objective',J,...
    'optimizationVariables',{x1,u},...
    'constraints',[ode,constraints],...
    'parameters',{p,k,Ts,r,umax,xmax,x0},...
    'outputExpressions',{J,x,u,x1warm,uwarm},...
    ...
    'simulinkLibrary','tmp_mpc_dcmotor_SL',...
    ...
    'compilerOptimization','-O0',...
    'solverVerboseLevel',2);
    
open_system('tmp_mpc_dcmotor_SL');



%% Simulate system

% set parameter values
p=2;k=1;
Ts=.1;

umax=1;
xmax=[.4;.3];

refmax=.35;
refomega=.5;

% set process initial condition
t0=0;
x0=[.2;.2];

mu0=1;
maxIter=50;
saveIter=false;

% cold start
u_cold=.1*randn(nu,T);
x1_cold=.1*randn(nx,T);

mu0=1;
maxIter=50;
saveIter=0;

open_system('mpc_dcmotor_simulink_S');
set_param('mpc_dcmotor_simulink_S','StartTime','0','StopTime','15');
out=sim('mpc_dcmotor_simulink_S');

fig=20;figure(fig);clf;
set(fig,'Name','Closed-loop');
plot(x.time,squeeze(x.signals.values),'.-',...
     u.time,squeeze(u.signals.values),'.-',...
     ref.time,squeeze(ref.signals.values),'.-');grid on;
legend('x1','x2','u','r');
xlabel('t');

fig=fig+1;figure(fig);clf;
set(fig,'Name','Solver');
yyaxis left
plot(iter.time,squeeze(iter.signals.values),'.-');grid on;
ylabel('# iterations');
yyaxis right
plot(stime.time,1000*squeeze(stime.signals.values),'.-');grid on;
ylabel('solver time [ms]');
xlabel('t');



