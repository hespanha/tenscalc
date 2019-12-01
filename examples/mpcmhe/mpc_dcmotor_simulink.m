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
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Generate solver

% Create symbolic optimization

nX=2;    % state size
nU=1;    % control size
T=30;    % forward horizon

% DC-motor-like transfer function
% [dot x1]=[0  1][x1]+[0]u
% [dot x2] [0  p][x2] [k]
%        y=x1     

Tvariable p  [1,1];
Tvariable k  [1,1];
dxFun=@(x,u,p,k)[0,1;0,p]*x+[0;k]*u;

[Ts,xMeas,xFut,~,uFut,ode]=TvariablesMPC(nX,nU,T,0,dxFun,p,k);

% parameters for state and input constraints
Tvariable uMax [];
Tvariable xMax [2,1];

% Criterium

Tvariable r [1,T];   % reference
J=norm2(xFut(1,:)-r)+norm2(uFut)/50;

% Warm start

xFutWarm=[xFut(:,2:end),zeros(nX,1)];
uFutWarm=[uFut(:,2:end),zeros(nU,1)];
xFutWarm=min(xFutWarm,.9*repmat(xMax,[1,T]));
xFutWarm=max(xFutWarm,-.9*repmat(xMax,[1,T]));
uFutWarm=min(uFutWarm,.9*uMax);
uFutWarm=max(uFutWarm,-.9*uMax);

% inequality constraints
constraints={uFut<=uMax, uFut>=-uMax,...
             xFut<=repmat(xMax,[1,T]),xFut>=-repmat(xMax,[1,T])};

[classname,code]=cmex2optimizeCS(...
    'classname','tmp_mpc_dcmotor',...
    ...
    'objective',J,...
    'optimizationVariables',{xFut,uFut},...
    'constraints',{ode,constraints{:}},...
    'parameters',{p,k,Ts,r,uMax,xMax,xMeas},...
    'outputExpressions',{J,[xMeas,xFut],uFut,xFutWarm,uFutWarm},...
    ...
    'simulinkLibrary','tmp_mpc_dcmotor_SL',...
    ...
    'compilerOptimization','-O0',...
    'solverVerboseLevel',2);
    
%open_system('tmp_mpc_dcmotor_SL');


%% Simulate system

% set parameter values
p=2;k=1;
Ts=.1;

uMax=1;
xMax=[.4;.3];

refmax=.35;
refomega=.5;

% set process initial condition
t0=0;
xMeas=[.2;.2];

% cold start
uFutCold=.1*randn(nU,T);
xFutCold=.1*randn(nX,T);

mu0=1;
maxIter=50;
saveIter=false;

open_system('mpc_dcmotor_simulink_S');
set_param('mpc_dcmotor_simulink_S',...
          'RelTol','1e-6',...
          'StartTime','0','StopTime','15');
out=sim('mpc_dcmotor_simulink_S');

fig=21;figure(fig);clf;
set(fig,'Name','Closed-loop (simulink)');
plot(x.time,squeeze(x.signals.values),'.-',...
     u.time,squeeze(u.signals.values),'.-',...
     ref.time,squeeze(ref.signals.values),'.-');grid on;
legend('x1','x2','u','r');
xlabel('t');

fig=fig+1;figure(fig);clf;
set(fig,'Name','Solver (simulink)');
subplot(2,1,1);
plot(cost.time,squeeze(cost.signals.values),'.-');grid on;
ylabel('MPC cost');
subplot(2,1,2);
yyaxis left
plot(iter.time,squeeze(iter.signals.values),'.-');grid on;
ylabel('# iterations');
yyaxis right
plot(stime.time,1000*squeeze(stime.signals.values),'.-');grid on;
ylabel('solver time [ms]');
xlabel('t');



