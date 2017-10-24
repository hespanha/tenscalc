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

clear all
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Generate solver

% Create symbolic optimization

nx=5;             % state size;
nu=1;             % control size
T=200;            % forward horizon;
controlDelay=0;   % delay in applying control, after state measurement

% create symbolic optimization
Tvariable Ts [];
Tvariable x [nx,T];  % [x(t+Ts), x(t+2*Ts), ..., x(t+T*Ts)]
Tvariable u [nu,T];  % [u(t), u(t+Ts), ..., u(t+(T-1)*Ts) ]

Tvariable v [];      % pursuer's velocity
Tvariable d [2,1];   % evader's velocity
Tvariable uMax [];   % max turning rate

% pursuer
% dot x1 = v cos x3
% dot x2 = v sin x3
% dot x3 = u            u in [-uMax,uMax]

% evader
% dot x4    = d1
% dot x5    = d2

% controlled output
% z1 = x1 - x4
% z2 = x2 - x5

dxFun=@(x,u,v,d,Ts,uMax)[v*cos(x(3,:));
                    v*sin(x(3,:));
                    u;
                    repmat(d,[1,size(x,2)])];

J=norm2(x(1:2,:)-x(4:5,:));

% create mpc object
mpc=Tmpc('reuseSolver',true,...
         'solverType','C',...
         'solverClassname','tmp2',...
         'sampleTime',Ts,...
         'controlVariable',u,...
         'controlDelay',controlDelay,...
         'stateVariable',x,...
         'stateDerivativeFunction',dxFun,....
         'objective',J,...
         'constraints',{u<=uMax, u>=-uMax},...
         'outputExpressions',{J,x,u},...
         'parameters',{v,d,Ts,uMax},...
         'solverParameters',{;
                    'compilerOptimization','-O0',...
                    'solverVerboseLevel',3 ...
                   });

%% Simulate system

% set parameter values
setParameter(mpc,'d',[1;0]);
setParameter(mpc,'v',2);
uMax=2;
setParameter(mpc,'uMax',uMax);
Ts=.01;
setParameter(mpc,'Ts',Ts);

% set process initial condition
t0=0;
x0=[-1;-1;0;0;0];
setInitialState(mpc,t0,x0);

fig=1;
figure(fig);clf;

mu0=1;
maxIter=20;
saveIter=false;

t=t0;
figure(fig);clf;

% cold start
u_warm=.1*randn(nu,T-controlDelay);               

for i=1:400
    % move warm start away from constraints
    u_warm=min(u_warm,.95*uMax);
    u_warm=max(u_warm,-.95*uMax);
    setSolverWarmStart(mpc,u_warm);
    
    [solution,J,x,u]=solve(mpc,mu0,maxIter,saveIter);
    
    fprintf('t=%g, J=%g computed in %g iterations & %g ms\n',t,J,solution.iter,1e3*solution.time);
    
    if mod(i,5)==0
        plot(x(1,:),x(2,:),'k.-',...
             x(4,:),x(5,:),'b.-');
        axis equal;
        grid on;
        legend('pursuer','evader','location','best');
        title(sprintf('t=%g, J=%g computed in %g iterations & %g ms\n',t,J,solution.iter,1e3*solution.time));
        %hold on;
        drawnow
        %pause
    end
    
    % apply 5 controls and get time, warm start for next iteration
    ufinal=zeros(nu,1);
    [t,u_warm]=applyControls(mpc,solution,5,ufinal); 
end

history=getHistory(mpc);

fig=fig+1;figure(fig);clf;
plot(history.x(1,:),history.x(2,:),'k.-',...
     history.x(4,:),history.x(5,:),'b.-');
axis equal;
grid on;
legend('pursuer','evader','location','best');
title('full trajectory');

fig=fig+1;figure(fig);clf;
yyaxis left
plot(history.t,history.iter,'.-');grid on;
ylabel('# iterations');
yyaxis right
plot(history.t(2:end),1000*history.stime(2:end),'.-');grid on;
ylabel('solver time [ms]');
xlabel('t');



