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

nx=2;  % state size;
nu=1;  % control size
T=30;  % forward horizon

% create symbolic optimization
Tvariable Ts [];
Tvariable x [nx,T];  % [x(t+Ts), x(t+2*Ts), ..., x(t+T*Ts)]
Tvariable u [nu,T];  % [u(t), u(t+Ts), ..., u(t+(T-1)*Ts) ]
Tvariable p [1,1];
Tvariable k [1,1];
Tvariable r [1,T];

% DC-motor-like transfer function
% [dot x1]=[0  1][x1]+[0]u
% [dot x2] [0  p][x2] [k]
%        y=x1     

dxFun=@(x,u,p,k,Ts,r)[0,1;0,p]*x+[0;k]*u;

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
         'constraints',{u<=1, u>=-1},...
         'outputExpressions',{J,x,u},...
         'parameters',{p,k,Ts,r},...
         'solverParameters',{;
                    'compilerOptimization','-O0',...
                    'solverVerboseLevel',3 ...
                   });

%% Simulate system

% set parameter values
setParameter(mpc,'p',2);
setParameter(mpc,'k',1);
Ts=.1;
setParameter(mpc,'Ts',Ts);

% set process initial condition
t0=0;
x0=[.2;.2];
setInitialState(mpc,t0,x0);

fig=1;
figure(fig);clf;

mu0=1;
maxIter=20;
saveIter=false;

ref=@(t).6*sin(t);

t=t0;
figure(fig);clf;

% cold start
u_warm=.1*randn(nu,T);               

for i=1:40
    % set reference signal
    r=ref(t+(0:T-1)*Ts);
    setParameter(mpc,'r',r);
    
    % move warm start away from constraints
    u_warm=min(u_warm,.95);
    u_warm=max(u_warm,-.95);
    setSolverWarmStart(mpc,u_warm);
    
    [solution,J,x,u]=solve(mpc,mu0,maxIter,saveIter);
    
    fprintf('t=%g, J=%g computed in %g iterations & %g secs\n',t,J,solution.iter,solution.time);
    
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
    
    % apply 3 controls and get time and warm start for next iteration
    ufinal=zeros(nu,1);
    [t,u_warm]=applyControls(mpc,solution,3,ufinal); 
end

history=getHistory(mpc);

fig=fig+1;figure(fig);clf;
plot(history.t,history.x,'.-',...
     history.t,history.u,'.-',history.t,ref(history.t),'.-');grid on;
legend('x1','x2','u','r');
xlabel('t');

fig=fig+1;figure(fig);clf;
yyaxis left
plot(history.t,history.iter,'.-');grid on;
ylabel('# iterations');
yyaxis right
plot(history.t(2:end),1000*history.stime(2:end),'.-');grid on;
ylabel('solver time [ms]');
xlabel('t');



