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

nx=2;  % state size
nu=1;  % control size
T=30;  % forward horizon

Tvariable Ts [];
Tvariable x0 [nx,1];  % [x(t)]
Tvariable x1 [nx,T]; % [x(t+Ts), x(t+2*Ts), ..., x(t+T*Ts)];
Tvariable u  [nu,T];  % [u(t), u(t+Ts), ..., u(t+(T-1)*wTs) ] -- held constant between sample times
Tvariable p  [1,1];
Tvariable k  [1,1];

x=[x0,x1];  % whole state;

% parameters for state and input constraints
Tvariable uMax [];
Tvariable xMax [2,1];

% DC-motor-like transfer function
% [dot x1]=[0  1][x1]+[0]u
% [dot x2] [0  p][x2] [k]
%        y=x1     

% ODE: (x(k+1)-x(k)) /Ts == [ f (x(k+1),u(k)) + f (x(k),u(k)) ]/2
ode={x1-x(:,1:end-1)==Ts*( .5*[0,1;0,p]*(x1+x(:,1:end-1))+[0;k]*u)};

% Criterium

Tvariable r [1,T];   % reference
J=norm2(x1(1,:)-r)+norm2(u)/50;

% Warm start

x1warm=[x1(:,2:end),zeros(nx,1)];
uwarm=[u(:,2:end),zeros(nu,1)];
x1warm=min(x1warm,.9*repmat(xMax,[1,T]));
x1warm=max(x1warm,-.9*repmat(xMax,[1,T]));
uwarm=min(uwarm,.9*uMax);
uwarm=max(uwarm,-.9*uMax);

% inequality constraints
constraints={u<=uMax, u>=-uMax,...
             x1<=repmat(xMax,[1,T]),x1>=-repmat(xMax,[1,T])};

[classname,code]=cmex2optimizeCS(...
    'pedigreeClass','tmp_mpc',...
    'executeScript','asneeded',...
    ...
    'objective',J,...
    'optimizationVariables',{x1,u},...
    'constraints',[ode,constraints],...
    'parameters',{p,k,Ts,r,uMax,xMax,x0},...
    'outputExpressions',{J,x,u,x1warm,uwarm},...
    ...
    'compilerOptimization','-O0',...
    'solverVerboseLevel',2);
    
classname=getValue(classname);
code=getValue(code);

obj=feval(classname);


%% Simulate system

% set parameter values
p=2;k=1;
setP_p(obj,p);
setP_k(obj,k);
Ts=.1;
setP_Ts(obj,Ts);

uMax=1;
xMax=[.4;.3];
setP_uMax(obj,uMax);
setP_xMax(obj,xMax);

refmax=.35;
refomega=.5;

% set process initial condition
t0=0;
x0=[.2;.2];

% reference signals
ref=@(t)-refmax*sign(sin(refomega*t));

% cold start
u_cold=.1*randn(nu,T);
x1_cold=.1*randn(nx,T);
setV_u(obj,u_cold);
setV_x1(obj,x1_cold);

mu0=1;
maxIter=50;
saveIter=false;

nSteps=150;

closedloop.t=nan(nSteps+1,1);
closedloop.J=nan(nSteps+1,1);
closedloop.iter=nan(nSteps+1,1);
closedloop.stime=nan(nSteps+1,1);
closedloop.x=nan(nx,nSteps+1);
closedloop.u=nan(nu,nSteps+1);
closedloop.r=nan(1,nSteps+1);

closedloop.tt=nan(0,1);
closedloop.xx=nan(0,nx);

closedloop.t(1)=t0;
closedloop.x(:,1)=x0;

fig=1;
figure(fig);clf;
set(fig,'Name','MPC solutions');

for i=1:nSteps
    
    % current state
    t=closedloop.t(i);
    setP_x0(obj,closedloop.x(:,i));

    % set reference signal
    r=ref(t+(0:T-1)*Ts);
    setP_r(obj,r);
    
    [closedloop.status,closedloop.iter(i),closedloop.stime(i)]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    [closedloop.J(i),x,u,x1warm,uwarm]=getOutputs(obj);
    
    fprintf('t=%g, J=%g computed in %g iterations & %g secs\n',...
            t,closedloop.J(i),closedloop.iter(i),closedloop.stime(i));
    if closedloop.status
        warning('solver failed')
    end
    
        
    if true %mod(t/Ts,5)==0
        if 0
            plot(t:Ts:t+Ts*T,x,'.-',t:Ts:t+Ts*(T-1),u,'.-',t:Ts:t+Ts*(T-1),r,'.-');grid on;
            legend('x1','x2','u','r');
        else
            plot(t:Ts:t+Ts*(T-1),r,'k.-',t:Ts:t+Ts*T,x(1,:),'.-');grid on;
            legend('x1','r');
        end
        hold on;
    end
    
    % apply control
    closedloop.r(i)=r(1);r(1)=[];
    closedloop.u(:,i)=u(:,1);u(:,1)=[];
    [tout,yout]=ode23(@(t,x)[0,1;0,p]*x+[0;k]*closedloop.u(:,i),closedloop.t(i)+[0,Ts],closedloop.x(:,i));
    closedloop.tt(end+1:end+length(tout))=tout;
    closedloop.xx(end+1:end+length(tout),:)=yout;    
    closedloop.x(:,i+1)=yout(end,:)';
    closedloop.t(i+1)=tout(end,:)';
    
    % apply warm start for next iteration
    setV_u(obj,uwarm);
    setV_x1(obj,x1warm);
end

fig=fig+1;figure(fig);clf;
set(fig,'Name','Closed-loop');
plot(closedloop.tt,closedloop.xx,'.-',...
     closedloop.t,closedloop.u,'.-',...
     closedloop.t,ref(closedloop.t),'.-');grid on;
legend('x1','x2','u','r');
xlabel('t');

fig=fig+1;figure(fig);clf;
set(fig,'Name','Solver');
subplot(2,1,1);
plot(closedloop.t,closedloop.J,'.-');grid on;
ylabel('MPC cost');
subplot(2,1,2);
yyaxis left
plot(closedloop.t,closedloop.iter,'.-');grid on;
ylabel('# iterations');
yyaxis right
plot(closedloop.t,1000*closedloop.stime,'.-');grid on;
ylabel('solver time [ms]');
xlabel('t');



