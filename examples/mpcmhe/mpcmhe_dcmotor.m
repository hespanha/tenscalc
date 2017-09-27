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

nx=2;    % state size;
nu=1;    % control size;
nd=1;    % disturbance size;
ny=1;    % output sizes
T=30;    % forward horizon
L=8;    % backward horizon
delay=1;

% create symbolic optimization
Tvariable Ts [];
Tvariable x        [nx,L+T+1];    % [x(t-L*Ts),...,x(t),...,x(t+T*Ts)]
Tvariable y_past   [ny,L+1];      % [y(t-L*Ts),...,y(t)]
Tvariable u_past   [nu,L+delay];  % [u(t-L*Ts),...,u(t+(delay-1)*Ts)]
Tvariable u_future [nu,T-delay];  % [u(t+delay*Ts), ...,u(t+(T-1)*Ts)]
Tvariable d        [nd,L+T];      % [d(t-L*Ts,...,x(t),...,x(t+(T-1)*Ts)]

Tvariable p  [1,1];
Tvariable k  [1,1];
Tvariable r  [1,T-delay];         % [r(t+(delay+1)*Ts), ...,u(t+T*Ts)]
Tvariable umax [];

Tvariable cc [4];                 % optimization weights

% DC-motor-like transfer function
% [dot x1]=[0  1][x1]+[0]u
% [dot x2] [0  p][x2] [k]
%        y=x1     

dxFun=@(x,u,d,p,k,Ts,r,cc,umax)[0,1;0,p]*x+[0;k]*(u+d);
yFun=@(x,u,d,p,k,Ts,r,cc,umax)x(1,:);

u=[u_past,u_future];
JJ=[norm2(yFun(x(:,L+delay+2:L+T+1))-r);
    norm2(u_future);
    norm2(d);
    norm2(yFun(x(:,1:L+1))-y_past)];
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
               'controlConstraints',{u_future<=umax, u_future>=-umax},...
               'disturbanceConstraints',{d<=1, d>=-1},...
               'outputExpressions',{J,JJ,cc,Jcc,x,u_past,u_future,d},...
               'parameters',{p,k,Ts,r,cc,umax},...
               'solverParameters',{;
                    'alphaMin',1e-4,...
                    'solverVerboseLevel',3 ...
                   });

% set parameter values
setParameter(mpcmhe,'p',2);
setParameter(mpcmhe,'k',5);
setParameter(mpcmhe,'umax',3);
Ts=.1;
setParameter(mpcmhe,'Ts',Ts);

cc=[30;.01;-1e4;-1e4];
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
maxIter=30;
saveIter=false;

ref=@(t).6*sign(sin(3*t));

% cold start
x_warm=.1*randn(nx,1);
u_warm=.1*randn(nu,T-delay);
d_warm=0*.1*randn(nd,L+T);

%x_warm=state(:,1); % cheating

for i=1:50
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
    
    % apply 3 optimal controls/disturbances and get time,
    % (noiseless) measurements, and warm start for the next iteration
    u_final=zeros(nu,3);
    d_final=zeros(nd,3);
    noise=zeros(ny,3);
    [t,t_y_past,y_past,t_u_past,u_past,x_warm,u_warm,d_warm,t_state,state]=...
        applyControl(mpcmhe,u_past,y_past,solution,3,noise,u_final,d_final); 
    
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

    if true %mod(i,5)==0
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
        break
    end

end

history=getHistory(mpcmhe);
fig=fig+1;figure(fig);clf;
plot(history.t,history.x,'.-',...
     history.t,history.u,'.-',history.t,history.d,'.-',...
     history.t,ref(history.t),'.-',history.t,history.y,'.-');grid on;
legend('x1','x2','u','d','r','y');
xlabel('t');

fig=fig+1;figure(fig);clf;
yyaxis right
plot(history.t(2:end),1000*history.stime(2:end),'.-');grid on;
ylabel('solver time [ms]');
yyaxis left
plot(history.t,history.iter,'.-');grid on;
ylabel('# iterations');
xlabel('t');



