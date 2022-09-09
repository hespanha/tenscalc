% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

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

% Criterion

Tvariable r [1,T];   % reference
J=norm2(xFut(1,:)-r)+norm2(uFut)/50;

% inequality constraints
constraints={uFut<=uMax, uFut>=-uMax,...
             xFut<=repmat(xMax,[1,T]),xFut>=-repmat(xMax,[1,T])};

% Warm start

xFutWarm=[xFut(:,2:end),zeros(nX,1)];
uFutWarm=[uFut(:,2:end),zeros(nU,1)];
xFutWarm=min(xFutWarm,.9*repmat(xMax,[1,T]));
xFutWarm=max(xFutWarm,-.9*repmat(xMax,[1,T]));
uFutWarm=min(uFutWarm,.9*uMax);
uFutWarm=max(uFutWarm,-.9*uMax);

[classname,code]=cmex2optimizeCS(...
    'pedigreeClass','tmp_mpc_dcmotor',...
    'executeScript','asneeded',...
    ...
    'objective',J,...
    'optimizationVariables',{xFut,uFut},...
    'constraints',{ode,constraints{:}},...
    'parameters',{p,k,Ts,r,uMax,xMax,xMeas},...
    'outputExpressions',{J,[xMeas,xFut],uFut,xFutWarm,uFutWarm},...
    ...
    'adjustAddEye2Hessian',false,...
    'compilerOptimization','-O0',...
    'solverVerboseLevel',2);
    
classname=getValue(classname);
code=getValue(code);
obj=feval(classname);


%% Set parameter values

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

%% Simulate system

% set process initial condition
t0=0;
xMeas=[.2;.2];

% reference signals
ref=@(t)-refmax*sign(sin(refomega*t));

% cold start
uFutCold=.1*randn(nU,T);
xFutCold=.1*randn(nX,T);
setV_uFut(obj,uFutCold);
setV_xFut(obj,xFutCold);

mu0=1;
maxIter=50;
saveIter=false;

nSteps=150;

closedloop.t=nan(nSteps+1,1);
closedloop.J=nan(nSteps+1,1);
closedloop.iter=nan(nSteps+1,1);
closedloop.stime=nan(nSteps+1,1);
closedloop.x=nan(nX,nSteps+1);
closedloop.u=nan(nU,nSteps+1);
closedloop.r=nan(1,nSteps+1);

closedloop.tt=nan(0,1);
closedloop.xx=nan(0,nX);

closedloop.t(1)=t0;
closedloop.x(:,1)=xMeas;

fig=1;
figure(fig);clf;
set(fig,'Name','MPC solutions');

for i=1:nSteps
    
    % current state
    t=closedloop.t(i);
    setP_xMeas(obj,closedloop.x(:,i));

    % set reference signal
    r=ref(t+(0:T-1)*Ts);
    setP_r(obj,r);
    
    [closedloop.status,closedloop.iter(i),closedloop.stime(i)]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    [closedloop.J(i),x,uFut,xFutWarm,uFutWarm]=getOutputs(obj);
    
    fprintf('t=%g, J=%g computed in %g iterations & %g secs\n',...
            t,closedloop.J(i),closedloop.iter(i),closedloop.stime(i));
    if closedloop.status
        warning('solver failed')
    end
    
        
    if true %mod(t/Ts,5)==0
        if 0
            plot(t:Ts:t+Ts*T,x,'.-',t:Ts:t+Ts*(T-1),u,'.-',t:Ts:t+Ts*(T-1),r,'.-');grid on;
            legend('xFut','x2','u','r');
        else
            plot(t+(0:T-1)*Ts,r,'k.-',t:Ts:t+Ts*T,x(1,:),'.-');grid on;
            legend('xFut','r');
        end
        hold on;
    end
    
    % apply control
    closedloop.r(i)=r(1);
    closedloop.u(:,i)=uFut(:,1);
    [tout,yout]=ode23(@(t,x)dxFun(x,closedloop.u(:,i),p,k),closedloop.t(i)+[0,Ts],closedloop.x(:,i));
    closedloop.tt(end+1:end+length(tout))=tout;
    closedloop.xx(end+1:end+length(tout),:)=yout;    
    closedloop.x(:,i+1)=yout(end,:)';
    closedloop.t(i+1)=tout(end,:)';
    
    % apply warm start for next iteration
    setV_uFut(obj,uFutWarm);
    setV_xFut(obj,xFutWarm);
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



