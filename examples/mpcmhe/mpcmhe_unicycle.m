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

clear all;
% remove previous solvers
% ATTENTION: for cmex version must not erase 
%          @tmp_mpcmhe_uni/tmp_mpcmhe_uni_WW.subscripts & @tmp_mpcmhe_uni/tmp_mpcmhe_uni_WW.values
%!rm -rf *toremove* tmp* % @tmp*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=8;    % backward horizon
T=12;   % forward horizon

% System dynamics:
% pursuer
% x1(k+1)    = x1(k) + v cos theta(k) 
% x2(k+1)    = x2(k) + v sin theta(k) 
% theta(k+1) = theta(k) + u(k)         u(k) in [-uMax,uMax]

% evader
% y1(k+1)    = y1(k) + d1(k)           d1(k), d2(k) in [-dMax,dMax]
% y2(k+1)    = y2(k) + d2(k)  

% controlled output
% z1 = y1 - x1
% z2 = y2 - x2

nX=5;
nU=1;
nD=2;
nY=4;

Tvariable v [];  % pursuer's velocity

dXfun=@(x,u,d,v)[x(1,:)+v*cos(x(3,:));
               x(2,:)+v*sin(x(3,:));
               x(3,:)+u;
               x(4:5,:)+d];
measuredOutput=@(x)x([1,2,4,5],:);

Tvariable x0    [nX,1];     % x(t-L)
Tvariable x1    [nX,L+T];   % x(t-L+1), ... , x(t+T)
Tvariable d     [nD,L+T];   % d(t-L), ... , d(t+T-1)
Tvariable uPast [nU,L];     % u(t-L), ... , u(t-1)
Tvariable uFut  [nU,T];     % u(t),   ... , u(t+T-1);
Tvariable yPast [nY,L+1];   % y(t-L), ..., y(t)

x=[x0,x1];                  % x(t-L), ... , x(t+T) 
uAll=[uPast,uFut];          % u(t-L), ... , u(t+T-1)
xFut=x1(:,L:end);           % x(t), ... , x(t+T)

dynamics = x1==dXfun(x(:,1:end-1),uAll,d,v);

n=measuredOutput(x(:,1:L+1))-yPast;  % noise


% parameters for state and input constraints
Tvariable uMax [];
Tvariable dMax [];

% Criterion

JJ=[norm2(xFut(1:2,:)-xFut(4:5,:)); % distance from pursuer to evader
    norm2(uFut);              % pursuer control
    -norm2(n);                % noise
    -norm2(d)];               % disturbance
Tvariable cc size(JJ);

J=cc*JJ;

% Warm start
tol=.9;
xWarm=[x(:,2:end),zeros(nX,1)];
uFutWarm=[uFut(:,2:end),zeros(nU,1)];
uFutWarm=min(uFutWarm,tol*uMax);
uFutWarm=max(uFutWarm,-tol*uMax);
dWarm=[d(:,2:end),zeros(nD,1)];
dWarm=min(dWarm,tol*dMax);
dWarm=max(dWarm,-tol*dMax);

if strcmp(computer,'PCWIN64');
    coder=@class2equilibriumLatentCS;
else
    coder=@cmex2equilibriumLatentCS;
end    
classname=coder(...
    ...%'pedigreeClass','tmp_mm_uni',...
    ...%'executeScript','yes',...
    'classname','tmp_mpcmhe_uni',...
    ...
    'P1objective',J,...
    'P2objective',-J,...
    'P1optimizationVariables',{uFut},...
    'P2optimizationVariables',{d,x0},...
    'latentVariables',{x1},...
    'P1constraints',{uFut<=uMax,uFut>=-uMax},...
    'P2constraints',{sum(d.*d,1)<=dMax*dMax},...
    'latentConstraints',{dynamics},...
    'parameters',{uMax,dMax,cc,v,uPast,yPast},...
    'outputExpressions',{J,JJ,uFut,d,x,full(uFutWarm),full(dWarm),full(xWarm)},...
    ...
    'muFactorAggressive',.5,...
    'muFactorConservative',.95,...
    ...;
    'allowSave',true,...
    ...%'smallerNewtonMatrix',true,...
    ...%'umfpack',true,...
    ...
    'compilerOptimization','-O0',...
    'solverVerboseLevel',2);
    
obj=feval(classname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v=.1;
setP_v(obj,v);
uMax=.5;
setP_uMax(obj,uMax);
dMax=.06;
setP_dMax(obj,dMax);

uWeight=8;
dWeight=100;
nWeight=1000;
cc=[1;uWeight;nWeight;dWeight];
setP_cc(obj,cc);

nSigma=.005;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% set initial condition

t0=0;
x0=[0;0;0;.5;.5];

mu0=1;
maxIter=200;
saveIter=-1; % on save on error

closedloop.t=t0;
closedloop.x=x0;
closedloop.y=zeros(nY,0);
closedloop.status=[];
closedloop.iter=[];
closedloop.stime=[];
closedloop.J=[];
closedloop.JJ=zeros(length(JJ),0);
closedloop.u=[];
closedloop.d=[];

figure(1);clf;
set(1,'Name','Trajectory');
figure(2);clf;
set(2,'Name','MPC-MHE solutions');

for k=1:200
    % get measurement
    n=nSigma*randn(nY,1);
    closedloop.y(:,end+1)=measuredOutput(closedloop.x(:,end))+n;
    
    if size(closedloop.u,2)<=L
        % open loop
        [closedloop.status(end+1,1),...
         closedloop.iter(end+1,1),...
         closedloop.stime(end+1,1)]=deal(nan,nan,nan);
        [closedloop.J(end+1,1),...
         closedloop.JJ(:,end+1),...
         uFut,hatd,hatx,uFutWarm,dWarm,xWarm]=deal(...
             nan,nan(size(JJ,1),1),...
             zeros(nU,T),zeros(nD,L+T),zeros(nX,L+T+1),...
             zeros(nU,T),zeros(nD,L+T),zeros(nX,L+T+1));        
    else
        % MPC-MHE optimization
        setP_uPast(obj,closedloop.u(:,end-L+1:end));
        setP_yPast(obj,closedloop.y(:,end-L:end));
        
        [closedloop.status(end+1,1),...
         closedloop.iter(end+1,1),...
         closedloop.stime(end+1,1)]=solve(obj,mu0,int32(maxIter),int32(saveIter));

        [closedloop.J(end+1,1),...
         closedloop.JJ(:,end+1),...
         uFut,hatd,hatx,uFutWarm,dWarm,xWarm]=getOutputs(obj);
    
        fprintf('t=%g, J=%g computed in %g iterations & %g ms\n',...
                closedloop.t(end),closedloop.J(end),closedloop.iter(end),1e3*closedloop.stime(end));

        if closedloop.status(end)==4
            error('failed to invert hessian, regenerate code to compile with saved values for hessian')
        end
        if closedloop.status(end)>0
            uFut=lastUfut(:,2:end);
            lastUfut=uFut;
            fprintf('solver failed at time %g, reusing result from time %g (paused)\n',...
                    closedloop.t(end),lastUfutT);
            pause;
        else
            lastUfut=uFut;
            lastUfutT=closedloop.t(end);
        end
    
    end
    
    
    % Apply control
     closedloop.u(:,end+1)=uFut(:,1);     
     % Apply disturbance 
     if closedloop.t(end)<55
         closedloop.d(:,end+1)=[.05;0];
     else
         closedloop.d(:,end+1)=hatd(:,L+1);
     end
     closedloop.x(:,end+1)=dXfun(closedloop.x(:,end),closedloop.u(:,end),closedloop.d(:,end),v);
     closedloop.t(end+1)=closedloop.t(end)+1;
     
     % apply warm start for next iteration     
     setV_uFut(obj,uFutWarm);
     setV_d(obj,dWarm);
     if size(closedloop.u,2)>L
         % fixe xWarm
         uPast=closedloop.u(:,end-L+1:end);
         uAll=[uPast,uFutWarm];
         for i=1:L+T
             xWarm(:,i+1)=dXfun(xWarm(:,i),uAll(:,i),dWarm(:,i),v);
         end
     end
     setV_x0(obj,xWarm(:,1));
     setV_x1(obj,xWarm(:,2:end));

     
     %subplot(4,4,mod(k,16)+1);
     figure(1);
     plot(closedloop.x(1,:),closedloop.x(2,:),'g.-',...
          closedloop.x(4,:),closedloop.x(5,:),'b.-',...
          hatx(1,L+1:end),hatx(2,L+1:end),'g:',...
          hatx(4,L+1:end),hatx(5,L+1:end),'b:');grid on
     legend('pursuer','evader','pursuer prediction','evader prediction');
     axis equal
     figure(2);
     subplot(3,1,1);
     plot(closedloop.t(1:end-1),closedloop.JJ,'*');grid on
     legend('x','u','n','d');
     subplot(3,1,2);
     plot(closedloop.t(1:end-1),closedloop.u,'g.-');grid on
     subplot(3,1,3);
     plot(closedloop.t(1:end-1),closedloop.d,'b.-');grid on
     
     title(sprintf('t=%d J=%.2f (status=%d, niter=%d, %.1fms)\n',...
                   closedloop.t(end),closedloop.J(end),...
                   closedloop.status(end),closedloop.iter(end),1e3*closedloop.stime(end)));
     drawnow
end

clear obj;

