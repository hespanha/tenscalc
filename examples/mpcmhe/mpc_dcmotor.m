clear all;
clear global;
!rm -rf toremove.m tmp* @tmp*

% Create dc-motor-like transfer function
% [dot x1]=[0  1][x1]+[0]u
% [dot x2] [0  p][x2] [k]
%        y=x1     

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
    % move warm start away from constraints
    u_warm=min(u_warm,.95);
    u_warm=max(u_warm,-.95);
    setSolverWarmStart(mpc,u_warm);
    
    % set reference signal
    r=ref(t+(0:T-1)*Ts);
    setParameter(mpc,'r',r);
    
    [solution,J,x,u]=solve(mpc,mu0,maxIter,saveIter);
    
    fprintf('J=%g computed in %g iterations & %g secs\n',J,solution.iter,solution.time);
    
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
figure(fig);clf;fig=fig+1;
plot(history.t,history.x,'.-',history.t,history.u,'.-',history.t,ref(history.t),'.-');grid on;
legend('x1','x2','u','r');
xlabel('t');

figure(fig);clf;fig=fig+1;
plot(history.t,history.iter,'.-',history.t(2:end),1000*history.stime(2:end),'.-');grid on;
legend('# iterations','solver time [ms]');
xlabel('t');



