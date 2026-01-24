function constraint=tsODE(x,uZOH,uC,ts,fun,method)
% constraint=tsODE(x,uZOH,uC,ts,fun,method)
%
% computes a Tcalculus constraint that encodes the ODE
%     \dot x = f(x,uZOH,uC,t)
% with the understanding that the input 'uC' is assumes continuous
% and the input 'uZOH' is assumed piecewise constant (zero-order
% hold).
%
% Inputs:
%   x [n x N]  - values of the state x at sampling times
%                (one time per column)
%   uZOH [n x N]  - values of the input u at sampling times
%                (one time per column)
%                uZOH(t) is assumed constant until the next
%                sampling time
%                Can be an emptuy vector if no such input is needed.
%   uC [n x N]  - values of the input d at sampling times
%                (one time per column)
%                uC(t) is assumed continuous until the next
%                sampling time
%                Can be an emptuy vector if no such input is needed.
%   ts [N x 1] - vector of times
%    or
%   ts [1 x 1] - (constant) sample interval
%   fun        - handle to a function that computes
%                    f(x,uZOH,uC,t)
%   method     - ODE integration method:
%                  'forwardEuler'  - explicit forward Euler method
%                  'backwardEuler' - implicit forward Euler method
%                  'midPoint'      - implicit mid-point method
%                The mid-point method typically leads to the best
%                accuracy (for a similar integration step).
%                HOWEVER, it often leads to numerical problems
%                when optimizing for the control input.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

[n,N]=size(x);
switch method
    case 'forwardEuler'
        if length(ts)>1
            constraint=(x(:,2:end)==x(:,1:end-1)...
                +tprod(ts(2:end)-ts(1:end-1),2,...
                fun(x(:,1:end-1),uZOH(:,1:end-1),uC(:,1:end-1),ts(1:end-1)),[1,2]));
        else
            constraint=(x(:,2:end)==x(:,1:end-1)...
                +ts*fun(x(:,1:end-1),uZOH(:,1:end-1),uC(:,1:end-1),ts*(0:N-1)'));
        end
    case 'backwardEuler'
        if length(ts)>1
            constraint=(x(:,2:end)==x(:,1:end-1)...
                +tprod(ts(2:end)-ts(1:end-1),2,...
                fun(x(:,2:end),uZOH(:,1:end-1),uC(:,2:end),ts(2:end)),[1,2]));
        else
            constraint=(x(:,2:end)==x(:,1:end-1)...
                +ts*fun(x(:,2:end),uZOH(:,1:end-1),uC(:,2:end),ts*(1:N)'));
        end
    case 'midPoint'
        % with N data points, one can only have N-1 independent
        % constraints so one constraint should be removed
        lhs=tsDerivative(x,ts);
        lhs=lhs(:,1:end-1);
        if ~isempty(uZOH)
            if 1
                uu=.5*(uZOH(:,1:end-1)+uZOH(:,2:end));
                uu=[uu(:,1),uu];
            else
                fprintf('ATTENTION: tsODE(...,''midpoint'') is currently not doing zero-order hold for input ''uZOH''\n');
                uu=u;
            end
        else
            uu=uZOH;
        end
        if length(ts)>1
            rhs=fun(x,uZOH,uC,ts);
        else
            rhs=fun(x,uZOH,uC,ts*(1:N)');
        end
        rhs=rhs(:,1:end-1);
        constraint=(lhs==rhs);
    otherwise
        error('tsODE: method ''%s'' not implemented\n',method);
end

end

function test()

% dot x = x

clear all;delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');
clear all
T=1;
ts=1/3/2;
N=round(T/ts);

Tvariable x [1,N];
u=Tzeros([1,N]);
d=Tzeros([1,N]);
a=1;
fun=@(x,u,d,t)a*x;
ode=tsODE(x,u,d,ts,fun,'forwardEuler');
ode=tsODE(x,u,d,ts,fun,'backwardEuler');
ode=tsODE(x,u,d,ts,fun,'midPoint');
Tvariable y [];
Tvariable Hess_ (1+N+N)*[1,1];
class=class2optimizeCS('pedigreeClass','tmp1','executeScript','asneeded',...
    'objective',norm2(y-x(1,end))+norm2(x),...
    'optimizationVariables',{y,x},...
    'constraints',{ode,x(1,1)==1},...
    'outputExpressions',{y,x,Hess_},...
    'equalTolerance',1e-8,...
    'addEye2Hessian',0e-12,...
    'debugConvergence',false,...
    'solverVerboseLevel',3);
obj=feval(getValue(class));
setV_y(obj,1);
setV_x(obj,randn(1,N));
[status,iter,time]=solve(obj,1,int32(10),int32(false));
[y1,x1,Hess]=getOutputs(obj);
t=(0:N-1)*ts;
clf;
subplot(2,1,1)
plot(t,x1,'.-',t,exp(a*t),'.-');grid on
legend('x from solver','x exact');
subplot(2,1,2)
plot(t,x1-exp(a*t),'.-');grid on
legend('error');
end
