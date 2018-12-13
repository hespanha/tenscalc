function constraint=tsODE(x,u,d,ts,fun,method)
% constraint=tsODE(x,u,d,ts,fun,method)
%
% computes a Tcalculus constraint that encodes the ODE
%     \dot x = f(x,u,d,t)
% with the understanding that the input 'd' is assumes continuous
% and the input 'u' is assumed piecewise constant (zero-oerder
% hold).
%
% Inputs:
%   x [n x N]  - values of the state x at sampling times
%                (one time per column)
%   u [n x N]  - values of the input u at sampling times
%                (one time per column)
%                u(t) is assumed constant until the next
%                sampling time
%                Can be an emptuy vector if no such input is needed.
%   d [n x N]  - values of the input d at sampling times
%                (one time per column)
%                d(t) is assumed continuous until the next
%                sampling time
%                Can be an emptuy vector if no such input is needed.
%   ts [N x 1] - vector of times 
%    or
%   ts [1 x 1] - (constant) sample interval
%   fun        - handle to a function that computes 
%                    f(x,u,d,t)
%   method     - ODE integration method: 
%                  'forwardEuler'  - explicit forward Euler method
%                  'backwardEuler' - implicit forward Euler method
%                  'midPoint'      - implicit mid-point method
%                The mid-point method typically leads to the best
%                accuracy (for a similar integration step).
%
%                The mid-point methid currently ignores the
%                difference between the inputs u & d (and treats
%                them both as continuous) because it appears to
%                lead to large errors for a zoh control input.
[n,N]=size(x);
switch method
  case 'forwardEuler'
    if length(ts)>1 
        constraint=(x(:,2:end)==x(:,1:end-1)...
                    +tprod(ts(2:end)-ts(1:end-1),2,...
                           fun(x(:,1:end-1),u(:,1:end-1),d(:,1:end-1),ts(1:end-1)),[1,2]));
    else
        constraint=(x(:,2:end)==x(:,1:end-1)...
                    +ts*fun(x(:,1:end-1),u(:,1:end-1),d(:,1:end-1),ts*(0:N-1)'));
    end
  case 'backwardEuler'
    if length(ts)>1 
        constraint=(x(:,2:end)==x(:,1:end-1)...
                    +tprod(ts(2:end)-ts(1:end-1),2,...
                           fun(x(:,2:end),u(:,1:end-1),d(:,2:end),ts(2:end)),[1,2]));
    else
        constraint=(x(:,2:end)==x(:,1:end-1)...
                    +ts*fun(x(:,2:end),u(:,1:end-1),d(:,2:end),ts*(1:N)'));
    end
  case 'midPoint'
    % with N data points, one can only have N-1 independent
    % constraints so one constraint should be removed
    lhs=tsDerivative(x,ts);
    lhs=lhs(:,1:end-1);
    if ~isempty(u)
        if 0
            uu=.5*(u(:,1:end-1)+u(:,2:end));
            uu=[uu(:,1),uu];
        else
            fprintf('ATTENTION: tsODE(...,''midpoint'') is currently not doing zero-order hold for input ''u''\n');
            uu=u;
        end
    else
        uu=u;
    end
    if length(ts)>1 
        rhs=fun(x,uu,d,ts);
    else
        rhs=fun(x,uu,d,ts*(1:N)');
    end
    rhs=rhs(:,1:end-1);
    constraint=(lhs==rhs);
  otherwise
    error('tsODE: method ''%s'' not implemented\n',method);
end

return

%% Test

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
