% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all;
% remove previous solvers
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% LASSO-like regression -- uses TClasso tool to generate the optimization.
%
% min || y - theta0 ones(size(y)) - H \theta || + lambda \|\theta\|_1
%
% w.r.t. theta0 in R, theta in R^n

if 0
    n=15;
    m=1000;
elseif 0
    n=36;
    m=2200;
else
    load RegressionTestData/regErrData.mat
    [m,n]=size(H);
    n=n-1;
    loadH=H(:,2:end);
    loady=y;
end


Tvariable lambda [];

Tvariable theta0   [];
Tvariable theta    [n];
Tvariable absTheta [n];
Tvariable v    [];

Tvariable y2   [];
Tvariable sumY []; 
Tvariable yH   [n];
Tvariable sumH [n];
Tvariable HH   [n,n];

v2=y2-2*sumY*theta0-2*yH*theta+2*theta0*sumH*theta+theta0*theta0*m+tprod(theta,[-1],HH,[-1,-2],theta,[-2]);

J =v+lambda*sum(absTheta,1);

% dual variables
Tvariable lambda1_ [n];
Tvariable lambda2_ [n];
Tvariable lambda3_ [];
Tvariable lambda4_ [];
Tvariable Hess_ [1,1]*(2*n+2+2*n+2);

preMultiplied=true;
useSqrt=true;      % true is slightly slower, but can be more robust
smoothSqrt=false;  % true generally leads to slower convergence

classname=TClasso('classname','tmp_lasso',...
                  'dimension',n,...
                  'nTrainingPoints',m,...
                  'addConstant',true,...
                  'preMultiplied',preMultiplied,....
                  'smoothSqrt',smoothSqrt,....
                  'useSqrt',useSqrt,...
                  'solverType','C',...
                  'solverVerboseLevel',2);
help(classname)
% Create object
obj=feval(classname);

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

scale_th=1;
if useSqrt
    lambda=10/scale_th;
else
    lambda=10/scale_th;
end

setP_l1weight(obj,lambda);

clear status iter time

for i=1:100
    % create problem instance
    if 0 && exist('loadH','var')  %% very bad example
        thisy=loady;
        thisH=loadH;
    else
        thisTheta=randn(n,1);
        thisTheta(rand(n,1)<.5)=0;
        thisTheta0=randn(1,1);
        thisH=randn(m,n);
        thisy=thisTheta0+thisH*thisTheta+.2*randn(m,1); 
    end
    if 1
        % rescale theta
        thisH=thisH/scale_th;
    end
    if 1
        % rescale residuals
        s=mean(abs(thisy));
        thisy=thisy/s;
        thisH=thisH/s;
    end
    if 1
        % check conditining
        fprintf('svd(H''H) = %s\n',mat2str(svd(thisH'*thisH),2));
    end
    
    t0=clock();

    % initial value for solver
    if 1
        theta_0=.01*randn(n,1);
        theta0_0=mean(thisy-thisH*theta_0);
    else
        % cheat
        theta_0=thisTheta;
        theta0_0=thisTheta0;
    end


    if useSqrt
        fprintf('Initial cost= %g+%g x %g = %g\n',...
                norm(theta0_0+thisH*theta_0-thisy,2),...
                lambda,...
                norm(theta_0,1),...
                norm(theta0_0+thisH*theta_0-thisy,2)+lambda*norm(theta_0,1));
    else
        fprintf('Initial cost= %g+%g x %g = %g\n',...
                norm(theta0_0+thisH*theta_0-thisy,2)^2,...
                lambda,...
                norm(theta_0,1),...
                norm(theta0_0+thisH*theta_0-thisy,2)^2+lambda*norm(theta_0,1));
    end
    
    absTheta_0=abs(theta_0)+1;
    e2_0=norm(thisy-theta0_0-thisH*theta_0,2)^2;
    e_0=sqrt(e2_0+1);

    if preMultiplied
        % Set parameters
        setP_XX(obj,thisH'*thisH);
        setP_y2(obj,thisy'*thisy);
        setP_yX(obj,thisH'*thisy);
        setP_sumY(obj,sum(thisy));
        setP_sumX(obj,sum(thisH,1)');
    else
        setP_y(obj,thisy);
        setP_X(obj,thisH);
    end
    % Initialize primal variables
    setV_W(obj,theta_0);
    setV_absW(obj,absTheta_0);
    setV_c(obj,theta0_0);
    if useSqrt && smoothSqrt
        setV_sqrte(obj,e_0);
    end
    
    % Solve optimization
    mu0=1;
    maxIter=30;
    saveIter=-1;
    [status(i),iter(i),time(i)]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    [theta,theta0,J]=getOutputs(obj);
    %theta,in2,in3,lambda1_,lambda2_,lambda3_,lambda4_
    if useSqrt
        fprintf('TC  J=%9.4f (reported) %9.4f (actual) (%3d iterations, %8.1f ms)\n',...
                J,norm(thisy-theta0-thisH*theta,2)+lambda*norm(theta,1),...
                iter(i),etime(clock(),t0)*1e3); 
    else
        fprintf('TC  J=%9.4f (reported) %9.4f (actual) (%3d iterations, %8.1f ms)\n',...
                J,norm(thisy-theta0-thisH*theta,2)^2+lambda*norm(theta,1),...
                iter(i),etime(clock(),t0)*1e3); 
    end
    disp([theta0,theta']);

    %regressionCheckHessian(lambda,thisy,thisH,J,theta0,theta,absTheta,v,lambda1_,lambda2_,lambda3_,lambda4_,Hess_);

    if 0
        %% Compare with CVX
        t0=clock();
        cvx_quiet(true)
        cvx_begin
        variable th(n)
        variable th0;
        if useSqrt
            J=norm(thisy-th0-thisH*th,2)+lambda*norm(th,1);
        else
            J=sum_square(thisy-th0-thisH*th)+lambda*norm(th,1);
        end
            minimize J
            cvx_end;
            if useSqrt
                fprintf('CVX J=%9.4f (%8.1f ms)\n',...
                        norm(thisy-th0-thisH*th,2)+lambda*norm(th,1),etime(clock(),t0)*1e3);
            else
                fprintf('CVX J=%9.4f (%8.1f ms)\n',...
                        norm(thisy-th0-thisH*th,2)^2+lambda*norm(th,1),etime(clock(),t0)*1e3);
            end
                disp([th0,th']);
    end

    if status(i)
        error()
    end
end

fprintf('iter = mean=%.2f [%d,%d], time = median=%.1f us [%.1f,%.1f]\n',...
        mean(iter),min(iter),max(iter),...
        1e6*median(time),1e6*min(time),1e6*max(time));
