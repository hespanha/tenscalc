% MPC state-feedback control for a quadcopter with continuous-time model:
%
%    ddot p = -b_{drag} dot p + g + R_{ib} u_thrust m_thrust
%
% where
%    p        = center of mass position
%    b_{drag} = drag coefficient
%               (assumed a scalar constant)
%    u_thrust = unit-vector along the thrust direction in body frame
%               (assumed constant and equal to [0;0;1] - all propellers pointing down)
%    m_thrust = magnitude
%    g        = gravity vector (in inertial frame)
%               (assumed constant)
%    R_{ib}   = rotation matrix from body to intertial frame
%
% All vectors assume a North-East-Down coordinate system (so z<0 means xabove ground!)
%
% We assume that a low level controller regulates 4DOFs:
%
%    1) the total thrust magnitude m_thrust (scalar),
%
%       subject to the following constraints
%           m_thrust >= min_thrust   [> 0, propellers only rotate in one direction]
%           m_thrust <= max_thrust
%
%    2) the orientation R_{ib} of the quadcopter (3-vector)
%
% To faciliate the optimization, we optimize with respect to the variable
%
%      u = R_{ib} u_thrust m_thrust
%
% which leads to the linear dynamics
%
%    ddot p = -b_{drag} dot p + g + u
%
% from which we recover the desired R_{ib} and m_thrust as follows:
%
%   a) u = R_{ib} [0;0;1] m_thrust
%           => m_thrust =  \|u\|       -- should not become zero to avoid singularity
%           => R_ib(:,3) = u/\|u\|
%
%   b)  y-axis of body selected to have no vertical component
%       (i.e., pointing sideway horizontally)
%           => R_ib(:,2) = R_ib(:,3) x [0,0,1]
%
%   c) R is a rotation matrix
%           => R_ib(:,1) = R_ib(:,2) x R_ib(:,3)
%
% The goal of MPC is to minimize the distance to a desired final
% position/velocity, penalizing total thrust
%
%         J = \|p(T)-pfinal\|^2 
%                     + lambda_v \|dot p(T)-vfinal\|^2 
%                     + lambda_thrust \int_0^T m_thrust(t) dt
%
% subject to the constraint mentioned on 1) above as well as a minimum
% altitude constraint:
%
%        p(3,:) <= -min_altitude    (recall that z-xis points down)
%
% The solver can be generated once and then called multiple
% times. Creating a solver typically takes much longer than solving
% the optimization by creating the solver.
%
% This script gives you the option of skipping creating the solver,
% but you can only skip generating a solver if you have done it once.

% This file is part of Tencalc.
%
% Copyright (C) 2012-22 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all;
% remove previous solvers
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

help mpc_quadcopter

createSolvers=input('Do you want to create the solvers (Y/n)? ','s');
createSolvers=~isequal(createSolvers,'n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=100;                        % horizon length in time-steps

if createSolvers

    %%%%%%%%%%%%%%%%%%%
    %% Generate solver
    %%%%%%%%%%%%%%%%%%%

    Tcalculus.clear();

    Tvariable Ts [];              % sampling time
    Tvariable p [3,T];            % position:              [p(1) p(2) ... p(T)]
    Tvariable u [3,T];            % optimization variable: [u(1) u(2) ... u(T)]

    Tvariable pinit [3,1];        % initial position:       p(1)
    Tvariable vinit [3,1];        % initial velocity:       v(1)

    Tvariable b_drag [];          % drag coefficient
    g=Tconstant([0;0;9.8],[3,1]); % gravity vector (in inertial frame)

    Tvariable min_thrust [];      % minimum thrust
    Tvariable max_thrust [];      % maximum thrust
    Tvariable max_omega  [];      % maximum angular velocity
    Tvariable min_altitude [];    % minimum altitude

    v=tsDerivative(p,Ts);         % velocity
    a=tsDerivative2(p,Ts);        % acceleration

    % System dynamics
    dynamics={
        a==-b_drag*v+repmat(g,[1,T])+u;
        p(:,1)==pinit;            % initial position
        v(:,1)==vinit;            % initial velocity
             };

    m_thrust = sqrt(sum(u.^2,1));      % m_thurst = \|u\| from a) above -- sqrt is dangerous because not diff
    R_ib3    = u./repmat(reshape(m_thrust,[1,T]),[3,1]); % reshape needed by repmat
    R_ib2    = tsCross(R_ib3,[Tzeros(2,T);Tones(1,T)]);
    R_ib1    = tsCross(R_ib2,R_ib3);

    % Constraints
    Tvariable positive2 [T];

    constraints={
        m_thrust >= min_thrust;
        ...%m_thrust <= max_thrust;
        max_thrust-m_thrust==positive2;positive2>0;  % this will work even if initial condition violates
        p(3,:) < -min_altitude;
                };

    % Criterion
    Tvariable pdesired [3,1];     % desired final position: p_{desired}(T)
    Tvariable vdesired [3,1];     % desired final velocity: v_{desired}(T)

    Tvariable lambda_v [];            % weigth for velocity error
    Tvariable lambda_thrust [];       % weigth for thrust

    Jp2=norm2(p(:,T)-pdesired);       % position error^2
    Jv2=norm2(v(:,T)-vdesired);       % velocity error^2
    Jp2=tsIntegral(sum((p-repmat(pdesired,[1,T])).^2,1),Ts);    % integral of position error^2
    Jv2=tsIntegral(sum((v-repmat(vdesired,[1,T])).^2,1),Ts);    % integral of velocity error^2
    Jthrust=tsIntegral(m_thrust,Ts);  % integral of m_thrust

    J = Jp2+lambda_v*Jv2+lambda_thrust*Jthrust;

    % Solver output
    outputExpressions=struct('J',J,...
                             'Jp2',Jp2,...
                             'Jv2',Jv2,...
                             'Jthrust',Jthrust,...
                             't',(0:T-1)*Ts,...
                             'u',u,...
                             'm_thrust',m_thrust,...
                             'p',p,...
                             'v',v,...
                             'a',a,...
                             'R_ib1',R_ib1,...
                             'R_ib2',R_ib2,...
                             'R_ib3',R_ib3,...
                             'positive2',positive2);

    classname=cmex2optimizeCS('classname','tmp_quadcopter',...
                              'objective',J,...
                              'optimizationVariables',{p,u,positive2},...
                              'constraints',[dynamics;constraints],...
                              'outputExpressions',outputExpressions,...
                              'parameters',{...
                                  Ts;...
                                  pinit;vinit;...
                                  pdesired;vdesired;...
                                  b_drag;...
                                  min_thrust;max_thrust;min_altitude;...
                                  lambda_v;lambda_thrust},...
                              'addEye2Hessian',true,...
                              'adjustAddEye2Hessian',true,...
                              'scaleInequalities',true,...
                              'solverVerboseLevel',2);  % use 4 to see solver iterations

end

%%%%%%%%%%%%%%%
%% Use solver
%%%%%%%%%%%%%%%

%% Create object
obj=tmp_quadcopter();
mu0=1e-3;       % low value will speed-up convergence (up to a point where it causes numerical issues)
maxIter=1000;
saveIter=-1;

s = RandStream('mt19937ar','Seed',0); % fix seed for reproducibility
RandStream.setGlobalStream(s);

%% Set parameters
Ts=.02;               % 20 ms
b_drag=.1;
min_altitude=-.1;     % must be a little lower than initial altitude
min_thrust=5;         % must be positive to avoid singularities in computing R_{ib}
max_thrust=20;
lambda_v=0.05;        % low weight just to regularize solution
lambda_thrust=.05;    % low weight just to regularize solution

setP_Ts(obj,Ts);
setP_b_drag(obj,b_drag);
setP_min_altitude(obj,min_altitude);
setP_min_thrust(obj,min_thrust);
setP_max_thrust(obj,max_thrust);
setP_lambda_v(obj,lambda_v);
setP_lambda_thrust(obj,lambda_thrust);


pinit=[0;0;0];            % initial position
vinit=[0;0;0];            % initial velocity
pdesired=[0;10;-5]/2;     % desired final position
vdesired=[0;0;0];         % desired final velocity

setP_pdesired(obj,pdesired);
setP_vdesired(obj,vdesired);

clear history;   % structure where system's state, control, etc. will be stored
nSteps=150;
for k=1:nSteps

    %% Set parameters specific to this optimization
    setP_pinit(obj,pinit);  
    setP_vinit(obj,vinit);

    %% Initialize primal variables using a spline that respects initial position and velocity
    pfinal=pdesired;                  % try to go all the way (may not respect constraints)
    pfinal=pinit+(pdesired-pinit)/10; % try to go only part of way (to respect constraints)
    t=(1:T)*Ts;
    cs=spline([t(1),t(2),t(3),t(end)],...
              [pinit,pinit+Ts*vinit,pinit+2*Ts*vinit,pfinal]);
    p=ppval(cs,t);
    setV_p(obj,p);

    v=tsDerivative(p,Ts);         % velocity
    a=tsDerivative2(p,Ts);        % acceleration
    g=[0;0;9.8];
    u=a+b_drag.*v-g; % ddot p = -b_{drag} dot p + g + u
    setV_u(obj,u);
    setV_positive2(obj,ones(T,1));

    %% Solve optimization
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    out=getOutputs(obj);
    k
    disp(out)

    %% Save applied sate, control, etc. in `history` structure
    if k==1
        history.t(k)=0;
    else
        history.t(k)=history.t(k-1)+Ts;
    end
    history.p(1:3,k)=pinit;
    history.v(1:3,k)=vinit;
    history.u(1:3,k)=u(:,1);
    history.R_ib1(1:3,k)=out.R_ib1(:,1);
    history.R_ib2(1:3,k)=out.R_ib2(:,1);
    history.R_ib3(1:3,k)=out.R_ib3(:,1);
    history.m_thrust(k)=out.m_thrust(1);
    history.J(k)=out.J;
    history.Jp2(k)=out.Jp2;
    history.Jv2(k)=out.Jv2;
    history.iter(k)=iter;
    history.time(k)=time;

    %% Apply (constant) control to update initial position and velocity
    % ddot p = -b_{drag} dot p + g + R_{ib} u_thrust m_thrust
    gRu=g+out.R_ib3(:,1)*out.m_thrust(1);  % 1st time
    [tout,yout]=ode23(@(t,pv)[pv(4:6,:);-b_drag*pv(4:6,:)+gRu],[0,Ts],[pinit;vinit]);
    pinit=yout(end,1:3)';
    vinit=yout(end,4:6)';

    if mod(k,10)==1 || k==nSteps || status
        %% Plot MPC solution
        fig=clearFigure('figureNumber',1,'figureName','Latest MPC solution');
        plotData(out,Ts)
        %% Plot actual control
        fig=clearFigure('figureNumber',fig+1,'figureName','Actual control');
        plotData(history,Ts)
        %% Plot solver statistics
        fig=clearFigure('figureNumber',fig+1,'figureName','Solver statistics');
        subplot(2,1,1);
        plot(history.t,history.J,'.-',...
             history.t,history.Jp2,'.-',...
             history.t,history.Jv2,'.-')
        legend('J','Jp2','Jv2','location','best');grid on;
        subplot(2,1,2);
        yyaxis left;plot(history.t,history.iter,'.-');ylabel('# solver iter');grid on;
        yyaxis right;plot(history.t,1000*history.time,'.-');ylabel('solver time [ms]');grid on;
        drawnow;
    end

    if status
        warning('optimization failed\n');
    end

end

% Delete object
clear obj

function plotData(data,Ts)
    % 1: x-y    2: px-t  3: vx-t 4: wx-t
    % 5: u-t    6: py-t  7: vx-t 8: wx-t
    % 9: th-t  10: py-t 11: vx-t 12: wx-t
    subplot(3,4,1);
    plot(data.p(1,:),data.p(2,:),'.-');
    xlabel('x'),ylabel('y');
    axis square;grid on;
    subplot(3,4,5);
    plot(data.t,data.u,'.-');grid on
    legend('u_x','u_y','u_z','location','best');
    subplot(3,4,9);
    plot(data.t,data.m_thrust,'.-');ylabel('m_{thrust}');grid on
    subplot(3,4,2);
    plot(data.t,data.p(1,:),'.-');ylabel('p_x');grid on
    subplot(3,4,6);
    plot(data.t,data.p(2,:),'.-');ylabel('p_y');grid on
    subplot(3,4,10);
    plot(data.t,data.p(3,:),'.-');ylabel('p_z');grid on
    subplot(3,4,3);
    plot(data.t,data.v(1,:),'.-');ylabel('v_x');grid on
    subplot(3,4,7);
    plot(data.t,data.v(2,:),'.-');ylabel('v_y');grid on
    subplot(3,4,11);
    plot(data.t,data.v(3,:),'.-');ylabel('v_z');grid on
    subplot(3,4,4);
    plot(data.t,data.R_ib3(1,:),'.-');ylabel('R_{ib}(x,down)');grid on
    subplot(3,4,8);
    plot(data.t,data.R_ib3(2,:),'.-');ylabel('R_{ib}(y,down)');grid on
    subplot(3,4,12);
    plot(data.t,data.R_ib3(3,:),'.-');ylabel('R_{ib}(z,down)');grid on


end
