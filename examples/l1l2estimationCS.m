
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%optimizeCS=@class2optimizeCS; % matlab
optimizeCS=@cmex2optimizeCS;  % c
minInstructions4loop=100;
N=200

codeType='C';
%compilerOptimization='-O0';
compilerOptimization='-O1';

allowSave=false;
debugConvergence=false;

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

%% l2 noise
noise=1;
t=cumsum(ceil(1*rand(N,1)));
truePosition=5*sin(abs(t-100)/5);
thisMeasurement=round(truePosition+noise*randn(N,1));

%% create outliers
p_outlier=.1;
k_outlier=find(rand(N,1)<p_outlier);
thisMeasurement(k_outlier)=round(10*randn(length(k_outlier),1));

thisDt1=1./(t(2:end)-t(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L2 optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% declare (symbolic) variables
Tvariable measurement N;          % measurements
Tvariable dt1 N-1;                % 1/time differences 
Tvariable weight2acceleration []; % cost-weight for acceleration's l2-norm
Tvariable position N;             % positions (to be optimized)

velocity=(position(2:end)-position(1:end-1)).*dt1;              % velocity
acceleration=(velocity(2:end)-velocity(1:end-1)).*dt1(1:end-1); % total acceleration

%% construct optimization criteria
noise2=measurement-position;                                    % l2 measurement noise
acceleration2=acceleration;                                     % l2 acceleration
J=norm2(noise2)+weight2acceleration*norm2(acceleration2);

%% parameter values
thisWeight2acceleration=10;

if 1
    [classname,code]=optimizeCS(...
        'classname','tmpC_minl2',...
        'objective',J,...
        'optimizationVariables',{position},...
        'outputExpressions',{J,position},...
        'parameters',{
            measurement;dt1;
            weight2acceleration},...
        'skipAffine',false,...
        'codeType',codeType,...
        'debugConvergence',debugConvergence,...
        'minInstructions4loop',minInstructions4loop,...
        'compilerOptimization',compilerOptimization,...
       ...%'callType','client-server','serverProgramName','tmpC_minl2_server','serverAddress','localhost','port',1968,...
        ...%'callType','client-server','serverComputer','glnxa64','compileStandalones',false,'serverProgramName','tmpC_minl2_server','serverAddress','tamina.ece.ucsb.edu','port',1968,...
        'allowSave',allowSave,...
        'profiling',false,...
        'solverVerboseLevel',2);

    obj=feval(classname);

    if 0
        if system('killall tmpC_minl2_server')==0
            pause(3);
        end
        if system('tmpC_minl2_server &')==0
            pause(3);
        end
    end

    if 0
        upload(obj)
    end
    
    %% initial condition
    position_0=zeros(position.size,1);

    % set parameters
    setP_measurement(obj,thisMeasurement);
    setP_dt1(obj,thisDt1);
    setP_weight2acceleration(obj,thisWeight2acceleration);
    % initialize primal
    setV_position(obj,position_0);
    
    % call solver
    mu=10*N;
    maxIter=50;
    saveIter=1;
    for i=1:4
        [status,iter,time]=solve(obj,mu,int32(maxIter),int32(saveIter));
        [J_star,position_l2_star]=getOutputs(obj);
        saveIter=-1;
    end
    clear obj;

    %% Plots
    figure(1),clf
    subplot(6,1,1)
    plot(t,truePosition,'-.',t,thisMeasurement,'-+',t,position_l2_star,'-o',...
         t(k_outlier),zeros(length(k_outlier),1),'*')
    grid on
    legend('true position','measurements','l2 estimates','outliers');
    title('l2 estimate')
    subplot(6,1,2)
    plot(t,thisMeasurement-position_l2_star,'-+')
    grid on
    legend('l2 estimated noise')
    subplot(6,1,3)
    plot(t,position_l2_star-truePosition,'-+')
    grid on
    legend('l2 estimation error')
    drawnow
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L1-L2 optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% declare (symbolic) variables
Tvariable weight1acceleration []; % cost-weight for acceleration's l1-norm
Tvariable weight1noise [];        % cost-weight for noise's l1-norm
Tvariable noise1   N;             % l1 measurement noise (to be optimized)
Tvariable acceleration1 N-2;      % l1 acceleration (to be optimized)
Tvariable noise1abs   N;          % noise1 absolute values (to be optimized)
Tvariable acceleration1abs N-2;   % acceleration1 absolute values (to be optimized)

%% construct optimization criteria
noise2=measurement-position-noise1;                             % l2 measurement noise
acceleration2=acceleration-acceleration1;                       % l2 acceleration
J=norm2(noise2) ...
  +weight2acceleration*norm2(acceleration2)...
  +weight1noise*sum(noise1abs,1) ...
  +weight1acceleration*sum(acceleration1abs,1);

%profile on
[classname,code]=optimizeCS(...
    'classname','tmpC_minl1l2',...
    'objective',J,...
    'optimizationVariables',{...
        position,noise1,acceleration1,noise1abs,acceleration1abs},...
    'constraints',{
        noise1<=noise1abs;
        noise1>=-noise1abs;
        acceleration1<=acceleration1abs;
        acceleration1>=-acceleration1abs;
                  },...
    'outputExpressions',{J,position,noise1,acceleration1,noise1abs,acceleration1abs},...
    'parameters',{
        measurement;dt1;
        weight1acceleration;weight2acceleration;weight1noise},...
    'skipAffine',false,...
    'profiling',false,...
    'codeType',codeType,...
    'debugConvergence',debugConvergence,...
    'minInstructions4loop',minInstructions4loop,...
    'compilerOptimization',compilerOptimization,...
...%'callType','client-server','serverProgramName','tmpC_minl1l2_server','serverAddress','localhost','port',1968,...
...%'callType','client-server','serverComputer','glnxa64','compileStandalones',false,'serverProgramName','tmpC_minl1l2_server','serverAddress','tamina.ece.ucsb.edu','port',1968,...
    'allowSave',allowSave,...
    'solverVerboseLevel',2);
%profile viewer

obj=feval(classname);

if 0
    if system('killall tmpC_minl1l2_server')==0
        pause(3);
    end
    if system('tmpC_minl1l2_server &')==0
        pause(3);
    end
end

if 0
    upload(obj)
end

%% parameter values
thisWeight1acceleration=1;
thisWeight1noise=1;
%% initial condition
position_0=rand(position.size,1);
noise1_0=rand(noise1.size,1);
acceleration1_0=rand(acceleration1.size,1);
noise1abs_0=10*ones(noise1abs.size,1);
acceleration1abs_0=10*ones(acceleration1abs.size,1);
%% solve


% set parameters
setP_measurement(obj,thisMeasurement);
setP_dt1(obj,thisDt1);
setP_weight1acceleration(obj,thisWeight1acceleration);
setP_weight2acceleration(obj,thisWeight2acceleration);
setP_weight1noise(obj,thisWeight1noise);
% initialize primal
setV_position(obj,position_0);
setV_noise1(obj,noise1_0);
setV_acceleration1(obj,acceleration1_0);
setV_noise1abs(obj,noise1abs_0);
setV_acceleration1abs(obj,acceleration1abs_0);
% call solver
mu=10*N;
maxIter=50;
saveIter=1;
for i=1:4
    [status,iter,time]=solve(obj,mu,int32(maxIter),int32(saveIter));
    [J_star,position_l1l2_star,noise1_star,acceleration1_star,noise1abs_star,acceleration1abs_star]=...
        getOutputs(obj);
    saveIter=-1;
end
clear obj


%% plots
figure(1)
subplot(6,1,4)
plot(t,truePosition,'-.',t,thisMeasurement,'-+',t,position_l2_star,'-o',t,position_l1l2_star,'-*',...
     t(k_outlier),zeros(length(k_outlier),1),'*')
grid on
legend('true position','measurements','l2 estimates','l1-l2 estimates','outliers');
title('l1-l2 estimation')
subplot(6,1,5)
plot(t,thisMeasurement-position_l1l2_star-noise1_star,'-.',t,noise1_star,'-+',t(1:end-2),acceleration1_star,'-o')
grid on
legend('l2 noise','l1 noise','l1 acceleration')
subplot(6,1,6)
plot(t,position_l1l2_star-truePosition,'-+')
grid on
legend('l1-l2 estimation error')
drawnow

if allowSave && isequal(optimizeCS,@cmex2optimizeCS)

    [osize,subscripts,values,WW]=...
        loadCSparse('@tmpC_minl1l2/tmpC_minl1l2_WW.subscripts',...
                    '@tmpC_minl1l2/tmpC_minl1l2_WW.values');
    WWs=sparse(double(subscripts(1,:)),double(subscripts(2,:)),...
               ones(size(values)),size(WW,1),size(WW,2));
    
    figure(2);clf
    subplot(1,2,1)
    spy(WW);
    title('WW from loadCSparse')
    subplot(1,2,2)
    spy(WWs);
    title('subscripts from loadCSparse')
end

return

%% Disassemble

nm -g tmpC_minl1l2.dylib

otool -vt -p _tmpC_minl1l2_compute tmpC_minl1l2.dylib | more
