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

%   sparsifyLU(89): size=[996,996]  values from "tmpC_minl1l2_WW.values", pivoting from lu, nnz=4174 (  0.42%) <- LU   (size=[996,996] , nnz=4174) (2.46 sec)
%   done ipmPD_CS symbolic computations (13.318 sec)
%   computeScalarInstructions... done (5.531 sec)
%   dependencyGroups... done (0.442 sec)
%   compile2C: 1499 added vectorized Operations, 141 unique vectorized Operations
%   compile2C: 45988 added instructions, 41810 unique instructions
%   compile2C: sizeScratchbook=41810 double, nGroups=1255, nSets=16, nGets=9, nCopies=3

%%  write codeType=C+asmLB... done (1.157 sec)
%   Compiling 16 gateways... done compiling 16 gateways (10.50 sec)
%%  Compiling library...  optimization = -O1 tmpC_minl1l2.c = 6692.808kB, tmpC_minl1l2.dylib = 1820.230kB  (3.462 sec)
% done cmex2optimizeCS (34.785 sec)
% ipmPD_CS (skipAffine=0,delta=3,allowSave=0): 996 primal variable, 0 equality constraints, 796 inequality constraints
%  14:status=0, cost=  5.84853e+01, |grad|=  4.99e-10, ineq=  6.61e-09, dual=  3.83e-11, gap=  8.02e-07, last alpha=  9.90e-01 (1931.1us,137.94us/iter)

%%  write codeType=C+asmLB... done (1.176 sec)
%   Compiling 16 gateways... done compiling 16 gateways (10.40 sec)
%%  Compiling library...  optimization = -O0 tmpC_minl1l2.c = 6692.808kB, tmpC_minl1l2.dylib = 1992.191kB  (1.342 sec)
% done cmex2optimizeCS (32.310 sec)
% ipmPD_CS (skipAffine=0,delta=3,allowSave=0): 996 primal variable, 0 equality constraints, 796 inequality constraints
%  14:status=0, cost=  5.84853e+01, |grad|=  4.99e-10, ineq=  6.61e-09, dual=  3.83e-11, gap=  8.02e-07, last alpha=  9.90e-01 (2014.7us,143.91us/iter)

%%  write codeType=C+asmLB... done (1.162 sec)
%   Compiling 16 gateways... done compiling 16 gateways (10.44 sec)
%%  Compiling library...  optimization = -Ofast tmpC_minl1l2.c = 6692.808kB, tmpC_minl1l2.dylib = 7532.285kB  (15.057 sec)
% done cmex2optimizeCS (46.098 sec)
% ipmPD_CS (skipAffine=0,delta=3,allowSave=0): 996 primal variable, 0 equality constraints, 796 inequality constraints
%  15:status=0, cost=  5.84853e+01, |grad|=  1.71e-10, ineq=  2.04e-09, dual=  1.14e-11, gap=  2.21e-07, last alpha=  9.90e-01 (2340.5us,156.04us/iter)

%%  write codeType=C... done (1.061 sec)
%   Compiling 16 gateways... done compiling 16 gateways (10.40 sec)
%%  Compiling library...  optimization = -Ofast tmpC_minl1l2.c = 2963.523kB, tmpC_minl1l2.dylib = 1852.141kB  (54.696 sec)
% done cmex2optimizeCS (85.566 sec)
% ipmPD_CS (skipAffine=0,delta=3,allowSave=0): 996 primal variable, 0 equality constraints, 796 inequality constraints
%  14:status=0, cost=  5.84853e+01, |grad|=  1.77e-10, ineq=  1.92e-09, dual=  1.19e-11, gap=  2.76e-07, last alpha=  9.90e-01 (1868.7us,133.48us/iter)

%%  write codeType=C... done (1.107 sec)
%   Compiling 16 gateways... done compiling 16 gateways (10.42 sec)
%%  Compiling library...  optimization = -O1 tmpC_minl1l2.c = 2963.523kB, tmpC_minl1l2.dylib = 1576.094kB  (22.734 sec)
% done cmex2optimizeCS (53.687 sec)
% ipmPD_CS (skipAffine=0,delta=3,allowSave=0): 996 primal variable, 0 equality constraints, 796 inequality constraints
%  15:status=0, cost=  5.84853e+01, |grad|=  1.17e-09, ineq=  1.39e-08, dual=  7.80e-11, gap=  1.50e-06, last alpha=  9.31e-01 (2178.0us,145.20us/iter)

%%  write codeType=C... done (1.062 sec)
%   Compiling 16 gateways... done compiling 16 gateways (10.37 sec)
%%  Compiling library...  optimization = -O0 tmpC_minl1l2.c = 2963.523kB, tmpC_minl1l2.dylib = 3172.312kB  (1.986 sec)
% done cmex2optimizeCS (32.838 sec)
% ipmPD_CS (skipAffine=0,delta=3,allowSave=0): 996 primal variable, 0 equality constraints, 796 inequality constraints
%  15:status=0, cost=  5.84853e+01, |grad|=  1.17e-09, ineq=  1.39e-08, dual=  7.80e-11, gap=  1.50e-06, last alpha=  9.31e-01 (3649.3us,243.29us/iter)


%% https://llvm.org/bugs/show_bug.cgi?id=23366
%optimizeCS=@class2optimizeCS; % matlab
optimizeCS=@cmex2optimizeCS;  % c
N=200
%codeType='C+asmLB';
%compilerOptimization='-O1';

codeType='C';
compilerOptimization='-O0';

allowSave=true;
saveIter=1;

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
        'compilerOptimization',compilerOptimization,...
       ...%'callType','client-server','serverProgramName','tmpC_minl2_server','serverAddress','localhost','port',1968,...
        ...%'callType','client-server','serverComputer','glnxa64','compileStandalones',false,'serverProgramName','tmpC_minl2_server','serverAddress','tamina.ece.ucsb.edu','port',1968,...
        'allowSave',allowSave,...
        'profiling',false,...
        'verboseLevel',1);

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
    for i=1:4
        [status,iter,time]=solve(obj,mu,int32(maxIter),int32(saveIter));
        [J_star,position_l2_star]=getOutputs(obj);
    end
    clear obj;

    %% Plots
    figure(1),clf
    subplot(4,1,1)
    plot(t,truePosition,'-.',t,thisMeasurement,'-+',t,position_l2_star,'-o',...
         t(k_outlier),zeros(length(k_outlier),1),'*')
    grid on
    legend('true position','measurements','l2 estimates','outliers');
    title('l2 estimate')
    subplot(4,1,2)
    plot(t,thisMeasurement-position_l2_star,'-+')
    grid on
    legend('l2 estimated noise')
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
    'compilerOptimization',compilerOptimization,...
...%'callType','client-server','serverProgramName','tmpC_minl1l2_server','serverAddress','localhost','port',1968,...
...%'callType','client-server','serverComputer','glnxa64','compileStandalones',false,'serverProgramName','tmpC_minl1l2_server','serverAddress','tamina.ece.ucsb.edu','port',1968,...
    'allowSave',allowSave,...
    'verboseLevel',1);
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
for i=1:4
    [status,iter,time]=solve(obj,mu,int32(maxIter),int32(saveIter));
    [J_star,position_l1l2_star,noise1_star,acceleration1_star,noise1abs_star,acceleration1abs_star]=...
        getOutputs(obj);
    saveIter=-1;
end
clear obj


%% plots
figure(1)
subplot(4,1,3)
plot(t,truePosition,'-.',t,thisMeasurement,'-+',t,position_l2_star,'-o',t,position_l1l2_star,'-*',...
     t(k_outlier),zeros(length(k_outlier),1),'*')
grid on
legend('true position','measurements','l2 estimates','l1-l2 estimates','outliers');
title('l1-l2 estimation')
subplot(4,1,4)
plot(t,thisMeasurement-position_l1l2_star-noise1_star,'-.',t,noise1_star,'-+',t(1:end-2),acceleration1_star,'-o')
grid on
legend('l2 noise','l1 noise','l1 acceleration')
title('l1-l2 estimate')
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
