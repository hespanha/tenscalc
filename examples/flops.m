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

% Many ways to skin a cat...
%
% This example is about finding the x that minimizes
%       J(x)= x' A x - b x,    with A positive definite
% which can be solved exactly in one Newton step:
%    x=inv(A+A')*b;
%
% It is really about testing the speed of the LDL factorization for
% a full matrix

clear all
!rm -rf toremove.m tmp* @tmp*;

solverGeneration=@cmex2optimizeCS;
solverVerboseLevel=1;

N=150;

Tvariable x [N];
Tvariable A [N,N];
Tvariable b [N];

J=x*A*x-b*x;

[classname,code]=solverGeneration('classname','tmp_flops',...
                                  'objective',J,...
                                  'optimizationVariables',{x},...
                                  'constraints',{},...
                                  'outputExpressions',{J,x},...
                                  'parameters',{A,b},...
                                  'compilerOptimization','-O1',...
                                  'profiling',true,...
                                  'solverVerboseLevel',solverVerboseLevel);

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

A=randn(N,N);
b=randn(N,1);

obj=feval(classname);

% Set parameters
setP_A(obj,A);
setP_b(obj,b);
% Initialize primal variables
x0=ones(N,1);
setV_x(obj,x0);

% Solve optimization
mu0=1;
maxIter=30;
saveIter=-1;
clear iter time;
nRepeats=1000;
for i=1:nRepeats
    [status,iter(i),time(i)]=solve(obj,mu0,int32(maxIter),int32(saveIter));
end

if 0
    [J,x]=getOutputs(obj);
    fprintf('error = %e\n',norm(inv(A+A')*b-x));
end

clear obj;


fprintf('N=%d, object code size = %d\n',N,code.statistics.createGateway.dFileSize);
fprintf('   iter     :   mean=%.2f, min=%.2f, max=%.2f\n',...
        mean(iter),min(iter),max(iter));
fprintf('   time [us]:   median=%.2f, min=%.2f, max=%.2f\n',...
        1e6*median(time),1e6*min(time),1e6*max(time));

% matlab's factorization;
clear mtime;
for i=1:nRepeats
    t0=clock();C=ldl(A);dt=etime(clock(),t0);
    mtime(i)=dt;
end
dt=min(dt);
fprintf('   matlab LDL time [us] = %.2f\n',1e6*dt);

fprintf('N, flops, time [s], dylib file size [b], matlab ldl time [s]\n');
fprintf('%3d, ,%.3e,%10d,%.3e; %% solverVerboseLevel=%d\n',...
        N,min(time),code.statistics.createGateway.dFileSize,dt,solverVerboseLevel);

stats=[% N, flops, time [s], dylib file size [b]
       10,   1086486/1000,9.640e-07,     30896,7.000e-06; % solverVerboseLevel=1       
       30,   9276356/1000,4.960e-06,    198848,2.800e-05; % solverVerboseLevel=1
       60,  36613411/1000,1.979e-05,   1091776,8.400e-05; % solverVerboseLevel=1
       100,101244351/1000,6.615e-05,   4331728,1.990e-04; % solverVerboseLevel=1
      ];

% compute exponents by fit
fprintf('flops scaling:');
sf=fit(stats(:,1),stats(:,2),'power2');
disp(sf)
fprintf('time scaling:');
st=fit(stats(:,1),stats(:,3),'power2');
disp(st)
fprintf('dylib scaling:');
sc=fit(stats(:,1),stats(:,4),'power2');
disp(sc)
fprintf('matlab time scaling:');
stm=fit(stats(:,1),stats(:,5),'power2');
disp(stm)

figure(1);clf;
title('22nm Ivy Bridge 2.6 GHz Intel Core i7, L2 Cache 256 KB/core, L3 Cache 6 MB');
loglog(stats(:,1),1e-9*stats(:,2)./stats(:,3),'x-',...
       stats(:,1),1e-6*stats(:,2),'+-',...
       stats(:,1),1e3*stats(:,3),'d-',...
       stats(:,1),1e-6*stats(:,4),'o-',...
       stats(:,1),1e3*stats(:,5),'^-');
grid on
xlabel('N');
legend(sprintf('Giga flops (mean=%.2f)',1e-9*mean(stats(:,2)./stats(:,3))),...
       sprintf('Mflops \\propto N^{%.1f}',sf.b),...
       sprintf('time [ms] \\propto N^{%.1f}',st.b),...
       sprintf('dylib code size [Mb] \\propto N^{%.1f}',sc.b),...
       sprintf('matlab LDL time [ms] \\propto N^{%.1f}',stm.b),...
       'location','best');
axis tight;
