% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

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
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

solverGeneration=@cmex2optimizeCS;
solverVerboseLevel=1;

N=60;


% -O0
compilerOptimization='-O0';
stats=[% N,        flops,  time [s], dylib [b], data [dbls]
        30, 1.962100e+04, 1.640e-05,    280824,      2138,2.100e-05; % solverVerboseLevel=1       
        60, 1.139460e+05, 8.797e-05,   1550584,      7868,8.100e-05; % solverVerboseLevel=1
       100, 4.492460e+05, 5.053e-04,   6109448,     21108,1.360e-04; % solverVerboseLevel=1
       200, 3.128496e+06, 4.625e-03,  42973448,     82208,7.230e-04; % solverVerboseLevel=1
       300, 1.003775e+07, 2.359e-02, 138610952,    183308,1.984e-03; % solverVerboseLevel=1
      ];

% -O1
compilerOptimization='-O1';
stats=[% N,        flops,  time [s], dylib [b], data [dbls]
        30, 1.962100e+04, 1.189e-05,    199032,      2138,2.000e-05; % solverVerboseLevel=1
        60, 1.139460e+05, 6.395e-05,   1096056,      7868,8.400e-05; % solverVerboseLevel=1      
       100, 4.492460e+05, 2.811e-04,   4331912,     21108,1.350e-04; % solverVerboseLevel=1
       150, 1.385121e+06, 1.382e-03,  13461896,     46658,5.840e-04; % solverVerboseLevel=1
      ];


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
                                  'compilerOptimization',compilerOptimization,...
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
    setP_A(obj,A);
    setP_b(obj,b);
    setV_x(obj,x0);
    [status,iter(i),time(i)]=solve(obj,mu0,int32(maxIter),int32(saveIter));
end

if 0
    [J,x]=getOutputs(obj);
    fprintf('error = %e\n',norm(inv(A+A')*b-x));
end

%% get profile data
clear obj;
prof=textread(sprintf('@%s/%s.c.profile',classname,classname),'%s');
k=min(find(strcmp(prof,'Total')));
totalFlops=str2double(prof{k+1})/nRepeats;

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

fprintf('N, flops, time [s], dylib file size [b], data size [doubles], matlab ldl time [s]\n');
fprintf('%3d, %10d, %.3e,%10d,%10d,%.3e; %% solverVerboseLevel=%d\n',...
        N,totalFlops,min(time),code.statistics.createGateway.dFileSize,code.statistics.sizeScratchbook,dt,solverVerboseLevel);

%% Scalability plots

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
fprintf('data scaling:');
sd=fit(stats(:,1),stats(:,5),'power2');
disp(sd)
fprintf('matlab time scaling:');
stm=fit(stats(:,1),stats(:,6),'power2');
disp(stm)

fig=clearFigure('figureNumber',1,'figureName',sprintf('flops scalability (%s)',compilerOptimization));
title('22nm Ivy Bridge 2.6 GHz Intel Core i7, L2 Cache 256 KB/core, L3 Cache 6 MB');
loglog(stats(:,1),1e-6*(stats(:,1).^3/3+11.5*stats(:,1).^2),'o-',...
       stats(:,1),1e-6*stats(:,2),'+-',...
       stats(:,1),1e-6*stats(:,4),'v-',...
       stats(:,1),1e-6*8*stats(:,5),'p-',...
       stats(:,1),1e-9*stats(:,2)./stats(:,3),'x-',...
       stats(:,1),1e3*stats(:,3),'d-'...%,...
       ...%stats(:,1),1e3*stats(:,6),'^-'
       );
grid on
xlabel('N');
legend(sprintf('N^3/3+11.5*N^2'),...
       sprintf('Mflops \\propto N^{%.1f}',sf.b),...
       sprintf('dylib code size [Mb] \\propto N^{%.1f}',sc.b),...
       sprintf('data size [Mb] \\propto N^{%.1f}',sd.b),...
       sprintf('Giga flops (max=%.2f)',1e-9*max(stats(:,2)./stats(:,3))),...
       sprintf('time [ms] \\propto N^{%.1f}',st.b),...
       ...%sprintf('matlab LDL time [ms] \\propto N^{%.1f}',stm.b),...
       'location','northwest');
axis tight;
saveFigure('figureNumber',fig,'folder','plots','createFolder',true);
