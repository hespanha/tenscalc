% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all

!rm -f tmpC*

N=6;

profile on
t0=clock;
fprintf('Creating code... ');
createGateway('template','testSubsref_raw.c',...
              'callType','include',...
              'compileGateways',true,...
              'compilerOptimization','-O0',...
              'preprocessParameters',{N},...
              'verboseLevel',2);
fprintf('done creating code (%.2f sec)\n',etime(clock,t0));
profile viewer

t=cumsum(rand(N,1));
dt1=1./(t(2:end)-t(1:end-1));
X=rand(N,1);

fprintf('Matlab:\n');
t0=clock;
DX1=(X(2:end)-X(1:end-1)).*dt1;             % velocity
DDX1=(DX1(2:end)-DX1(1:end-1)).*dt1(1:end-1); % acceleration
fprintf('  %.1f us\n',1e6*etime(clock,t0))
DX1,
DDX1,

fprintf('C code:\n');
t0=clock;
[DX,DDX]=tmpC_testSubsref(X,dt1);
fprintf('  csparse: %.1f us\n',1e6*etime(clock,t0))
DX,
DDX,

if ~isequal(DX,DX1) 
    error('mismatch %f\n',norm(DX(:)-DX1(:)))
end

if ~isequal(DDX,DDX1)
    error('mismatch %f\n',norm(DDX(:)-DDX1(:)))
end

