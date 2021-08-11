% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

!rm -f tmpC*

clear all

n1=30
n2=20
n3=20

profile on
t0=clock;
fprintf('Creating code... ');
createGateway('template','testTprod_raw.c',...
              'callType','include',...
              'compileGateways',true,...
              'preprocessParameters',{n1,n2,n3},...
              'verboseLevel',2);
fprintf('done creating code (%.2f sec)\n',etime(clock,t0));
profile viewer

A=rand(n1,n2,n3);
B=rand(n1,n2,n3);
C=rand(n1,n2,n3);

fprintf('Matlab:\n');
ind=[-1,1,2];
DD=mytprod(A,ind,B,ind,C,ind);
t0=clock;
DD=mytprod(A,ind,B,ind,C,ind);
fprintf('  mytprod: %.1f us\n',1e6*etime(clock,t0))
t0=clock;
DD=mytprod(A,ind,B,ind,C,ind);
fprintf('  mytprod: %.1f us\n',1e6*etime(clock,t0))

fprintf('C code:\n');
D=tmpC_testTprod(A,B,C);
t0=clock;
D=tmpC_testTprod(A,B,C);
fprintf('  csparse: %.1f us\n',1e6*etime(clock,t0))

if ~isequal(DD,D)
    error('mismatch %f\n',norm(DD(:)-D(:)))
end

