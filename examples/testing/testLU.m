% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all

streamSelect=RandStream.create('mt19937ar','seed',0);
RandStream.setGlobalStream(streamSelect);

n=6
m=2

profile on
t0=clock;
fprintf('Creating code... ');
createGateway('template','testLU_raw.c',...
              'callType','include',...
              'compileGateways',true,...
              'compilerOptimization','-O0',...
              'preprocessParameters',{n,m},...
              'verboseLevel',2);
fprintf('done creating code (%.2f sec)\n',etime(clock,t0));
profile viewer

WW=rand(n);
WW(abs(WW(:))<.6)=0;
WW=WW*WW';

B=rand(n,m);

fprintf('Matlab:\n');
WW

% matlab's sparse LU
sWW=sparse(WW);
[l,u,p,q]=lu(sWW,'vector');
full(l),full(u),p,q
full(l*u-sWW(p,q))


t0=clock;
[l,u,p,q]=lu(sWW,'vector');
fprintf('  lu: %.1f us\n',1e6*etime(clock,t0))
t0=clock;
[l,u,p,q]=lu(sWW,'vector');
fprintf('  lu: %.1f us\n',1e6*etime(clock,t0))

if 0
    X=WW\B;
    t0=clock;
    X=WW\B;
    fprintf('  mldivide: %.1f us\n',1e6*etime(clock,t0))
    t0=clock;
    X=WW\B;
    fprintf('  mldivide: %.1f us\n',1e6*etime(clock,t0))
end

fprintf('C code:\n');
[LU,X]=tmpC_testLU(WW,B);
t0=clock;
[LU,X]=tmpC_testLU(WW,B);
fprintf('  csparse: %.1f us\n',1e6*etime(clock,t0))

[i,j,v]=find(LU);
k=i>j;
L=full(sparse(i(k),j(k),v(k),n,n)+speye(n));
k=i<=j;
U=full(sparse(i(k),j(k),v(k),n,n));

if 1
    LU
    L
    U
    
    WW
    L*U
    
    X
    WW\B
end

if norm(WW-L*U)>eps
    fprintf('mismatch WW~=LU: %e\n',norm(WW-L*U))
end

if norm(X-WW\B)>eps
    fprintf('mismatch X~=WW\\B: %e\n',norm(X-WW\B))
end

