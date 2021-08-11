% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all

!rm -f tmpC*

streamSelect=RandStream.create('mt19937ar','seed',0);
RandStream.setGlobalStream(streamSelect);

n=6
m=2

profile on
t0=clock;
fprintf('Creating code... ');
createGateway('template','testLDL_raw.c',...
              'callType','include',...
              'compileGateways',true,...
              'compilerOptimization','-O0',...
              'preprocessParameters',{n,m},...
              'verboseLevel',2);
fprintf('done creating code (%.2f sec)\n',etime(clock,t0));
profile viewer

WW=rand(n);
WW=WW*diag(randn(n,1))*WW';
WW(abs(WW(:))<.1)=0;

B=rand(n,m);

fprintf('Matlab:\n');
WW

% matlab's sparse LDL
sWW=sparse(WW);
[l,d,p]=ldl(sWW,'vector');
full(l),full(d),p
full(l*d*l'-sWW(p,p))


t0=clock;
[l,u,p]=ldl(sWW,'vector');
fprintf('  ldl: %.1f us\n',1e6*etime(clock,t0))
t0=clock;
[l,u,p]=ldl(sWW,'vector');
fprintf('  ldl: %.1f us\n',1e6*etime(clock,t0))

if 1
    X=WW\B;
    t0=clock;
    X=WW\B;
    fprintf('  mldivide: %.1f us\n',1e6*etime(clock,t0))
    t0=clock;
    X=WW\B;
    fprintf('  mldivide: %.1f us\n',1e6*etime(clock,t0))
end

fprintf('C code:\n');
[LDL,X]=tmpC_testLDL(WW,B);
t0=clock;
[LDL,X]=tmpC_testLDL(WW,B);
fprintf('  csparse: %.1f us\n',1e6*etime(clock,t0))
L=tril(LDL,-1)+eye(n);
U=triu(LDL,1)+eye(n);
D=diag(diag(LDL));
X

if 1
    LDL
    L
    D
    U
    
    WW
    L*D*U
    
    X
    WW\B
end

if norm(WW-L*D*U)>eps
    WW
    L*D*U
    fprintf('mismatch WW~=LDL: %e (perhaps due to permutation)\n',norm(WW-L*D*U))
end

if norm(X-WW\B)>eps
    fprintf('mismatch X~=WW\\B: %e (should match even with permutation)\n',norm(X-WW\B))
end

