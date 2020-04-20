clear all

streamSelect=RandStream.create('mt19937ar','seed',0);
RandStream.setGlobalStream(streamSelect);

n=3
m=1

%profile on
t0=clock;
fprintf('Creating code... ');
createGateway('template','testLUatomic_raw.c',...
              'callType','include',...
              'compileGateways',true,...
              'compilerOptimization','-O0',...
              'preprocessParameters',{n},...
              'verboseLevel',2);
fprintf('done creating code (%.2f sec)\n',etime(clock,t0));
%profile viewer

WW=rand(n);
WW(abs(WW(:))<.6)=0;
WW=WW*WW';

B=rand(n,m);

fprintf('Matlab:\n');
%WW

% matlab's sparse LU
sWW=sparse(WW);
if 0
    [l,u,p,q]=lu(sWW,'vector');
    full(l),full(u),p,q
    full(l*u-sWW(p,q))
end

t0=clock;
[l,u,p,q]=lu(sWW,'vector');
fprintf('  lu: %.1f us\n',1e6*etime(clock,t0))
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
t0=clock;
[X]=tmpC_testLUatomic(WW,B);
fprintf('  csparse: %.1f us\n',1e6*etime(clock,t0))
t0=clock;
[X]=tmpC_testLUatomic(WW,B);
fprintf('  csparse: %.1f us\n',1e6*etime(clock,t0))
t0=clock;
[X]=tmpC_testLUatomic(WW,B);
fprintf('  csparse: %.1f us\n',1e6*etime(clock,t0))

if 0
    X
    WW\B
end

if norm(X-WW\B)>eps
    fprintf('mismatch X~=WW\\B: %e\n',norm(X-WW\B))
end

