!rm -f tmpC*

clear all

n1=30
n2=40
n3=2


profile on
t0=clock;
fprintf('Creating code... ');
createGateway('template','testSum_raw.c',...
              'callType','include',...
              'compileGateways',true,...
              'preprocessParameters',{n1,n2,n3},...
              'verboseLevel',2);
fprintf('done creating code (%.2f sec)\n',etime(clock,t0));
profile viewer

for i=1:3
    A=rand(n1,n2,n3);
    B=rand(n1,n2,n3);
    C=rand(n1,n2,n3);

    t0=clock;
    DD=-A+B+C;
    fprintf('  Matlab plus(%d): %.1f us\n',i,1e6*etime(clock,t0))

    t0=clock;
    D=tmpC_testSum(A,B,C);
    fprintf('  csparse(%d)    : %.1f us\n',i,1e6*etime(clock,t0))
end

if ~isequal(DD,D)
    error('mismatch %f\n',norm(DD(:)-D(:)))
end

