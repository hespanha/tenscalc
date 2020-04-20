!rm -f tmpC*

clear all

N=4

profile on
t0=clock;
fprintf('Creating code... ');
createGateway('template','testTprod2_raw.c',...
              'callType','include',...
              'compileGateways',true,...
              'preprocessParameters',{N},...
              'verboseLevel',2);
fprintf('done creating code (%.2f sec)\n',etime(clock,t0));
profile viewer

A=rand(N,N);

fprintf('Matlab:\n');
t0=clock;
eye_21=eye(N);
plus_26=eye_21(2:end,:)-eye_21(1:end-1,:)
plus_30=plus_26(2:end,:)-plus_26(1:end-1,:)
tprod_37=tprod(plus_30,[-1,1],plus_30,[-1,2])
fprintf('  mytprod: %.1f us\n',1e6*etime(clock,t0))

fprintf('C code:\n');
[p26,p30,t37]=tmpC_testTprod2(A)
t0=clock;

[p26,p30,t37]=tmpC_testTprod2(A)
fprintf('  csparse: %.1f us\n',1e6*etime(clock,t0))

if ~isequal(tprod_37,t37)
    error('mismatch %f\n',norm(tprod_37(:)-t37(:)))
end

