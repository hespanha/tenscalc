% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all

!rm -rf @tmp* tmp*

n1=3%0
n2=4%0
n3=2


%profile on
t0=clock;
fprintf('Creating code... ');
createGateway('template','testSet_raw.c',...
              'callType','include',...
              'compileGateways',true,...
              'preprocessParameters',{n1,n2,n3},...
              'verboseLevel',1);
fprintf('done creating code (%.2f sec)\n',etime(clock,t0));
%profile viewer

for i=1:3
    A=rand(n1,n2,n3);
    B=rand(n1,n2,n3);
    C=rand(n1,n2,n3);
    
    ABC=cat(1,A,B,C);
    
    AB=ABC(1:2*n1,:,:);
    BC=ABC(n1+1:end,:,:);
    

    t0=clock;
    DD=-AB+BC;
    fprintf('  Matlab plus(%d): %.1f us\n',i,1e6*etime(clock,t0))

    t0=clock;
    D=tmpC_testSet(A,B,C);
    fprintf('  csparse(%d)    : %.1f us\n',i,1e6*etime(clock,t0))
end

if ~isequal(DD,D)
    error('mismatch %f\n',norm(DD(:)-D(:)))
end

