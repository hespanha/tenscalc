% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all;
% remove previous classes
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Function to approximate


N=1500000;
batchSize=5;
lambda0=1e-3;          % regularization weigth
alpha0=1e-2/batchSize; % step size

%% Symbolic variables and computations

n=[1;30;10;10;1]; % # of neurons from 1st (input) to last (output) layer
Tvariable u [1];  % input
Tvariable y [1];  % desired output;

x{1}=u;
for i=1:length(n)-1
    % parameters
    W{i}=Tvariable(sprintf('W%d',i),[n(i+1),n(i)]);
    b{i}=Tvariable(sprintf('b%d',i),n(i+1));
end
[Wb,packCmd,unpackCmd,W{:},b{:}]=packVariables({W{:},b{:}},'x',W{:},b{:});

Tvariable gWb size(Wb);

for i=1:length(n)-1
    % node values
    if i<length(n)-1
        x{i+1}=relu(W{i}*x{i}+b{i});
    else
        x{i+1}=W{i}*x{i}+b{i};
    end
end

% loss function
Tvariable lambda []; % regularization weight

Jreg=norm2(Wb);
Jloss=norm2(x{end}-y);
J=Jloss+lambda*Jreg;

Tvariable alpha []; % step size

% gradient update (withiun batch)
new_gWb=gWb+gradient(J,Wb);
% parameter update (at batch end)
new_Wb=Wb-alpha*gWb;

%% Generate code
code=csparse();

declareSet(code,u,'set_u');
declareSet(code,y,'set_y');
declareGet(code,x{end},'get_output');
declareSet(code,Wb,sprintf('set_Wb',i));
declareSet(code,alpha,'set_alpha');
declareSet(code,lambda,'set_lambda');

declareCopy(code,gWb,full(Tzeros(size(gWb))),'resetGradient');
declareCopy(code,gWb,new_gWb,'updateGradient');
declareCopy(code,Wb,new_Wb,'updateParameters');

declareGet(code,[J;Jloss;Jreg],'get_J');

declareGet(code,{W{:},b{:}},sprintf('get_Wb',i));

classname='tmpNN';
classname=cmex2compute('csparseObject',code,'pedigreeClass',classname,'executeScript','asneeded','compilerOptimization','-O0');
%classname=class2compute('classname',classname,'executeScript','asneeded','csparseObject',code,'compilerOptimization','-O0');

help(classname);
obj=feval(classname);

%% Train

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

% initialize parameters and gradients
set_Wb(obj,(rand(size(Wb),1)-.5)/5);

training.J=nan(N,3);
set_alpha(obj,alpha0);
set_lambda(obj,lambda0);
resetGradient(obj);
sumJ=zeros(1,3);
batch=1;
for k=1:N
    %fprintf('k=%d\n',k);
    u=4*pi*rand(1,1)-2*pi;
    y=sin(u);
    set_u(obj,u);    
    set_y(obj,y);
    
    updateGradient(obj);
    training.J(k,:)=get_J(obj)';

    %% End-Of-Batch
    if mod(k,batchSize)==0
        updateParameters(obj);

        resetGradient(obj);
        batch=batch+1;
        sumJ=zeros(1,3);
    end

    %% Output - plot
    if mod(k,100000)==0 || k==N
        % low-pass filter with 100 time constants over total # of points
        figure(1);clf
        tau=k/100;
        Jlow=lsim(ss(-1/tau*eye(size(training.J,2)),1/tau*eye(size(training.J,2)),...
                     eye(size(training.J,2)),0),...
                  training.J(1:k,:),0:k-1,training.J(1,:));
        subplot(3,1,1)
        semilogy(Jlow(:,1),'r');grid on;
        ylim([min(Jlow(:,1)),1.01*max(Jlow(min(k,100):end,1))])
        legend(sprintf('J     -> %10.3e',Jlow(k,1)));
        subplot(3,1,2)
        semilogy(1:k,training.J(1:k,:),'g',1:k,Jlow(:,2),'r');grid on;
        ylim([min(Jlow(:,2)),1.01*max(Jlow(min(batch-1,100):end,2))])
        legend(sprintf('Jloss -> %10.3e',Jlow(k,2)));
        subplot(3,1,3)
        semilogy(Jlow(:,3),'r');grid on;
        ylim([min(Jlow(:,3)),1.01*max(Jlow(min(batch-1,100):end,3))])
        xlabel(sprintf('batch # out of %.0e',N/batchSize));
        legend(sprintf('Jreg  -> %10.3e',Jlow(k,3)));
        drawnow;
    end

end

[W{:},b{:}]=get_Wb(obj);
for i=1:length(n)-1
    fprintf('norm2(W{%d})=%10.3e, norm2(b{%d})=%10.3e\n',i,norm2(W{i}),i,norm2(b{i}));
end

%% Test

% always same batch
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);
K=400;
test.u=nan(K,1);
test.y=nan(K,1);
test.yhat=nan(K,1);
for k=1:K
    u=4*pi*rand(1,1)-2*pi;
    y=sin(u);
    set_u(obj,u);    
    set_y(obj,y);
    yhat=get_output(obj);
    J=get_J(obj);
    fprintf('%10.4f',[u',NaN,y',NaN,yhat',NaN,sqrt(J(2))]);
    fprintf('\n');
    test.u(k,:)=u;
    test.y(k,:)=y;
    test.yhat(k,:)=yhat;
end

figure(2);clf;
plot(test.u(:,1),test.y(:,1),'*g',test.u(:,1),test.yhat(:,1),'or');grid on;
xlabel('u(1)')
legend('y','yhat');
title(sprintf('mean err = %.2e',norm(test.y(:)-test.yhat(:))/prod(size(test.y))))
