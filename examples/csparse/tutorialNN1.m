% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all;
% remove previous classes
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Function to approximate


N=1000000;
batchSize=5;


%% Symbolic variables and computations

n=[1;30;5;10;1]; % # of neurons from 1st (input) to last (output) layer
Tvariable u [1];  % input
Tvariable y [1];  % desired output;

x{1}=u;
for i=1:length(n)-1
    % parameters
    W{i}=Tvariable(sprintf('W%d',i),[n(i+1),n(i)]);
    b{i}=Tvariable(sprintf('b%d',i),n(i+1));
    % gradients
    gW{i}=Tvariable(sprintf('gW%d',i),[n(i+1),n(i)]);
    gb{i}=Tvariable(sprintf('gb%d',i),n(i+1));
    % node values
    if i<length(n)-1
        x{i+1}=relu(W{i}*x{i}+b{i});
    else
        x{i+1}=W{i}*x{i}+b{i};
    end
end

% loss function
Jreg=0;
%JJreg=[];
for i=1:length(n)-1
    Jreg=Jreg+norm2(W{i})+norm2(b{i});
end
Jloss=norm2(x{end}-y);

lambda=1e-3;
J=Jloss+lambda*Jreg;

% gradient step
Tvariable alpha [];

for i=1:length(n)-1
    % gradient update (withiun batch)
    new_gW{i}=gW{i}+gradient(J,W{i});
    new_gb{i}=gb{i}+gradient(J,b{i});
    % parameter update (at batch end)
    new_W{i}=W{i}-alpha*gW{i};
    new_b{i}=b{i}-alpha*gb{i};
    % reset values for gradient
    zero_gW{i}=full(Tzeros(size(gW{i})));
    zero_gb{i}=full(Tzeros(size(gb{i})));
end

%% Generate code
code=csparse();

declareSet(code,u,'set_u');
declareSet(code,y,'set_y');
declareGet(code,x{end},'get_output');
for i=1:length(n)-1
    declareSet(code,W{i},sprintf('set_W%d',i));
    declareSet(code,b{i},sprintf('set_b%d',i));
end    
declareSet(code,alpha,'set_alpha');

declareCopy(code,{gW{:},gb{:}},{zero_gW{:},zero_gb{:}},'resetGradient');
declareCopy(code,{gW{:},gb{:}},{new_gW{:},new_gb{:}},'updateGradient');
declareCopy(code,{W{:},b{:}},{new_W{:},new_b{:}},'updateParameters');

declareGet(code,[J;Jloss;Jreg],'get_J');

declareGet(code,{x{:}},sprintf('get_x',i));
declareGet(code,{W{:},b{:}},sprintf('get_Wb',i));
declareGet(code,{gW{:},gb{:}},sprintf('get_gWb',i));


classname='tmpNN';
classname=cmex2compute('csparseObject',code,'pedigreeClass',classname,'executeScript','asneeded','compilerOptimization','-O0');
%classname=class2compute('classname',classname,'executeScript','asneeded','csparseObject',code,'compilerOptimization','-O0');

help(classname);
obj=feval(classname);

%% Train

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

% initialize parameters and gradients
for i=1:length(n)-1
    % initialize parameters
    W{i}=(rand(n(i+1),n(i))-.5)/n(i);     %  both signs
    feval(sprintf('set_W%d',i),obj,W{i});
    b{i}=rand(n(i+1),1);                  % positive to get neurons activated
    feval(sprintf('set_b%d',i),obj,b{i});
end

training.meanJ=nan(N/batchSize,3);
set_alpha(obj,1e-2/batchSize);
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
    sumJ=sumJ+get_J(obj)';

    %% End-Of-Batch
    if mod(k,batchSize)==0
        training.meanJ(batch,:)=sumJ/batchSize;
        
        updateParameters(obj);

        resetGradient(obj);
        batch=batch+1;
        sumJ=zeros(1,3);
    end

    %% Output - plot
    if mod(k,100000)==0 || k==N
        % low-pass filter with 100 time constants over total # of points
        figure(1);clf
        tau=batch/100;
        Jlow=lsim(ss(-1/tau*eye(size(training.meanJ,2)),1/tau*eye(size(training.meanJ,2)),...
                     eye(size(training.meanJ,2)),0),...
                  training.meanJ(1:batch-1,:),0:batch-2,training.meanJ(1,:));
        subplot(3,1,1)
        semilogy(Jlow(:,1));grid on;
        ylim([min(Jlow(:,1)),1.01*max(Jlow(min(batch-1,100):end,1))])
        legend(sprintf('J     -> %10.3e',Jlow(batch-1,1)));
        subplot(3,1,2)
        semilogy(Jlow(:,2));grid on;
        ylim([min(Jlow(:,2)),1.01*max(Jlow(min(batch-1,100):end,2))])
        legend(sprintf('Jloss -> %10.3e',Jlow(batch-1,2)));
        subplot(3,1,3)
        semilogy(Jlow(:,3));grid on;
        ylim([min(Jlow(:,3)),1.01*max(Jlow(min(batch-1,100):end,3))])
        xlabel(sprintf('batch # out of %d',N/batchSize));
        legend(sprintf('Jreg  -> %10.3e',Jlow(batch-1,3)));
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
