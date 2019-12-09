clear all;
% remove previous classes
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Function to approximate


ex='log-det';
ex='sin';
ex='inv';
switch (ex)
  case 'inv'
    ninput=4;
    noutput=4;
    n=[ninput;30;15;15;noutput]; % # of neurons from 1st (input) to last (output) layer
    fun=@(u)reshape(inv(reshape(u,[2,2])),[4,1]);
    randData=@(u)reshape(eye(2)+5*reshape(u,[2,2])'*reshape(u,[2,2]),[4,1]);
    N=2000000;
    batchSize=200;
    alpha0=1e-3/batchSize;
    lambda=1e-3;
  case 'sin'
    ninput=1;
    noutput=1;
    n=[ninput;20;10;10;noutput]; % # of neurons from 1st (input) to last (output) layer
    fun=@(u)sin(u);
    randData=@(u)2*pi*u;    
    N=2000000;
    batchSize=5;
    alpha0=1e-2/batchSize;
    lambda=1e-2;
  case 'log-det'
    ninput=4;
    noutput=1;
    n=[ninput;30;10;10;noutput]; % # of neurons from 1st (input) to last (output) layer
    fun=@(u)log(det(reshape(u,[2,2])));
    randData=@(u)reshape(eye(2)+5*reshape(u,[2,2])'*reshape(u,[2,2]),[4,1]);
    N=5000000;
    batchSize=50;
    alpha0=5e-3/batchSize;
    lambda=1e-3;
end


%% Symbolic variables and computations

Tvariable u [n(1)];   % input layer
Tvariable y [n(end)]; % desired output;

x{1}=u;
totalalive=[];
whichalive=[];
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
    if i<length(n)-1
        totalalive=[totalalive;sum(heaviside(x{i+1}-eps),1)];
        whichalive=[whichalive;heaviside(x{i+1}-eps)];
    end
end

% loss function
Jreg=0;
%JJreg=[];
Jreg=Jreg+norm2(Wb);

Jloss=norm2(x{end}-y);

J=Jloss+lambda*Jreg;

% gradient step
Tvariable alpha [];

% gradient update (withiun batch)
new_gWb=gWb+gradient(J,Wb);
% parameter update (at batch end)
new_Wb=Wb-alpha*gWb;

% derivative of x{i} in the direction of the gradient (for relu's outputs)
for i=2:length(n)-1
    dx{i}=gradient(x{i},Wb)*gWb;
end

%% Generate code
code=csparse();

declareSet(code,u,'set_u');
declareSet(code,y,'set_y');
declareGet(code,x{end},'get_output');
declareSet(code,Wb,sprintf('set_Wb',i));
declareSet(code,alpha,'set_alpha');

declareCopy(code,gWb,full(Tzeros(size(gWb))),'resetGradient');
declareCopy(code,gWb,new_gWb,'updateGradient');
declareCopy(code,Wb,new_Wb,'updateParameters');

declareGet(code,[J;Jloss;Jreg],'get_J');
%declareGet(code,JJreg,'get_JJreg');

declareGet(code,cat(1,x{2:end-1}),sprintf('get_x',i));
declareGet(code,cat(1,dx{2:end}),sprintf('get_dx',i));
declareGet(code,{totalalive,whichalive},sprintf('get_alive',i));
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
set_Wb(obj,rand(size(Wb),1)-.5);

figure(1);clf;
figure(2);clf;
figure(3);clf;
training.meanJ=nan(N/batchSize,3);
training.totalalive=nan(N/batchSize,length(n)-2);
training.whichalive=nan(N/batchSize,sum(n(2:end-1)));
training.alpha=nan(N/batchSize,1);
set_alpha(obj,alpha0);
resetGradient(obj);
totalalive=zeros(length(n)-2,1);
whichalive=zeros(sum(n(2:end-1)),1);
sumJ=zeros(1,3);
batch=1;
minx=+inf(sum(n(2:end-1)),1);
for k=1:N
    %fprintf('k=%d\n',k);
    u=randData(2*rand(ninput,1)-1);%randn(2);u=eye(2)+u'*u;u=u(:);
    y=fun(u);
    set_u(obj,u);    
    set_y(obj,y);
    [a1,a2]=get_alive(obj);
    totalalive=totalalive+a1;
    whichalive=whichalive+a2;

    x=get_x(obj);
    j=x>eps;
    minx(j)=min(minx(j),x(j));

    updateGradient(obj);
    sumJ=sumJ+get_J(obj)';

    %% End-Of-Batch
    if mod(k,batchSize)==0
        training.meanJ(batch,:)=sumJ/batchSize;
        training.totalalive(batch,:)=totalalive/batchSize;
        training.whichalive(batch,:)=whichalive/batchSize;
        
        dx=get_dx(obj);
        j=minx>0 & dx>0;
        alpha=minx(j)./dx(j);
        alpha=min(alpha);
        if isempty(alpha) || alpha<alpha0
            %alpha0,
            training.alpha(batch)=alpha0;
        else
            %alpha,
            training.alpha(batch)=alpha;
        end
        set_alpha(obj,training.alpha(batch));
        %[minx,dx/1e5,minx-training.alpha(batch)*dx],;pause
        updateParameters(obj);
        
        if any(totalalive==0)
            error('netowrk is dead\n');
        end
    
        resetGradient(obj);
        batch=batch+1;
        totalalive=zeros(length(n)-2,1);
        whichalive=zeros(sum(n(2:end-1)),1);
        sumJ=zeros(1,3);
        minx=+inf(sum(n(2:end-1)),1);
    end

    %% Output - print
    if mod(k,50000)==0 || k==N
        fprintf('batch %4d/%4d: J=%10.3e Jloss=%10.3e Jreg=%10.3e\n',batch-1,N/batchSize,training.meanJ(batch-1,:));
        fprintf('   alive = ');
        fprintf(' %4.2f/%2d',[training.totalalive(batch-1,:);n(2:end-1)']);
        fprintf('\n');
    end

    %% Output - plot
    if mod(k,100000)==0 || k==N
        % low-pass filter with 100 time constants over total # of points
        tau=batch/100;
        Alow=lsim(ss(-1/tau*eye(size(training.totalalive,2)),1/tau*eye(size(training.totalalive,2)),...
                     eye(size(training.totalalive,2)),0),...
                  training.totalalive(1:batch-1,:),0:batch-2,training.totalalive(1,:));
        figure(1);
        subplot(2,1,1);
        plot(Alow);grid on;
        leg=cellstr(num2str((2:length(n)-1)','layer %d'));legend(leg{:},'location','best')
        xlabel('batch #');ylabel('# neurons alive');
        subplot(2,1,2);
        semilogy(training.alpha);grid on
        xlabel('batch #');ylabel('\alpha');
        figure(2);
        Alow=lsim(ss(-1/tau*eye(size(training.whichalive,2)),1/tau*eye(size(training.whichalive,2)),...
                     eye(size(training.whichalive,2)),0),...
                  training.whichalive(1:batch-1,:),0:batch-2,training.whichalive(1,:));
        xx=union(1:ceil(batch/100):batch-1,batch-1);
        surf(xx,1:size(Alow,2),...
             Alow(xx,:)');shading flat;
        caxis([0,1]);%set(gca,'colorscale','log');
        view(2);
        colorbar;colormap(gray);
        xlabel('batch #');ylabel('neuron #');
        axis tight;
        Jlow=lsim(ss(-1/tau*eye(size(training.meanJ,2)),1/tau*eye(size(training.meanJ,2)),...
                     eye(size(training.meanJ,2)),0),...
                  training.meanJ(1:batch-1,:),0:batch-2,training.meanJ(1,:));
        figure(3);
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
        xlabel('batch #');
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
test.u=nan(K,ninput);
test.y=nan(K,noutput);
test.yhat=nan(K,noutput);
for k=1:K
    u=randData(2*rand(ninput,1)-1);%randn(2);u=eye(2)+u'*u;u=u(:);
    y=fun(u);
    set_u(obj,u);    
    set_y(obj,y);
    yhat=get_output(obj);
    J=get_J(obj);
    %[x{:}]=get_x(obj);
    %alive=cellfun(@(x)sum(abs(x)>eps),x(2:end-1));
    alive=get_alive(obj);
    fprintf('%10.4f',[u',NaN,y',NaN,yhat',NaN,sqrt(J(2))]);
    fprintf('   alive = ');
    fprintf(' %2d/%2d',[alive';n(2:end-1)']);
    fprintf('\n');
    test.u(k,:)=u;
    test.y(k,:)=y;
    test.yhat(k,:)=yhat;
end

figure(4);clf;
plot(test.u(:,1),test.y(:,1),'*g',test.u(:,1),test.yhat(:,1),'or');grid on;
xlabel('u(1)')
legend('y','yhat');

