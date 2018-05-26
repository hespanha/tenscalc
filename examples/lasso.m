clear all;
%!rm -fr @tmp* tmp*

% debluring N*N image
N=256;
X=spdiags(ones(N*N,3),-1:1,N*N,N*N);
X=X+spdiags(ones(N*N,3),[-N,0,N],N*N,N*N);

classname=TClasso('classname','tmpLasso',...
                  'x',X,...
                  'addConstant',false,...
                  'solverType','matlab',...
                  'addEye2Hessian',1e-10,...
                  'useLDL',false,...
                  'smallerNewtonMatrix',true,...
                  'solverVerboseLevel',3);
                  
%help(classname);

%% Create data

% initialize seed for repeatability
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

% original image
im=zeros(N,N);
for i=1:ceil(N/5):N-4
    im(i:i+4,:)=1;
    im(:,i:i+4)=1;
end
w=reshape(sparse(im),N*N,1);
% bluring matrix
y=X*w+.1*randn(N*N,1);

% draw
figure(1);clf;
colormap(gray);
subplot(1,3,1);
imagesc(reshape(w,N,N));
axis image
subplot(1,3,2);
imagesc(reshape(y,N,N));
axis image


lambda=1;

%% Create object
obj=feval(classname);
% Set parameters
setP_y(obj,y);
setP_l1weight(obj,lambda);
% Initialize primal variables
w0=y;
setV_W(obj,w0);
setV_absW(obj,abs(w0)+1);
% Solve optimization
mu0=1;
maxIter=50;
saveIter=-1;
%profile on
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
%profile viewer;

% Get outputs
[w1,J1]=getOutputs(obj);


subplot(1,3,3);
imagesc(reshape(w1,N,N));
axis image




