clear all;
% remove previous classes
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Create symbolic ecxpressions

Tvariable theta [6];  % initial position/velocity
Tvariable t [];       % time
Tvariable M [3,3];    % camera matrix parameter
Tvariable p [3];      % camera position
Tvariable invS [2,2]; % error covariance;

q=theta(1:3)+theta(4:6)*t;
mu=(M(1:2,:)*(p-q))./(M([3,3],:)*(p-q));


g=gradient(mu,theta);

FIM=tprod(g,[-1,1],invS,[-1,-2],g,[-2,2]);
%FIM=g'*invS*g;

%% Generate code
code=csparse();

declareSet(code,theta,'set_theta');
declareSet(code,t,'set_t');
declareSet(code,M,'set_M');
declareSet(code,p,'set_p');
declareSet(code,invS,'set_invS');

declareGet(code,FIM,'get_FIM');

cmex2compute('csparseObject',code,'classname','tmpFIM');
%class2compute('csparseObject',code,'classname','tmpFIM');

obj=tmpFIM();

%% Call code

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

theta=rand(6,1);
set_theta(obj,theta);
M=eye(3)+rand(3,3);
set_M(obj,M);
invS=rand(2,2);invS=invS'*invS;
set_invS(obj,invS);

FIM=zeros(6,6);
fprintf('Comoputing FIMs... ');
t0=clock();
for i=1:100000
    t=rand(1);
    set_t(obj,t);
    p=5+rand(3,1); % 5+ to keep p away from q
    set_p(obj,p);
    FIM=FIM+get_FIM(obj);
end
fprintf('done %.3f sec\n',etime(clock(),t0));
FIM,;
