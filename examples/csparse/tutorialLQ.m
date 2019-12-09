clear all;
% remove previous classes
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Defines sizes

N=100; % length of x
n=2;   % length of x0
k=10;  % length of u

%% Create symbolic ecxpressions

Tvariable A  [N,n];
Tvariable x0 [n];
Tvariable B  [N,k];
Tvariable u  [k];

x=A*x0+B*u;

J=norm2(x)+norm2(u);

g=gradient(J,u);
h=gradient(g,u);

factor=ldl(h); % faster but less robust and needs to be symmetric
               %factor=lu(h);

ustar=-(factor\g);

%% Generate code
code=csparse();

declareSet(code,A, 'set_A','set variable A');
declareSet(code,x0,'set_x0','set variable x0');
declareSet(code,B, 'set_B','set variable B');
declareSet(code,u, 'set_u','set variable u');

declareGet(code,{J,g,h},'get_Jgh','get cost, gradient and hessian');
declareGet(code,{ustar},'get_ustar','optimal u');
declareCopy(code,u,ustar,'copy_ustar2u','copy ustar to u');

cmex2compute('csparseObject',code,'classname','tmpLQ');
%class2compute('csparseObject',code,'classname','tmpLQ');

obj=tmpLQ();

%% Call code

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

A=rand(N,n);
x0=rand(n,1);
B=rand(N,k);

set_A(obj,A);
set_B(obj,B);
set_x0(obj,x0);
set_u(obj,zeros(k,1));

[J,g,h]=get_Jgh(obj); % cost, grad, and hessian for zero 0
fprintf('cost for u=0 is %f\n',J);

ustar=get_ustar(obj);
copy_ustar2u(obj); % be careful: after copy ustar changes!

[J,g,h]=get_Jgh(obj); % cost, grad, and hessian for optimal u
fprintf('cost for optimal u is %f\n',J);
