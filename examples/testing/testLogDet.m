clear all;
% remove previous classes
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

%% Create symbolic ecxpressions
n=5;

sparsity=ones(n);
if 0
    sparsity([1,5],[1,5])=0;
    sparsity(3,3)=0;
elseif 1
    sparsity(1,n)=0;
    sparsity(n,1)=0;
    sparsity(1,3)=0;
    sparsity(3,1)=0;
else
    sparsity(1,1)=0;
end

Tvariable A [n,n];

AS=A.*sparsity;

%% Generate code
code=csparse();

umfpack=false;
if umfpack
    % not yet implemented
    luAS=lu_sym(AS);
    luAS=declareAlias(code,lu(AS),'luAS',true)    ;
else
    luAS=lu(AS);    
end
declareSet(code,A,'set_A');
declareGet(code,full(AS),'get_AS');
declareGet(code,trace(AS),'get_traceAS');
declareGet(code,full(gradient(trace(AS),A)),'get_gtraceAS');
declareGet(code,logdet(lu(AS)),'get_logdetLUAS');
declareGet(code,logdet(ldl(AS)),'get_logdetLDLAS');
declareGet(code,lu(AS)\eye(n),'get_mldLUAS');
declareGet(code,full(gradient(logdet(ldl(AS)),A)),'get_glogdetLDLAS');
declareGet(code,full(gradient(logdet(lu(AS)),A)),'get_glogdetLUAS');
declareGet(code,traceinv(lu(AS)),'get_traceinvLUAS');
declareGet(code,traceinv(ldl(AS)),'get_traceinvLDLAS');
declareGet(code,full(gradient(traceinv(ldl(AS)),A)),'get_gtraceinvLDLAS');
declareGet(code,full(gradient(traceinv(lu(AS)),A)),'get_gtraceinvLUAS');
declareGet(code,inv(lu(AS)),'get_invLUAS');
declareGet(code,inv(ldl(AS)),'get_invLDLAS');
declareGet(code,full(gradient(inv(ldl(AS)),A)),'get_ginvLDLAS');

cmex2compute('csparseObject',code,'classname','tmpTest');
%class2compute('csparseObject',code,'classname','tmpTest');

obj=tmpTest();

%% Call code

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

A=rand(n);
if 1
    % symmetric case
    A=A*A';
end
AS=A.*sparsity;
if det(AS)<0
    A=-A;
    AS=-AS;
end
set_A(obj,A);
fprintf('AS =\n');
disp(AS);
AS1=get_AS(obj);
if norm2(AS-AS1)>sqrt(eps)
    AS1
    error('mismatch');
end

% trace
fprintf('trace(AS) = %f\n',trace(AS));
traceAS=get_traceAS(obj);
fprintf('trace(AS) = %f\n',full(traceAS));
if abs(trace(AS)-traceAS)>sqrt(eps)
    error('mismatch');
end

% grad trace
gtraceDiff=numericalGradient(@(A)trace(A.*sparsity),A,sqrt(eps));
fprintf('grad-diff trace(AS) =\n');
disp(full(gtraceDiff))
gtraceAS=get_gtraceAS(obj);
fprintf('grad trace(AS) =\n');
disp(full(gtraceAS))
if norm(gtraceDiff-gtraceAS)>sqrt(eps)
    norm(gtraceDiff-gtraceAS)
    error('mismatch');
end


% log-det
fprintf('log(det(AS)) = %f\n',log(det(AS)));
logdetLUAS=get_logdetLUAS(obj);
fprintf('logdetLU(AS) = %f\n',full(logdetLUAS));
logdetLDLAS=get_logdetLDLAS(obj);
fprintf('logdetLDL(AS) = %f\n',full(logdetLDLAS));
if abs(log(det(AS))-logdetLUAS)>sqrt(eps)
    abs(log(det(AS))-logdetLUAS)
    error('mismatch');
end
if abs(log(det(AS))-logdetLDLAS)>sqrt(eps)
    abs(log(det(AS))-logdetLDLAS)
    error('mismatch');
end

% grad log-det
glogdetDiff=numericalGradient(@(A)log(det(A.*sparsity)),A,1e-5);
fprintf('grad-diff logdet(AS) =\n');
disp(full(glogdetDiff))
glogdetLUAS=get_glogdetLUAS(obj);
fprintf('grad logdetLU(AS) =\n');
disp(full(glogdetLUAS))
glogdetLDLAS=get_glogdetLDLAS(obj);
fprintf('grad logdetLDL(AS) =\n');
disp(full(glogdetLDLAS))
if norm(glogdetDiff-glogdetLUAS)>sqrt(eps)
    norm(glogdetDiff-glogdetLUAS)
    error('mismatch');
end
if norm(A-A')<eps && norm(glogdetDiff-glogdetLDLAS)>sqrt(eps)
    norm(glogdetDiff-glogdetLDLAS)
    error('mismatch');
end


% trace-inv
fprintf('trace(inv(AS)) = %f\n',trace(inv(AS)));
traceinvLUAS=get_traceinvLUAS(obj);
fprintf('traceinvLU(AS) = %f\n',full(traceinvLUAS));
traceinvLDLAS=get_traceinvLDLAS(obj);
fprintf('traceinvLDL(AS) = %f\n',traceinvLDLAS);
if abs(trace(inv(AS))-traceinvLUAS)>sqrt(eps)
    error('mismatch');
end
if norm(A-A')<eps && abs(trace(inv(AS))-traceinvLDLAS)>sqrt(eps)
    error('mismatch');
end

% grad trace-inv
gtraceinvDiff=numericalGradient(@(A)trace(inv(A.*sparsity)),A,1e-5);
fprintf('grad-diff traceinv(AS) =\n');
disp(full(gtraceinvDiff))
gtraceinvLUAS=get_gtraceinvLUAS(obj);
fprintf('grad traceinvLU(AS) =\n');
disp(full(gtraceinvLUAS))
gtraceinvLDLAS=get_gtraceinvLDLAS(obj);
fprintf('grad traceinvLDL(AS) =\n');
disp(full(gtraceinvLDLAS))
if norm(gtraceinvDiff-gtraceinvLUAS)>10*sqrt(eps)
    norm(gtraceinvDiff-gtraceinvLUAS)
    error('mismatch');
end
if norm(A-A')<eps && norm(gtraceinvDiff-gtraceinvLDLAS)>sqrt(eps)
    norm(gtraceinvDiff-gtraceinvLDLAS)
    error('mismatch');
end

% inv
fprintf('inv(A) =\n');
disp(inv(AS));
mldLUAS=get_mldLUAS(obj);
invLUAS=get_invLUAS(obj);
invLDLAS=get_invLDLAS(obj);
% disp(full(mldLUAS));
% disp(full(invLUAS));
% disp(full(invLDLAS));
if norm(inv(AS)-mldLUAS)>sqrt(eps)
    error('mismatch');
end
if norm(inv(AS)-invLUAS)>sqrt(eps)
    error('mismatch');
end
if norm(A-A')<eps && norm(inv(AS)-invLDLAS)>sqrt(eps)
    error('mismatch');
end

% grad-inv
ginvLDLAS=get_ginvLDLAS(obj)

