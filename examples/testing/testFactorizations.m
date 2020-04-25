clear all;
!rm -fr tmp* @tmp*

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

gen={@cmex2compute,@class2compute};
    
N=5;

% Make a matrix symmetric as ldl does it.
sym=@(A)tril(A)+tril(A)'-diag(diag(A));

close=@(A,B,tol)norm(A(:)-B(:),inf)<=tol;

for g=1:length(gen)
    Tcalculus.clear();
    
    fprintf('*****trying %s\n',char(gen{g}));
    code=csparse();
    tc=struct();
        
    % non symmetric
    [i,j]=find(ones(N,N));
    k=randi(N*N,floor(.25*N*N),1);
    i(k)=[];
    j(k)=[];
    %[i,j]
    tc.a=Tvariable('a',length(i));
    tc.A=vec2tensor(tc.a,[N,N],[i,j]);
    
    % symmetric
    [i,j]=find(ones(N,N));
    k=randi(N*N,floor(.25*N*N),1);
    k(i(k)==j(k))=[];   % do not remove diagonal elements to not mess up ldl fact.
    i(k)=[];
    j(k)=[];
    k0=find(i==j);
    k1=find(i>j);
    is=[i(k0);i(k1);j(k1)];
    js=[j(k0);j(k1);i(k1)];
    %[is,js]
    tc.as=Tvariable('as',length(k0)+length(k1));
    as=[tc.as;tc.as(end-length(k1)+1:end)];
    tc.As=vec2tensor(as,[N,N],[is,js]);
    
    tc.ldl_d_As=ldl_d(ldl(tc.As));
    tc.ldl_l_As=ldl_l(ldl(tc.As));

    tc.lu_l_A=lu_l(lu(tc.A));
    tc.lu_u_A=lu_u(lu(tc.A));
    tc.lu_d_A=lu_d(lu(tc.A));
    
    
    declareSet(code,tc.a,'set_a');
    declareSet(code,tc.as,'set_as');
    declareGet(code,tc,'getAll');
    
    classname=gen{g}('classname',sprintf('tmp%d',g),'csparseObject',code);
        
    obj=feval(classname);

    in.a=rand(size(tc.a),1);
    set_a(obj,in.a);
    in.as=rand(size(tc.as),1);
    set_as(obj,in.as);
    out=getAll_struct(obj),
    
    fprintf('LDL factorization\n');
    full(out.As)
    full(out.ldl_d_As)
    full(out.ldl_l_As)
    
    full(out.ldl_l_As*diag(out.ldl_d_As)*out.ldl_l_As')
    if ~close(out.As,out.ldl_l_As*diag(out.ldl_d_As)*out.ldl_l_As',1e3*eps)
        
        error('error in LDL computation of %s',char(gen{g}));
    end
    
    
    fprintf('LU factorization\n');
    full(out.A),
    full(out.lu_l_A),
    full(out.lu_u_A),
    full(out.lu_d_A),
    full(out.lu_l_A*out.lu_u_A),
    if ~close(out.A,out.lu_l_A*out.lu_u_A,1e3*eps),
        
        error('error in LU computation of %s',char(gen{g}));
    end
    
    clear obj
end

