clear all;
!rm -fr tmp* @tmp*

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

gen={@cmex2compute,@class2compute};
    
N=2;
n=1;

N=12;

close=@(A,B,tol)norm(A(:)-B(:),inf)<=tol;

for g=1:length(gen)
    fprintf('*****trying with %s\n',char(gen{g}));
    
    Tcalculus.clear();
    code=csparse();
    tc=struct();

    subs=[1:3:N;
          2:3:N]';
    n=size(subs,2);
    N=max(subs,[],'all');
    
    tc.x=Tvariable('x',[size(subs,1)-1]);
    tc.x0=[0;tc.x];

    tc.diagx=diag(tc.x0);
    tc.ddiagx=full(gradient(tc.diagx,tc.x));

    tc.xx=vec2tensor(tc.x0,[N,N],subs);
    tc.dxx=full(gradient(tc.xx,tc.x));
    tc.dxxTimes=full(gradient(tc.xx.*tc.xx,tc.x));
    tc.dxxMtimes=full(gradient(tc.xx*tc.xx,tc.x));
    tc.dxxPlus=full(gradient(tc.xx+tc.xx,tc.x));

    % Creates NxN diagonal matrix
    tc.v1=Tvariable('v1',[N]);
    tc.diagv=diag(tc.v1,0);
    tc.diagvp1=diag(tc.v1,1);
    tc.diagvm1=diag(tc.v1,-1);
    
    % Creates NxN lower triangular matrix
    [i,j]=find(ones(N));
    k=find(i>=j);
    tc.v2=Tvariable('v2',length(k));
    tc.trilv=vec2tensor(tc.v2,[N,N],[i(k),j(k)]);
    
    % Creates NxN lower triangular matrix
    [i,j]=find(ones(N));
    kl=find(i>j);
    k0=find(i==j);
    tc.v3=Tvariable('v3',length(k0)+length(kl));
    tc.symv=vec2tensor([tc.v3;tc.v3(1:length(kl))],[N,N],[i(kl),j(kl);i(k0),j(k0);j(kl),i(kl)]);

    declareSet(code,tc.x,'set_x');
    declareSet(code,tc.v1,'set_v1');
    declareSet(code,tc.v2,'set_v2');
    declareSet(code,tc.v3,'set_v3');
    f=fields(tc);
    for i=1:length(fields(tc))
        disp(getfield(tc,f{i}))
        declareGet(code,getfield(tc,f{i}),sprintf('get_%s',f{i}));
    end
    
    classname=gen{g}('classname',sprintf('tmp%d',g),'csparseObject',code,'profiling',true);

    obj=feval(classname);
    
    in.x=randn(length(tc.x),1);
    in.v1=randn(length(tc.v1),1);
    in.v2=randn(length(tc.v2),1);
    set_x(obj,in.x);
    set_v1(obj,in.v1);
    set_v2(obj,in.v2);
    set_v3(obj,in.v2);
    out=struct();
    for i=1:length(fields(tc))
        out=setfield(out,f{i},feval(sprintf('get_%s',f{i}),obj));
        %disp(getfield(out,f{i}))        
    end
    out,;
    
    if ~close(out.diagv,diag(out.v1),eps)
        full(out.diagv),full(diag(out.v1))
        error('mismatch error in diagv');
    end
    
    if ~close(out.diagvp1,diag(out.v1,1),eps)
        full(out.diagvp1),full(diag(out.v1,1))
        error('mismatch error in diagv');
    end
    
    if ~close(out.diagvm1,diag(out.v1,-1),eps)
        full(out.diagvm1),full(diag(out.v1,-1))
        error('mismatch error in diagv');
    end
    
    if ~close(triu(out.trilv,1),zeros(size(out.trilv)),eps)
        full(out.trilv),
        error('mismatch error in trilv');
    end
    
    if ~close(out.symv,out.symv',eps)
        full(out.symv),
        error('mismatch error in symv');
    end
    
    fun.diagx=@(x)diag([0;x]);
    num.ddiagx=numericalGradient(fun.diagx,in.x,eps.^.5);
    if ~close(out.ddiagx,num.ddiagx,sqrt(eps))
        out.ddiagx,num.ddiagx
        error('mismatch error in ddiagx');
    end
    
    fun.xx=@(x)sparse(subs(2:end,1),subs(2:end,2),x,N,N);
    num.dxx=numericalGradient(fun.xx,in.x,eps.^.5);
    if ~close(out.dxx,num.dxx,sqrt(eps))
        out.dxx,num.dxx
        error('mismatch error in dxx');
    end
    
    num.dxxTimes=numericalGradient(@(x)fun.xx(x).*fun.xx(x),in.x,eps^.25);
    if ~close(out.dxxTimes,num.dxxTimes,sqrt(eps))
        out.dxxTimes,num.dxxTimes
        error('mismatch error in dxxTimes');
    end
    
    num.dxxMtimes=numericalGradient(@(x)fun.xx(x)*fun.xx(x),in.x,eps^.25);
    if ~close(out.dxxMtimes,num.dxxMtimes,sqrt(eps))
        out.dxxMtimes,num.dxxMtimes
        error('mismatch error in dxxMtimes');
    end
    
    num.dxxPlus=numericalGradient(@(x)fun.xx(x)+fun.xx(x),in.x,eps^.25);
    if ~close(out.dxxPlus,num.dxxPlus,sqrt(eps))
        out.dxxPlus,num.dxxPlus
        error('mismatch error in dxxPlus');
    end
    type profileView.profile
end

