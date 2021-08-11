% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all;
!rm -fr tmp* @tmp*

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

gen={@class2compute,@cmex2compute};
    
N=2;
n=1;

N=5;
n=3;

% Make a matrix symmetric as ldl does it.
sym=@(A)tril(A)+tril(A)'-diag(diag(A));

in.A=randn(N,N);
if det(in.A)<0
    in.A=-in.A;
end
%in.A=diag(2:3);
in.As=sym(in.A);
in.B=randn(N,n);


% matrices needed for the derivatives
dDiag=zeros(N,N,N);
dDiag(1:N*N+N+1:end)=1;
dTranspose=zeros(N,N,N,N);
dInverse=zeros(N,N,N,N);
dInverses=zeros(N,N,N,N);
for i=1:N
    for j=1:N
        dTranspose(i,j,j,i)=1;
        eij=zeros(N,N);
        eij(i,j)=1;
        dInverse(:,:,i,j)=-inv(in.A)*eij*inv(in.A);
        dInverses(:,:,i,j)=-inv(in.As)*eij*inv(in.As);
    end
end


% function, derivative

close=@(A,B,tol)norm(A(:)-B(:),inf)<=tol;
samesign=@(A,B,tol)norm(sign(sort(A(:)))-sign(sort(B(:))),inf)<=tol;

fixdiff=@(A,B,tol)norm2(tril(A)-tril(B,-1)/2-diag(diag(B)))<=tol;
fixdiff1=@(A,B,tol)tril(A)-tril(B,-1)/2-diag(diag(B))

fun={;
     close,    @(A,As,B)mldivide(lu(A),B),   @(A,As,B)mldivide(A,B),  @(A,As,B)tprod(dInverse,[1,-1,3,4],B,[-1,2]),     @(A,As,B)zeros(N,n,N,N);
     close,    @(A,As,B)mldivide(ldl(As),B), @(A,As,B)mldivide(As,B), @(A,As,B)zeros(N,n,N,N),                     @(A,As,B)tprod(dInverses,[1,-1,3,4],B,[-1,2]);
     close,    @(A,As,B)mldivide(ldl(A),B),  @(A,As,B)mldivide(sym(A),B), NaN,                     NaN;
     close,    @(A,As,B)det(lu(A)),          @(A,As,B)det(A),         @(A,As,B)det(A)*inv(A'), @(A,As,B)zeros(N,N);
     close,    @(A,As,B)det(ldl(As)),        @(A,As,B)det(As),        @(A,As,B)zeros(N,N),     @(A,As,B)det(As)*inv(As');
     close,    @(A,As,B)traceinv(lu(A)),     @(A,As,B)traceinv(A),    @(A,As,B)-inv(A')^2,     @(A,As,B)zeros(N,N);
     close,    @(A,As,B)traceinv(ldl(As)),   @(A,As,B)traceinv(As),   @(A,As,B)zeros(N,N),     @(A,As,B)-inv(As)^2;
     fixdiff,    @(A,As,B)traceinv(ldl(A)),    @(A,As,B)traceinv(sym(A)),@(A,As,B)-inv(As)^2,    @(A,As,B)zeros(N,N);
     close,    @(A,As,B)logdet(lu(A)),       @(A,As,B)logdet(A),      @(A,As,B)inv(A'),        @(A,As,B)zeros(N,N);
     close,    @(A,As,B)logdet(ldl(As)),     @(A,As,B)logdet(As),     @(A,As,B)zeros(N,N),     @(A,As,B)inv(As);
     fixdiff,  @(A,As,B)logdet(ldl(A)),      @(A,As,B)logdet(sym(A)), @(A,As,B)inv(As),        @(A,As,B)zeros(N,N);
     close,    @(A,As,B)diag(A),             [],                      @(A,As,B)dDiag,          @(A,As,B)zeros(N,N,N);
...%     close,    @(A,As,B)diag(A,1),           [],                      @(A,As,B)dDiag,          @(A,As,B)zeros(N,N,N);
...%     close,    @(A,As,B)diag(A,-1),          [],                      @(A,As,B)dDiag,          @(A,As,B)zeros(N,N,N);
     close,    @(A,As,B)trace(A),            [],                      @(A,As,B)eye(N),         @(A,As,B)zeros(N,N);
     close,    @(A,As,B)transpose(A),        [],                      @(A,As,B)dTranspose,     @(A,As,B)zeros(N,N,N,N);
     close,    @(A,As,B)ctranspose(A),       [],                      @(A,As,B)dTranspose,     @(A,As,B)zeros(N,N,N,N);
     close,    @(A,As,B)inv(lu(A)),          @(A,As,B)inv(A),         @(A,As,B)dInverse,       @(A,As,B)zeros(N,N,N,N);
     close,    @(A,As,B)inv(ldl(As)),        @(A,As,B)inv(As),        @(A,As,B)zeros(N,N,N,N), @(A,As,B)dInverses;
    };

for f=1:size(fun,1)
    for g=1:length(gen)
        Tcalculus.clear();
        
        fprintf('*****trying %s with %s\n',char(fun{f,2}),char(gen{g}));
        code=csparse();
        tc=struct();
        
        tc.A=Tvariable('A',[N,N]);
        tc.As=Tvariable('As',[N,N]);
        tc.B=Tvariable('B',[N,n]);
        
        tc.y=fun{f,2}(tc.A,tc.As,tc.B);
        if isa(fun{f,4},'function_handle')
            tc.dy=gradient(tc.y,tc.A);
            tc.dys=gradient(tc.y,tc.As);
        else
            tc.dy=Tzeros([]);
            tc.dys=Tzeros([]);
        end
        declareSet(code,tc.A,'set_A');
        declareSet(code,tc.As,'set_As');
        declareSet(code,tc.B,'set_B');
        declareGet(code,{tc.y,full(tc.dy),full(tc.dys)},'get');
    
        classname=gen{g}('classname',sprintf('tmp%d%d',f,i),'csparseObject',code);
        
        obj=feval(classname);
          
        set_A(obj,in.A);
        set_As(obj,in.As);
        set_B(obj,in.B);
        [out.y,out.dy,out.dys]=get(obj);
        
        if isempty(fun{f,3})
            fun{f,3}=fun{f,2};
        end
        y=fun{f,3}(in.A,in.As,in.B);
               
        if ~close(out.y,y,1e3*eps)
            y,out.y,
            error('error in computation of %s',char(fun{f,2}));
        end
        
        
        
        if isa(fun{f,4},'function_handle')
            %% numerical
            h=sqrt(eps);
            if numel(y)==1
                dy=zeros(N,N);
                dys=zeros(N,N);
                for i=1:N
                    for j=1:N
                        eij=zeros(N,N);
                        eij(i,j)=h;
                        dy(i,j)=(fun{f,3}(in.A+eij,in.As,in.B)-y)/h;
                        dys(i,j)=(fun{f,3}(in.A,in.As+eij,in.B)-y)/h;
                    end
                end
            elseif numel(y)==N
                dy=zeros(N,N,N);
                dys=zeros(N,N,N);
                for i=1:N
                    for j=1:N
                        eij=zeros(N,N);
                        eij(i,j)=h;
                        dy(:,i,j)=(fun{f,3}(in.A+eij,in.As,in.B)-y)/h;
                        dys(:,i,j)=(fun{f,3}(in.A,in.As+eij,in.B)-y)/h;
                    end
                end
            elseif numel(y)==N*N
                dy=zeros(N,N,N,N);
                dys=zeros(N,N,N,N);
                for i=1:N
                    for j=1:N
                        eij=zeros(N,N);
                        eij(i,j)=h;
                        dy(:,:,i,j)=(fun{f,3}(in.A+eij,in.As,in.B)-y)/h;
                        dys(:,:,i,j)=(fun{f,3}(in.A,in.As+eij,in.B)-y)/h;
                    end
                end
            elseif length(size(y))==2
                [m1,m2]=size(y);
                dy=zeros(m1,m2,N,N);
                dys=zeros(m1,m2,N,N);
                for i=1:N
                    for j=1:N
                        eij=zeros(N,N);
                        eij(i,j)=h;
                        dy(:,:,i,j)=(fun{f,3}(in.A+eij,in.As,in.B)-y)/h;
                        dys(:,:,i,j)=(fun{f,3}(in.A,in.As+eij,in.B)-y)/h;
                    end
                end
            end
            if ~fun{f,1}(out.dy,dy,1e3*sqrt(eps))
                dy,out.dy,
                error('error in computation of %s',char(fun{f,2}));
            end
            if ~fun{f,1}(out.dys,dys,1e3*sqrt(eps))
                dys,out.dys,
                error('error in computation of %s',char(fun{f,2}));
            end
            
            %% exact
            dy=fun{f,4}(in.A,in.As,in.B);
            if ~close(out.dy,dy,1e3*eps)
                dy,out.dy,
                error('error in computation of %s',char(fun{f,2}));
            end
            dys=fun{f,5}(in.A,in.As,in.B);
            if ~close(out.dys,dys,1e3*eps)
                dys,out.dys,
                error('error in computation of %s',char(fun{f,2}));
            end
        end
        clear obj
    end
end

