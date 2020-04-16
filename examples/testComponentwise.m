clear all;
!rm -fr tmp* @tmp*

gen={@cmex2compute,@class2compute};
    
%% Test componentwise

if 0
    N=2;
    
    Tvariable a [N];
    %A=diag(a);
    A=a;
    
    cA=componentwise(A,@cos,'cos(%s)');
    sA=componentwise(A,@sin,'sin(%s)',@(x)componentwise(x,@cos,'cos(%s)',@(x)-componentwise(x,@sin,'sin(%s)')));
    dA=componentwise(A,@(x)2*x,'2*%s',@(x)2*Tones(size(x)),@(x)Tzeros(size(x)));
    
    dsA=gradient(sA,a);
    d2sA=gradient(dsA,a);
    
    ddA=gradient(dA,a);
    d2dA=gradient(ddA,a);
    
    for i=1:length(gen)
        code=csparse();
        declareSet(code,a,'set_a');
        declareGet(code,{A,cA,sA,dA,full(dsA),full(d2sA),full(ddA),full(d2dA)},'getAll');
        
        classname=gen{i}('classname','tmp','csparseObject',code);
        
        str(code,1);
        
        obj=feval(classname);
        set_a(obj,(2:N+1)');
        [A_,cA_,sA_,dA_,dsA_,d2sA_,ddA_,d2dA_]=getAll(obj)
        
    end
end

%% Test functions derived from componentwise

% function, derivative
fun={@normpdf,@(x)-x.*exp(-x.^2/2)/sqrt(2*pi),@(x)(x.^2-1).*exp(-x.^2/2)/sqrt(2*pi);
     @atan,@(x)1./(1+x.^2),@(x)-2*x./(1+x.^2).^2;
     @sqr,@(x)2*x,@(x)2*ones(size(x));
     @cube,@(x)3*x.^2,@(x)6*x;
     @sqrt,@(x).5./sqrt(x),@(x)-.25.*x^(-1.5);
     @exp,@exp,@exp;
     @log,@(x)1./x,@(x)-1./x.^2;
     @sin,@cos,@(x)-sin(x);
     @cos,@(x)-sin(x),@(x)-cos(x);
     @tan,@(x)sec(x).^2,@(x)2*sec(x).^2.*tan(x);
    };

for f=1:size(fun,1)
    for i=1:length(gen)
        fprintf('*****trying %s with %s\n',char(fun{f,1}),char(gen{i}));
        code=csparse();
        
        Tvariable x [2];
        y=fun{f,1}(x);
        dy=gradient(y,x);
        d2y=gradient(dy,x);
        declareSet(code,x,'set_x');
        declareGet(code,{y,dy,full(d2y)},'get');
    
        classname=gen{i}('classname',sprintf('tmp%d%d',f,i),'csparseObject',code);
        
        obj=feval(classname);
        x=rand(2,1);
        set_x(obj,x);
        [y_,dy_,d2y_]=get(obj);
        
        y=fun{f,1}(x);
        dy=diag(fun{f,2}(x));
        d2y=zeros(2,2,2);
        d2y(1,1,1)=fun{f,3}(x(1));
        d2y(2,2,2)=fun{f,3}(x(2));
               
        if norm(y-y_,inf)>10*eps
            y,y_,
            error('error in computation of %s',char(fun{f,1}));
        end
        if norm(dy-dy_,inf)>10*eps
            dy,dy_,
            error('error in computation of %s''',char(fun{f,1}));
        end
        if any(abs(d2y-d2y_)>10*eps,'all')
            d2y,d2y_,
            error('error in computation of %s''''',char(fun{f,1}));
        end
        clear obj
    end
end

