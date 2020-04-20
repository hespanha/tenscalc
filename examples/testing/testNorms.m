clear all;
%!rm -fr tmp* @tmp*

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

gen={@cmex2compute,@class2compute};

m=3;
n=5;

in.x=randn(n,1);
in.xx=randn(n,m);

close=@(A,B,tol)norm(A(:)-B(:),inf)<=tol;

%% vectors norms
%   apply to mat, gradient for vect, gradient for mat
fun={;
     false, true,  true, @(x)norm(x,.5), [];
     true, true,  false, @(x)norm(x,1), [];
     true, false, false, @(x)norm(x,inf), [];
     true, true,  true, @(x)norm2(x), [];
     false, true,  true, @(x)norm(x), [];
     false, true,  true, @(x)norm(x,2), [];
    };

for f=1:size(fun,1)
    for g=1:length(gen)
        fprintf('*****trying %s with %s\n',char(fun{f,4}),char(gen{g}));

        Tcalculus.clear();
        code=csparse();
        tc=struct();
        
        tc.x=Tvariable('x',[n]);
        tc.xx=Tvariable('xx',[n,m]);
        
        tc.y=fun{f,4}(tc.x);
        if fun{f,1}
            tc.yy=fun{f,4}(tc.xx);
        else
            tc.yy=Tzeros(0);
        end
        if fun{f,2}
            tc.dy=gradient(tc.y,tc.x);
            tc.d2y=gradient(tc.dy,tc.x);
        else
            tc.dy=Tzeros(n);
            tc.d2y=Tzeros(n,n);
        end
        if fun{f,3}
            tc.dyy=gradient(tc.yy,tc.xx);
            tc.d2yy=gradient(tc.dyy,tc.xx);
        else
            tc.dyy=Tzeros(n);
            tc.d2yy=Tzeros(n,n);
        end

        declareSet(code,tc.x,'set_x');
        declareSet(code,tc.xx,'set_xx');
        declareGet(code,{tc.y,full(tc.dy),full(tc.d2y),tc.yy,full(tc.dyy),full(tc.d2yy)},'get');
    
        classname=gen{g}('classname',sprintf('tmp%d%d',f,i),'csparseObject',code);
        
        obj=feval(classname);
          
        set_x(obj,in.x);
        set_xx(obj,in.xx);
        [out.y,out.dy,out.d2y,out.yy,out.dyy,out.d2yy]=get(obj);
        
        if isempty(fun{f,5})
            fun{f,5}=fun{f,4};
        end
        y=fun{f,5}(in.x);
        if fun{f,1}
            yy=fun{f,5}(in.xx);
        end
        
        if ~close(out.y,y,1e3*eps)
            y,out.y,
            error('error in computation of %s for vector',char(fun{f,4}));
        end
        
        if fun{f,1} && ~close(out.yy,yy,1e3*eps)
            yy,out.yy,
            error('error in computation of %s for matrix',char(fun{f,4}));
        end
        
        if fun{f,2}
            %% numerical differentiation
            h1=eps^.25;
            h2=eps^.25;
            
            dy=zeros(n,1);
            d2y=zeros(n,n);
            for i=1:n
                ei=zeros(n,1);
                ei(i)=h1;
            dy(i)=(fun{f,5}(in.x+ei)-fun{f,5}(in.x))/h1;                
            for j=1:n
                ej=zeros(n,1);
                ej(j)=h2;
                dyj=(fun{f,5}(in.x+ei+ej)-fun{f,5}(in.x+ej))/h1;                
                d2y(i,j)=(dyj-dy(i))/h2;
            end
            end
            if ~close(out.dy,dy,20*h1)
                dy,out.dy,
                error('error in computation of 1st derivative of %s',char(fun{f,4}));
            end
            if ~close(out.d2y,d2y,100*h2)
                d2y,out.d2y,
                error('error in computation of 2nd derivative of %s',char(fun{f,4}));
            end            
        end
        clear obj
    end
end

