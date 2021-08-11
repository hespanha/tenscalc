function grad=numericalGradient(f,x,h)
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.
    
    sizef=size(f(x));
    sizex=size(x);
    sizegrad=[sizef,sizex];
    if isequal(sizef,[1,1]) && length(sizex)<=2
        grad=nan(sizex);
        for j1=1:size(x,1)
            for j2=1:size(x,2)
                hj=zeros(sizex);
                hj(j1,j2)=h/2;
                grad(j1,j2)=(f(x+hj)-f(x-hj))/h;
            end
        end
    elseif length(sizef)<=2 && sizex(2)<=1
        grad=nan([sizef,sizex(1)]);
        for j1=1:size(x,1)
            hj=zeros(sizex);
            hj(j1)=h/2;
            grad(:,:,j1)=(f(x+hj)-f(x-hj))/h;
        end
        if sizef(2)==1
            sizegrad(2)=[]
            reshape(grad,sizegrad);
        end
    else
        sizef,sizex,
        error('numerical Gradient not fully implemented');
    end
end

