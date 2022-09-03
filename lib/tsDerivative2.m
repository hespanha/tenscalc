function [ddx,ts]=tsDerivative2(x,ts,invDts,invD2ts)
% [ddx,ts]=tsDerivative2(x,ts);
%
% [ddx,ts]=tsDerivative2(x,ts,invDts,invD2ts)
%
% Differentiates twice a vector-valued time series (x,ts). The time
% derivative is computed assuming that the input time-series is
% piecewise quadratic.
%
% Inputs:
%   x [n x N]   - values of the function at the given times
%                 (one time per column)
%   ts [N x 1]  - vector of times
%    or
%   ts [1 x 1] - (constant) sample interval
%
% the remaining parameters are needed if ts is of class Tcalculus
%   invDts  [N-1 x 1] = 1./(ts(2:end)-ts(1:end-1))
%   invD2ts [N-2 x 1] = 1./(ts(3:end)-ts(1:end-2))
%    or
%   invD2ts   [1 x 1] = 1/ts
%
% Output
%   ddx [n x N] - 2nd derivative of the function at the given times
%                  (one time per column)
%   ts [N x 1]  - vector of times (equal to the corresponding input)
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    if isa(ts,'Tcalculus')
        if isempty(size(ts))
            scalarTs=true;
        elseif length(size(ts))==1
            scalarTs=false;
        else
            error('tsDerivative2: requires Tcalculus times vector to be a scalars or a one-dimensional (not size=[%s])\n',...
                  index2str(size(ts)));
        end
    else
        if numel(ts)==1
            scalarTs=true;
        elseif size(ts,2)==numel(ts)
            scalarTs=false;
        else
            error('tsDerivative2: requires times vector to be a column vector (not size=[%s])\n',...
                  index2str(size(ts)));
        end
    end

    xsize=size(x);
    if length(xsize)~=2
        error('tsDerivative2: only implemented for time series of vectors ([%s])\n',...
              index2str(xsize));
    end

    if ~scalarTs && length(ts)~=xsize(2)
        error('tsDerivative2: length of sample times does not match size of input (%d,[%s])\n',...
              length(ts),index2str(xsize));
    end

    if scalarTs

        if nargin<4
            invD2t=1/(ts*ts);
        end

        if isa(ts,'Tcalculus')
            whichtprod=@tprod;
            if ~isequal(class(invD2t),'Tcalculus')
                invD2t=Tconstant(invD2t,[]);
            end
        else
            whichtprod=@mytprod;            
        end

        ddx=invD2t*[x(:,1)-2*x(:,2)+x(:,3),...
                    x(:,1:end-2)-2*x(:,2:end-1)+x(:,3:end),...
                    x(:,end-2)-2*x(:,end-1)+x(:,end)];
    else % if scalarTs

        if nargin<4
            if isa(ts,'Tcalculus')
                two=Tconstant(2,1);
                ones2=Tones(length(ts)-2);
            else
                two=2;
                ones2=ones(length(ts)-2,1);
            end

            t1=[two./(ts(1)-ts(2))./(ts(1)-ts(3));
                2*ones2./(ts(1:end-2)-ts(2:end-1))./(ts(1:end-2)-ts(3:end));
                two./(ts(end-2)-ts(end-1))./(ts(end-2)-ts(end))];

            t2=[two./(ts(2)-ts(1))./(ts(2)-ts(3));
                2*ones2./(ts(2:end-1)-ts(1:end-2))./(ts(2:end-1)-ts(3:end));
                two./(ts(end-1)-ts(end-2))./(ts(end-1)-ts(end))];

            t3=[two./(ts(1)-ts(3))./(ts(2)-ts(3));
                2*ones2./(ts(1:end-2)-ts(3:end))./(ts(2:end-1)-ts(3:end));
                two./(ts(end-2)-ts(end))./(ts(end-1)-ts(end))];
        else
            t1=[2*invDts(1).*invD2ts(1);
                2*invDts(1:end-1).*invD2ts;
                2*invDts(end-1).*invD2ts(end)];

            t2=[-2*invDts(1).*invDts(2);
                -2*invDts(1:end-1).*invDts(2:end);
                -2*invDts(end-1).*invDts(end)];

            t3=[2*invD2ts(1).*invDts(2);
                2*invD2ts.*invDts(2:end);
                2*invD2ts(end).*invDts(end)];
        end

        if isa(ts,'Tcalculus')
            whichtprod=@tprod;
            if ~isequal(class(ts),'Tcalculus')
                t1=Tconstant(t1,length(t1));
                t2=Tconstant(t2,length(t2));
                t3=Tconstant(t3,length(t3));
            end
        else
            whichtprod=@mytprod;            
        end

        ddx=[
            ...%x(:,1)*t1(1)+x(:,2)*t2(1)+x(:,3)*t3(1),...
            whichtprod(t1(1),[2],x(:,1),[1,2])+...
            whichtprod(t2(1),[2],x(:,2),[1,2])+...
            whichtprod(t3(1),[2],x(:,3),[1,2]),...
            whichtprod(t1(2:end-1),[2],x(:,1:end-2),[1,2])+...
            whichtprod(t2(2:end-1),[2],x(:,2:end-1),[1,2])+...
            whichtprod(t3(2:end-1),[2],x(:,3:end),[1,2]),...
            ...%x(:,end-2)*t1(end)+x(:,end-1)*t2(end)+x(:,end)*t3(end)...
            whichtprod(t1(end),[2],x(:,end-2),[1,2])+...
            whichtprod(t2(end),[2],x(:,end-1),[1,2])+...
            whichtprod(t3(end),[2],x(:,end),[1,2])...
            ];
    end

    % turn into linear operation to remove number of times that x appears
    if false && isa(ts,'Tcalculus')
        g=gradient(ddx,x);
        %g=eval(str(g))
        ddx=whichtprod(g,[1,2,-1,-2],x,[-1,-2]);
    end
end

function test
    % Numeric
    ts=pi/2:pi/10:4*pi+pi/2;
    x=[sin(ts);cos(ts)];
    dx=tsDerivative2(x,ts);

    plot(ts,x','-x',ts,dx','-+');
    legend('sin','cos','d2 sin','d2 cos')

    % Symbolic
    Tvariable x [2,length(ts)]
    dx=tsDerivative2(x,ts)
end
