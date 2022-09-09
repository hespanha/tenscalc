function [dx,ts]=tsDerivative(x,ts,invDts,invD2ts);
% [dx,ts]=tsDerivative(x,ts);
%
% [dx,ts]=tsDerivative(x,ts,invDts,invD2ts);
%
% Differentiates a vector-valued time series (x,ts). The time
% derivative is computed assuming that the input time-series is
% piecewise quadratic.
%
% Inputs:
%   x [n x N]  - values of the function at the given times
%                (one time per column)
%   ts [N x 1] - vector of times
%    or
%   ts [1 x 1] - (constant) sample interval
%
% the remaining parameters are needed if ts is of class Tcalculus
%   invDts  [N-1 x 1] = 1./(ts(2:end)-ts(1:end-1))
%   invD2ts [N-2 x 1] = 1./(ts(3:end)-ts(1:end-2))
%    or
%   invDts    [1 x 1] = 1/ts
%
% Output
%   dx [n x N] - derivative of the function at the given times
%                (one time per column)
%   ts [N x 1] - vector of times (equal to the corresponding input)
%
% Attention: Using this function to set constraints for MPC dynamics,
%    as in
%               tsDerivative(x,ts) == f(x,u)
%    for the system
%               dot x = f(x,u)
%    ignores that u is piecewise-constant (ZOH) and often leads to
%       1) optimal controls with high-frequency components
%       2) make the output of tsDerivative a bad approximation for
%          the derivative
%       3) overall bad MPC performance
%
%    In  view of this, it is often better to constraints of the form
%         x^+ = x + Ts * f(x,u)                [forward Euler]
%    or
%         x^+ = x + Ts * ( f(x^+,u)+f(x,u) )   [trapesiodal with ZOH for u]
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
            error('tsDerivative: requires Tcalculus times vector to be a scalars or a one-dimensional (not size=[%s])\n',...
                  index2str(size(ts)));
        end
    else
        if numel(ts)==1
            scalarTs=true;
        elseif size(ts,2)==numel(ts)
            scalarTs=false;
        else
            error('tsDerivative: requires times vector to be a column vector (not size=[%s])\n',...
                  index2str(size(ts)));
        end
    end

    xsize=size(x);
    if length(xsize)~=2
        error('tsDerivative: only implemented for time series of vectors ([%s])\n',...
              index2str(xsize));
    end

    if ~scalarTs && length(ts)~=xsize(2)
        error('tsDerivative: length of sample times does not match size of input (%d,[%s])\n',...
              length(ts),index2str(xsize));
    end

    if scalarTs

        if nargin<3
            t1=[-1.5;2;-.5]/ts;
            t2=-.5/ts;
            t3=[.5;-2;1.5]/ts;
        else
            t1=[-1.5;2;-.5]*invDts;
            t2=-.5*invDts;
            t3=[.5;-2;1.5]*invDts;
        end

        if isa(ts,'Tcalculus')
            whichtprod=@tprod;
            if ~isa(ts,'Tcalculus')
                t1=Tconstant(t1,[3]);
                t2=Tconstant(t2,[]);
                t3=Tconstant(t3,[3]);
            end
        else
            whichtprod=@mytprod;
        end

        dx=[
            reshape(whichtprod(t1,-1,x(:,1:3),[1,-1]),[xsize(1),1]),...
            t2*(x(:,1:end-2)-x(:,3:end)),...
            reshape(whichtprod(t3,-1,x(:,end-2:end),[1,-1]),[xsize(1),1])
           ];

    else % if scalarTs

        if nargin<4
            t1=[(2*ts(1)-ts(2)-ts(3))./(ts(1)-ts(3))./(ts(1)-ts(2));
                (ts(2:end-1)-ts(3:end))./(ts(1:end-2)-ts(3:end))./(ts(1:end-2)-ts(2:end-1));
                (ts(end)-ts(end-1))./(ts(end-2)-ts(end-1))./(ts(end-2)-ts(end))];

            t2=[(ts(3)-ts(1))./(ts(2)-ts(3))./(ts(1)-ts(2));
                (ts(1:end-2)+ts(3:end)-2*ts(2:end-1))./(ts(2:end-1)-ts(3:end))./(ts(1:end-2)-ts(2:end-1));
                (ts(end-2)-ts(end))./(ts(end-2)-ts(end-1))./(ts(end-1)-ts(end))];

            t3=[(ts(1)-ts(2))./(ts(1)-ts(3))./(ts(2)-ts(3));
                (ts(1:end-2)-ts(2:end-1))./(ts(1:end-2)-ts(3:end))./(ts(3:end)-ts(2:end-1));
                (ts(end-2)+ts(end-1)-2*ts(end))./(ts(end-2)-ts(end))./(ts(end)-ts(end-1))];
        else
            t1=[(2*ts(1)-ts(2)-ts(3)).*invD2ts(1)*invDts(1);
                (ts(2:end-1)-ts(3:end)).*invD2ts.*invDts(1:end-1);
                (ts(end)-ts(end-1)).*invD2ts(end).*invDts(end-1)];

            t2=[(ts(3)-ts(1)).*invDts(2).*invDts(1);
                (ts(1:end-2)+ts(3:end)-2*ts(2:end-1)).*invDts(2:end).*invDts(1:end-1);
                (ts(end-2)-ts(end)).*invDts(end-1).*invDts(end)];

            t3=[(ts(1)-ts(2)).*invD2ts(1).*invDts(2);
                (ts(2:end-1)-ts(1:end-2)).*invD2ts.*invDts(2:end);
                (2*ts(end)-ts(end-2)-ts(end-1)).*invD2ts(end).*invDts(end)];
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

        dx=[
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
        g=gradient(dx,x);
        %g=eval(str(g))
        dx=whichtprod(g,[1,2,-1,-2],x,[-1,-2]);
    end
end

function test
    % Numeric
    ts=pi/2:pi/10:4*pi+pi/2;
    x=[sin(ts);cos(ts)];
    dx1=tsDerivative(x,ts');
    dx2=tsDerivative(x,ts(2)-ts(1));

    plot(ts,x','-x',ts,dx1','-+',ts,dx2','-x');
    legend('sin','cos','d sin','d cos','d sin','d cos')

    % Symbolic
    Tvariable x [2,length(ts)]
    dx=tsDerivative(x,ts)
end