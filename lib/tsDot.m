function [y,ts]=tsDot(x1,x2,ts);
% [y,ts]=tsDot(x1,x2,ts);
%
% Computes a scalar-valued time-series (y,ts) that represents the
% dot product of two n-vector time-series (x1,ts) and (x2,ts)
%
% Inputs:
%   x1 [n x N] - values of the 1st time series at the given times
%                (one time per column)
%   x2 [n x N] - values of the 2nd time series at the given times
%                (one time per column)
%   ts [N x 1] - vector of times
%                (for TC variables size(ts)=N)
%
% Output
%   y [1 x N]  - values of the dot product at the given times
%                (for TC variables size(y)=N)
%                (one time per column)
%   ts [N x 1] - vector of times (equal to the corresponding input)
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    if nargin<3
        ts=(1:size(x1,2))';
    else
        ts=ts(:);
    end

    if length(size(x1))~=2 || length(size(x2))~=2
        error('tsDot: inputs must be time series of vectors ([%s],[%s])\n',...
              index2str(size(x1)),index2str(size(x2)));
    end

    if length(ts)~=size(x1,2) || length(ts)~=size(x2,2)
        error('tsDot: length of sample times does not match size of inputs (%d,[%s],[%s])\n',...
              length(ts),index2str(size(x1)),index2str(size(x2)));
    end

    if isequal(class(x1),'Tcalculus') || isequal(class(x2),'Tcalculus')
        y=tprod(x1,[-1,1],x2,[-1,1]);
    else
        %y=sum(x1.*x2,1);
        y=mytprod(x1,[-1,1],x2,[-1,1])';
    end
end

function test
    % Numeric
    ts=pi/2:pi/10:4*pi+pi/2;
    x=[sin(ts).*cos(ts);cos(ts).*cos(ts);sin(ts)];
    y=tsDot(x,x,ts);

    plot(ts,x','-x',ts,y','-+');
    legend('x1','x2','x3','y')

    % Symbolic
    Tvariable x [2,length(ts)]
    y=tsDot(x,x,ts)
end