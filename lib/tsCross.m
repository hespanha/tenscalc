function [y,ts]=tsCross(x1,x2,ts);
% [y,ts]=tsCross(x1,x2,ts);
%
% Computes a scalar-valued time-series (y,ts) that represents the
% cross product of two n-vector time-series (x1,ts) and (x2,ts)
%
% Inputs:
%   x1 [3 x N] - values of the 1st time series at the given times
%                (one time per column)
%   x2 [3 x N] - values of the 2nd time series at the given times
%                (one time per column)
%   ts [N x 1] - vector of times
%                (for TC variables size(ts)=N)
%
% Output
%   y  [3 x N] - values of the dot product at the given times
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

if length(size(x1))~=2 || length(size(x2))~=2 || size(x1,1)~=3 || size(x2,1)~=3
    error('tsCross: inputs must be time series of 3-vectors ([%s],[%s])\n',...
        index2str(size(x1)),index2str(size(x2)));
end

if length(ts)~=size(x1,2) || length(ts)~=size(x2,2)
    error('tsCross: length of sample times does not match size of inputs (%d,[%s],[%s])\n',...
        length(ts),index2str(size(x1)),index2str(size(x2)));
end

if isequal(class(x1),'Tcalculus') || isequal(class(x2),'Tcalculus')
    % y=[x1(2,:).*x2(3,:)-x1(3,:).*x2(2,:);
    %    x1(3,:).*x2(1,:)-x1(1,:).*x2(3,:);
    %    x1(1,:).*x2(2,:)-x1(2,:).*x2(1,:);];
    y=x1([2,3,1],:).*x2([3,1,2],:)-x1([3,1,2],:).*x2([2,3,1],:);
else
    y=cross(x1,x2,1);
end
end

function test
% Numeric
ts=pi/2:pi/10:4*pi+pi/2;
x=[sin(ts).*cos(ts);cos(ts).*cos(ts);sin(ts)];
y=tsCross(x,x,ts);

plot(ts,x','-x',ts,y','-+');
legend('x1','x2','x3','y1','y2','y3')

% Symbolic
Tvariable x [3,length(ts)]
y=tsCross(x,x,ts)
end