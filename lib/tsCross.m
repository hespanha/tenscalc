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
% Copyright 2012-2017 Joao Hespanha

% This file is part of Tencalc.
%
% TensCalc is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% TensCalc is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.

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