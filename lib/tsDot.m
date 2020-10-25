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