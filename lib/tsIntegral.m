function y=tsIntegral(x,ts);
% y=tsIntegral(x,ts)
%
% Integrates a vector-valued time series using the trapesoidal method.
%
% Inputs:
%   x [n x N]  - values of the function at the given times
%                (one time per column)
%   ts [N x 1] - vector of times
%
% Output
%   y  [n x 1] - integral of the function from ts(1) to ts(N)
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

    ts=ts(:);

    osize=size(x);
    if length(osize)<1 || length(osize)>2 
        error('tsIntegral: only implemented for time series of vectors ([%s])\n',...
              index2str(osize));
    end
    
    if length(ts)>1 && length(ts)~=osize(end)
        error('tsIntegral: length of sample times does not match size of input (%d,[%s])\n',...
              length(ts),index2str(osize));
    end
    
    if ~isequal(class(x),'Tcalculus') 
        whichtprod=@mytprod;
    else
        whichtprod=@tprod;
    end

    if length(ts)>1
        dt=.5*[ts(2)-ts(1);ts(3:end)-ts(1:end-2);ts(end)-ts(end-1)];
        y=whichtprod(dt,-1,x,[1:length(osize)-1,-1]);
    else
        % not very efficient to create this vector of mostly equal entries
        dt=[.5*ts;ts*ones(osize(end)-2,1);.5*ts];
        y=whichtprod(dt,-1,x,[1:length(osize)-1,-1]);
    end
    
end

function test
    % Numeric
    ts=pi/2:pi/10:4*pi+pi/2;
    x=[sin(ts);cos(ts)];

    plot(ts,x','-x');
    legend('sin','cos')

    y=tsIntegral(x,ts)
    
    % Symbolic
    Tvariable x [2,length(ts)]
    y=tsIntegral(x,ts)

end
