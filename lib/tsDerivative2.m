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

    if isequal(class(ts),'Tcalculus') && length(size(ts))~=1 && ~isempty(size(ts))
        error('tsDerivative: requires times vector to be one-dimensional ([%s])\n',...
              index2str(size(ts)));
    end

    if ~isequal(class(ts),'Tcalculus') && size(ts,2)~=1
        error('tsDerivative: requires times vector to be a column vector ([%s])\n',...
              index2str(size(ts)));
    end

    if length(size(x))~=2
        error('tsDerivative2: only implemented for time series of vectors ([%s])\n',...
              index2str(size(x)));
    end

    if length(ts)~=1 && length(ts)~=size(x,2)
        error('tsDerivative2: length of sample times does not match size of input (%d,[%s])\n',...
              length(ts),index2str(size(x)));
    end

    if length(ts)>1
        if nargin<4
            if isequal(class(x),'Tcalculus')
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

        if ~isequal(class(x),'Tcalculus')
            whichtprod=@mytprod;
        else
            whichtprod=@tprod;
            if ~isequal(class(ts),'Tcalculus')
                t1=Tconstant(t1,length(t1));
                t2=Tconstant(t2,length(t2));
                t3=Tconstant(t3,length(t3));
            end
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
    else % if length(ts)>1
        if nargin<4
            invD2t=1/(ts*ts);
        end

        if ~isequal(class(x),'Tcalculus')
            whichtprod=@mytprod;
        else
            whichtprod=@tprod;
            if ~isequal(class(invD2t),'Tcalculus')
                invD2t=Tconstant(invD2t,[]);
            end
        end

        ddx=invD2t*[x(:,1)-2*x(:,2)+x(:,3),...
                    x(:,1:end-2)-2*x(:,2:end-1)+x(:,3:end),...
                    x(:,end-2)-2*x(:,end-1)+x(:,end)];
    end

    % turn into linear operation to remove number of times that x appears
    if false && isequal(class(ddx),'Tcalculus')
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
