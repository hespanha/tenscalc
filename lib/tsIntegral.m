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
            error('tsIntegral: requires Tcalculus times vector to be a scalars or a one-dimensional (not size=[%s])\n',...
                  index2str(size(ts)));
        end
    else
        if numel(ts)==1
            scalarTs=true;
        elseif size(ts,2)==numel(ts)
            scalarTs=false;
        else
            error('tsIntegral: requires times vector to be a column vector (not size=[%s])\n',...
                  index2str(size(ts)));
        end
    end

    xsize=size(x);
    if length(xsize)<1 || length(xsize)>2
        error('tsIntegral: only implemented for time series of vectors ([%s])\n',...
              index2str(xsize));
    end

    if ~scalarTs && length(ts)~=xsize(end)
        error('tsIntegral: length of sample times does not match size of input (%d,[%s])\n',...
              length(ts),index2str(xsize));
    end

    if isa(ts,'Tcalculus')
        whichtprod=@tprod;
    else
        whichtprod=@mytprod;
    end

    if scalarTs
        % not very efficient to create this vector of mostly equal entries
        dt=[.5*ts;ts*ones(xsize(end)-2,1);.5*ts];
        y=whichtprod(dt,-1,x,[1:length(xsize)-1,-1]);
    else
        dt=.5*[ts(2)-ts(1);ts(3:end)-ts(1:end-2);ts(end)-ts(end-1)];
        y=whichtprod(dt,-1,x,[1:length(xsize)-1,-1]);
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
