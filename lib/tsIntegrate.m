function [intX,ts]=tsIntegrate(x,x0,ts,method);
% [intX,ts]=tsIntegrate(x,x0,ts,method);
%
% Integrates a vector-valued time series (x,ts). The time
% derivative is computed assuming that the input time-series is
% piecewise quadratic.
%
% Inputs:
%   x [n x N]  - values of the function at the given times
%                (one time per column)
%   x0 [n x 1] - initial value
%   ts [N x 1] - vector of times
%    or
%   ts [1 x 1] - (constant) sample interval
%
%   method     - integration method (euler, trapesoidal, invtsDerivative)
%
%
% Output
%   intX [n x N] - integral of the function at the given times
%                (one time per column)
%   ts [N x 1] - vector of times (equal to the corresponding input)
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

if nargin<4
    method='trapesoidal';
    method='euler';
end

if isequal(class(ts),'Tcalculus') && length(size(ts))~=1 && ~isempty(size(ts))
    error('tsIntegrate: requires times vector to be scalars or one-dimensional ([%s])\n',...
        index2str(size(ts)));
end

if ~isequal(class(ts),'Tcalculus') && size(ts,2)~=1
    error('tsIntegrate: requires times vector to be a column vector ([%s])\n',...
        index2str(size(ts)));
end

if length(size(x))~=2
    error('tsIntegrate: only implemented for time series of vectors ([%s])\n',...
        index2str(size(x)));
end

if length(ts)~=1 && length(ts)~=size(x,2)
    error('tsIntegrate: length of sample times does not match size of input (%d,[%s])\n',...
        length(ts),index2str(size(x)));
end

if length(ts)>1
    switch (method)
        case 'euler'
            intX=[x0,x0+cumsum(repmat(diff(ts'),size(x,1),1).*x(:,1:end-1),2)];
        otherwise
            error('not implemented');
    end
else % if length(ts)>1
    switch (method)
        case 'euler'
            % causal but not quite symmetric
            intX=[x0,x0+ts'*cumsum(x(:,1:end-1),2)];
        case 'trapesoidal'
            % not quite causal, but more accurate
            intX=[x0,x0+(ts/2)*cumsum(x(:,1:end-1)+x(:,2:end),2)];
        case 'invtsDerivative'
            % inverse of tsDerivative (up to an errorin sample
            % before last) -- introduces high frequency
            %x
            intXodd=cumsum([x0,2*ts*x(:,2:2:end-1)],2);
            intXeven=cumsum([x0+ts/2*(x(:,1)+x(:,2)),2*ts*x(:,3:2:end-1)],2);
            intX=zeros(size(x));
            intX(:,1:2:end)=intXodd;
            intX(:,2:2:end)=intXeven;
            intX(:,end)=(2*intX(:,end-1)-.5*intX(:,end-2)+ts*x(:,end))/1.5;
            %intX
        otherwise
            error('unknown mehtod %s',method);
    end
end
end

function test
% Numeric
ts=pi/2:pi/10:4*pi+pi/2;
x=[sin(ts);cos(ts)];
intX1=tsIntegrate(x,[0;0],ts');
intX2=tsIntegrate(x,[0;0],ts(2)-ts(1));

plot(ts,x','-x',ts,intX1','-+',ts,intX2','-x');
legend('sin','cos','d sin','d cos','d sin','d cos')

% Symbolic
Tvariable x [2,length(ts)]
intX=tsIntegrate(x,[0;0],ts)
end