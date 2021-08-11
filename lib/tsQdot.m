function [y,ts]=tsQdot(q1,q2,ts);
% [y,ts]=tsQdot(q1,q2,ts);
%
% Computes a 4-vector time-series (y,ts) that represents the product
% of two quaternions (q1,ts), (q2,ts). Specifically,
%    y = q1 x q2
% If any of the input quaternions is a 3-vector, it is assumed to be a
% pure quaternion.

% Inputs:
%   q1 [4 x N]  - values of the quaternion at the given times
%                 (one time per column)
%      or [3 x N], in which case q1 is a pure quaternion
%   q1 [4 x N]  - values of the quaternion at the given times
%                 (one time per column)
%      or [3 x N], in which case q2 is a pure quaternion
%   ts [N x 1] - vector of times
%
% Output
%   y [4 x N]  - values of the rotated function at the given times
%                (one time per column)
%   ts [N x 1] - vector of times (equal to the corresponding input)
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    if nargin<3
        ts=(1:size(q1,2))';
    else
        ts=ts(:);
    end

    if length(size(q1))~=2 || size(q1,1)<3 ||  size(q1,1)>4
        error('tsQdot: first input must be a time series of 3-vectors or 4-vectors ([%s])\n',...
              index2str(size(q1)));
    end

    if length(size(q2))~=2 || size(q2,1)<3 ||  size(q2,1)>4
        error('tsQdot: first input must be a time series of 3-vectors or 4-vectors ([%s])\n',...
              index2str(size(q2)));
    end

    if length(ts)~=size(q1,2) || length(ts)~=size(q2,2)
        error('tsQdot: length of sample times does not match size of inputs (%d,[%s],[%s])\n',...
              length(ts),index2str(size(q1)),index2str(size(q2)));
    end

    if ~isequal(class(q1),'Tcalculus') && ~isequal(class(q2),'Tcalculus')
        whichtprod=@mytprod;
    else
        whichtprod=@tprod;
    end

    % pq=[p0*q0-p'*q;
    %     q0*p+p0*q+cross(p,q)];

    if size(q1,1)==4 && size(q2,1)==3
        % full x pure quaternion
        if isequal(class(q1),'Tcalculus')
            q10=reshape(q1(1,:),length(ts));
        else
            q10=q1(1,:)';
        end
        y=[reshape(-tsDot(q1(2:4,:),q2,ts),1,size(q1,2));
           whichtprod(q10,[2],q2,[1,2])+tsCross(q1(2:4,:),q2,ts)];
    elseif size(q1,1)==4 && size(q2,1)==4
        % full x full quaternion
        if isequal(class(q1),'Tcalculus')
            q10=reshape(q1(1,:),length(ts));
        else
            q10=q1(1,:)';
        end
        if isequal(class(q2),'Tcalculus')
            q20=reshape(q2(1,:),length(ts));
        else
            q20=q2(1,:)';
        end
        y00=q10.*q20;
        if ~isequal(class(y00),'Tcalculus')
            y00=y00';
        end
        y=[reshape(y00-tsDot(q1(2:4,:),q2(2:4,:),ts),1,size(q1,2));
           whichtprod(q10,[2],q2(2:4,:),[1,2])+whichtprod(q20,[2],q1(2:4,:),[1,2])+tsCross(q1(2:4,:),q2(2:4,:),ts)];
    else
        error('tsQdot: not implemented for sizes [%s], [%s]\n', ...
              index2str(size(q1)),index2str(size(q2)));
    end
end

function test
    % Numeric
    ts=pi/2:pi/10:4*pi+pi/2;
    x=ones(3,length(ts))
    omega=rand(3,1);
    omega=omega/norm(omega,2);
    theta=ts;
    q=[cos(theta/2);omega*sin(theta/2)];
    y=tsQdot(q,x,ts)

    plot(ts,q','-x',ts,y','-+');
    legend('q1','q2','q3','q4','y1','y2','y3','y4')

    % Symbolic
    Tvariable x [3,length(ts)]
    Tvariable q [4,length(ts)]
    y=tsQdot(q,x,ts)
end