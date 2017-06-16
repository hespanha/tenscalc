function [y,ts]=tsRotation(q,x,ts);
% [y,ts]=tsRotation(q,x,ts);
%
% Computes a 3-vector time-series (y,ts) that represents the rotation
% of a 3-vector time-series (x,ts) by a 4-vector time-series (q,ts)
% representing a quaternion. Specifically,
%    y = q x q^*
% In terms of the rotation matrix
%    R=expm(theta*skew(omega)),
%                 skew(w1,w2,w3)=[0,-w3,w2;w3,0,-w1;-w2,w1,0]
% with    
%    q=(cos(theta/2),omega*sin(theta/2))
% this corresponds to
%    y=R x
%
% Inputs:
%   x [3 x N]  - values of the function at the given times
%                (one time per column)
%   q [4 x N]  - values of the quaternion at the given times
%                (one time per column)
%   ts [N x 1] - vector of times
%
% Output
%   y [3 x N]  - values of the rotated function at the given times 
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
        ts=(1:size(q,2))';
    else
        ts=ts(:);
    end

    if ~isequal(class(q),'Tcalculus') && ~isequal(class(x),'Tcalculus') 
        whichtprod=@mytprod;
    else
        whichtprod=@tprod;
    end

    if length(size(q))~=2 || size(q,1)~=4
        error('tsRotation: first input must be a time series of 4-vectors ([%s])\n',...
              index2str(size(q)));
    end
    
    if length(size(x))~=2 || size(x,1)~=3
        error('tsRotation: second input must be a time series of 3-vectors ([%s])\n',...
              index2str(size(x)));
    end
    
    if length(ts)~=size(x,2) || length(ts)~=size(q,2)
        error('tsRotation: length of sample times does not match size of inputs (%d,[%s],[%s])\n',...
              length(ts),index2str(size(q)),index2str(size(x)));
    end
    
    % pq=[p0*q0-p'*q;
    %     q0*p+p0*q+cross(p,q)];
    if isequal(class(q),'Tcalculus')
        q0=reshape(q(1,:),length(ts));
    else
        q0=q(1,:)';
    end
    minus_qx0=tsDot(q(2:4,:),x,ts);
    if ~isequal(class(minus_qx0),'Tcalculus')
        minus_qx0=minus_qx0';
    end
    qx1=whichtprod(q0,[2],x,[1,2])+tsCross(q(2:4,:),x,ts);
    y=whichtprod(minus_qx0,[2],q(2:4,:),[1,2])...
      +whichtprod(q0,[2],qx1,[1,2])...
      +tsCross(q(2:4,:),qx1,ts);
        
    % y1=whichtprod(minus_qx0,[2],q(2:4,:),[1,2]);
    % y2=whichtprod(q0,[2],qx1,[1,2]);
    % y3=tsCross(q(2:4,:),qx1,ts);
    % whos
end

function test
    % Numeric
    ts=pi/2:pi/10:4*pi+pi/2;
    x=ones(3,length(ts))
    omega=rand(3,1);

    %x=[0;0;1];x=repmat(x,1,length(ts));
    %omega=[0;1;0]

    omega=omega/norm(omega,2);
    theta=ts;
    q=[cos(theta/2);omega*sin(theta/2)];
    y=tsRotation(q,x,ts)

    % Compare with rotation matrix
    yy=y;
    for i=1:length(ts)
        th=theta(i);
        % exponential formula
        R=expm(th*[0,-omega(3),omega(2);omega(3),0,-omega(1);-omega(2),omega(1),0]);
        % explicit formula (http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation)
        
        qi=q(:,i);
        R=[qi(1)^2+qi(2)^2-qi(3)^2-qi(4)^2,2*qi(2)*qi(3)-2*qi(1)*qi(4),2*qi(2)*qi(4)+2*qi(1)*qi(3);
           2*qi(2)*qi(3)+2*qi(1)*qi(4),qi(1)^2-qi(2)^2+qi(3)^2-qi(4)^2,2*qi(3)*qi(4)-2*qi(1)*qi(2);
           2*qi(2)*qi(4)-2*qi(1)*qi(3),2*qi(3)*qi(4)+2*qi(1)*qi(2),qi(1)^2-qi(2)^2-qi(3)^2+qi(4)^2];
           
        yy(:,i)=R*x(:,i);
    end
    y
    yy
    if norm(y-yy)>sqrt(eps)
        error('tsRotation error %g\n',norm(y-yy))
    end
    
    plot(ts,x','-x',ts,y','-+');
    legend('x1','x2','x3','y1','y2','y3')

    % Symbolic
    Tvariable x [3,length(ts)]
    Tvariable q [4,length(ts)]
    y=tsRotation(q,x,ts)

end