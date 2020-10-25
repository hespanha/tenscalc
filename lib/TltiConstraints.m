function [stateConstraints,y,z]=TltiConstraints(A,B,C,D,G,H,x0,x,u,Ty,Tz);
% [stateConstraints,y,z]=TltiConstraints(A,B,C,D,G,H,x0,x,u,Ty,Tz);
%
% These commands generate constraints for optimization involving a discrete-time LTI system
%
% Inputs:
% A,B,C,D,G,H - state-space model of an LTI system with nx states, nu inputs,
%               ny measured outputs (C,D), and nz controlled
%               outputs (G,H)
%               C and/or G may be empty if the outputs y and/ot z
%               are not needed
% x0 - nx column vector with initial state   x(0)
% x  - nx by Tu Tcalculus tensor variable with state vectors x(1), x(2), ..., x(Tu)
% u  - nu by Tu matrix with control inputs   u(0), u(1), ..., u(Tu-1)
% Ty - desired time length for the measured output vector (see output y)
% Tz - desired time length for the controlle output vector (see output z)
%
% Output:
% stateConstraints  - constraints of the form
%                     x(t+1) = A x(t) + B u(t) for times 0,1, ... Tu-1
% y  - ny by Ty matrix with outputs          y(0), y(1), ..., y(Ty-1)
%      (only returned is C is not empty)
% z  - nz by Tz matrix with outputs          y(0), y(1), ..., y(Tz-1)
%      (only returned is G is not empty)
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

if ismember(class(B),{'Tvariable','Tcalculus'})
    nx=msize(B,1);
    nu=msize(B,2);
else
    nx=size(B,1);
    nu=size(B,2);
end
if ismember(class(u),{'Tvariable','Tcalculus'})
    Tu=msize(u,2);
else
    Tu=size(u,2);
end

u=reshape(u,[nu,Tu]);
A=reshape(A,[nx,nx]);
B=reshape(B,[nx,nu]);
stateConstraints=(x==A*[x0,x(:,1:Tu-1)]+B*u);


if ~isempty(C)
    if ismember(class(C),{'Tvariable','Tcalculus'})
        ny=msize(C,1);
    else
        ny=size(C,1);
    end
    C=reshape(C,[ny,nx]);
    D=reshape(D,[ny,nu]);
    y=C*[x0,x(:,1:Ty-1)]+D*u(:,1:Ty);
end


if ~isempty(G)
    if ismember(class(G),{'Tvariable','Tcalculus'})
        nz=msize(G,1);
    else
        nz=size(G,1);
    end
    G=reshape(G,[nz,nx]);
    H=reshape(H,[nz,nu]);
    z=G*[x0,x(:,1:Tz-1)]+H*u(:,1:Tz);
end
