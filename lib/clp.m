function [alpha,k]=clp(x,dx)
% [alpha,k]=clp(x,dx)
%
% Canonical LP: For a given x with all entries >=0 ,
% computes the scalar
%     max { alpha >0 : x + alpha dx >= 0 }
% x and dx must have the same size.
% Used by ipm Newton solvers to determine the step size.
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


    k=find(dx<0);
    if isempty(k)
        alpha=inf;
    else
        [alpha,kk]=min(-x(k)./dx(k));
        k=k(kk);
    end
end