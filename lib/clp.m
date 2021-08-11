function [alpha,k]=clp(x,dx)
% [alpha,k]=clp(x,dx)
%
% Canonical linear proram (LP): For a given x with all entries >=0 ,
% computes the scalar
%     max { alpha >0 : x + alpha dx >= 0 }
% x and dx must have the same size.
%
% Used by ipm Newton solvers to determine the step size.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    k=find(dx<0);
    if isempty(k)
        alpha=inf;
    else
        [alpha,kk]=min(-x(k)./dx(k));
        k=k(kk);
    end
end