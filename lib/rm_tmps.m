% Removes from the current directory all solvers whose names start with 'tmp'
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

% erase all objects before erasing classes
clear all;

% remove previous solvers from current directory
delete('toremove.m','tmp*');
rc=rmdir('@tmp*','s');


