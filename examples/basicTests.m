% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

%% Optimizations in tenscalc/examples
examples=which('cmex2optimizeCS');
examples=fileparts(examples);
cd(fullfile(examples,'../examples'));

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

mls
sls
l1l2estimationCS

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Optimizations in tenscalc/examples/mpcmhe
examples=which('cmex2optimizeCS');
examples=fileparts(examples);
cd(fullfile(examples,'../examples/mpcmhe'));

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

mpc_dcmotor
mpcmhe_dcmotor
mpc_unicycle_pursuit

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Back to examples
examples=which('cmex2optimizeCS');
examples=fileparts(examples);
cd(fullfile(examples,'../examples'));
