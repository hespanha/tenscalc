% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

%% Optimizations in tenscalc/examples
examplesFolder=which('cmex2optimizeCS');
examplesFolder=fileparts(examplesFolder);
cd(fullfile(examplesFolder,'../examples'));

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
examplesFolder=which('cmex2optimizeCS');
examplesFolder=fileparts(examplesFolder);
cd(fullfile(examplesFolder,'../examples/mpcmhe'));

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

mpc_dcmotor
mpcmhe_dcmotor
mpc_unicycle

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%% Back to examples
examplesFolder=which('cmex2optimizeCS');
examplesFolder=fileparts(examplesFolder);
cd(fullfile(examplesFolder,'../examples'));
