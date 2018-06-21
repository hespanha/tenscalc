mls
sls
l1l2estimationCS

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

cd mpcmhe
mpc_dcmotor
mpcmhe_dcmotor
mpc_unicycle_pursuit

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');
