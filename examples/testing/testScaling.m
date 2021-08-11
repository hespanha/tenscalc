% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all;
!rm -fr tmp* @tmp*

%% min_u max_d (x-1)^2-2*d^2  with x=u+d;
% without latent variables

Tvariable u;
Tvariable d;
x=u+d;

J=(x-1).^2-2*d.^2;

pars.P1optimizationVariables={u};
pars.P2optimizationVariables={d};
pars.P1objective=J;
pars.P2objective=-J;
pars.P1constraints={d>-.1};
pars.P2constraints={d>-.1};
pars.outputExpressions=struct('u',u,'d',d,'x',x);

pars.scaleCost=1e-3;

%pars.desiredDualityGap=1e-10;
pars.classname='tmp';
pars.debugConvergence=false;
pars.solverVerboseLevel=3;

classname=class2equilibriumLatentCS(pars);

obj=feval(classname);

setV_u(obj,1);
setV_d(obj,1);
sol=solve(obj),
out=getOutputs(obj),


clear all;

%% min_u max_d (x-1)^2-2*d^2  with x=u+d;
% with equality but no latent variables

Tvariable u;
Tvariable d;
Tvariable xu;
Tvariable xd;

Ju=(xu-1).^2-2*d.^2;
Jd=(xd-1).^2-2*d.^2;

pars.P1optimizationVariables={u,xu};
pars.P2optimizationVariables={d,xd};
pars.P1objective=Ju;
pars.P2objective=-Jd;
pars.P1constraints={xu==u+d,d>-.1};
pars.P2constraints={xd==u+d,d>-.1};
pars.outputExpressions=struct('u',u,'d',d,'xu',xu,'xd',xd);

pars.scaleCost=1e-3;

%pars.desiredDualityGap=1e-10;
pars.classname='tmp';
pars.debugConvergence=false;
pars.solverVerboseLevel=3;

classname=class2equilibriumLatentCS(pars);

obj=feval(classname);

setV_u(obj,1);
setV_d(obj,1);
setV_xu(obj,1);
setV_xd(obj,-1);
sol=solve(obj),
out=getOutputs(obj),

clear all
%% min_u max_d (x-1)^2-2*d^2  with x=u+d;
% with latent variables

Tvariable u;
Tvariable d;
Tvariable x;

pars.latentVariables={x};
pars.latentConstraints={x==u+d};

J=(x-1).^2-2*d.^2;

pars.P1optimizationVariables={u};
pars.P2optimizationVariables={d};
pars.P1objective=J;
pars.P2objective=-J;
pars.P1constraints={d>-.1};
pars.P2constraints={d>-.1};
pars.outputExpressions=struct('u',u,'d',d,'x',x);

pars.scaleCost=1e-3;

%pars.desiredDualityGap=1e-10;
pars.classname='tmp';
pars.debugConvergence=false;
pars.solverVerboseLevel=3;

classname=class2equilibriumLatentCS(pars);

obj=feval(classname);

setV_u(obj,1);
setV_d(obj,1);
setV_x(obj,1);
sol=solve(obj),
out=getOutputs(obj),

