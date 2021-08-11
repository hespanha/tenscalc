function [Ts,xMeas,xFut,uPast,uFut,dynamics]=TvariablesMPC(nX,nU,T,delay,fun,varargin);
% [Ts,xMeas,xFut,uPast,uFut,dynamics]=TvariablesMPC(nX,nU,T,delay,fun,varargin);
%
% Create the key Tvariables needed to create an MPC solver
% and also the Tcalculus constraints corresponding to the dynamics
%
% Inputs (all scalars):
%   nX    - state dimension
%   nU    - input dimension
%   T     - horizon length
%   delay - input delay
%   fun   - function handle for the ODE dynamics, of the form
%           @(x,u,parameter1,parameter2,...)
%   parameter1, parameter2, .... - parameters for fun
%
% Outputs:
%   Ts []              - Tcalculus variable with sampling interval
%   xMeas  [nX,1]      - Tcalculus variable with (current) measured state
%                          [ x(t) ]
%   xFut   [nX,T]      - Tcalculus variable with future state vector
%                          [ x(t+Ts), x(t+2*Ts), ..., x(t+T*Ts) ];
%   uPast [nU,delay]   - Tcalculus variable with previously computed inputs
%                          [u(t), u(t+Ts), ..., u(t+(delay-1)*Ts) ]
%                        (only needed if delay>1)
%   uFut [nU,T-delay]  - Tcalculus variable with previously computed inputs
%                          [u(t+delay*Ts), ... , u(t+(T-1)*Ts) ]
%   dynamics           - Tcalculus constraint encoding trapezoidal integration,
%                        but with u assumed held constant (ZOH)
%                        between sampling times.
%
% Attention: The variables created will have precisely the names
%                Ts,xMeas,xFut,uPast,uFut
%            which will be important to when calling the corresponding
%                setV_... and setP_...
%            functions prior to calling the solver
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    Tvariable Ts     [];

    Tvariable xMeas  [nX,1];
    Tvariable xFut   [nX,T];

    Tvariable uPast  [nU,delay];
    Tvariable uFut   [nU,T-delay];

    xPast=[xMeas,xFut(:,1:end-1)];
    uAll=[uPast,uFut];
    dynamics = xFut-xPast==.5*Ts*(fun(xFut,uAll,varargin{:})+fun(xPast,uAll,varargin{:}));

end