function [G,F,nus,lambdas,outputExpressions,template]=...
    parseConstraints(code,classname,constraints,outputExpressions,prefix)
% [G,F,nus,lambdas,template,outputExpressions]=...
%    parseConstraints(code,classname,constraints,outputExpressions,prefix)
%
% Parses an array of constraints into different types in order to
% 1) separates contraints into different type (equality and inequality)
% 2) create the corresponding dual variables
% 3) adds appropriate sets for the dual variables to a csparse
%    object
% 4) generates createGateway templates to generate cmex functions
%    to set values for the dual variable
% 5) replaces appearences of the dual variables in
%    outputExpressions to the given names.
%
% 
% Input: 
%
% code              - csparse object 
% classname         - classname to be used in the names of the functions
% constraints       - cell array of constrainsts as TC symbolic variables
% outputExpressions - original output expressions
% prefix            - string with prefix to be used for the name of
%                     the dual variables
%
% Outputs:
%
% Gcells - cell array of TC symbolic variables corresponding to
%          equality constraints (value should be made equal to
%          zero)
% nus    - cell array of dual variables corresponding to the 
%          equality constraints in Gcells
% Fcells - cell array of TC symbolic variables corresponding to
%          inequality constraints (value should be made positive)
% lambdas - cell array of dual variables corresponding to the 
%          inequality constraints in Fcells
% outputExpressions - output expressions with the dual variables
%          replaced by the ones that appear in nus & lambdas
% template - struct array with templates to be passed to
%          createGateway() to generate cmex functions to
%          initialize all dual variables (nus and lambdas)
%          (this output argument is optional)
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
 
createSets4Duals=false;

if nargin<5
    prefix='';
end

if nargout>=6
    % template for createGateway
    template=struct('MEXfunction',{},...% string
                    'Sfunction',{},...  % string
                    'Cfunction',{},...  % string
                    'method',{},...     % string
                    'inputs',struct(...  
                        'type',{},...   % string
                        'name',{},...   % cell-array of strings (one per dimension)
                        'sizes',{}),... % cell-array of strings (one per dimension)
                    'outputs',struct(... % string
                        'type',{},...   % string
                        'name',{},...   % cell-array of strings (one per dimension)
                        'sizes',{}),... % cell-array of strings (one per dimension)
                    'preprocess',{},... % strings (starting with parameters in parenthesis)'
                    'includes',{});     % cell-array of strings (one per file)
end

Gcells={};
Fcells={};
nus={};
lambdas={};
for k=1:length(constraints)
    switch (type(constraints{k}))
      case {'iszero'}
        % remove iszero
        op1=Tcalculus(operands(constraints{k}));
        Gcells{end+1}=op1;
        % create appropriate nu
        vname=sprintf('%snu%d_',prefix,length(nus)+1);
        nus{end+1}=Tvariable([vname,'_'],size(op1));
        if nargout>=6 && createSets4Duals
            template(end+1,1).MEXfunction=sprintf('%s_set_%s',classname,name(nus{end}));
            template(end).Cfunction=sprintf('%s_set_%s',classname,name(nus{end}));
            template(end).method=sprintf('setD_%s',name(nus{end}));
            template(end).inputs =struct('type','double',...
                                         'name',name(nus{end}),...
                                         'sizes',size(nus{end}));
        end
        declareSet(code,nus{end},sprintf('%s_set_%s',classname,name(nus{end})));
        outputExpressions=substitute(outputExpressions,...
                                     Tvariable(vname,size(op1)),...
                                     nus{end});
      case {'ispositive'}
        % remove ispositive
        op1=Tcalculus(operands(constraints{k}));
        Fcells{end+1}=op1;
        % create appropriate lambda;
        vname=sprintf('%slambda%d_',prefix,length(lambdas)+1);
        lambdas{end+1}=Tvariable([vname,'_'],size(op1));
        if nargout>=6 && createSets4Duals
            template(end+1,1).MEXfunction=sprintf('%s_set_%s',classname,name(lambdas{end}));
            template(end).Cfunction=sprintf('%s_set_%s',classname,name(lambdas{end}));
            template(end).method=sprintf('setD_%s',name(lambdas{end}));
            template(end).inputs =struct('type','double',...
                                         'name',name(lambdas{end}),...
                                         'sizes',size(lambdas{end}));
        end
        declareSet(code,lambdas{end},sprintf('%s_set_%s',classname,name(lambdas{end})));
        outputExpressions=substitute(outputExpressions,...
                                     Tvariable(vname,size(op1)),...
                                     lambdas{end});
      otherwise
        error('constraint of type ''%s'' not implemented\n',type(constraints{k}));
    end
end

%% Pack constraints 
G=packExpressions(Gcells);
F=packExpressions(Fcells);
    

