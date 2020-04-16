function varargout=class2compute(varargin)
% To get help, type class2compute('help')
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
    
    %% Function global help
    declareParameter(...
        'Help', {
            'Creates a set of matlab functions for performing a csparse computation.'
            ' '
            'The computation is performed through a matlab class with methods'
            'for the set, get, and copy operations.'
                });
    
    localVariables_=parameters4compute(localVariables_);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Retrieve parameters and inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [stopNow,params]=setParameters(nargout,varargin);
    if stopNow
        return 
    end
    

    %% transfer any folder in classname into folder
    [folder,classname]=fileparts(fsfullfile(folder,classname));
    
    %% create folder if it does not exist
    if ~strcmp(folder,'.') && ~exist(folder,'dir')
        fprintf('class2compute: outputs folder ''%s'' does not exist, creating it\n',folder);
        if mkdir(folder)==0
            error('Unable to create folder ''%s''\n',folder)
        end
    end
    
    rmpath(folder);
    addpath(folder);
    
    %% Fix class when gotten from pedigree
    classname=regexprep(classname,'+TS=','_TS_');
    classname=regexprep(classname,'-','_');
    classname=regexprep(classname,'+','_');

    classHelp=helpFromTemplate(classname,csparseObject.template);

    fprintf(' creating matlab code... ');
    t_compile2matlab=clock();
    compile2matlab(csparseObject,...
                   fsfullfile(folder,sprintf('%s.m',classname)),...
                   fsfullfile(folder,sprintf('%s.log',classname)),...
                   classHelp,profiling);
    statistics.time.compile2matlab=etime(clock,t_compile2matlab);
    
    fprintf(' done creating matlab code (%.3f sec)\n',etime(clock,t_compile2matlab));
    
    rehash path;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    varargout=setOutputs(nargout,params);
end
        
