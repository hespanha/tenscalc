function varargout=cmex2compute(varargin)
% To get help, type cmex2compute('help')
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
            'Creates a set of cmex C functions for performing a csparse computation.'
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
        fprintf('cmex2compute: outputs folder ''%s'' does not exist, creating it\n',folder);
        if mkdir(folder)==0
            error('Unable to create folder ''%s''\n',folder)
        end
    end

    rmpath(folder);
    addpath(folder);

    classFolder=fsfullfile(folder,sprintf('@%s',classname));
    if ~exist(classFolder,'dir')
        fprintf('class classFolder @%s does not exist, creating it... ',classname);
        mkdir(classFolder);
    end


    %% Fix class when gotten from pedigree
    classname=regexprep(classname,'+TS=','_TS_');
    classname=regexprep(classname,'-','_');
    classname=regexprep(classname,'+','_');

    %% Compile code
    fprintf(' Creating C code... ');
    t_compile2C=clock();
    compile2C(csparseObject,minInstructions4loop,maxInstructionsPerFunction,...
              sprintf('%s.c',classname),...
              sprintf('%s.h',classname),...
              sprintf('%s.log',classname),...
              classFolder,profiling);
    csparseObject.statistics.time.compile2C=etime(clock,t_compile2C);

    fprintf('  done creating C code (%.3f sec)\n',etime(clock,t_compile2C));

    tpl=csparseObject.template;
    % add classname to cmex & S-functions to make names unique
    for i=1:length(tpl)
        tpl(i).MEXfunction=sprintf('%s_%s',classname,tpl(i).MEXfunction);
        tpl(i).Sfunction=sprintf('%s_%s',classname,tpl(i).Sfunction);
        % no simulink functions
        tpl(i).Sfunction='';
    end

    t_createGateway=clock();

    statistics=createGateway('template',tpl,...
                             'CfunctionsSource',fsfullfile(classFolder,sprintf('%s.c',classname)),...
                             'callType','dynamicLibrary',...
                             'dynamicLibrary',classname,...
                             'folder',folder,...
                             'className',classname,...
                             'absolutePath',absolutePath,...
                             'compilerOptimization',compilerOptimization);

    statistics.time.createGateway=etime(clock,t_createGateway);
    statistics.createGateway=statistics;

    fprintf('  done creating & compiling gateways & library (%.2f sec)\n',...
            etime(clock,t_createGateway));

    rehash path;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    varargout=setOutputs(nargout,params);
end
