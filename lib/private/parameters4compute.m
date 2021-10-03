function localVariables_=parameters4compute(localVariables_)
% Declare input parameters common to the 4 tenscalc functions:
%   cmex2compute.m
%   class2compute.m
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Code generation parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    declareParameter(...
        'VariableName','csparseObject',...
        'Description', {
            'csparse object with computations to be performed.'
                       });

    declareParameter(...
        'VariableName','classname',...
        'DefaultValue',getFromPedigree(),...
        'Description', {
            'Name of the class to be created.'
            'A matlab class will be created with this name plus a |.m| extension.';
            ' '
            'One can look "inside" this class to find the name of the cmex functions.'
                       });

    declareParameter(...
        'VariableName','folder',...
        'DefaultValue','.',...
        'Description', {
            'Path to the folder where the files will be created.';
            'Needs to be in the Matlab path.'
                       });

    declareParameter(...
        'VariableName','simulinkLibrary',...
        'DefaultValue','',...
        'Description', {
            'Name of a simulink library to be created with Simulink blocks that can be used'
            'to call the different solver functions.'
            'The blocks are created with direct feedthrough.'
            'No library and no simulink blocks are created if |simulinkLibrary| is an empty string.'
                       });

    declareParameter(...
        'VariableName','scratchbookType',...
        'DefaultValue','double',...
        'AdmissibleValues',{'double','float'},...
        'Description',{
            'C variable type used for the scratchbook. See csparse documentation.';
            ' ';
            'This parameter is only used for C-code solvers.';
                      });
    declareParameter(...
        'VariableName','fastRedundancyCheck',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'when true, very intensive operations (like the lu factorization) '
            'do not check if it is possible to reuse computations.'
            'See csparse documentation.'
            ' ';
            'This parameter is only used for C-code solvers.';
                      });

    declareParameter(...
        'VariableName','minInstructions4loop',...
        'DefaultValue',50,...
        'Description',{
            'Minimum number of similar instruction that will be execute as part of a'
            'for(;;) loop, rather than being executed as independent C commands.'
            'When equal to ''inf'', instructions will never be grouped into foor loops.'
            ' ';
            'This parameter is only used for C-code solvers.';
                      });

    declareParameter(...
        'VariableName','maxInstructionsPerFunction',...
        'DefaultValue',100,...
        'Description',{
            'Maximum number of instructions to be included into a single function.'
            'When equal to ''inf'', there is no limit on the size of a single function.';
            ' ';
            'Large values of |maxInstructionsPerFunction| and therefore large functions';
            'give more opportunities for compiler optimization, but can resul in very';
            'slow compilation (especially with compiler optimization turned on).';
            ' ';
            'This parameter is only used for C-code.';
                      });

    declareParameter(...
        'VariableName','useUmfpack',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true| the system of linear equations needed for the Newton step';
            'is solved using the umfpack package.';
            'Using this package tyically increases the solve time, but results in';
            'code that is more robust.'
            ' ';
            'This parameter is only used for C-code solvers.';
                      });

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compilation/linking parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    declareParameter(...
        'VariableName','compilerOptimization',...
        'DefaultValue','-O0',...
        'AdmissibleValues',{'-O0','-O1','-O2','-O3','-Os','-Ofast'},...
        'Description',{
            'Optimization parameters passed to the C compiler.'
            '* |-Ofast| often generates the fastest code'
            '* |-Os| often generates the smallest code'
            '* |-O1| includes the most important optimizations, while still compiling fast'
            '* |-O0| compiles the fastest'
            'Only used when either |compileGateways|, |compileLibrary|, or '
            '|compileStandalones| are set to |true|.'
            ' ';
            'This parameter is only used for C-code solvers.';
        });

    declareParameter(...
        'VariableName','compileGateways',...
        'DefaultValue',false,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'When |true| the gateway functions are compiled using |cmex|.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });
    declareParameter(...
        'VariableName','compileLibrary',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'When |true| the dynamicLibrary is compiled using |gcc|.'
            'This parameter is used only when |callType=''dynamicLibrary''|.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });
    declareParameter(...
        'VariableName','compileStandalones',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'When |true| the standalone/server executable is compiled using |gcc|.'
            'This parameter is used only when |callType| has one of the values:';
            '   |''standalone''| or |''client-server''|'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','callType',...
        'DefaultValue','dynamicLibrary',...
        'AdmissibleValues',{'dynamicLibrary','client-server'},...
        'Description', {
            'Method used to interact with the solver:';
            '* |dynamicLibrary| - the solver is linked with matlab using a dynamic library';
            '* |client-server|   - the solver runs as a server in independent process,';
            '                     and a socket is used to exchange data.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','targetComputer',...
        'DefaultValue',lower(computer),...
        'AdmissibleValues',{'maci64','glnxa64','pcwin64'},...
        'Description', {
            'OS where the C files will be compiled.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','serverComputer',...
        'DefaultValue',lower(computer),...
        'AdmissibleValues',{'maci64','glnxa64','pcwin64'},...
        'Description', {
            'OS where the server will be compiled.'
            'This parameter is used only when |callType=''client-server''|.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','absolutePath',...
        'DefaultValue',true,...
        'AdmissibleValues',{true,false},...
        'Description', {
            'When ''true'' the the cmex functions use an absolute path to open';
            'the dynamic library, which means that the dynamic library cannot';
            'be moved away from the folder where it was created.';
            ' '
            'When ''false'' no path information about the dynamic library is'
            'included in the cmex function, which must then rely on the OS-specific';
            'method used to find dynamic libraries. See documentation of ''dlopen'''
            'for linux and OSX or ''LoadLibrary'' for Microsoft Windows.'
            ' '
            'This parameter is used only when ''callType''=''dynamicLibrary''.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','serverProgramName',...
        'DefaultValue','',...
        'Description', {
            'Name of the executable file for the server executable.'
            'This parameter is used only when |callType=''client-server''|.'
                       });

    declareParameter(...
        'VariableName','serverAddress',...
        'DefaultValue','localhost',...
        'Description', {
            'IP address (or name) of the server.'
            'This parameter is used only when |callType=''client-server''|.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });

    declareParameter(...
        'VariableName','port',...
        'DefaultValue',1968,...
        'Description', {
            'Port number for the socket that connects client and server.'
            'This parameter is used only when |callType=''client-server''|.'
            ' ';
            'This parameter is only used for C-code solvers.';
                       });


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Debugging parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    declareParameter(...
        'VariableName','profiling',...
        'DefaultValue',false,...
        'AdmissibleValues',{true,false},...
        'Description',{
            'When |true|, adds profiling to the C code.';
            ' '
            'Accumulated profiling information is diplayed on the screen when the'
            'dynamic library is unloaded.'
            ' ';
            'This parameter is only used for C-code solvers.';
                      });

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Output parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    declareOutput(...
        'VariableName','classname',...
        'Description', {
            'Name of the class created.'
                       });

    declareOutput(...
        'VariableName','statistics',...
        'Description', {
            'Structure with various statistics, including the file sizes and compilations times'
                       });

end
