function localVariables_=parameters4compute(localVariables_)
% Declare input parameters common to the 4 tenscalc functions:
%   cmex2compute.m
%   class2compute.m

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CreateGateway parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                      });

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CreateGateway parameters
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
        'VariableName','compilerOptimization',...
        'DefaultValue','-Ofast',...
        'AdmissibleValues',{'-O0','-O1','-O2','-O3','-Ofast'},...
        'Description', {
            'Optimization flag used for compilation.'
            'Only used when compileGateways, compileLibrary, or compileStandalones'
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
