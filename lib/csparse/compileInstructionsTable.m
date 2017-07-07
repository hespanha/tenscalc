
%% Unload previous library
if libisloaded('instructionsTable')
    unloadlibrary('instructionsTable')
end

%% Go to csparse folder
s=which('compileInstructionsTable');
p=fileparts(s);
thisDir=pwd;
cd(p);

%% Create hesder file with enum & defines
fid=fopen('instructionsTableTypes.h','w');
fprintf(fid,'// This file is automatically generated by compileInstructionsTable.m\n');

[iTypes,pTypes]=instructionTypes();

% enum with instruction types
f=fields(iTypes);
fprintf(fid,'typedef enum instructionType_e {\n');
for i=1:length(f)
    fprintf(fid,'    %s=%d,\n',f{i},getfield(iTypes,f{i}));
end
fprintf(fid,'} instructionType_t;\n');

% defines with profile types
f=fields(pTypes);
fprintf(fid,'// profile defines\n');
for i=1:length(f)
    fprintf(fid,'#define %s %d\n',f{i},getfield(pTypes,f{i}));
end
fprintf(fid,'#define P_nCountFlops %d\n',length(f));
fclose(fid);

%% Call cmextools
try
    if 1
        createGateway('template','instructionsTableUTHash.c',...
                      'callType','dynamicLibrary',...
                      'CfunctionsSource','instructionsTableUTHash.c',...
                      'dynamicLibrary','instructionsTable',...
                      'verboseLevel',1);
    else
        createGateway('template','instructionsTable.c',...
                      'callType','dynamicLibrary',...
                      'CfunctionsSource','instructionsTable.c',...
                      'dynamicLibrary','instructionsTable',...
                      'verboseLevel',1);
    end
    %!nm -a instructionsTable.dylib
    %!otool -Tvt instructionsTable.dylib 
    
    cd(thisDir)
catch me
    cd(thisDir)
    rethrow(me)
end

%% Load new library
[notfound,warnings]=loadlibrary('instructionsTable');



