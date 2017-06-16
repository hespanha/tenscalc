thisDir=pwd;

if libisloaded('instructionsTable')
    unloadlibrary('instructionsTable')
end

s=which('compileInstructionsTable');
p=fileparts(s);
cd(p);

types=instructionTypes();
f=fields(types);

fid=fopen('instructionsTableTypes.h','w');
fprintf(fid,'typedef enum instructionType_e {\n');
for i=1:length(f)
    fprintf(fid,'    %s=%d,\n',f{i},getfield(types,f{i}));
end
fprintf(fid,'} instructionType_t;\n');
fclose(fid);

try
    createGateway('template','instructionsTable.c',...
                  'callType','dynamicLibrary',...
                  'CfunctionsSource','instructionsTable.c',...
                  'dynamicLibrary','instructionsTable',...
                  'verboseLevel',1);
    %!nm -a instructionsTable.dylib
    %!otool -Tvt instructionsTable.dylib 
    
    cd(thisDir)
catch me
    cd(thisDir)
    rethrow(me)
end

[notfound,warnings]=loadlibrary('instructionsTable');



