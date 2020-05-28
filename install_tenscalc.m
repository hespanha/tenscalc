fprintf('Seeting up path:\n');
home=[fileparts(which('install_tenscalc')),'/lib'];
folders={home;[home,'/csparse']};

s=path;
if ispc
    old=regexp(s,'[^;]*(tenscalc.lib.csparse|tenscalc.lib)[^/;]*','match');
else
    old=regexp(s,'[^:]*(tenscalc.lib.csparse|tenscalc.lib)[^/:]*','match');
end
if ~isempty(old)
    fprintf('  removing from path:\n');
    disp(old')
    rmpath(old{:})
end

fprintf('  adding to path:\n');
addpath(folders{:});
disp(folders)

fprintf('  saving path...');
try
    savepath;
catch me
    fprintf('ATTENTION: unable to save path. This was probably caused because of insufficient permissions. Either change the permissions of your ''matlabroot'' folder or add following strings to the matlab path:');
    disp(folders)
    rethrow(me)
end
fprintf('done with path!\n');


fprintf('Compiling...\n');
compileInstructionsTable;
fprintf('done!\n');

pth=findSuiteSparse();
if isempty(findSuiteSparse)
    fprintf('The umfpack not installed or not in the path, ''umfpack'' option will not be available. This should not be a problem. However, if you have umfpack installed add ''SuiteSparse/UMFPACK/MATLAB/'' to your matlab path\n');
end