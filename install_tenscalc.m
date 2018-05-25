fprintf('Seeting up path...\n');
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

fprintf('Compiling...\n');
compileInstructionsTable;

fprintf('Saving path...');
try
    savepath;
catch me
    fprintf('ATTENTION: unable to save path, add following strings to the matlab path:');
    disp(folders)
    rethrow(me)
end

fprintf('done!\n');
