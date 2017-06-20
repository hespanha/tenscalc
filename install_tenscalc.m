fprintf('Seeting up path...');
home=[fileparts(which('install_tenscalc')),'/lib'];
folders={home;[home,'/csparse']};

s=path;
old=regexp(s,'[^:]*tenscalc[^:]*','match');
if ~isempty(old)
    fprintf('removing from path:\n');
    disp(old')
    rmpath(old{:})
end

fprintf('adding to path:\n');
addpath(folders{:});
disp(folders)

fprintf('Compiling...\n');
compileInstructionsTable;

fprintf('saving path...');
try
    savepath;
catch me
    fprintf('ATTENTION: unable to save path, add following strings to the matlab path:');
    disp(folders)
    rethrow
end

fprintf('done!\n');