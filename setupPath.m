function setupPath()
%SETUPPATH Add only the runtime folders belonging to this distribution.
rootDir = fileparts(mfilename('fullpath'));
folders = {'core', 'association', 'ggiw', 'scenario', 'metrics', 'math', 'scripts'};
addpath(rootDir);
for k = 1:numel(folders)
    addpath(fullfile(rootDir, folders{k}));
end
end
