function add_BED_path()
directory = fileparts(fileparts(mfilename('fullpath')));
cd(directory)
addpath(genpath('.'))
end