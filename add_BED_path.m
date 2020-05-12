function add_BED_path()
directory = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('.'))
end
