function dispersioneering(params_handle,varargin)

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));

%% Check for valid params function handle
if ~exist('params_handle','var')
    if ~isa(params_handle,'functionhandle')
        error('Invalid parameters function handle, exiting.');
    end
end

global phys;
phys = constants();

%% Create options
opts = OPTS(varargin);

%% Create parameters
params = params_handle();

%% Run dispersioneering with these parameters
[output] = dispersioneering_run(opts,params);

end


