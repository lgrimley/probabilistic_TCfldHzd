%==========================================================================
% Diagnostic Visualization of TCR Rainfall Grids
%
% Purpose:
%   Generates timestep-by-timestep visualizations of TCR rainfall output
%   for a single storm. Each rainfall snapshot is plotted on its native
%   dynamic grid and saved as an image file for quality control and
%   debugging.
%
% Use Case:
%   - Verify rainfall structure and magnitude
%   - Inspect storm evolution in space and time
%   - Confirm longitude conversion and grid alignment
%
% Inputs:
%   - TCR rainfall MAT file containing:
%       * PLAT_SAVE   latitude grid (dynamic)
%       * PLONG_SAVE  longitude grid (dynamic)
%       * rain        rainfall values
%
% Outputs:
%   - PNG images, one per timestep, saved to an output directory
%
% Assumptions:
%   - Rainfall grids are defined on a moving (storm-relative) grid
%   - Longitude values may be in 0–360 and are converted to −180–180
%
% Notes:
%   - Intended for diagnostic use, not production processing
%   - No spatial interpolation or regridding is performed
%
% Load TCR rainfall MAT file
% Convert longitude to −180–180
% 
% FOR each timestep
%     Extract dynamic latitude, longitude, and rainfall grids
%     Plot rainfall using native grid
%     Save plot as PNG
% END
%==========================================================================

% Directory containing TCR rainfall MAT files
rain_folder = 'Z:\Data-Expansion\users\lelise\Chapter3\NCEP_Reanalysis\rain\TCR_RainOutput\';

% Specify storm rainfall file to visualize
filename = [rain_folder,'0245.mat'];

% Load rainfall and grid variables
% Loads: PLAT_SAVE, PLONG_SAVE, rain
load(filename)

% Output directory for rainfall plots
outdir = 'Z:\Data-Expansion\users\lelise\Chapter3\NCEP_Reanalysis\rain\tc_245';

% Get dimensions of rainfall array
trk = size(rain);

% Convert longitude from 0–360 to −180–180 convention
PLONG_SAVE(PLONG_SAVE>180) = PLONG_SAVE(PLONG_SAVE>180)-360;


% Example timestep extraction (not used in loop below)
t = 5;
xt = PLONG_SAVE(:,:,t);
yt = PLAT_SAVE(:,:,t);
rt = rain(:,:,t);


% Loop through all timesteps in the storm
for t = 1:trk(3)

    % Extract longitude, latitude, and rainfall at timestep t
    xt = PLONG_SAVE(:,:,t);
    yt = PLAT_SAVE(:,:,t);
    rt = rain(:,:,t);

    % Construct output image filename
    stormfile = strcat(num2str(t), '.png');

    % Plot rainfall using native TCR grid
    pcolor(xt, yt, rt)
    shading flat
    colorbar()

    % Save figure to output directory
    saveas(gcf,sprintf('%s/%s',outdir,stormfile))

    % Close figure to avoid memory buildup
    close(gcf)
end
