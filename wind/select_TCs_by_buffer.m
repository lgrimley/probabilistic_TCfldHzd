%% Clear workspace and command window
clear
clc

%% ------------------------------------------------------------------------
% This script identifies tropical cyclones (TCs) whose tracks pass within
% a distance (200 km) of a specific ADCIRC gage location (in this example: Cape Fear outlet, Gage 194).
%
% The TC track data come from NCEP reanalysis, and the ADCIRC gage locations
% come from a low-resolution coastline file.
%
% The selected storms are then compared against storms available in ADCIRC
% water level (WL) simulation outputs to find common storm IDs.
% ------------------------------------------------------------------------

%% Load ADCIRC coastline (x,y) locations
% These vectors define the longitude (xvec) and latitude (yvec) of ADCIRC
% nodes or gage locations.
load 'C:\Users\lelise\Documents\globus_data\coastline_lores'

%% Load tropical cyclone (TC) track data
% This file contains latitude and longitude tracks for NCEP reanalysis storms.
% lat100 and lon100 are assumed to be arrays of size:
%   (number_of_storms x number_of_track_points)
load 'Z:\Data-Expansion\users\lelise\Chapter3\NCEP_Reanalysis\tracks\UScoast6_AL_ncep_reanal_roEst1rmEst1_trk100.mat'

%% Change directory to where the distance calculation function is stored
% m_lldist_L computes great-circle distances (km) between lat/lon points.
cd Z:\Data-Expansion\users\lelise\Chapter3\scripts\wind

%% Define ADCIRC gage location (Cape Fear outlet, Gage 194)
% Extract latitude and longitude for gage 194 from coastline vectors.
clat  = yvec(194);
clong = xvec(194);

%% Initialize arrays to store storm indices
select_tcs       = [];  % Storms that pass within 200 km of gage 194
not_selected_tcs = [];  % Storms farther than 200 km

%% Loop through all storms and calculate distance to gage
% Total number of storms is assumed to be 5018
for i = 1:5018

    % Extract latitude track for storm i and remove zero-padding
    tr_lat = lat100(i,:);
    tr_lat = tr_lat(tr_lat ~= 0);

    % Extract corresponding longitude track
    tr_lon = lon100(i,:);
    tr_lon = tr_lon(1:length(tr_lat));

    % Number of valid track points for this storm
    trk = length(tr_lat);

    % Compute distance from each track point to the ADCIRC gage
    d = [];
    for j = 1:trk
        % m_lldist_L returns distance in km between two lat/lon points
        [d(j), ~, ~] = m_lldist_L([tr_lon(j), clong], ...
                                 [tr_lat(j), clat]);
    end

    % Check if the minimum distance exceeds 200 km
    % If yes, storm is not selected
    if min(d) > 200
        not_selected_tcs(end+1) = i;
        continue
    end

    % Otherwise, storm is selected
    select_tcs(end+1) = i;

end

%% Load ADCIRC water level (WL) time series output
% This file contains:
%
% WL_inds_adj:
%   - Each row corresponds to an ADCIRC gage (e.g., 355 gages)
%   - Each row lists storm IDs relevant to that gage
%
% WL_mat:
%   - Dimensions: (gages, storms, time)
%   - Water level time series for each storm at each gage
load 'Z:\Data-Expansion\users\lelise\Chapter3\NCEP_Reanalysis\ncep_reanal_WLseries.mat'

%% Extract storm indices for ADCIRC gage 194
WL_inds = WL_inds_adj;
gage_194 = WL_inds(194,:);

%% Find storms common to:
% (1) storms passing within 200 km of gage 194
% (2) storms available in ADCIRC WL simulations at gage 194
[commonValues, idx1, idx2] = intersect(select_tcs, gage_194);
