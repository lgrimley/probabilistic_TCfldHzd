%==========================================================================
% TCRM Rainfall Post-Processing Script
%
% Purpose:
%   Converts dynamic TCRM rainfall output grids into static, user-defined
%   spatial grids clipped to a specified bounding box. Each storm is
%   processed independently and written to a NetCDF file suitable for
%   hydrodynamic or flood modeling applications (e.g., SFINCS).
%
% Key Features:
%   - Processes multiple CMIP6 GCMs under SSP585
%   - Embeds moving TCRM rainfall grids into a large static grid
%   - Clips rainfall to a fixed geographic domain
%   - Writes gridded rainfall time series to NetCDF (one file per storm)
%
% Inputs:
%   - TCRM rainfall .mat files containing:
%       * PLAT_SAVE   (latitude grid per timestep)
%       * PLONG_SAVE  (longitude grid per timestep)
%       * rain        (rainfall rate per timestep)
%   - CSV files listing selected storm IDs
%
% Outputs:
%   - NetCDF files containing:
%       * precip(x, y, time)  rainfall (mm)
%       * x, y                longitude and latitude coordinates
%       * time                hourly timestep index
%
% Assumptions:
%   - Rainfall grid spacing is uniform and matches user-defined increments
%   - Storm rainfall ends when total rainfall becomes zero
%   - Longitude convention is converted from 0–360 to −180–180
%
% Author:
%   Lauren Grimley and Avantika Gori
%
%==========================================================================

clear
clc
% Clear workspace variables and command window
% Script for post-processing TCRM rainfall into a static user-defined grid


% ======================
% USER INPUT DESCRIPTION
% ======================
% bb: [LonMin, LatMin;
%      LonMax, LatMax] bounding box for spatial clipping
% lat_inc : latitude grid increment (degrees)
% lon_inc : longitude grid increment (degrees)
% rain_folder: directory containing TCRM rainfall .mat files
% outfolder: directory where output NetCDF files will be written
% Ns: number of storms to process


% =========
% OUTPUTS
% =========
% Rainfall NetCDF files:
%   - One NetCDF per storm
%   - Contains gridded rainfall time series on a fixed spatial grid
%   - Rainfall is clipped to the user-defined bounding box


% Specify bounding box (SFINCS domain for NC/SC)
bb = [-83.67,31.95;-75.11,36.81];

% List of CMIP6 GCMs under SSP585 scenario
gcm = {'canesm_ssp585'; 'cnrm6_ssp585'; 'ecearth6_ssp585'; 'ipsl6_ssp585'; 'miroc6_ssp585'};

% Change working directory to project data location
cd 'Z:\Data-Expansion\users\lelise\Chapter3\CMIP6_585';


% Loop over each GCM
for kk = 1:length(gcm)

    % Define input folder containing TCRM rainfall outputs
    inputfolder = fullfile('.\rain', strcat(gcm{kk}, 'cal'));

    % Define output folder for gridded rainfall NetCDFs
    outputfolder = fullfile('.\rain', strcat('TCR_Gridded_' ,gcm{kk}, 'cal'));

    % Read table of selected storm IDs and associated metadata
    selected_tracks = readtable(fullfile('.\stormTide', strcat('stormTide_TCIDs_and_gageCounts_', gcm{kk}, '.csv')));

    % Create output directory if it does not exist
    if ~isfolder(outputfolder)
        mkdir(outputfolder);
    end


    % ==========================
    % Create global static grid
    % ==========================
    lat_inc = 0.05;
    lon_inc = 0.05;

    % Longitude and latitude vectors covering bounding box
    xvec = floor(bb(1,1)):lon_inc:ceil(bb(2,1));
    yvec = floor(bb(1,2)):lat_inc:ceil(bb(2,2));

    % Store bounding box limits
    min_lon = xvec(1);
    min_lat = yvec(1);
    max_lon = xvec(end);
    max_lat = yvec(end);


    % =====================
    % Loop through storms
    % =====================
    for ii = 1:length(selected_tracks.tc_id)

        % Extract storm ID
        stnum = selected_tracks.tc_id(ii);

        % Construct output NetCDF filename using zero-padded storm number
        if stnum < 10
            stormfile = strcat('000',num2str(stnum),'.nc');
        elseif stnum >= 10 && stnum < 100
            stormfile = strcat('00', num2str(stnum), '.nc');
        elseif stnum >= 100 && stnum < 1000
            stormfile = strcat('0', num2str(stnum), '.nc');
        else
            stormfile = strcat(num2str(stnum), '.nc');
        end

        % Full path to output NetCDF file
        f = sprintf('%s%s',outputfolder,stormfile);

        % Skip processing if this storm has already been processed
        if isfile(f) == 0
            disp(stnum)

            % Construct input MAT filename (zero-padded)
            if stnum < 10
                stormfile = strcat('000',num2str(stnum),'.mat');
            elseif stnum >= 10 && stnum < 100
                stormfile = strcat('00', num2str(stnum), '.mat');
            elseif stnum >= 100 && stnum < 1000
                stormfile = strcat('0', num2str(stnum), '.mat');
            else
                stormfile = strcat(num2str(stnum), '.mat');
            end

            % Load TCRM rainfall output
            % Loads: PLAT_SAVE, PLONG_SAVE, rain
            filename = fullfile(inputfolder,stormfile);
            load(filename)

            % Convert longitude convention from 0–360 to -180–180
            PLONG_SAVE(PLONG_SAVE>180) = PLONG_SAVE(PLONG_SAVE>180)-360;

            % Determine number of time steps
            [~, ~, trk] = size(rain);

            % Start time index (skip early spin-up timesteps)
            tstart = 4;


            % ===================================================
            % Embed dynamic TCRM grid into a large static grid
            % ===================================================
            Lgrid_lon = -100:lon_inc:-50;
            Lgrid_lat = 0:lat_inc:60;

            % Create meshgrid for plotting/reference
            % Longitude varies along rows, latitude along columns
            [lYn, lXn] = meshgrid(Lgrid_lat, Lgrid_lon);

            % Initialize large rainfall grid
            Lgrid = zeros(length(Lgrid_lon), length(Lgrid_lat));

            % Loop through each timestep
            for j = tstart:trk

                % Extract TCRM grid at timestep j
                tcr_lat = PLAT_SAVE(1,:,j);
                tcr_lon = PLONG_SAVE(:,1,j);
                tcr_rain = rain(:,:,j);

                % Stop processing if rainfall is zero everywhere
                if sum(tcr_rain, 'all') == 0
                    break
                else
                    % Find nearest indices in large grid
                    [~, J] = min(abs(Lgrid_lon - tcr_lon(1,1)));
                    [~, I] = min(abs(Lgrid_lat- tcr_lat(1,1)));

                    % Determine index ranges for inserting dynamic grid
                    lat_inds = I:I+size(tcr_lat,2)-1;
                    lon_inds = J:J+size(tcr_lon,1)-1;

                    % Insert rainfall into the larger static grid
                    Lgrid(lon_inds, lat_inds, j) = tcr_rain;
                end
            end


            % ===================================
            % Clip grid to user-defined domain
            % ===================================
            lat_idx = find(Lgrid_lat >= min_lat & Lgrid_lat <= max_lat);
            lon_idx = find(Lgrid_lon >= min_lon & Lgrid_lon <= max_lon);

            % Extract clipped rainfall and coordinates
            clipped_data = Lgrid(lon_idx, lat_idx, :);
            clipped_latitudes = lYn(lat_idx, lat_idx);
            clipped_longitudes = lXn(lon_idx, lon_idx);


            % ===================================
            % Write rainfall grid to NetCDF
            % ===================================
            x = clipped_longitudes(:,1);
            y = clipped_latitudes(1,:);
            [dim1, dim2, dim3] = size(clipped_data);

            % Time coordinate (hourly index)
            time = linspace(1, dim3, dim3);

            % Create NetCDF file
            [~, name, ~] = fileparts(filename);
            outFile = fullfile(outputfolder, [name, '.nc']);
            ncid = netcdf.create(outFile, 'NETCDF4');

            % Define dimensions
            dimid_x = netcdf.defDim(ncid, 'x', length(x));
            dimid_y = netcdf.defDim(ncid, 'y', length(y));
            dimid_time = netcdf.defDim(ncid, 'time', length(time));

            % Define variables
            varid_x = netcdf.defVar(ncid, 'x', 'double', dimid_x);
            varid_y = netcdf.defVar(ncid, 'y', 'double', dimid_y);
            varid_time = netcdf.defVar(ncid, 'time', 'double', dimid_time);
            varid_precip = netcdf.defVar(ncid, 'precip', 'double', [dimid_x, dimid_y, dimid_time]);

            % Add variable attributes
            netcdf.putAtt(ncid, varid_x, 'units', 'degrees_east');
            netcdf.putAtt(ncid, varid_y, 'units', 'degrees_north');
            netcdf.putAtt(ncid, varid_time, 'units', 'hours');
            netcdf.putAtt(ncid, varid_precip, 'units', 'mm');

            % Exit define mode
            netcdf.endDef(ncid);

            % Write data
            netcdf.putVar(ncid, varid_x, x);
            netcdf.putVar(ncid, varid_y, y);
            netcdf.putVar(ncid, varid_time, time);
            netcdf.putVar(ncid, varid_precip, clipped_data);

            % Close NetCDF file
            netcdf.close(ncid);

        else
            % Skip storm if output already exists
            continue
        end
    end
end
