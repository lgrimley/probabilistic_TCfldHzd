%% Clear workspace and set up paths
clear   % remove all variables from workspace
clc     % clear command window

% Add path for functions used in bias correction workflow
addpath 'C:\Users\lelise\Documents\GitHub\flood_model_carolinas\syntheticTCs_cmpdfld\tracks\bias_correction'

% Change directory to location of bias correction data
cd Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\BiasCorrection

%% Map numeric IDs to GCM names
% Used to loop through multiple GCMs
model_map = containers.Map({1, 2, 3, 4, 5}, {'canesm', 'cnrm6', ...
    'ecearth6', 'ipsl6', 'miroc6'});

% rps = 5:500;  % (commented out) could be used for sample sizes

%% Load historical NCEP landfall wind speeds
temp = readtable("vmax/NCEP/ncep_landfall_vmax_ZerosRemoved.csv");
xoh = temp.vstore100;  % observed historical Vmax
xoh_tot = length(xoh); % total number of observations

% Fit a statistical distribution (GPD/KDE) to observed Vmax
[xoh_fit, ~, oh_type] = fit_gpd(xoh, 15, 4); 
xcomp = min(xoh):1:max(xoh)*1.2;  % range for plotting or evaluation

%% Initialize variables for loop
iter = 1;
G_comp = [];              % placeholder (unused in this script)
total = 0;                % placeholder
st_weights_hist = {};     % cell array to store historical storm weights
st_weights_proj = {};     % cell array to store projected storm weights

%% Loop over each GCM model
% Currently set to run all five models, but only one may be used
for i = 1:length(model_map)
    m = i; 
    % name = strcat(model_map(i), '_20thcal_TCchars');  % optional naming

    %% Load historical modeled tracks from GCM
    fname = sprintf('vmax/CMIP6_20thcal/UScoast6_AL_%s_20thcal_roEst1rmEst1_trk100.mat_select_vstore100.csv', ...
        model_map(i));
    temp = readtable(fname);
    xmh = temp.vstore100;  % modeled historical Vmax
    inds = temp.tc_id;     % corresponding TC IDs
    xmh_tot = length(xmh);
    
    % Skip if too few storms
    if xmh_tot < 50
        continue
    end

    % Fit distribution to modeled historical data
    [xmh_fit, ~, mh_type] = fit_gpd(xmh, 15, 4); 

    %% Load future GCM tracks
    fname = sprintf('vmax/CMIP6_585/%s_ssp585_landfall_vmax_ZerosRemoved.csv', ...
        model_map(i));
    temp = readtable(fname);
    tc_ids = temp.tc_id;
    xmp = temp.vstore100;  % modeled projected (future) Vmax
    xmp_tot = length(xmp); 

    % Fit distribution to projected GCM data
    [xmp_fit, ~, mp_type] = fit_gpd(xmp, 15, 4); 

    %% Bias correction using Quantile Delta Mapping (QDM)
    % Correct future modeled Vmax based on observed vs historical modeled differences
    % xphat = bias-corrected future
    % xhhat = bias-corrected historical modeled
    [xphat, xhhat] = QDM(xoh, xmh, xmp, 4);

    % Fit distribution to bias-corrected future data
    [xphat_fit, ~, phat_type] = fit_gpd(xphat, 15, 4); 

    %% Calculate weights for historical storms
    % Adjust storm contributions so that modeled historical matches observed
    if strcmp(oh_type,'kde') || strcmp(mh_type,'kde')
        bw = (max(xmh)-min(xmh))/12;  % bin width for histogram
        edges = 0:round(bw):max(xmh)+round(bw); 
        c_mh = histcounts(xmh, 'BinEdges', edges, 'Normalization', 'probability');  % modeled histogram
        c_oh = histcounts(xoh, 'BinEdges', edges, 'Normalization', 'probability');  % observed histogram

        tau = 999*ones(length(c_oh),1);  % initialize array
        tau(c_mh~=0)  = c_oh(c_mh~=0)./c_mh(c_mh~=0);  % ratio of obs to modeled
        tau(tau==999) = mean(tau(c_mh~=0));  % fill gaps with mean
        weight = tau./sum(tau);  % normalize weights

        wi_h = zeros(xmh_tot, 1);
        % Loop through each storm and assign weight based on Vmax bin
        for j = 1:xmh_tot
            ind = find(xmh(j)<edges, 1); 
            wi_h(j) = weight(max(ind-1,1)); 
        end
        wi_h = wi_h./sum(wi_h);  % normalize again
    else
        % If GPD fit was used instead of KDE
        tau = pdf(xoh_fit, xmh)./pdf(xmh_fit, xmh);
        tau(tau==0) = min(tau(tau>0));
        wi_h = tau./sum(tau);
    end

    %% Calculate weights for projected (future) storms
    st_weights_hist{iter} = wi_h;  % save historical weights

    if strcmp(mp_type,'kde') || strcmp(phat_type,'kde')
        bw = (max(xmp)-min(xmp))/12; 
        edges = 0:round(bw):max(xmp)+round(bw); 
        c_mp = histcounts(xmp, 'BinEdges', edges, 'Normalization', 'probability'); 
        c_phat = histcounts(xphat, 'BinEdges', edges, 'Normalization', 'probability');
        tau = 999*ones(length(c_phat),1);
        tau(c_mp~=0)  = c_phat(c_mp~=0)./c_mp(c_mp~=0);
        tau(tau==999) = mean(tau(c_mp~=0));
        weight = tau./sum(tau);

        wi_p = zeros(xmp_tot, 1);
        for j = 1:xmp_tot
            ind = find(xmp(j)<edges, 1); 
            wi_p(j) = weight(max(ind-1,1)); 
        end
        wi_p(wi_p == 0) = 1*10^-5; 
        wi_p = wi_p./sum(wi_p);  % normalize weights
    else
        tau = pdf(xphat_fit, xmp)./pdf(xmp_fit, xmp); 
        tau(tau==0) = min(tau(tau>0));
        wi_p = tau./sum(tau);
    end

    st_weights_proj{iter} = wi_p;  % save projected weights

    %% Save bias-corrected Vmax and weights to CSV
    out = [tc_ids, xmp, wi_p];
    fname = sprintf('%s_ssp585_weighted.csv', model_map(i));
    writematrix(out, fname);

    iter = iter + 1;  % increment iterator for next model
end

%% Save all storm weights
% Historical and projected weights correspond to the order of TC IDs
save('storm_weights.mat', 'st_weights_hist', 'st_weights_proj')
