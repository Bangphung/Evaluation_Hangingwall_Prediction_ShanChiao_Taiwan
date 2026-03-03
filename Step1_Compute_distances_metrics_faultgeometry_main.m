clc; clear; close all;

%% Load all the sites that produces VS30 and Z1.0
df = readtable('MapGrid.csv'); 
df = df(df.StaLat>=24.46, :); % Consider the Northern Area to reduce the computation time.

% Plot to check the all site coordinates
figure;
geoplot(df.StaLat, df.StaLon, '.'); hold on;

%% Load fault Geometry of the Shanchiao Fault.
    addpath('D:\My_Study\NCU_Study_pvb\ASK14_HW_model\Shanchiao_Hangingwall_Evaluation\ShanChiaoFault');
    
    fault_filename = 'ShanChiaoFault_NCDR.xlsx';

%% All sites grid

    site_lat_lon_all = [df.StaLat, df.StaLon];

%% Compute all distances metrics and fault geometry
    Dist_alls = [];
for i = 1:length(df.StaLat)
    [Rrup, Rjb, Rx, W, Ztor, dip] = Step1_Compute_Dist_Metrics_fault_func(fault_filename, df.StaLat(i), df.StaLon(i));
     Dist_alls = cat(1, Dist_alls, [Rrup, Rjb, Rx, W, Ztor, dip]);
end
    Dist_alls_tab = array2table(Dist_alls, "VariableNames", {'Rrup','Rjb','Rx','W','Ztor','Dip'});

%% Merge two tables into One larger table.
    df_dist_alls = [Dist_alls_tab df];
    % Save the distance metrics to create the flatfile.
    writetable(df_dist_alls, 'df_Dists_hw_alls.csv');