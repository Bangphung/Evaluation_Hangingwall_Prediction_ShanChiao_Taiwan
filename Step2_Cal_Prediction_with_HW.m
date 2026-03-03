clc; clear; close all;
%% Path to the necessary programs
addpath("D:\My_Study\NCU_Study_pvb\ASK14_HW_model\Shanchiao_Hangingwall_Evaluation\gmm_lib");
addpath('D:\My_Study\NCU_Study_pvb\ASK14_HW_model\Shanchiao_Hangingwall_Evaluation\ShanChiaoFault');

M = 7.2; % Consider Maximum Magnitude of the Shanchiao Fault
%% Create a folder to store the figures
file_dir = fileparts(mfilename ('fullpath'));
if ~isempty(file_dir)
    file_dir = [file_dir,'\'];
end
dir_output = [file_dir, sprintf('Figures_/')];

if ~(exist(dir_output, 'dir')==7)
   mkdir(dir_output)
end
%% Load the scenario for prediction.
df = readtable("df_Dists_hw_alls.csv"); 
% Rx = nan(length(df.Rx), 1);
% id1 = df.Rx < 0; id2 = df.Rx > 0;
% Rx(id1) = abs(df.Rx(id1));
% Rx(id2) = -df.Rx(id2);
% df.Rx = [];
% % Remove Rx from the table df. 
% df.Rx = Rx; 

sitecond = 'soil';
%%: Rock site condition
if strcmp(sitecond,'rock')
    df.Vs30 = 760*ones(length(df.Rrup),1);
    df.Z1_0 = 999*ones(length(df.Rrup),1);
end

W = df.W(1);
Ztor = 0;
Dip = df.Dip(1);
Ry0 = 0; fas = 0; HW = 0;  region = 6; % input parameters for ASK14
lambda = -90; Fhw = 0; rgn = 1; flg_AS = 0; % input parameters for Ph20
FRO = 1; FSS = 0; FNO = 0; Finter = 0; Fintra = 0; Fas = 0; Fma = 0; FVS30 = 0; % input parameters for Ch20
model_id1 = 1; model_id2 = 2; model_id3 = 3; model_id4 = 4; model_id5 = 5;
Rjb = 0; 
Rx = 0;
% Period of consideration for the hanging wall model
periods = [-1, 0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, ...
               0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5];

for T = periods(2) 

    if T==-1
        Tp = 1;
    elseif T == 0
        Tp = 0.01;
    else
        Tp = T;
    end

    %% Cal GMMs Prediction
    SaASK14 = zeros(length(df.W),1);
    SaCh20 = zeros(length(df.W),1);
    SaP20 = zeros(length(df.W),1);
    f1_hw = zeros(length(df.W),1);
    f2_hw = zeros(length(df.W),1);
    f3_hw = zeros(length(df.W),1);
    f4_hw = zeros(length(df.W),1);
    f5_hw = zeros(length(df.W),1);
    
    for i = 1:length(df.W)
        f1_hw(i,1) = swus_hanging_wall_model(Tp, M, Dip, df.Rx(i), df.Rrup(i), df.Rjb(i), W, Ztor, model_id1);
        f2_hw(i,1) = swus_hanging_wall_model(Tp, M, Dip, df.Rx(i), df.Rrup(i), df.Rjb(i), W, Ztor, model_id2);
        f3_hw(i,1) = swus_hanging_wall_model(Tp, M, Dip, df.Rx(i), df.Rrup(i), df.Rjb(i), W, Ztor, model_id3);
        f4_hw(i,1) = swus_hanging_wall_model(Tp, M, Dip, df.Rx(i), df.Rrup(i), df.Rjb(i), W, Ztor, model_id4);
        f5_hw(i,1) = swus_hanging_wall_model(Tp, M, Dip, df.Rx(i), df.Rrup(i), df.Rjb(i), W, Ztor, model_id5);
    
        SaASK14(i,1) = ASK_2014_nga(M, T, df.Rrup(i), Rjb, Rx, Ry0, Ztor, Dip, lambda, fas, HW, 999, df.Z1_0(i), df.Vs30(i), FVS30, region);
        SaCh20(i,1) = Chao19H(T, M, df.Rrup(i), Ztor, df.Vs30(i), df.Z1_0(i), FRO, FSS,FNO,Finter,Fintra,Fas,Fma,FVS30);
        SaP20(i,1) = Phung_2019h_NGAw2_TW(M, T, df.Rrup(i), Rjb, Rx, Ztor, Dip, Fhw, lambda, df.Vs30(i), df.Z1_0(i), rgn, flg_AS);
    end
    
    Pred_Arr = [SaASK14, SaCh20, SaP20, f1_hw, f2_hw, f3_hw, f4_hw, f5_hw];
    tabPred_Arr = array2table(Pred_Arr,'VariableNames',...
           {'SA_ASK14','SA_Ch20','SA_Ph20','f1_hw','f2_hw','f3_hw','f4_hw','f5_hw'});
    df_tabPredArr = [tabPred_Arr df];
    
    % Save data. 
    if strcmp(sitecond,'rock')
        writetable(df_tabPredArr, ['dfhwPredtionGrid_rSA_M',num2str(M),'_',num2str(T),'s.csv']);
    else
        writetable(df_tabPredArr, ['dfhwPredtionGrid_vzSA_M',num2str(M),'_',num2str(T),'s.csv']);
    end
    
end
    

%% Load the coordinate of ShanChiao Faults
% Read fault data
fault_filename = 'ShanChiaoFault_NCDR.xlsx';
dT = readtable(fault_filename, 'Sheet', 'SC');
dT = dT(strcmp(dT.Entity, '3dFace'), :);  % Filter 3dFace

% Unique tri_id
tri_ids = unique(dT.tri_id);
% Preallocate faces
faces = cell(length(tri_ids), 1);
surface_points = [];
for i = 1:length(tri_ids)
    idx = dT.tri_id == tri_ids(i);
    verts = [dT.Latitude(idx), dT.Longitude(idx), dT.Z(idx)];
    % Triangulate quad: two triangles
    faces{i} = {verts(1:3, :), verts([1 3 4], :)};
    % Surface points (Z==0)
    surface_idx = verts(:,3) == 0;
    surface_points = [surface_points; verts(surface_idx, 1:2)];
end


%% Combine HW effects:
    model_collects = {'ASK14', 'Ch20', 'Ph20'};
    
    Sa_pred = [SaASK14, SaCh20, SaP20];

for k = 1

    model = model_collects{k};

    fig_name = [model,'_T=', num2str(T)];
    
    Sahw = Sa_pred(:,k).*exp(f1_hw);
    
    %%
    % Assuming 'faces' is precomputed from calculate_fault_distances or similar function
    fig1 = figure('Name', 'Shanchiao Fault with ASK14 Prediction', 'NumberTitle', 'off', ...
                  'Position', [100 100 600 600]);
    
    % Create a map axes with Mercator projection
    ax = axesm('mercator', 'MapLatLimit', [24.0 25.2], 'MapLonLimit', [120.0 122.0]);
    hold(ax, 'on');
    
    % Define projected CRS (TWD97 TM2, EPSG:3826)
    proj = projcrs(3826);
    
    % Convert ASK14 site coordinates to projected coordinates
    [site_x, site_y] = projfwd(proj, df.StaLat, df.StaLon);
    % Convert surface points
   

    % Convert fault faces to projected coordinates (already in TWD97 TM2, so use directly)
    fault_x = cell(length(faces), 1);
    fault_y = cell(length(faces), 1);
    for f = 1:length(faces)
        fault_x{f} = cell(2, 1);
        fault_y{f} = cell(2, 1);
        for t = 1:2
            tri = faces{f}{t};  % [X Y Z] in projected coordinates
            fault_x{f}{t} = tri(:, 1);  % X coordinates
            fault_y{f}{t} = tri(:, 2);  % Y coordinates
        end
    end
    
    % Overlay ASK14 prediction map
    scatter(ax,  df.StaLon, df.StaLat, 10 * ones(length(df.W), 1), 0*log(Sahw), 'o', 'filled'); hold on;
    scatter(ax,  df.StaLon(df.Rx<=0), df.StaLat(df.Rx<=0), 10 * ones(length(df.W(df.Rx<=0)), 1),...
                    df.Rx(df.Rx<=0), 'gr*'); hold on;
    % Plot fault geometry using patch with triangulation

    colors = lines(length(faces));  % Color each face differently
    for f = 1:length(faces)
        for t = 1:2
            % Create patch data for fault triangles
            x_data = fault_x{f}{t};
            y_data = fault_y{f}{t};
            patch(ax, y_data, x_data,  0 * x_data, 'blue', ...
                  'FaceAlpha', 0.01, 'EdgeColor', 'k', 'LineWidth', 0.3, ...
                  'DisplayName', sprintf('Face %d', f));
        end
    end
    
    % Set title and other properties
    % title(ax, [model, ': T=', num2str(T), 's, M=7.0, REV, Ztor = 0 km, Dip = 49°']);
    % c = colorbar(ax);
    % c.Label.String = "ln(IM), (g)";
    % c.Limits = [-6, 1];
    colormap(ax, 'jet');  % Use copper colormap
    % Add a basemap (approximate topographic effect with a grid)
    gridm('on');
    setm(ax, 'MLabelParallel', 'south');
    setm(ax, 'PLabelLocation', 1);
    
    % Legend
    % legend(ax, 'show', 'Location', 'best');
    
    % Save figure
    
    saveas(fig1, fullfile(pwd, 'Figures_', [fig_name, '.png']));
    
    % Turn hold off to clean up
    hold(ax, 'off');

end


%% 
% fig1 = figure('Name', 'Shanchiao Fault with ASK14 Prediction', 'NumberTitle', 'off', ...
%               'Position', [100 100 800 600]);
% geoscatter(df.StaLat, df.StaLon,  10*ones(length(df.W), 1), log(SaASK14hw), 'o','filled');
% 
% hold on;
% colors = lines(length(faces));  % Color each face differently
% for f = 1:length(faces)
%     for t = 1:2
%         tri = faces{f}{t};
%         lat_data = tri(:,1);
%         lon_data = tri(:,2);
%         % Use geoshow with FaceColor and transparency
%         geoshow(lat_data, lon_data, 'DisplayType', 'polygon', ...
%                 'FaceColor', colors(f,:), 'FaceAlpha', 0.7, ...
%                 'EdgeColor', 'k', 'LineWidth', 0.8, ...
%                 'DisplayName', sprintf('Face %d', f));
%     end
% end
% title(['ASK14: T=',num2str(T),'s, M=7.0, REV, Ztor = 0 km, and Dip = 49^o']);
% % colormap("copper");
% c = colorbar;
% c.Label.String = "ln(IM), (g)";
% c.Limits = [-6, 1];
% geobasemap topographic;
% geolimits([24.00 25.2],[120.0 122.0]);
% 
% % Legend
% legend('show', 'Location', 'best');
% 
% % Save figure
% fig_name = ['ASK14_T=', num2str(T)];
% saveas(fig1, fullfile(pwd, 'Figures_', [fig_name, '.png']));
% 
% % Hold off to prevent further modifications
% hold off;
% 
% 
% fig2=figure;
% geoscatter(df.StaLat, df.StaLon,  10*ones(length(df.W), 1), log(SaP20hw), 'o','filled');
% title(['Ch20: T=',num2str(T),'s, M=7.0, REV, Ztor = 0 km, and Dip = 49^o']);
% % colormap("copper");
% c = colorbar;
% c.Label.String = "ln(IM), (g)";
% c.Limits = [-6, 1];
% geobasemap topographic;
% geolimits([24.00 25.2],[120.0 122.0]);
% 
% fig_name = ['Ch20_T=',num2str(T)];
% saveas(fig2,[pwd,'\Figures_\',fig_name,'.png']);
% 
% 
% fig3=figure;
% geoscatter(df.StaLat, df.StaLon,  10*ones(length(df.W), 1), log(SaP20hw), 'o','filled');
% title(['Ph20: T=',num2str(T),'s, M=7.0, REV, Ztor = 0 km, and Dip = 49^o']);
% % colormap("copper");
% c = colorbar;
% c.Label.String = "ln(IM), (g)";
% c.Limits = [-6, 1];
% geobasemap topographic;
% geolimits([24.00 25.2],[120.0 122.0]);
% 
% fig_name = ['Ph20_T=',num2str(T)];
% saveas(fig3,[pwd,'\Figures_\',fig_name,'.png']);
