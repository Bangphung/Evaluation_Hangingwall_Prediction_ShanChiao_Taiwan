clc; clear; close all;
% MATLAB script to plot the 3D geometry and surface projection of the Shanchiao Fault
addpath('D:\My_Study\NCU_Study_pvb\ASK14_HW_model\Shanchiao_Hangingwall_Evaluation\ShanChiaoFault');
% Read the data from the Excel file
filename = 'ShanChiaoFault_NCDR.xlsx';
sheet = 'SC';
dT = readtable(filename, 'Sheet', sheet);

% Unique tri_id
tri_ids = unique(dT.tri_id);
% Preallocate faces
faces = cell(length(tri_ids), 1);
for i = 1:length(tri_ids)
    idx = dT.tri_id == tri_ids(i);
    verts = [dT.Latitude(idx), dT.Longitude(idx), dT.Z(idx)];
    % Triangulate quad: two triangles
    faces{i} = {verts(1:3, :), verts([1 3 4], :)};
    % Surface points (Z==0)
end

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
% Considered period, T
 T = -1;  
% Give the coordinate of NPP sites
 xyNPP1 = [25.2861, 121.5879];
 xyNPP2 = [25.2026, 121.6627];
 xyNPP4 = [25.0389, 121.9247];
 xyPoint = [24.96, 121.91];

% Load the prediction table
sitecond = 'soil';
%%: Rock site condition
if strcmp(sitecond,'rock')
    df = readtable(['dfhwPredtionGrid_rSA_M7.06_',num2str(T),'s.csv']); 
else
    df = readtable(['dfhwPredtionGrid_vzSA_M7.06_',num2str(T),'s.csv']); 
end

% Model collection
model_collects = {'ASK14', 'Ch20', 'Ph20'};
% The corresponding prediction without HW
Sa_pred_arr = [df.SA_ASK14, df.SA_Ch20, df.SA_Ph20];

pred_AB = []; pred_CD = [];
% Load 
load("TaipeiBasement.mat");

for k = 1:3
    model = model_collects{k};

    fig_name = [model,'_T=', num2str(T)];
  
    % Add the map of the prediction
    Sa_pred = Sa_pred_arr(:,k).*exp(df.f1_hw); % consider the first hanging wall model.   
    % Next, plot the surface projection (Z = 0)
    fig1 = figure('Position', [50, 50, 600, 500]);
    title('Surface Projection of Shanchiao Fault (Z = 0)');
    geoscatter(df.StaLat, df.StaLon,  5*ones(length(df.W), 1), log(Sa_pred), 'o','filled');
    hold on; 
    % Plot the counterlines of Taipei
    geoplot(TaipeiBasement.StaLat, TaipeiBasement.StaLon,'k','Color',[0.8,0.5,0.5]); hold on;
    for f = 1:length(faces)
            for t = 1:2
                tri = faces{f}{t};  % [X Y Z] in projected coordinates
                x_data = tri(:, 1);  % X coordinates
                y_data = tri(:, 2);  % Y coordinates
                geoplot(x_data, y_data, 'k','Color',[0.5,0.5,0.5]); hold on;
            end
    end
    % Add the location of NPPs sites. 
    geoplot(xyNPP1(1), xyNPP1(2),'k^','MarkerFaceColor','r'); hold on;
    text(xyNPP1(1), xyNPP1(2),'NPP1'); 
    geoplot(xyNPP2(1), xyNPP2(2),'k^','MarkerFaceColor','b');  hold on;
    text(xyNPP2(1), xyNPP2(2),'NPP2'); 
    geoplot(xyNPP4(1), xyNPP4(2),'k^','MarkerFaceColor','c'); hold on;
    text(xyNPP4(1), xyNPP4(2),'NPP4'); 
    % Set title and other properties.
    if T==0.01
        title([model, ':PGA, M=7.06, NM, Ztor = 0 km, Dip = 49°'], 'FontSize', 14);
    elseif T==-1
        title([model, ':PGV, M=7.06, NM, Ztor = 0 km, Dip = 49°'], 'FontSize', 14);
    else
        title([model, ':SA(T=', num2str(T), 's), M=7.06, NM, Ztor = 0 km, Dip = 49°'], 'FontSize', 14);
    end
    c = colorbar;
    c.Label.String = "ln(IM), (g)";
    % c.Limits = [-6, 1];
    colormap('jet');  % Use copper colormap
    % Add a line perpendicular to the fault strike.
    % Point A:  % Point B:
    
    Acoord = xyNPP1; % 
    Bcoord = xyPoint;% 
    Ccoord = [25.12, 121.34]; 
    Dcoord = [24.79, 121.82];

    [interp_SaAB, interp_RxAB, interp_Vs30, latlon_line] = Interpolate2D(Acoord, Bcoord, df, Sa_pred);
    [interp_SaCD, interp_RxCD, ~, latlon_line_CD] = Interpolate2D(Ccoord, Dcoord, df, Sa_pred);

    pred_AB = cat(2,pred_AB, interp_SaAB); 
    pred_CD = cat(2,pred_CD, interp_SaCD); 

    geoplot(latlon_line(:,1), latlon_line(:,2), 'gr--','LineWidth',1); hold on;
    geoplot(latlon_line_CD(:,1), latlon_line_CD(:,2), 'gr--','LineWidth',1); hold on;
    text(Acoord(1)-0.02, Acoord(2)-0.02,'A','FontSize',14); text(Bcoord(1), Bcoord(2)-0.02,'B','FontSize',14);
    text(Ccoord(1)-0.02, Ccoord(2)-0.02,'C','FontSize',14); text(Dcoord(1), Dcoord(2)-0.02,'D','FontSize',14);
    geolimits([24.7 25.5],[121.0 122.2]);
    c.Limits = [-3.3, 1];

    saveas(fig1, fullfile(pwd, 'Figures_', [fig_name, '.png']));
end

if strcmp(sitecond,'rock')
    predprofile.pred_AB = pred_AB;
    predprofile.pred_CD = pred_CD;
    save(['predprofile_T',num2str(T),'s.mat'],"predprofile");
end


%% Plot the hanging wall effects.
% Prediction with VS30 = 760 m/s;
filename = ['predprofile_T',num2str(T),'s.mat'];
load(filename);

fig2 = figure('Position', [50, 50, 1200, 500]);
fig_name = ['HW_pred_alongAB_T',num2str(T),'s.csv'];
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
h1=plot(interp_RxAB, pred_AB, '-','LineWidth',2); hold on;
plot(interp_RxAB, predprofile.pred_AB(:,1), '--','color',h1(1,1).Color,'LineWidth',1); hold on;
plot(interp_RxAB, predprofile.pred_AB(:,2), '--','color',h1(2,1).Color,'LineWidth',1); hold on;
plot(interp_RxAB, predprofile.pred_AB(:,3), '--','color',h1(3,1).Color,'LineWidth',1); 
legend('ASK14(V_{30}^S, Z_{1.0})','Ch20(V_{30}^S, Z_{1.0})','Ph20(V_{30}^S, Z_{1.0})',...
      'ASK14(V_{30}^S=760m/s)','Ch20(V_{30}^S=760m/s)','Ph20(V_{30}^S=760m/s)','Location','best');
xlabel('R_X (km)'); grid on;
ylabel('Pred SA (g) along AB profile');
if T==0.01
    title('PGA at AB profile');
elseif T ==-1
    title('PGV at AB profile');
else 
    title(['SA (T=', num2str(T), 's) at AB profile']);
end

nexttile;
h2=plot(interp_RxCD, pred_CD, '-','LineWidth',2); hold on;
plot(interp_RxCD, predprofile.pred_CD(:,1), '--','color',h2(1,1).Color,'LineWidth',1); hold on;
plot(interp_RxCD, predprofile.pred_CD(:,2), '--','color',h2(2,1).Color,'LineWidth',1); hold on;
plot(interp_RxCD, predprofile.pred_CD(:,3), '--','color',h2(3,1).Color,'LineWidth',1); 
legend('ASK14(V_{30}^S, Z_{1.0})','Ch20(V_{30}^S, Z_{1.0})','Ph20(V_{30}^S, Z_{1.0})',...
      'ASK14(V_{30}^S=760m/s)','Ch20(V_{30}^S=760m/s)','Ph20(V_{30}^S=760m/s)','Location','best');
ylabel('Pred SA (g) along CD profile'); grid on;
xlabel('R_X (km)');
if T==0.01
    title('PGA at CD profile');
elseif T ==-1
    title('PGV at CD profile');
else 
    title(['SA (T=', num2str(T), 's) at CD profile']);
end

if T==0.01
    title(t, 'PGA, M=7.06, NM, Ztor = 0 km, Dip = 49°', 'FontSize', 16);
elseif T==-1
    title(t, 'PGV, M=7.06, NM, Ztor = 0 km, Dip = 49°', 'FontSize', 16);
else
    title(t, ['SA(T=', num2str(T), 's), M=7.06, NM, Ztor = 0 km, Dip = 49°'], 'FontSize', 16);
end

saveas(fig2, fullfile(pwd, 'Figures_', [fig_name, '.png']));


%%
function [interp_Sa, interp_Rx, interp_Vs30, latlon_line] = Interpolate2D(Ecoord, Scoord, df, Sa_pred)
    % Find the prediction along a line determined by two points
    % Generate N interpolated points between the two coordinates
    N = 500;
    lat_line = linspace(Ecoord(1), Scoord(1), N);
    lon_line = linspace(Ecoord(2), Scoord(2), N);
    latlon_line = [lat_line', lon_line'];
    % Create interpolant function
    F = scatteredInterpolant(df.StaLat, df.StaLon, log(Sa_pred), 'linear', 'none');
    FRx = scatteredInterpolant(df.StaLat, df.StaLon, df.Rx, 'linear', 'none');
    FVS30 = scatteredInterpolant(df.StaLat, df.StaLon, df.Vs30, 'linear', 'none');
    % Interpolate values along the line
    interp_Sa = exp(F(lat_line', lon_line'));
    interp_Rx = FRx(lat_line', lon_line');
    interp_Vs30 = FVS30(lat_line', lon_line');
end