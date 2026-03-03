clc; clear; close all;

%% Create a folder to save figures
file_dir = fileparts(mfilename ('fullpath'));
if ~isempty(file_dir)
    file_dir = [file_dir,'\'];
end
dir_output = [file_dir, sprintf('Figures_/')];

if ~(exist(dir_output, 'dir')==7)
   mkdir(dir_output)
end

%% MATLAB script to plot the 3D geometry and surface projection of the Shanchiao Fault
addpath("D:\My_Study\NCU_Study_pvb\ASK14_HW_model\Shanchiao_Hangingwall_Evaluation\gmm_lib");
addpath("D:\My_Study\NCU_Study_pvb\ASK14_HW_model\Shanchiao_Hangingwall_Evaluation\ShanChiaoFault");
% Read the data from the Excel file
filename = 'ShanChiaoFault_NCDR.xlsx';
sheet = 'SC';
T = readtable(filename, 'Sheet', sheet);
% Get unique tri_id values
% Unique tri_id
tri_ids = unique(T.tri_id);

% Preallocate faces
faces = cell(length(tri_ids), 1);
all_points = [];
surface_points = [];

for i = 1:length(tri_ids)
    idx = T.tri_id == tri_ids(i);
    verts = [T.X(idx), T.Y(idx), T.Z(idx)];
    % Triangulate quad: two triangles
    faces{i} = {verts(1:3, :), verts([1 3 4], :)};
    % Collect all points for plane fitting
    all_points = [all_points; verts];
    % Surface points (Z==0)
    surface_idx = verts(:,3) == 0;
    surface_points = [surface_points; verts(surface_idx, 1:2)];
end
surface_points = unique(surface_points, 'rows');

% MATLAB Script to Estimate Magnitude and Plot Strike/Perpendicular Lines
% Parameters for Shan-chiao fault
nu = 49;          % Fault dip in degrees
Ztor = 0;         % Top of rupture depth (km)
ZBSZ = 25;        % Seismogenic basement depth (km)
region = 'TW';    % Region: Taiwan

% Step 1: Estimate A, L, W from fault geometry
all_xy = [];
for i = 1:length(tri_ids)
    idx = T.tri_id == tri_ids(i);
    verts = [T.X(idx), T.Y(idx)];  % Ignore Z for projection
    all_xy = [all_xy; verts];
end
all_xy = unique(all_xy, 'rows'); % Remove duplicates

% Convex hull for projection boundary
hull_idx = convhull(all_xy(:,1), all_xy(:,2));
hull_xy = all_xy(hull_idx, :);
A0 = polyarea(hull_xy(:,1), hull_xy(:,2)) / 1e6; % Convert m^2 to km^2

% Step 2: Compute TRUE RUPTURE PLANE AREA (3D surface area)
rupture_plane_areas = [];  % Individual triangle areas
for i = 1:length(tri_ids)
    verts = [T.X(tri_ids(i) == T.tri_id), T.Y(tri_ids(i) == T.tri_id), T.Z(tri_ids(i) == T.tri_id)];
    % Two triangles per quad
    tri1 = verts(1:3, :);  % Triangle 1
    tri2 = verts([1 3 4], :);  % Triangle 2
    
    % Area of tri1: 0.5 * || (B-A) x (C-A) ||
    u1 = tri1(2, :) - tri1(1, :);
    v1 = tri1(3, :) - tri1(1, :);
    area1 = 0.5 * norm(cross(u1, v1));
    
    % Area of tri2
    u2 = tri2(2, :) - tri2(1, :);
    v2 = tri2(3, :) - tri2(1, :);
    area2 = 0.5 * norm(cross(u2, v2));
    
    rupture_plane_areas = [rupture_plane_areas, area1 / 1e6, area2 / 1e6];  % km^2
end

rupture_plane_area = sum(rupture_plane_areas);
A = rupture_plane_area;
fprintf('Rupture Plane Area (3D surface): %.2f km^2\n', rupture_plane_area);
fprintf('Number of triangles: %d\n', length(rupture_plane_areas));
fprintf('Average triangle area: %.4f km^2\n', mean(rupture_plane_areas));

% Approximate rupture length (L) and width (W)
% L as max extent along strike (fit line direction)
% Fit line to surface points for strike
[p, ~] = polyfit(surface_points(:,2), surface_points(:,1), 1);  % X vs Y
 slope = p(1); intercept = p(2);

L = max(max(all_xy(:,1)) - min(all_xy(:,1)), max(all_xy(:,2)) - min(all_xy(:,2))) / 1000; % km
% W as max depth / sin(dip) (from earlier W calculation)
% max_depth = -min(cellfun(@(x) min(x(:,3)), faces)) / 1000; % km
max_depth = -min(all_points(:,3))/1000;
W = max_depth / sind(nu); % km (approximate)

% Step 2: Estimate Magnitude using the provided functions
[M_A, M_L, M_W] = Huang_and_Abrahamson_2024_calculate_M(A, L, W, nu, Ztor, ZBSZ, 'TW');

[Mw_A, Mw_L, Mw_W] = Ye_and_Ma_2011_calculate_M(A, L, W);

% Display results
fprintf('Estimated Magnitude from Area (M_A): %.2f\n', M_A);
fprintf('Estimated Magnitude from Length (M_L): %.2f\n', M_L);
fprintf('Estimated Magnitude from Width (M_W): %.2f\n', M_W);

fprintf('Estimated Ye and Ma (2011) Magnitude from Area (M_A): %.2f\n', Mw_A);
fprintf('Estimated Ye and Ma (2011) Magnitude from Length (M_L): %.2f\n', Mw_L);
fprintf('Estimated Ye and Ma (2011) Magnitude from Width (M_W): %.2f\n', Mw_W);





% PLOT_FAULT_LINES Plot fault surface projection with strike and perpendicular lines
%
%   plot_fault_lines(fault_filename)
%
% Inputs:
%   fault_filename: Path to the ShanChiaoFault_NCDR.xlsx file
%
% Outputs:
%   Plots the fault surface projection with:
%   - Strike line (along fault strike, 100 points)
%   - Perpendicular line (normal to strike, 100 points)
%
% Notes:
%   - Improved surface projection using union of valid projected triangles
%   - Handles degenerate triangles with < 3 unique points
%   - Strike line fitted to surface points, extended over fault extent
%   - Perpendicular line centered at fault midpoint
%   - Uses 100 points for each line
%   - Requires Mapping Toolbox


% Unique tri_id
tri_ids = unique(T.tri_id);

% Collect surface points (Z == 0) for lines
surface_points = [];
% Collect projected triangles for union
projected_polys = cell(length(tri_ids) * 2, 1);  % For union
valid_polys = 0;
all_xy = [];  % Fallback for all projected points
for i = 1:length(tri_ids)
    idx = T.tri_id == tri_ids(i);
    verts = [T.X(idx), T.Y(idx), T.Z(idx)];  % [X Y Z]
    % Triangulate quad: two triangles
    tri1 = verts(1:3, 1:2);  % Projected X,Y
    tri2 = verts([1 3 4], 1:2);  % Projected X,Y
    
    % Check for valid triangles (at least 3 unique points)
    if size(unique(tri1, 'rows'), 1) >= 3
        projected_polys{valid_polys + 1} = polyshape(tri1(:,1), tri1(:,2));
        valid_polys = valid_polys + 1;
    end
    if size(unique(tri2, 'rows'), 1) >= 3
        projected_polys{valid_polys + 1} = polyshape(tri2(:,1), tri2(:,2));
        valid_polys = valid_polys + 1;
    end
    % Collect surface points
    surface_idx = verts(:, 3) == 0;
    surface_points = [surface_points; verts(surface_idx, 1:2)];
    % Collect all projected points for fallback
    all_xy = [all_xy; tri1; tri2];
end
projected_polys = projected_polys(1:valid_polys);  % Trim to valid polys

% Compute union of valid projected triangles for complete projection
if ~isempty(projected_polys)
    full_projection = union(projected_polys{1}, projected_polys{2});
    for i = 3:length(projected_polys)
        full_projection = union(full_projection, projected_polys{i});
    end
else
    % Fallback: Use convex hull of all projected points if no valid polys
    all_xy = unique(all_xy, 'rows');
    full_projection = polyshape(all_xy(:,1), all_xy(:,2));
end

% Fit line to surface points for strike
[p, ~] = polyfit(surface_points(:,2), surface_points(:,1), 1);  % X vs Y
slope = p(1); intercept = p(2);

% Fault midpoint for centering lines
mid_point = mean(surface_points, 1);
mid_x = mid_point(1);
mid_y = mid_point(2);

% Strike direction vector (unit vector along line)
strike_dir = [1, slope];  % [dx, dy] along strike
strike_dir = strike_dir / norm(strike_dir);

% Perpendicular direction vector
perp_dir = [-slope, 1];  % Rotate by 90 degrees
perp_dir = perp_dir / norm(perp_dir);

% Generate 100 points along strike line
line_length = 100e3;  % 100 km total length (adjust as needed)
t_strike = linspace(-line_length/2, line_length/2, 100);  % 100 points
strike_x = mid_x + t_strike * strike_dir(1);
strike_y = mid_y + t_strike * strike_dir(2);

% Generate 100 points along perpendicular line
t_perp = linspace(-line_length/2, line_length/2, 100);  % 100 points
perp_x = mid_x + t_perp * perp_dir(1);
perp_y = mid_y + t_perp * perp_dir(2);

% Create figure with map axes (Mercator projection)
figure('Name', 'Fault Surface Projection with Strike and Perpendicular Lines', ...
       'Position', [100 100 800 600]);
ax = axesm('mercator', 'MapLatLimit', [24.0 25.2], 'MapLonLimit', [120.0 122.0]);
hold(ax, 'on');

% Define projected CRS (TWD97 TM2, EPSG:3826)
proj = projcrs(3826);

% Convert union projection to lat/lon for plotting
[proj_x, proj_y] = boundary(full_projection);
if isempty(proj_x)
    warning('Projection boundary is empty; using convex hull of all points.');
    [proj_x, proj_y] = boundary(polyshape(all_xy(:,1), all_xy(:,2)));
end
[proj_lat, proj_lon] = projinv(proj, proj_x, proj_y);
plotm(proj_lat, proj_lon, 'b-', 'LineWidth', 2, 'DisplayName', 'Fault Projection');

% Convert strike line to lat/lon
[strike_lat, strike_lon] = projinv(proj, strike_x, strike_y);
plotm(strike_lat, strike_lon, 'r-', 'LineWidth', 2, 'DisplayName', 'Strike Line');

% Convert perpendicular line to lat/lon
[perp_lat, perp_lon] = projinv(proj, perp_x, perp_y);
plotm(perp_lat, perp_lon, 'g-', 'LineWidth', 2, 'DisplayName', 'Perpendicular Line');

% Set map properties
gridm('on');
setm(ax, 'MLabelParallel', 'south');
setm(ax, 'PLabelLocation', 1);
title(ax, 'Fault Surface Projection with Strike and Perpendicular Lines');
legend(ax, 'show', 'Location', 'best');

% Save figure
saveas(gcf, fullfile(pwd, 'Figures_', 'Fault_Projection_Lines.png'));

% Hold off
hold(ax, 'off');


%% 