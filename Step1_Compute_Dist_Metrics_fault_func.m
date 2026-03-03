function [Rrup, Rjb, Rx, W, Ztor, dip] = Step1_Compute_Dist_Metrics_fault_func(fault_filename, site_lat, site_lon)
% PLOT_FAULT_GEOMETRY Plot 3D geometry of the fault mesh with an example site and distance metrics
%
%   plot_fault_geometry(fault_filename, site_lat, site_lon)
%
% Inputs:
%   fault_filename: Path to the ShanChiaoFault_NCDR.xlsx file
%   site_lat: Site latitude (degrees)
%   site_lon: Site longitude (degrees)
%
% Outputs:
%   Generates a 3D plot with fault mesh, site, and visualized distances (Rrup, Rjb, Rx).
%
% Notes:
%   - Requires Mapping Toolbox for coordinate projection.
%   - Fault in TWD97 TM2 (EPSG:3826) projected coordinates.
%   - Triangulates quad faces.
%   - Assumes hanging wall east (higher X).
%   - Z negative for depth.
%   - Distances in km, but plot in meters.
%   - Visualizes:
%     - Rrup: Red line to closest point on mesh.
%     - Rjb: Green line to closest point on surface projection (horizontal).
%     - Rx: Blue line perpendicular to fault trace.
%
% Example:
%   fault_file = 'ShanChiaoFault_NCDR.xlsx';
%   site_lat = 25.0338;  % Taipei 101
%   site_lon = 121.5645;
%   plot_fault_geometry(fault_file, site_lat, site_lon);

% Define projected CRS
proj = projcrs(3826);  % EPSG:3826 TWD97 TM2

% Convert site lat/lon to projected X,Y
[site_x, site_y] = projfwd(proj, site_lat, site_lon);
site_z = 0;  % Surface site
site_pos = [site_x, site_y, site_z];

% Read fault data
T = readtable(fault_filename, 'Sheet', 'SC');
T = T(strcmp(T.Entity, '3dFace'), :);  % Filter 3dFace

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
all_points = unique(all_points, 'rows');

% COMPUTE FAULT GEOMETRY PARAMETERS (dip, W, Ztor)
% 1. Fit plane to estimate average dip
mean_pt = mean(all_points, 1);
centered = all_points - mean_pt;
[~, ~, V] = svd(centered);
n = V(:, end);  % Normal vector (smallest singular value direction)
n = n / norm(n);

% Dip angle: angle between plane and horizontal = acos(|n_z|)
dip_rad = acos(abs(n(3)));
dip = rad2deg(dip_rad);  % degrees

% 2. Compute Ztor: depth to top of rupture (minimum |Z|)
Ztor = -min(all_points(:, 3)) / 1000;  % km (Z negative for depth)

% 3. Compute W: down-dip width
max_depth = -min(all_points(:, 3));  % meters
W = max_depth / 1000 / sin(dip_rad);  % km


%% Compute and visualize distances

% 1. Rrup: Closest to mesh
[Rrup, closest_rup_pt] = compute_rrup(site_pos, faces);

% 2. Rjb: Closest to surface projection
all_xy = [];
for f = 1:length(faces)
    for t = 1:2
        all_xy = [all_xy; faces{f}{t}(:,1:2)];
    end
end
all_xy = unique(all_xy, 'rows');
hull_idx = convhull(all_xy(:,1), all_xy(:,2));
hull_xy = all_xy(hull_idx, :);

[in, ~] = inpolygon(site_x, site_y, hull_xy(:,1), hull_xy(:,2));
if in
    Rjb = 0;
    closest_jb_pt = [site_x, site_y, 0];  % Directly below
    % disp('Site is above rupture projection; Rjb = 0');
else
    [Rjb, closest_jb_pt] = compute_rjb([site_x site_y], hull_xy);
    closest_jb_pt = [closest_jb_pt 0];  % At surface
end


% 3. Rx: Perp to trace
Rx = calculate_rx_multilinear(surface_points, site_x, site_y, dip);
[p, ~] = polyfit(surface_points(:,2), surface_points(:,1), 1);
slope = p(1); intercept = p(2);
% 
% dist = abs((slope * site_y - site_x + intercept) / sqrt(slope^2 + 1)) / 1000;
predicted_x = slope * site_y + intercept;
% 
% if site_x > predicted_x
%     Rx = dist;
% else
%     Rx = -dist;
% end

% Closest point on trace line
% Infinite line projection
dir_vec = [1, slope]; dir_vec = dir_vec / norm(dir_vec);
pt_on_line = [predicted_x, site_y, 0];
perp_vec = [site_x - predicted_x, site_y - site_y, 0];
proj_len = dot(perp_vec(1:2), dir_vec);
closest_rx_pt = pt_on_line + [proj_len * dir_vec, 0];

% Clamp to trace extent if needed (approximate)
min_y = min(surface_points(:,2)); max_y = max(surface_points(:,2));
if closest_rx_pt(2) < min_y || closest_rx_pt(2) > max_y
    % disp('Rx projection outside trace; clamping to endpoint.');
    if closest_rx_pt(2) < min_y
        closest_rx_pt(2) = min_y;
    else
        closest_rx_pt(2) = max_y;
    end
    closest_rx_pt(1) = slope * closest_rx_pt(2) + intercept;
end
closest_rx_pt(3) = 0;



end

function [Rrup, closest_pt] = compute_rrup(site_pos, faces)
Rrup = inf;
closest_pt = [0 0 0];
for f = 1:length(faces)
    for t = 1:2
        tri = faces{f}{t};
        [dist, pt] = point_to_triangle_distance(site_pos, tri);
        if dist < Rrup
            Rrup = dist;
            closest_pt = pt;
        end
    end
end
Rrup = Rrup / 1000;  % km
end

function [Rjb, closest_pt] = compute_rjb(site_xy, hull_xy)
Rjb = inf;
closest_pt = [0 0];
n = size(hull_xy, 1);
for k = 1:n
    p1 = hull_xy(k, :);
    p2 = hull_xy(mod(k, n)+1, :);
    [dist, pt] = point_to_line_segment_distance(site_xy, p1, p2);
    if dist < Rjb
        Rjb = dist;
        closest_pt = pt;
    end
end
Rjb = Rjb / 1000;  % km
end


function [dist, closest_pt] = point_to_triangle_distance(point, tri)
% Modified to return closest point
a = tri(1,:); b = tri(2,:); c = tri(3,:);
ab = b - a; ac = c - a; ap = point - a;
n = cross(ab, ac); n = n / norm(n);
plane_dist = abs(dot(ap, n));
proj = point - dot(ap, n) * n;
if is_point_in_triangle(proj, tri)
    dist = plane_dist;
    closest_pt = proj;
else
    [d1, p1] = point_to_line_segment_distance(point, a, b);
    [d2, p2] = point_to_line_segment_distance(point, a, c);
    [d3, p3] = point_to_line_segment_distance(point, b, c);
    [dist, idx] = min([d1 d2 d3]);
    pts = {p1, p2, p3};
    closest_pt = pts{idx};
end
end

% Improved Rx Calculation: Perpendicular Distance to Multi-Linear Fault Trace
% This replaces the straight-line polyfit with a polyline distance computation
% Assumes surface_points are the surface trace points (X, Y in projected coordinates)
% Outputs Rx with sign based on fault dip direction (assuming westward dip for Taiwan faults)

function Rx = calculate_rx_multilinear(surface_points, site_x, site_y, dip)
% CALCULATE_RX_MULTILINEAR Signed perpendicular distance to multi-linear fault trace (km)
%
% Inputs:
%   surface_points: Nx2 matrix [X Y] of surface trace points (meters)
%   site_x, site_y: Site coordinates (meters)
%   dip: Fault dip angle (degrees, for sign convention)
%
% Output:
%   Rx: Signed distance (km, positive on hanging wall)

% Step 1: Sort surface_points to form ordered polyline (along strike)
% Approximate strike direction from principal component (PCA)
centered = surface_points - mean(surface_points, 1);
[~, ~, V] = svd(centered);
strike_dir = V(:,1);  % First principal component (along strike, 2x1)

% Project points onto strike direction to sort
projs = sum(surface_points * strike_dir, 2);
[~, sort_idx] = sort(projs);
polyline = surface_points(sort_idx, :);  % Ordered polyline

% Step 2: Find closest segment to site
min_dist = inf;
closest_segment_idx = 1;
for seg = 1:size(polyline, 1)-1
    a = polyline(seg, :);
    b = polyline(seg+1, :);
    dist = point_to_line_segment_distance([site_x site_y], a, b);
    if dist < min_dist
        min_dist = dist;
        closest_segment_idx = seg;
    end
end
min_dist = min_dist / 1000;  % km

% Step 3: Compute perpendicular distance to closest segment
a = polyline(closest_segment_idx, :);
b = polyline(closest_segment_idx + 1, :);
vec = b - a;
proj = dot([site_x site_y] - a, vec) / dot(vec, vec);
proj = max(0, min(1, proj));  % Clamp to segment
closest_pt = a + proj * vec;
perp_dist = norm([site_x site_y] - closest_pt) / 1000;  % km

    % Step 4: Determine sign (hanging wall vs footwall)
    % FIXED: Use 2D cross product for side determination (avoids 3D concatenation)
    diff_vec = [site_x site_y] - closest_pt;
    side = strike_dir(1) * diff_vec(2) - strike_dir(2) * diff_vec(1);  % Z-component equivalent
    side = sign(side);
    
    % For westward dip (Taiwan faults), positive side (left of strike) is hanging wall
    % Adjust based on dip: if dip > 0, hanging wall is on the side where fault dips toward
    if dip < 0  % Typical reverse fault dip
        Rx = side * perp_dist;  % Positive for hanging wall (adjust sign if eastward dip)
    else
        Rx = -side * perp_dist;  % For normal faults
    end

end

function [dist, closest_pt] = point_to_line_segment_distance(p, a, b)
% Modified to return closest point (works for 2D/3D)
if size(a,2) == 2
    p = p(1:2); a = a(1:2); b = b(1:2);  % 2D
end
vec = b - a;
proj = dot(p - a, vec) / dot(vec, vec);
if proj < 0
    dist = norm(p - a);
    closest_pt = a;
elseif proj > 1
    dist = norm(p - b);
    closest_pt = b;
else
    closest_pt = a + proj * vec;
    dist = norm(p - closest_pt);
end
end

function in = is_point_in_triangle(p, tri)
    a = tri(1,:); b = tri(2,:); c = tri(3,:);
    v0 = c - a; v1 = b - a; v2 = p - a;
    dot00 = dot(v0, v0); dot01 = dot(v0, v1); dot02 = dot(v0, v2);
    dot11 = dot(v1, v1); dot12 = dot(v1, v2);
    inv_den = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * inv_den;
    v = (dot00 * dot12 - dot01 * dot02) * inv_den;
    in = (u >= 0) && (v >= 0) && (u + v <= 1);
end


% Improved Rx Calculation: Perpendicular Distance to Multi-Linear Fault Trace
% This replaces the straight-line polyfit with a polyline distance computation
% Assumes surface_points are the surface trace points (X, Y in projected coordinates)
% Outputs Rx with sign based on fault dip direction (assuming westward dip for Taiwan faults)

