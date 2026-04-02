clc; clear; close all;
% Load Bea20 and GC2
% Load the ShanChiaoFault geometry
%% ==========================================================
% 1. READ FAULT GEOMETRY (ShanChiao Fault)
% ===========================================================
filename = 'ShanChiaoFault_NCDR.xlsx';
sheet = 'SC';
dT = readtable(filename,'Sheet',sheet);
%% ==========================================================
% 2. EXTRACT SURFACE TRACE (Z ≈ 0)
% ===========================================================

tol = 1e-3;
surf_nodes = dT(abs(dT.Z) < tol, :);

% Unique surface coordinates (lon,lat)
trace_ll = unique([surf_nodes.Longitude surf_nodes.Latitude],'rows');

% Determine the xy coordinate of the top edges
xy_top_edge = [dT.Latitude(dT.Z==0), dT.Longitude(dT.Z==0)]; % this is equivalent to the trace.
xy_top_edge = unique(xy_top_edge,'rows');

% Determine the xy coordinate of the bottom edges
[Z_in_sort, izd] = sort(abs(unique(dT.Z)), 'descend');
Z_bottom = Z_in_sort(1:7); % define 7 bottom nodes corresponding to the 7 largest depth.

xy_bot_edge = [];
for j = 1:length(Z_bottom)
    idz = abs(dT.Z) == Z_bottom(j);
    xy_bot_edge = [xy_bot_edge; [dT.Latitude(idz), dT.Longitude(idz)]];
end
xy_bot_edge = unique(xy_bot_edge, 'rows');

strike_start = xy_top_edge(1,:);
strike_end   = xy_top_edge(end,:);
%% ==========================================================
% 3. CONVERT TO LOCAL CARTESIAN (km)
% ===========================================================

lat0 = mean(trace_ll(:,2));
lon0 = mean(trace_ll(:,1));
R = 6371;  % km

x = R*cosd(lat0).*deg2rad(trace_ll(:,1)-lon0);
y = R*deg2rad(trace_ll(:,2)-lat0);

trace_xy = [x y];
%% ==========================================================
% 4. DENSIFY TRACE (for accurate projection)
% ===========================================================

% Arc-length parameterization
ds = 1;  % 1000 m spacing

dxy = diff(trace_xy);
seglen = sqrt(sum(dxy.^2,2));
s = [0; cumsum(seglen)];
Ltotal = s(end);

s_dense = (0:ds:Ltotal)';
x_dense = interp1(s, trace_xy(:,1), s_dense);
y_dense = interp1(s, trace_xy(:,2), s_dense);

trace_xy = [x_dense y_dense];
%% ==========================================================
% 5. DEFINE HYPOCENTER and its Depth
% ===========================================================

hypo_lat = 25.15;
hypo_lon = 121.59;
hypo_depth = 10;   % km

hypo_x = R*cosd(lat0)*deg2rad(hypo_lon-lon0);
hypo_y = R*deg2rad(hypo_lat-lat0);

%% ==========================================================
% 6. COMPUTE MEAN STRIKE AND DIP PROJECTION
% ===========================================================

dxy = diff(trace_xy);
strike_all = atan2d(dxy(:,1), dxy(:,2));

strike_all(strike_all<0)=strike_all(strike_all<0)+360;
strike_mean = mean(strike_all);

dip_angle = 50; % Specify for ShanChiao Fault;
horizontal_shift = hypo_depth / tand(dip_angle);
dip_dir = mod(strike_mean + 90,360);

dx = horizontal_shift*sind(dip_dir);
dy = horizontal_shift*cosd(dip_dir);

hypo_surface = [hypo_x + dx , hypo_y + dy];
%% ==========================================================
% 7. EXACT PROJECTION ONTO FAULT TRACE
% ===========================================================
minDist = inf;

for i = 1:size(trace_xy,1)-1
    A = trace_xy(i,:);
    B = trace_xy(i+1,:);
    AB = B-A;
    AP = hypo_surface-A;

    t = dot(AP, AB)/dot(AB, AB);
    t = max(0, min(1, t));

    proj = A + t*AB;
    d = norm(proj-hypo_surface);

    if d < minDist
        minDist = d;
        origin = proj;
    end
end

%% Find projection of hypocenter surface point

% dist = vecnorm(trace_xy - hypo_surface, 2, 2);
% [~, idx] = min(dist);
% 
% origin = trace_xy(idx,:);

%% ==========================================================
% 8. SHIFT COORDINATE SYSTEM TO GC2 ORIGIN
% ===========================================================
trace_xy = trace_xy - origin;


%% ==========================================================
% 9. SEGMENTATION BASED ON STRIKE CHANGE
% ===========================================================

dxy = diff(trace_xy);
strike = atan2d(dxy(:,1), dxy(:,2));
strike(strike<0)=strike(strike<0)+360;

dstrike = abs(diff(strike));
break_idx = find(dstrike>5);

segment_nodes = unique([1; break_idx+1; size(trace_xy,1)]);

segments = {};
for i = 1:length(segment_nodes)-1
    segments{end+1} = trace_xy(segment_nodes(i):segment_nodes(i+1),:);
end

%% ==========================================================
% 10. COMPUTE SEGMENT LENGTH AND STRIKE
% ===========================================================

nseg = length(segments);
L = zeros(1,nseg);
Strike_seg = zeros(1,nseg);

for i = 1:nseg
    seg = segments{i};
    L(i) = sum(sqrt(sum(diff(seg).^2,2)));

    dx = seg(end,1)-seg(1,1);
    dy = seg(end,2)-seg(1,2);
    Strike_seg(i) = atan2d(dx,dy);
    if Strike_seg(i)<0
        Strike_seg(i)=Strike_seg(i)+360;
    end
end


%% ==========================================================
% 11. BUILD CORRECT GC2 STRUCTURE (MULTI-SEGMENT)
% ===========================================================

dx = diff(trace_xy(:,1));
dy = diff(trace_xy(:,2));

l = sqrt(dx.^2 + dy.^2);
strike = atan2d(dx,dy);
strike(strike<0) = strike(strike<0)+360;

ftraces(1).trace  = trace_xy;
ftraces(1).l      = l';
ftraces(1).strike = strike';



%% Define the Local Cartesian Coordinate
SiteX = -60:1:100;  % strike-parallel
SiteY = -80:1:80;   % strike-normal

nt = length(ftraces);
M = 7.2;
Version = 1;
Tdo = 3;

Rake = -90;
Ztor = 0;

type.epi = origin; %[0 , 0];%origin;
type.po  = origin; %[0 , 0];%origin;

%% Call the Spudich and Chiou (2015) GC2 function
type.str = 'JB'; 
discordant = false;
gridflag = true;
[T, U, W, reference_axis, p_origin, nominal_strike, Upo] = GC2(ftraces,SiteX,SiteY,type,discordant,gridflag);

% Calculate the maximum value of S in each direction for this hypocenter; it is U calculated at the nominal strike ends
[~, Uend,~,~,~,~,~,~] = GC2(ftraces,nominal_strike.a(1,1), nominal_strike.a(1,2),type,discordant,gridflag);
[~, Uend2,~,~,~,~,~,~] = GC2(ftraces,nominal_strike.a(2,1), nominal_strike.a(2,2),type,discordant,gridflag);

Smax1 = min(Uend, Uend2);
Smax2 = max(Uend, Uend2); 


dx = diff(trace_xy(:,1));
dy = diff(trace_xy(:,2));
seglen = hypot(dx,dy);
Ucum = [0; cumsum(seglen)];

% Find node closest to origin (hypocenter surface projection)
[~,idx_origin] = min(vecnorm(trace_xy,2,2));
U_origin = Ucum(idx_origin);

% Total rupture length
L_total = Ucum(end);


% ----- CHOOSE RUPTURE DIRECTION -----
% Example: rupture propagates toward positive U direction

% Smax1 = 0;
% Smax2 = L_total - U_origin;

% If rupture propagates toward negative direction instead:
% Smax1 = -U_origin;
% Smax2 = 0;
%% Call the directivity model
fDi = zeros(size(U));
for ii=1:size(U,2)
    [fD,fDi(:,ii),PhiRed,PhiRedi,PredicFuncs,Other]=Bea24(M,U(:,ii),T(:,ii),Smax1,Smax2,Ztor,Rake,Tdo,Version);
    
    S2(:,ii)=Other.S2;
    fs2(:,ii)=PredicFuncs.fs2;
    ftheta(:,ii)=PredicFuncs.ftheta;
    fdist(:,ii)=PredicFuncs.fdist;
    fGprime(:,ii)=PredicFuncs.fGprime;   
end  



%% Plot U T contours
figure;  set(gcf,'position',[311   188    747 391 ]); 
subplot(1,2,1)
    Z=[fliplr(0:-5:round(min(min(T)))) 5:5:round(max(max(T)))]; % contour interval
    V=[fliplr(0:-20:round(min(min(T)))) 20:20:round(max(max(T)))]; % label interval
    [c,h]=contour(SiteX,SiteY,T,Z); hold on
    clabel(c,h,V)
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    axis square
    title('GC2, T Coordinate')
    xlabel('Easting (km)')
    ylabel('Northing (km)')

subplot(1,2,2)
    Z=[fliplr(0:-5:round(min(min(U)))) 5:5:round(max(max(U)))]; % contour interval
    V=[fliplr(0:-20:round(min(min(U)))) 20:20:round(max(max(U)))]; % label interval
    [c,h]=contour(SiteX,SiteY,U,Z); hold on
    clabel(c,h,V)
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    axis square
    title('GC2, U Coordinate')
    xlabel('Easting (km)')
    ylabel('Northing (km)')

    
 %%  plot the directivity model
 
figure; set(gcf,'position',[54 16 996 758]);
subplot(2,3,1)
    contourf(SiteX,SiteY,S2,'linestyle','none'); hold on
    colorbar; 
    title('\itS2')
    colormap(othercolor('BuDRd_18', 256));
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    %axis equal
    ylabel('Northing (km)')
    axis([-80 80 -80 170])
subplot(2,3,2)
    contourf(SiteX,SiteY,fs2,'linestyle','none'); hold on
    colorbar; 
    title('\itf_{S2}')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    %axis equal
    axis([-80 80 -80 170])
subplot(2,3,3)
    contourf(SiteX,SiteY,ftheta,'linestyle','none'); hold on
    colorbar;
    title('\itf_{\theta}')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    %axis equal
    axis([-80 80 -80 170])
subplot(2,3,4)
    contourf(SiteX,SiteY,fGprime,'linestyle','none'); hold on
    colorbar; 
    title('\itf_G''')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    xlabel('Easting (km)')
    ylabel('Northing (km)')
    %axis equal
    axis([-80 80 -80 170])
subplot(2,3,5)
    contourf(SiteX,SiteY,fdist,'linestyle','none'); hold on
    colorbar; 
    title('\itf_{dist}')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    xlabel('Easting (km)')
    %axis equal
    axis([-80 80 -80 170])
subplot(2,3,6)
    contourf(SiteX,SiteY,exp(fDi),'linestyle','none'); hold on
    colorbar; 
    title('Amplification, T=3 sec')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    xlabel('Easting (km)')
    %axis equal
    axis([-80 80 -80 170])
    % clim([.6 1.8])

    
 disp('Finished with example script');

 [XX, YY] = meshgrid(SiteX, SiteY);

 X_global = XX + origin(1);
 Y_global = YY + origin(2);

 lat_site = lat0 + rad2deg(Y_global / R);
 lon_site = lon0 + rad2deg(X_global ./ (R * cosd(lat0)));

 lat_origin_test = lat0 + rad2deg(origin(2)/R);
 lon_origin_test = lon0 + rad2deg(origin(1)/(R*cosd(lat0)));



%% =====================================================
% Convert Fault Trace to Geographic Coordinates
% =====================================================

fault_lat = cell(nt,1);
fault_lon = cell(nt,1);

for ii = 1:nt
    
    X_fault = ftraces(ii).trace(:,1) + origin(1);
    Y_fault = ftraces(ii).trace(:,2) + origin(2);
    
    fault_lat{ii} = lat0 + rad2deg(Y_fault / R);
    fault_lon{ii} = lon0 + rad2deg(X_fault ./ (R*cosd(lat0)));
    
end

% Convert epicenter
epi_lat = lat0 + rad2deg(type.epi(2)/R);
epi_lon = lon0 + rad2deg(type.epi(1)/(R*cosd(lat0)));


figure; set(gcf,'position',[54 16 1200 800]);
subplot(2,3,1)
contourf(lon_site, lat_site, S2,'linestyle','none'); hold on
colorbar
title('\itS2')
colormap(othercolor('BuDRd_18',256))
for ii=1:nt
    plot(fault_lon{ii}, fault_lat{ii},'k','linewidth',2)
end
plot(epi_lon, epi_lat,'kp','markerfacecolor','r','markersize',12)

axis equal
xlabel('Longitude')
ylabel('Latitude')

subplot(2,3,2)
    contourf(lon_site, lat_site, fs2, 'linestyle','none'); hold on
    colorbar; 
    title('\itf_{S2}')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(fault_lon{ii}, fault_lat{ii},'k','linewidth',2)
    end
    plot(epi_lon, epi_lat,'kp','markerfacecolor','r','markersize',12)
    axis equal
    xlabel('Longitude')
    ylabel('Latitude')
    
subplot(2,3,3)
    contourf(lon_site, lat_site, ftheta,'linestyle','none'); hold on
    colorbar;
    title('\itf_{\theta}')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(fault_lon{ii}, fault_lat{ii},'k','linewidth',2)
    end
    plot(epi_lon, epi_lat,'kp','markerfacecolor','r','markersize',12)
    axis equal
    xlabel('Longitude')
    ylabel('Latitude')

subplot(2,3,4)
    contourf(lon_site, lat_site, fGprime,'linestyle','none'); hold on
    colorbar; 
    title('\itf_G''')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(fault_lon{ii}, fault_lat{ii},'k','linewidth',2)
    end
    plot(epi_lon, epi_lat,'kp','markerfacecolor','r','markersize',12)
    xlabel('Longitude')
    ylabel('Latitude')
    axis equal

subplot(2,3,5)
    contourf(lon_site, lat_site, fdist,'linestyle','none'); hold on
    colorbar; 
    title('\itf_{dist}')
     colormap(othercolor('BuDRd_18',256))
    
    for ii=1:nt
        plot(fault_lon{ii}, fault_lat{ii},'k','linewidth',2)
    end
    plot(epi_lon, epi_lat,'kp','markerfacecolor','r','markersize',12)
    xlabel('Longitude')
    ylabel('Latitude')
    axis equal


subplot(2,3,6)
contourf(lon_site, lat_site, exp(fDi),'linestyle','none'); hold on
colorbar
title('Amplification, T=3 sec')
colormap(othercolor('BuDRd_18',256))

for ii=1:nt
    plot(fault_lon{ii}, fault_lat{ii},'k','linewidth',2)
end
plot(epi_lon, epi_lat,'kp','markerfacecolor','r','markersize',12)

axis equal
xlabel('Longitude')
ylabel('Latitude')




%% ==========================================================
% 12. VISUAL CHECK
% ===========================================================

% figure; hold on
% plot(trace_xy(:,1),trace_xy(:,2),'k')
% plot(0,0,'ro','LineWidth',2);
% axis equal
% title('GC2 Origin at Hypocenter Surface Projection');
