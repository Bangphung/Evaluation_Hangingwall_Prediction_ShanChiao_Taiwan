function f_hw = swus_hanging_wall_model(T, M, dip, Rx, Rrup, Rjb, W, Ztor, model_id)
% HANGING_WALL_MODEL Computes the hanging-wall effect factor f_HW
%   f_hw = hanging_wall_model(T, M, dip, Rx, Rrup, Ztor, model_id)
%
% Inputs:
%   T: Spectral period (sec), scalar
%   M: Magnitude
%   dip: Fault dip angle (degrees)
%   Rx: Horizontal distance from top of rupture perpendicular to strike (km), positive on hanging wall
%   Rrup: Closest distance to rupture plane (km)
%   Ztor: Depth to top of rupture (km)
%   model_id: Model number (1 to 5) for alternative HW models
%
% Output:
%   f_hw: Hanging wall factor added to ln(ground motion)
%
% Notes:
%   - Assumes distances in km, dip in degrees.
%   - Interpolates coefficients in log-period space.
%   - Clips T to [0.01, 10] sec.
%   - Sets f_hw = 0 if Rx < 0 (footwall).
%   - Based on SWUS GMC hanging wall model common form (Eq. 2-3a, b, c).

    if Rx < 0
        f_hw = 0;
        return;
    end

    % Period grid
    periods = [0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, ...
               0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10];

% C1 coefficients for models 1 to 5 (rows: periods, columns: models)
% C1 = [ ...
%     0.868, 0.982, 1.038, 1.095, 1.209; ...
%     0.867, 0.987, 1.046, 1.106, 1.226; ...
%     0.856, 0.997, 1.067, 1.138, 1.278; ...
%     0.840, 1.027, 1.121, 1.215, 1.402; ...
%     0.857, 1.041, 1.133, 1.226, 1.410; ...
%     0.848, 1.040, 1.135, 1.231, 1.422; ...
%     0.868, 1.009, 1.080, 1.150, 1.292; ...
%     0.850, 1.005, 1.082, 1.160, 1.315; ...
%     0.868, 0.985, 1.044, 1.102, 1.219; ...
%     0.839, 0.974, 1.041, 1.108, 1.242; ...
%     0.780, 0.934, 1.011, 1.089, 1.243; ...
%     0.741, 0.902, 0.982, 1.063, 1.223; ...
%     0.613, 0.869, 0.997, 1.125, 1.380; ...
%     0.621, 0.788, 0.872, 0.955, 1.123; ...
%     0.506, 0.662, 0.740, 0.818, 0.974; ...
%     0.391, 0.537, 0.609, 0.682, 0.828; ...
%     0.128, 0.245, 0.304, 0.362, 0.480; ...
%     0.000, 0.034, 0.088, 0.138, 0.231; ...
%     0.000, 0.000, 0.000, 0.000, 0.040; ...
%     0.000, 0.000, 0.000, 0.000, 0.000; ...
%     0.000, 0.000, 0.000, 0.000, 0.000  ...
% ];

% % C2 coefficients
% C2 = [0.2160, 0.2172, 0.2178, 0.2199, 0.2218, 0.2223, 0.2214, 0.2212, ...
%       0.2198, 0.2019, 0.2090, 0.2063, 0.1713, 0.1571, 0.1559, 0.1559, ...
%       0.1616, 0.1616, 0.1616, 0.1616, 0.1616];
% 
% % C3 coefficients
% C3 = [2.0289, 2.0260, 2.0163, 1.9870, 1.9906, 1.9974, 2.0162, 1.9746, ...
%       1.9931, 2.0179, 2.0249, 2.0041, 1.8687, 1.8526, 1.8336, 1.7996, ...
%       1.6740, 1.6740, 1.6740, 1.6740, 1.6740];
% 
% % C4 coefficients
% C4 = [0.1675, 0.1666, 0.1670, 0.1699, 0.1817, 0.1717, 0.1814, 0.1834, ...
%       0.1767, 0.1658, 0.1624, 0.1749, 0.1866, 0.3143, 0.3195, 0.3246, ...
%       0.3314, 0.3314, 0.3314, 0.3314, 0.3314];

    C1_HW1	= [0.868 	0.867 	0.856 	0.840 	0.857 	0.848 	0.868 	0.850 	0.868 	0.839 	0.780 	0.741 	0.613 	0.621 	0.506 	0.391 	0.128 	0.000 	0.000 	0.000 	0.000]; 
    C1_HW2	= [0.982 	0.987 	0.997 	1.027 	1.041 	1.040 	1.009 	1.005 	0.985 	0.974 	0.934 	0.902 	0.869 	0.788 	0.662 	0.537 	0.245 	0.034 	0.000 	0.000 	0.000];
    C1_HW3	= [1.038 	1.046 	1.067 	1.121 	1.133 	1.135 	1.080 	1.082 	1.044 	1.041 	1.011 	0.982 	0.997 	0.872 	0.740 	0.609 	0.304 	0.088 	0.000 	0.000 	0.000]; 
    C1_HW4	= [1.095 	1.106 	1.138 	1.215 	1.226 	1.231 	1.150 	1.160 	1.102 	1.108 	1.089 	1.063 	1.125 	0.955 	0.818 	0.682 	0.362 	0.138 	0.000 	0.000 	0.000]; 
    C1_HW5	= [1.209 	1.226 	1.278 	1.402 	1.410 	1.422 	1.292 	1.315 	1.219 	1.242 	1.243 	1.223 	1.380 	1.123 	0.974 	0.828 	0.480 	0.231 	0.040 	0.000 	0.000]; 
    C1 = [C1_HW1; C1_HW2; C1_HW3; C1_HW4; C1_HW5]';
    
    C2	= [0.2160 	0.2172 	0.2178 	0.2199 	0.2218 	0.2213 	0.2169 	0.2131 	0.1988 	0.2019 	0.2090 	0.2053 	0.1713 	0.1571 	0.1559 	0.1559 	0.1616 	0.1616 	0.1616 	0.1616 	0.1616]; 
    C3	= [2.0289 	2.0260 	2.0163 	1.9870 	1.9906 	1.9974 	2.0162 	1.9746 	1.9931 	2.0179 	2.0249 	2.0041 	1.8697 	1.8526 	1.8336 	1.7996 	1.6740 	1.6740 	1.6740 	1.6740 	1.6740]; 
    C4	= [0.1675 	0.1666 	0.1670 	0.1699 	0.1817 	0.1717 	0.1814 	0.1834 	0.1767 	0.1658 	0.1624 	0.1719 	0.1866 	0.3143 	0.3195 	0.3246 	0.3314 	0.3314 	0.3314 	0.3314 	0.3314]; 
    
    % Clip T to range
    T = max(0.01, min(10, T));
    
    % Interpolate in log-period (linear interpolation in log space)
    logT = log10(T);
    log_periods = log10(periods);
    
    c1 = interp1(log_periods, C1(:, model_id), logT, 'linear', 'extrap');
    c2 = interp1(log_periods, C2, logT, 'linear', 'extrap');
    c3 = interp1(log_periods, C3, logT, 'linear', 'extrap');
    c4 = interp1(log_periods, C4, logT, 'linear', 'extrap');
    
    % Compute terms
    cos_dip = cos(dip * pi / 180);
    rx_term = c2 + (1 - c2) * tanh(c3*Rx / (W*cos_dip));
    mag_term = 1 + c4 * (M - 7);
    r_taper = 1 - Rjb / (Rrup + 0.1);
    z_taper = max(0, 1 - Ztor / 12);
    
    % Hanging wall factor
    f_hw = c1 * cos_dip * rx_term * mag_term * r_taper * z_taper;

end