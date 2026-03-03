%% Ye and Ma (2011)
% Matlab Script to calculate Moment Magnitude (M) from Rupture Area (A), Rupture Length (L), and Rupture Width (W)
% Based on Yen and Ma (2011) empirical source-scaling laws for Taiwan earthquakes
% Uses piecewise scaling with break at M0 = 10^20 N m
% Units: A in km^2, L and W in km
% Mw calculated using Mw = (2/3) log10(M0) - 6.07, with M0 in N m
% Example usage:
% M_A = calculate_M_from_area(1000);
% M_L = calculate_M_from_length(50);
% M_W = calculate_M_from_width(20);
function [M_A, M_L, M_W]=Ye_and_Ma_2011_calculate_M(A, L, W)
    M_A = calculate_M_from_area(A);
    M_L = calculate_M_from_length(L);
    M_W = calculate_M_from_width(W);
end
function Mw = calculate_M_from_area(A)
    logA = log10(A);
    % Assume small regime
    logM0 = logA + 16.16;
    if 10^logM0 <= 10^20
        % Small
    else
        % Large regime
        logM0 = 1.5 * logA + 15.165;
    end
    Mw = (2/3) * logM0 - 6.07;
end

function Mw = calculate_M_from_length(L)
    logL = log10(L);
    % Assume small
    logM0 = 2 * logL + 16.16;
    if 10^logM0 <= 10^20
        % Small
    else
        % Large
        logM0 = 3 * logL + 14.52;
    end
    Mw = (2/3) * logM0 - 6.07;
end

function Mw = calculate_M_from_width(W)
    logW = log10(W);
    % Assume small
    logM0 = 2 * logW + 16.16;
    if 10^logM0 <= 10^20
        % Small
    else
        % Large
        logM0 = 3 * logW + 15.81;
    end
    Mw = (2/3) * logM0 - 6.07;
end