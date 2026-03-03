% Matlab Script to calculate Moment Magnitude (M) from Rupture Area (A), Rupture Length (L), and Rupture Width (W)
% Based on Huang et al. (2024) empirical source-scaling laws
% Requires inputs: A or L or W (in km or km^2), fault dip nu (degrees), Ztor (km), ZBSZ (km), region ('CA', 'TW', 'JP', 'Other')
% Uses iteration to handle dependency on NZbor
% Example usage:
% M_A = calculate_M_from_area(1000, 60, 0.5, 25, 'CA');
% M_L = calculate_M_from_length(50, 60, 0.5, 25, 'CA');
% M_W = calculate_M_from_width(20, 60, 0.5, 25, 'CA');
% Include the original functions (copy-pasted as provided)

function [M_A, M_L, M_W] = Huang_and_Abrahamson_2024_calculate_M(A, L, W, nu, Ztor, ZBSZ, region)
    if isempty(region)
        region = 'CA';
    end
    M_A = calculate_M_from_area(A, nu, Ztor, ZBSZ, region);
    M_L = calculate_M_from_length(L, nu, Ztor, ZBSZ, region);
    M_W = calculate_M_from_width(W, nu, Ztor, ZBSZ, region);
end

function M = calculate_M_from_area(A, nu, Ztor, ZBSZ, region)
    logA = log10(A);
    M = logA + 4;
    tol = 1e-3;
    max_iter = 50;
    for iter = 1:max_iter
        NZbor_guess = 0.5;
        asp_log = compute_log_asp(M, nu, NZbor_guess, region);
        asp = 10^asp_log;
        W_est = sqrt(A / asp);
        Zbor = Ztor + W_est * sind(nu);
        NZbor = min(1, max(0, Zbor / ZBSZ));
        logA_model = compute_log_area(M, nu, NZbor, Ztor < 1.0);
        delta = logA - logA_model;
        M = M + delta;
        if abs(delta) < tol
            break;
        end
    end
end

function M = calculate_M_from_length(L, nu, Ztor, ZBSZ, region)
    logL = log10(L);
    M = 2 * (logL - 2);
    tol = 1e-3;
    max_iter = 50;
    for iter = 1:max_iter
        NZbor_guess = 0.5;
        asp_log = compute_log_asp(M, nu, NZbor_guess, region);
        asp = 10^asp_log;
        A_est = L^2 / asp;
        W_est = L / asp;
        Zbor = Ztor + W_est * sind(nu);
        NZbor = min(1, max(0, Zbor / ZBSZ));
        logA_model = compute_log_area(M, nu, NZbor, Ztor < 1);
        logL_model = 0.5 * (logA_model + asp_log);
        delta = logL - logL_model;
        M = M + 2 * delta;
        if abs(delta) < tol
            break;
        end
    end
end

function M = calculate_M_from_width(W, nu, Ztor, ZBSZ, region)
    logW = log10(W);
    M = 2 * (logW - 2);
    tol = 1e-3;
    max_iter = 50;
    for iter = 1:max_iter
        Zbor = Ztor + W * sind(nu);
        NZbor = min(1, max(0, Zbor / ZBSZ));
        asp_log = compute_log_asp(M, nu, NZbor, region);
        logA_model = compute_log_area(M, nu, NZbor, Ztor < 1);
        logW_model = 0.5 * (logA_model - asp_log);
        delta = logW - logW_model;
        M = M + 2 * delta;
        if abs(delta) < tol
            break;
        end
    end
end

function log_asp = compute_log_asp(M, nu, NZbor, region)
    c0 = log10(1.3);
    c1 = 0.532;
    c2 = 7.17;
    c3 = -0.0105;
    c4 = 0.126;
    c5 = -0.196;
    switch region
        case 'CA'
            c6 = 0.21;
        otherwise
            c6 = 0;
    end
    MBP = c2 + c3 * nu;
    if M <= MBP
        S2 = c0;
    else
        S2 = c0 + c1 * (M - MBP);
    end
    S3 = S2 + c4 + c5 * NZbor;
    if M <= MBP - 0.5
        RT = 0;
    elseif M < MBP + 0.5
        RT = c6 * (M - MBP + 0.5);
    else
        RT = c6;
    end
    log_asp = S3 + RT;
end

function log_area = compute_log_area(M, nu, NZbor, is_surface)
    d0 = -0.215;
    d1 = -0.068;
    d2 = 1.371;
    if is_surface
        d3 = -0.248;
        d4 = 0.408;
    else
        d3 = -0.269;
        d4 = 0.296;
    end
    c2 = 7.17;
    c3 = -0.0105;
    MBP = c2 + c3 * nu;
    MBPL = MBP + d2;
    S1A = M - 4;
    term = d0 * (max(M, MBP) - MBPL);
    adjustment = max(term, d1);
    S2A = S1A + adjustment;
    log_area = S2A + d3 + d4 * NZbor;
end