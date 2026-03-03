clc; clear; close all;

addpath('D:\My_Study\NCU_Study_pvb\ASK14_HW_model\Shanchiao_Hangingwall_Evaluation\gmm_lib');

% This script is used to compute the Site Amplification Factor presented in
% Table 1.

% Period of consideration.
T = 1;

% NPP1: 
M = 7.06; 
npp1.Rrup = 7; npp1.Rjb = 7; npp1.Rx = -10.4; npp1.VS30 = 447.4; npp1.Z1 = 138.7;
% NPP2: 
npp2.Rrup = 3.9; npp2.Rjb = 0; npp2.Rx = 1.5; npp2.VS30 = 707.9; npp2.Z1 = 22.7;
% NPP4: 
npp4.Rrup = 26.5; npp4.Rjb = 10; npp4.Rx = 32.2; npp4.VS30 = 1364.1; npp4.Z1 = 2.3;


% Cal prediction: 
Dip = 49; Ztor = 0;
Ry0 = 0; fas = 0; HW = 0;  region = 6; % input parameters for ASK14
lambda = 90; Fhw = 0; rgn = 1; flg_AS = 0; % input parameters for Ph20
FRO = 1; FSS = 0; FNO = 0; Finter = 0; Fintra = 0; Fas = 0; Fma = 0; FVS30 = 0; % input parameters for Ch20

VS30r = 760;
Z_1r1 = exp(-3.73/2*log((VS30r.^2+290.53^2)/(1750^2+290.53^2))); % Phung et al.
Z_1r2 = exp(-4.08/2*log( (VS30r.^2+355.4^2)/(1750^2+355.4^2) )); % Chao et al.
Z_1r3 = exp(-7.67/4 * log((VS30r^4 + 610^4)/(1360^4 + 610^4))); % ASK14
% Z_1r = min([Z_1r1, Z_1r2, Z_1r3]); % take the min value of Z_1r.
Z_1r = 999;

VS760 = 760; 
Z_1er1 = exp(-3.73/2*log((VS760.^2+290.53^2)/(1750^2+290.53^2))); % Phung et al.
Z_1er2 = exp(-4.08/2*log( (VS760.^2+355.4^2)/(1750^2+355.4^2) )); % Chao et al.
Z_1er3 = exp(-7.67/4 * log((VS760^4 + 610^4)/(1360^4 + 610^4))); % ASK14
Z_1er = min([Z_1er1, Z_1er2, Z_1er3]);

s = 1; Af760 = cell(1,3); Af = cell(1,3);

for npp = [npp1, npp2, npp4]

    SaASK14VS1000 = ASK_2014_nga(M, T, npp.Rrup, npp.Rjb, npp.Rx, Ry0, Ztor, Dip, lambda, fas, HW, 999, 999, VS30r, FVS30, region);
    SaCh20VS1000 = Chao19H(T, M, npp.Rrup, Ztor, VS30r, 999, FRO, FSS, FNO, Finter, Fintra, Fas, Fma, FVS30);
    SaP20VS1000 = Phung_2019h_NGAw2_TW(M, T, npp.Rrup, npp.Rjb, npp.Rx, Ztor, Dip, Fhw, lambda, VS30r, 999, rgn, flg_AS);

    SaASK14VS760 = ASK_2014_nga(M, T, npp.Rrup, npp.Rjb, npp.Rx, Ry0, Ztor, Dip, lambda, fas, HW, 999, 999, VS760, FVS30, region);
    SaCh20VS760 = Chao19H(T, M, npp.Rrup, Ztor, VS760, 999, FRO, FSS, FNO, Finter, Fintra, Fas, Fma, FVS30);
    SaP20VS760 = Phung_2019h_NGAw2_TW(M, T, npp.Rrup, npp.Rjb, npp.Rx, Ztor, Dip, Fhw, lambda, VS760, 999, rgn, flg_AS);
    

    SaASK14VS760_v = ASK_2014_nga(M, T, npp.Rrup, npp.Rjb, npp.Rx, Ry0, Ztor, Dip, lambda, fas, HW, 999, -999, VS760, FVS30, region);
    SaCh20VS760_v = Chao19H(T, M, npp.Rrup, Ztor, VS760, -999, FRO, FSS, FNO, Finter, Fintra, Fas, Fma, FVS30);
    SaP20VS760_v = Phung_2019h_NGAw2_TW(M, T, npp.Rrup, npp.Rjb, npp.Rx, Ztor, Dip, Fhw, lambda, VS760, -999, rgn, flg_AS);
     
    
    SaASK14 = ASK_2014_nga(M, T, npp.Rrup, npp.Rjb, npp.Rx, Ry0, Ztor, Dip, lambda, fas, HW, 999, npp.Z1, npp.VS30, FVS30, region);
    SaCh20 = Chao19H(T, M, npp.Rrup, Ztor, npp.VS30, npp.Z1, FRO, FSS, FNO, Finter, Fintra, Fas, Fma, FVS30);
    SaP20 = Phung_2019h_NGAw2_TW(M, T, npp.Rrup, npp.Rjb, npp.Rx, Ztor, Dip, Fhw, lambda, npp.VS30, npp.Z1, rgn, flg_AS);
    
    SaASK14_v = ASK_2014_nga(M, T, npp.Rrup, npp.Rjb, npp.Rx, Ry0, Ztor, Dip, lambda, fas, HW, 999, -999, npp.VS30, FVS30, region);
    SaCh20_v = Chao19H(T, M, npp.Rrup, Ztor, npp.VS30, -999, FRO, FSS, FNO, Finter, Fintra, Fas, Fma, FVS30);
    SaP20_v = Phung_2019h_NGAw2_TW(M, T, npp.Rrup, npp.Rjb, npp.Rx, Ztor, Dip, Fhw, lambda, npp.VS30, -999, rgn, flg_AS);


%% Calculate the site amp
    % total VS760 siteampl to 1000.
    Sf760ASK14 = SaASK14VS760./SaASK14VS1000;
    Sf760Ch20 = SaCh20VS760./SaCh20VS1000;
    Sf760Ph20 = SaP20VS760./SaP20VS1000;
    % VS30760 amp rel to 1000.
    Sf760ASK14_v = SaASK14VS760_v./SaASK14VS1000; 
    Sf760Ch20_v = SaCh20VS760_v./SaCh20VS1000;
    Sf760Ph20_v = SaP20VS760_v./SaP20VS1000;

    % total site amplification
    SfSpeASK14 = SaASK14./SaASK14VS1000;
    SfSpeCh20 = SaCh20./SaCh20VS1000;
    SfSpePh20 = SaP20./SaP20VS1000;
    % actual VS30 amp rel to 1000.
    SfSpeASK14_v = SaASK14_v./SaASK14VS1000;
    SfSpeCh20_v = SaCh20_v./SaCh20VS1000;
    SfSpePh20_v = SaP20_v./SaP20VS1000;
    
    % Ratio for engineering bedrock. 
    Sf760 = [Sf760ASK14, Sf760Ch20, Sf760Ph20]';
    Sf760_v = [Sf760ASK14_v, Sf760Ch20_v, Sf760Ph20_v]';
    
    Sf760_r = 999*ones(6,1);
    Sf760_r([1, 3, 5],:) = Sf760;
    Sf760_r([2, 4, 6],:) = Sf760_v;
    
    % Ratio for actual Site.
    Sf = [SfSpeASK14, SfSpeCh20, SfSpePh20]';
    Sf_v = [SfSpeASK14_v, SfSpeCh20_v, SfSpePh20_v]';
    
    Sf_r = 999*ones(6,1);
    Sf_r([1, 3, 5],:) = Sf;
    Sf_r([2, 4, 6],:) = Sf_v;

    Af760{s} = Sf760_r;
    Af{s} = Sf_r;
    s = s + 1;
end

Af760mat = cell2mat(Af760);
Afmat = cell2mat(Af);
Amp = [Af760mat Afmat];