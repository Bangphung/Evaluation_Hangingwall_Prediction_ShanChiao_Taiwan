% coded by Yue Hua
%               Stanford University
%               yuehua@stanford.edu
%
% CY_14 attenuation equation - Update of the Chiou and Youngs NGA Model for
% the Average Horizontal Component of Peak Ground Motion and Response
% Spectra, Chiou and Youngs (2014)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
% M     = Moment Magnitude
% T     = Period (sec); Use Period = -1 for PGV computation
%                 Use 1000 for output the array of Sa with original period
%                 (no interpolation)
% Rup   = Closest distance (km) to the ruptured plane
% Rjb   = Joyner-Boore distance (km); closest distance (km) to surface
%       projection of rupture plane
% Rx    =    Horizontal distance from top of rupture measured perpendicular 
%       to fault strike (km). See Figures a, b and c for illustation% Ztor  = Depth(km) to the top of ruptured plane
% delta = Fault dip angle (in degree)
% lamda = Rake angle      (in degree)
% region        = 0 for global (incl. Taiwan)
%               = 1 for California
%               = 2 for Japan
%               = 3 for China 
%               = 4 for Italy 
%               = 5 for Turkey
% Z10            = Basin depth (km); depth from the groundsurface to the
%                   1km/s shear-wave horizon.
%               = 999 if unknown
%               = 'na' if unknow
% Vs30          = shear wave velocity averaged over top 30 m in m/s
%               = ref: 1130
% FVS30         = 1 for Vs30 is inferred from geology
%               = 0 for measured  Vs30
% d_DPP         = DPP centered on the site and earthquake specific average
%               DPP, = 0 for median calc (not included as a variable here)
%
% Output Variables
% Sa: Median spectral acceleration prediction
% sigma: logarithmic standard deviation of spectral acceleration
%          prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sa, yrefij, sigma] = Ph20_PGV_GMM(M, Rup, Rjb, Rx, Ztor, delta, lambda, Z10, Vs30, Fhw, FVS30, region)

delta=delta*pi()/180.0;
frv = lambda >= 30 & lambda <= 150; % frv: 1 for lambda between 30 and 150, 0 otherwise
fnm = lambda >= -120 & lambda <= -60; % fnm: 1 for lambda between -120 and -60, 0 otherwise

if Fhw == 1
    HW = 1;
elseif Fhw == 0
    HW = 0;
else
    HW = Rx>=0;
end

d_DPP=0; % for median calculatio, d_DPP=0.


[Sa, yrefij, sigma] = CY_2014_sub (M, 1, Rup, Rjb, Rx, Ztor, delta, frv, fnm, HW, Z10, Vs30, FVS30, region, d_DPP);

function [Sa, yrefij, sigma] = CY_2014_sub (M, ip, R_RUP, R_JB, Rx, Ztor, delta, F_RV, F_NM, HW, Z10, Vs30, FVS30, region, d_DPP)

%% Coefficients
c2 = 1.06;
c4 = -2.1;
c4_a = -0.5;
c_RB	= 50;
c8	= 0.2154;
c8_a	= 0.2695;
c1	= 2.2686;
c1_a	= 0.1650;
c1_b	= -0.0626;  
c1_c	= -0.1650;
c1_d	= 0.0626;
c_n	= 3.3024;
c_m	= 5.4230;	
c3	= -1.7691;
   
c5	= 5.8096;
c_HM	= 3.0514;
c6	= 0.4407;
c7	= 0.0324;

c7_b	= 0.0097;
c8_b	= 5.0000;

c9	= 0.3079;
c9_a	= 0.1000;
c9_b	= 6.5000;
c11	=  0.0000;

c11_b=	-0.3834;
c_g1	= -0.005633;
c_g2	= -0.007403;
c_g3	= 4.3439;
phi1	= -0.7161; 
phi2	= -0.0699;
   
phi3	= -0.008444;
phi4	= 5.410000;
phi5	=  0.1239;
phi6	= 300;

tau1	= 0.3882;
tau2	= 0.2578;
sigma1	= 0.4785;

sigma2	= 0.3629;
sigma3	= 0.7504;
sigma2_JP = 0.3918;
gamma_JP_IT = 2.2306;
gamma_Wn	=0.3350;
phi1_JP=	-0.7966;
phi5_JP=	0.9488;
phi6_JP=	800;


if region==2
    sigma2=sigma2_JP;
     phi1=phi1_JP; 
     phi5=phi5_JP; 
     phi6=phi6_JP;
end

%% fmag
term6=c2(ip)*(M-6);
term7=(c2(ip)-c3(ip))/c_n(ip)*log(1+exp(c_n(ip)*(c_m(ip)-M)));

%% Distance Scaling and attenuation term
term8 = c4(ip)*log(R_RUP+c5(ip)*cosh(c6(ip)*max(M-c_HM(ip),0)));
term9 = (c4_a(ip)-c4(ip))*log(sqrt(R_RUP^2+c_RB(ip)^2));
term10 = (c_g1(ip)+c_g2(ip)/(cosh(max(M-c_g3(ip),0))))*R_RUP;

if region == 2 || region == 4
    if M>6 && M<6.9
    term10= gamma_JP_IT (ip)*term10;
    end
end
if region == 3
    term10 = gamma_Wn(ip)* term10;
end

%% Style of faulting term
term2=(c1_a(ip)+c1_c(ip)/(cosh(2*max(M-4.5,0))))*F_RV;
term3=(c1_b(ip)+c1_d(ip)/cosh(2*max(M-4.5,0)))*F_NM;


%% Ztor term
    if F_RV==1
       E_Ztor = (max(3.5384-2.60*max(M-5.8530,0),0)).^2;
    else
       E_Ztor = (max(2.7482-1.7639*max(M-5.5210,0),0)).^2;
    end 
      
    if Ztor == -999
        Ztor = E_Ztor;
    end
    delta_ZTOR=Ztor-E_Ztor;
    
    term4=(c7(ip)+c7_b(ip)/cosh(2*max(M-4.5,0)))*delta_ZTOR;
%% Hanging wall term
term12=c9(ip)*HW*cos(delta)*(c9_a(ip)+(1-c9_a(ip))*tanh(Rx/c9_b(ip)))...
    *(1-sqrt(R_JB^2+Ztor^2)/(R_RUP+1));

%% Basin Depth term
% Z1.0 (m) ~ Vs30 (m/s) relationship
z_1 = exp(-3.73/2*log((Vs30.^2+290.53^2)/(1750^2+290.53^2)));

if Z10 ==-999 || Z10 == 999
    d_Z1 = 0;
else
    d_Z1 = Z10 -z_1;
end
%% Dip term
term5=(c11(ip)+c11_b(ip)/cosh(2*max(M-4.5,0)))*(cos(delta)^2);

%% Directivity
term11=c8(ip)*max(1-max(R_RUP-40,0)/30,0)*min(max(M-5.5,0)/0.8,1)...
     *exp(-c8_a(ip)*(M-c8_b(ip))^2)*d_DPP;

term1 = c1(ip);

ln_yrefij=term1+term2+term3+term4+term5+term6+term7+term8+term9+term10+...
           term11+term12;
        
yrefij=exp(ln_yrefij);

%% Site response
term14=phi1(ip)*min(log(Vs30/1130),0);
term15=phi2(ip)*(exp(phi3(ip)*(min(Vs30,1130)-360))-exp(phi3(ip)*(1130-360)))*log((yrefij+phi4(ip))/phi4(ip));
term16=phi5(ip)*(1-exp(-d_Z1/phi6(ip)));

Sa= yrefij*exp(term14+term15+term16);
%% Compute standard deviation
Finferred=(FVS30==0); % 1: Vs30 is inferred from geology.
Fmeasured=(FVS30==1); % 1: Vs30 is measured.


NL0=phi2(ip)*(exp(phi3(ip)*(min(Vs30,1130)-360))-exp(phi3(ip)*(1130-360)))*(yrefij/(yrefij+phi4(ip)));

sigmaNL0 = (sigma1(ip)+(sigma2(ip) - sigma1(ip))/1.5*(min(max(M,5),6.5)-5))*sqrt((sigma3(ip)*Finferred + 0.7* Fmeasured) + (1+NL0)^2);

tau = tau1(ip) + (tau2(ip)-tau1(ip))/1.5 * (min(max(M,5),6.5)-5);
sigma=sqrt((1+NL0)^2*tau^2+sigmaNL0^2);





