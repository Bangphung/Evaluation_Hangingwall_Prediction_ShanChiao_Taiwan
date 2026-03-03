clc; clear; close all;

addpath("D:\My_Study\Gaussian Process&RVT Study\ASK14_HW_model\gmm");
% Considered period, T
 T = 0.01;    
% Load the prediction table
df = readtable(['dfhwPredtionGrid_SA_M7.06_',num2str(T),'s.csv']);
df = df(df.Rrup <= 30.0,:); [Rx, id] = sort(df.Rx);
df = df(id, :);

model_id1 = 1; model_id2 = 2; model_id3 = 3;
model_id4 = 4; model_id5 = 5;

M = 7.06;
W = df.W(1);
Rx = [-20:1:30]';
Rrup = 5; 
Rjb = 0;
Ztor = 0;
Dip = df.Dip(1);

f1_hw = zeros(length(df.W),1);
f2_hw = zeros(length(df.W),1);
f3_hw = zeros(length(df.W),1);
f4_hw = zeros(length(df.W),1);
f5_hw = zeros(length(df.W),1);

for i = 1:length(Rx)
    f1_hw(i,1) = swus_hanging_wall_model(T, M, Dip, Rx(i), Rrup, Rjb, W, Ztor, model_id1);
    f2_hw(i,1) = swus_hanging_wall_model(T, M, Dip, Rx(i), Rrup, Rjb, W, Ztor, model_id2);
    f3_hw(i,1) = swus_hanging_wall_model(T, M, Dip, Rx(i), Rrup, Rjb, W, Ztor, model_id3);
    f4_hw(i,1) = swus_hanging_wall_model(T, M, Dip, Rx(i), Rrup, Rjb, W, Ztor, model_id4);
    f5_hw(i,1) = swus_hanging_wall_model(T, M, Dip, Rx(i), Rrup, Rjb, W, Ztor, model_id5);
end

Rx_new = [-30:0.5:30]; [Rx, ip] = unique(Rx);
f1new_hw = interp1(Rx, f1_hw(ip), Rx_new);
f2new_hw = interp1(Rx, f2_hw(ip), Rx_new);
f3new_hw = interp1(Rx, f3_hw(ip), Rx_new);
f4new_hw = interp1(Rx, f4_hw(ip), Rx_new);
f5new_hw = interp1(Rx, f5_hw(ip), Rx_new);

fig=figure;
plot(Rx_new,f1new_hw,'r','LineWidth',2); hold on;
plot(Rx_new,f2new_hw,'b','LineWidth',2); hold on;
plot(Rx_new,f3new_hw,'c','LineWidth',2); hold on;
plot(Rx_new,f4new_hw,'gr','LineWidth',2); hold on;
plot(Rx_new,f5new_hw,'m','LineWidth',2); 
legend('model_1','model_2','model_3','model_4','model_5',...
    'Location','best');
xlabel('R_X (km)'); ylabel('HW factor'); grid on;
title(['T=',num2str(T),'s, M=7.0, REV, Ztor = 0 km, and Dip = 49^o']);

fig_name = ['Fhw_T=', num2str(T)];
saveas(fig, fullfile(pwd, 'Figures_', [fig_name, '.png']));