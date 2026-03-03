clc; clear; close all;

%% Save data
file_dir = fileparts(mfilename ('fullpath'));
if ~isempty(file_dir)
    file_dir = [file_dir,'\'];
end
dir_output = [file_dir, sprintf('Figures_/')];

if ~(exist(dir_output, 'dir')==7)
   mkdir(dir_output)
end

% Considered period, T
%     periods = [0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, ...
%                0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5];
    periods = [-1, 0];
    pred_models_ask14 = [];
    pred_models_ch20 = [];
    pred_models_ph20 = [];
% Give the coordinate of NPP sites
 xyNPP1 = [25.2861, 121.5879];
 xyNPP2 = [25.2026, 121.6627];
 xyNPP4 = [25.0389, 121.9247];
 % Select site condition
 sitecond = 'soil';

 for T = periods
    if strcmp(sitecond, 'rock')
        df = readtable(['dfhwPredtionGrid_rSA_M7.06_',num2str(T),'s.csv']);
        fig_name = 'NPPs_Site_VS760_SA';
    else
    % Load the prediction table
        df = readtable(['dfhwPredtionGrid_vzSA_M7.06_',num2str(T),'s.csv']);
        fig_name = 'NPPs_Site_Specific_SA';
    end
    % Model collection
    model_collects = {'ASK14', 'Ch20', 'Ph20'};
    % The corresponding prediction without HW
    Sa_pred_arr = [df.SA_ASK14, df.SA_Ch20, df.SA_Ph20];
    
    pred_models = [];
    % Load 
    for k = 1:3
        model = model_collects{k};
        % Add the map of the prediction
        Sa_pred = Sa_pred_arr(:,k).*exp(df.f1_hw);    
       
        [Sa_npp1, Sa_npp2, Sa_npp3, npp1, npp2, npp4] = Interpolate2D(xyNPP1, xyNPP2, xyNPP4, df, Sa_pred);
        pred_models = cat(1, pred_models, [Sa_npp1, Sa_npp2, Sa_npp3]); 
    end

    pred_models_ask14 = cat(1, pred_models_ask14, pred_models(1,:));
    pred_models_ch20 = cat(1, pred_models_ch20, pred_models(2,:));
    pred_models_ph20 = cat(1, pred_models_ph20, pred_models(3,:));
 end
%% Plot the hanging wall effects.

if length(periods)>2
    
    periods_set = [ 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, ...
                   0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5];
    
    fig = figure('Position', [50, 50, 1200, 400]);
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    loglog(periods_set, pred_models_ask14(:,1),'LineWidth',2); hold on;
    loglog(periods_set, pred_models_ch20(:,1),'LineWidth',2); hold on;
    loglog(periods_set, pred_models_ph20(:,1),'LineWidth',2); 
    legend('ASK14','Ch20','Ph20','Location','best');
    xlabel('Period (s)'); ylabel('SA Spectra (g)'); grid on;
    % set(gca,'Xtick',[0.001, 0.004, 0.01, 0.1, 1, 5]);
    axis([0.01, 5, 1e-2, 3]);
    if strcmp(sitecond, 'rock')
        title({['NPP1: M= 7.06, Rrup=',num2str(round(npp1.Rup,1)),...
                ', Rjb=',num2str(round(npp1.Rjb,1))],...
                 [', RX=',num2str(round(npp1.Rx,1)),...
                 ', V^S_{30}=760m/s',...
                  ', Z_{1.0}=default']});
    else
        title({['NPP1: M= 7.06, Rrup=',num2str(round(npp1.Rup,1)),...
                ', Rjb=',num2str(round(npp1.Rjb,1))],...
                 [', RX=',num2str(round(npp1.Rx,1)),...
                 ', VS30=',num2str(round(npp1.Vs30,1)),'m/s',...
                  ', Z_{1.0}=',num2str(round(npp1.z1,1)),'m']});
    end
    
    nexttile;
    loglog(periods_set, pred_models_ask14(:,2),'LineWidth',2); hold on;
    loglog(periods_set, pred_models_ch20(:,2),'LineWidth',2); hold on;
    loglog(periods_set, pred_models_ph20(:,2),'LineWidth',2); 
    legend('ASK14','Ch20','Ph20','Location','best');
    xlabel('Period (s)'); ylabel('SA Spectra (g)'); grid on;
    axis([0.01, 5, 1e-2, 3]);
    
    if strcmp(sitecond, 'rock')
        title({['NPP2: M= 7.06, Rrup=',num2str(round(npp2.Rup,1)),...
                ', Rjb=',num2str(round(npp2.Rjb,1))],...
                 [', RX=',num2str(round(npp2.Rx,1)),...
                 ', V^S_{30}=760m/s',...
                  ', Z_{1.0}=default']});
    else
        title({['NPP2: M= 7.06, Rrup=',num2str(round(npp2.Rup,1)),...
                ', Rjb=',num2str(round(npp2.Rjb,1))],...
                 [', RX=',num2str(round(npp2.Rx,1)),...
                 ', VS30=',num2str(round(npp2.Vs30,1)),'m/s',...
                  ', Z_{1.0}=',num2str(round(npp2.z1,1)),'m']});
    end
    
    nexttile;
    loglog(periods_set, pred_models_ask14(:,3),'LineWidth',2); hold on;
    loglog(periods_set, pred_models_ch20(:,3),'LineWidth',2); hold on;
    loglog(periods_set, pred_models_ph20(:,3),'LineWidth',2); 
    legend('ASK14','Ch20','Ph20','Location','best');
    xlabel('Period (s)'); ylabel('SA Spectra (g)'); grid on;
    axis([0.01, 5, 1e-2, 3]);
    if strcmp(sitecond, 'rock')
        title({['NPP4: M= 7.06, Rrup=',num2str(round(npp4.Rup,1)),...
                ', Rjb=',num2str(round(npp4.Rjb,1))],...
                 [', RX=',num2str(round(npp4.Rx,1)),...
                 ', V^S_{30}=760m/s',...
                  ', Z_{1.0}=default']});
    else
        title({['NPP4: M= 7.06, Rrup=',num2str(round(npp4.Rup,1)),...
                ', Rjb=',num2str(round(npp4.Rjb,1))],...
                 [', RX=',num2str(round(npp4.Rx,1)),...
                 ', VS30=',num2str(round(npp4.Vs30,1)),'m/s',...
                  ', Z_{1.0}=',num2str(round(npp4.z1,1)),'m']});
    end
    
      if strcmp(sitecond, 'rock')
        title(t, 'NPPs with V^S_{30} = 760 m/s Spectral Acceleration (g)', 'FontSize', 16);
      else
        title(t, 'NPPs Site-Specific Spectral Acceleration (g)', 'FontSize', 16);
      end
    
    saveas(fig, fullfile(pwd, 'Figures_', [fig_name, '.png']));

end
%%
function [Sa_npp1, Sa_npp2, Sa_npp3, npp1, npp2, npp4] = Interpolate2D(xyNPP1, xyNPP2, xyNPP4, df, Sa_pred)
    % Find the prediction along a line determined by two points
    % Generate N interpolated points between the two coordinates
    % Create interpolant function
    F = scatteredInterpolant(df.StaLat, df.StaLon, log(Sa_pred), 'linear', 'none');
    Fup = scatteredInterpolant(df.StaLat, df.StaLon, df.Rrup, 'linear', 'none');
    Fjb = scatteredInterpolant(df.StaLat, df.StaLon, df.Rjb, 'linear', 'none');
    Fx = scatteredInterpolant(df.StaLat, df.StaLon, df.Rx, 'linear', 'none');
    FVS30 = scatteredInterpolant(df.StaLat, df.StaLon, df.Vs30, 'linear', 'none');
    Fz1 = scatteredInterpolant(df.StaLat, df.StaLon, df.Z1_0, 'linear', 'none');
    % Interpolate values along the line
    Sa_npp1 = exp(F(xyNPP1(1), xyNPP1(2)));
    Sa_npp2 = exp(F(xyNPP2(1), xyNPP2(2)));
    Sa_npp3 = exp(F(xyNPP4(1), xyNPP4(2)));
    % Interpoltate the distance
    npp1.Rup = Fup(xyNPP1(1), xyNPP1(2));
    npp1.Rjb = Fjb(xyNPP1(1), xyNPP1(2));
    npp1.Rx = Fx(xyNPP1(1), xyNPP1(2));
    npp1.Vs30 = FVS30(xyNPP1(1), xyNPP1(2));
    npp1.z1 = Fz1(xyNPP1(1), xyNPP1(2));

    npp2.Rup = Fup(xyNPP2(1), xyNPP2(2));
    npp2.Rjb = Fjb(xyNPP2(1), xyNPP2(2));
    npp2.Rx = Fx(xyNPP2(1), xyNPP2(2));
    npp2.Vs30 = FVS30(xyNPP2(1), xyNPP2(2));
    npp2.z1 = Fz1(xyNPP2(1), xyNPP2(2));

    npp4.Rup = Fup(xyNPP4(1), xyNPP4(2));
    npp4.Rjb = Fjb(xyNPP4(1), xyNPP4(2));
    npp4.Rx = Fx(xyNPP4(1), xyNPP4(2));
    npp4.Vs30 = FVS30(xyNPP4(1), xyNPP4(2));
    npp4.z1 = Fz1(xyNPP4(1), xyNPP4(2));

end