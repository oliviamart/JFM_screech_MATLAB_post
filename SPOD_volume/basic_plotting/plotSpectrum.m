% Olivia Martin
% plot the SPOD spectrum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultAxesFontWeight','default');
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultAxesLineWidth',1.5);
set(0,'DefaultAxesTickDir','in');
set(0,'DefaultAxesTickLength',[0.025 0.08]);
set(0,'DefaultLineLineWidth',2.5);
set(0,'defaultTextInterpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spod info
mode_idx = 29;
filename1 = '../SPOD_data/Mj_169/spod_energy.mat';

% flow conditions
M = 1.69; 
dt = 0.2;
numT = 200;

% saving info
doSave = 0;
figName = 'figures/Mj_169_spectra.jpg';

% geometry
h = 0.010; b = 0.040;
c = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L1 = load(filename1).L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET FREQUENCY SCALING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = dt:dt:T;
fs = 1/dt;
numT = length(f);
f = fs*(0:(numT/2))/numT;

De = 2 * sqrt(h*b/pi);
U = M*(1+0.2*M^2)^(-1/2);
f = f * c / h; % Convert from acoustic time units back to regular
St = f * De / U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPECTRUM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col = colormap('lines');

g = figure
subplot(1, 1, 1)
for i = 1:size(L1,2)
    colorval = col(i,:);
    lin_width = 2;
    if i == 1
        lin_width = 3;
    end
    semilogy(St, L1(:,i), '-', 'color', colorval, 'lineWidth', lin_width)

    xlabel('S$t = fD_{e}/U_{j}$', 'fontSize', 24, 'fontweight','normal')
    ylabel('$\lambda$', 'fontSize', 26, 'fontweight','normal')
    xlim([min(St) max(St)])
    ylim([10^-1.5 10^2.5])
    set(gca,'defaultTextInterpreter','tex');
    hold on

end
legend('mode 1', 'mode 2', 'mode 3','mode 4', 'mode 5', 'mode 6', 'location', 'southwest', 'FontSize', 24)

if doSave
    saveas(gcf, figName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRA PLOTTING THINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xlim([0 0.3])
%% Print energy contribution from first mode
%energy = L1(mode_idx,1)/sum(L1(mode_idx,:))*100;
%fprintf(strcat('\n\nMode 1 contains: \n', num2str(energy),' percent of total energy\n'))
%%
% %%
% subplot(1, 3, 2)
% for i = 1:size(L2,2)
%     colorval = col(i*40,:);
%     lin_width = 1;
%     if i == 1
%         %colorval = 'red';
%         lin_width = 2;
%     end
%     loglog(L2(:,i), 'color', colorval, 'lineWidth', lin_width)
%     %set(gca,'fontsize', 20, 'fontweight','normal');
%     xlabel('S$t = fD_{e}/U_{j}$', 'fontSize', 24, 'fontweight','normal')
%     ylabel('$\lambda$', 'fontSize', 26, 'fontweight','normal')
%     %xlim([min(St) max(St)])
%     %ylim([10^-1 10^3.5])
%     set(gca,'defaultTextInterpreter','tex');
%     hold on
% 
% end
% ylim([10^-3 10^3])
% %legend('mode 1', 'mode 2', 'mode 3','mode 4', 'mode 5', 'mode 6', 'location', 'southwest')
% 
% 
% subplot(1, 3, 3)
% for i = 1:size(L2,2)
%     colorval = col(i*40,:);
%     lin_width = 1;
%     if i == 1
%         %colorval = 'red';
%         lin_width = 2;
%     end
%     loglog((L2(:,i) - L1(:,i))./L2(:,i)*100, 'color', colorval, 'lineWidth', lin_width)
%     %set(gca,'fontsize', 20, 'fontweight','normal');
%     xlabel('S$t = fD_{e}/U_{j}$', 'fontSize', 24, 'fontweight','normal')
%     ylabel('$\lambda$', 'fontSize', 26, 'fontweight','normal')
%     %xlim([min(St) max(St)])
%     %ylim([10^-1 10^3.5])
%     set(gca,'defaultTextInterpreter','tex');
%     hold on
% 
% end
% ylim([10^-3 10^3])
% legend('mode 1', 'mode 2', 'mode 3','mode 4', 'mode 5', 'mode 6', 'location', 'eastoutside')
% 
