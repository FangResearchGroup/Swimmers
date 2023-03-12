clear;
mainfolder = 'Q_sets_window';
addpath(genpath(mainfolder));
subfolder = {'S','+','-','=','F'};
Q_save_all = cell(1,5);
filename = 'Q_*.mat';

load('c_string.mat');

for kk = 1:5
    filelist = dir(fullfile(mainfolder,subfolder{kk},filename));
    Q_all = [];
    for i = 1:length(filelist)
        load(fullfile(filelist(i).folder,filelist(i).name));
        Q_all = [Q_all;Q_s_save];
    end
    Q_save_all{kk} = Q_all;
end

%
Q_mean_shear = mean(Q_save_all{1});
Q_mean_inv = mean(Q_save_all{2});
Q_mean_for = mean(Q_save_all{3});
Q_mean_neu = mean(Q_save_all{4});
Q_mean_fluc = mean(Q_save_all{5});

% figure properties
LineWidth = 2;

fig = figure(); hold on
ax = gca;

plot(L/ppcm,Q_mean_shear,'LineWidth',LineWidth);
plot(L/ppcm,Q_mean_inv,'LineWidth',LineWidth);
plot(L/ppcm,Q_mean_for,'LineWidth',LineWidth);
plot(L/ppcm,Q_mean_neu,'LineWidth',LineWidth);
plot(L/ppcm,Q_mean_fluc,'-','LineWidth',LineWidth);

colororder({c_string{2}{2},c_string{5}{1},c_string{51}{4},c_string{6}{3},c_string{2}{3}});

xlb = xlabel('L (cm)');
ylb = ylabel('Q_{s} (m^2/s^3)');

lgd = legend('Pure shear','\theta_s < \pi/4', '\theta_s > \pi/4', '\theta_s \approx \pi/4', 'No shear');
lgd.NumColumns = 1;
lgd.Units = 'centimeters';
lgd.Position = [1.04,0.9,1.73,1.48];
lgd.ItemTokenSize = [10,1];
legend boxoff

axis([0, 7, -4e-3, 3e-3]);

set(gcf,'MenuBar','figure',...
    'Units','centimeters',...
    'Position',[5,5,6,6],...
    'Resize',0,...
    'PaperUnits','centimeters',...
    'PaperSize',[6 6]);

set(gca,'FontName','Arial',...
    'FontSize',7,... 
    'FontUnits','Points',...       
    'Xtick',0:1:7,...
    'Box','on',...
    'LineWidth',0.3,...
    'Units','centimeters',...
    'looseInset',[0,0,0,0]...
    );
ax.PlotBoxAspectRatio = [1,1,1];