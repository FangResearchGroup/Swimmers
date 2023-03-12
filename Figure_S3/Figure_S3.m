%% Low Re case
clear;
addpath(genpath('C:\Users\xinyu\Documents\MATLAB_master\MatLab_Toolbox'));
% addpath(genpath('X:\Xinyu_Si_X\MATLAB_master\MatLab_Toolbox'));
crtfd = pwd;
% Add data folder to path and load data
mainfolder = 'low';
addpath(genpath(mainfolder));

load('Q_mode.mat');
load('c_string.mat');

fig = figure(); hold on
ax = gca;
LineWidth = 2;


plot(L/ppcm,mean(Q_save,1),'-','LineWidth',LineWidth,'Color',c_string{5}{1});
plot(L/ppcm,mean(Q_L_save,1),'-','LineWidth',LineWidth,'Color',c_string{51}{4});
plot(L/ppcm,mean(Q_c_save,1),'-','LineWidth',LineWidth,'Color',c_string{2}{3});
plot(L/ppcm,mean(Q_s_save,1),'-','LineWidth',LineWidth,'Color',c_string{6}{3});


axis([0, 16, -9e-6, 6e-6]);
% axis([0 20 -2e-3 2e-3]);



xlb = xlabel('L (cm)');
ylb = ylabel('{Q}^{(L)} (cm^2/s^3)');
lgd = legend('Q','Q_{Leo}', 'Q_c', 'Q_s');
lgd.NumColumns = 1;
lgd.Units = 'centimeters';
lgd.Position = [3.58,0.89,1.73,1.48];
lgd.ItemTokenSize = [10,1];
legend boxoff

set(gcf,'MenuBar','figure',...
    'Units','centimeters',...
    'Position',[5,5,6,6],...
    'Resize',0,...
    'PaperUnits','centimeters',...
    'PaperSize',[6 6]);

set(gca,'FontName','Arial',...
    'FontSize',7,... 
    'FontUnits','Points',...       
    'Xtick',0:4:16,...
    'Box','on',...
    'LineWidth',0.3,...
    'Units','centimeters',...
    'looseInset',[0,0,0,0]...
    );
ax.PlotBoxAspectRatio = [1,1,1];

%% High Re case

clear;
addpath(genpath('C:\Users\xinyu\Documents\MATLAB_master\MatLab_Toolbox'));
% addpath(genpath('X:\Xinyu_Si_X\MATLAB_master\MatLab_Toolbox'));
crtfd = pwd;
% Add data folder to path and load data
mainfolder = 'high';
addpath(genpath(mainfolder));

load('Q_mode.mat');
load('c_string.mat');

fig = figure(); hold on
ax = gca;
LineWidth = 2;


plot(L/ppcm,mean(Q_save,1),'-','LineWidth',LineWidth,'Color',c_string{5}{1});
plot(L/ppcm,mean(Q_L_save,1),'-','LineWidth',LineWidth,'Color',c_string{51}{4});
plot(L/ppcm,mean(Q_c_save,1),'-','LineWidth',LineWidth,'Color',c_string{2}{3});
plot(L/ppcm,mean(Q_s_save,1),'-','LineWidth',LineWidth,'Color',c_string{6}{3});


axis([0, 16, -1.4e-2, 6e-3]);
% axis([0 20 -0.2 0.15]);


xlb = xlabel('L (cm)');
ylb = ylabel('{Q}^{(L)} (cm^2/s^3)');
lgd = legend('Q','Q_{Leo}', 'Q_c', 'Q_s');
lgd.NumColumns = 1;
lgd.Units = 'centimeters';
lgd.Position = [3.58,0.89,1.73,1.48];
lgd.ItemTokenSize = [10,1];
legend boxoff

set(gcf,'MenuBar','figure',...
    'Units','centimeters',...
    'Position',[5,5,6,6],...
    'Resize',0,...
    'PaperUnits','centimeters',...
    'PaperSize',[6 6]);

set(gca,'FontName','Arial',...
    'FontSize',7,... 
    'FontUnits','Points',...       
    'Xtick',0:4:16,...
    'Box','on',...
    'LineWidth',0.3,...
    'Units','centimeters',...
    'looseInset',[0,0,0,0]...
    );
ax.PlotBoxAspectRatio = [1,1,1];
