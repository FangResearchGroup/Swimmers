
clear;

load('P1.mat');
load('set1.mat');
load('c_string.mat');
load('redblue.mat');

color1 = c_string{50}{1};

% Change here
% load('Inverse.mat');
% load('Forward.mat');
load('Neutral.mat');

ppcm = 1024/8.4; %140.3; % pixel/cm
fps = 60;
nv = 1.25e-2;
d0 = 10;

L = 4*ppcm;


region = [50 50 900 900]; %[left bottom width height]

RO = 1; % remove velocity outliers or not
sigma = 0.5;
vel_outlier_ratio = 3;


region_window = [50 950 50 950];
x_select = x_grid>=region_window(1) & x_grid<=region_window(2);
y_select = y_grid>=region_window(3) & y_grid<=region_window(4);
mean_v_cut = mean_v(y_select',x_select);
mean_u_cut = mean_u(y_select',x_select);
x_plot = x_grid(x_select);
y_plot = y_grid(y_select);

Down_Samp = 5:10:85;

% calculate energy flux
[Q,Q_L,Q_c,Q_s,~,~] = fst_extrapolation_demo(pvx,pvy,px,py,L,region,d0,RO,...
    ppcm,fps,nv,sigma,vel_outlier_ratio,P1,'gaussian');

disp(['Q_mean: ', num2str(mean(Q,'all')), '; Q_L: ', num2str(mean(Q_L,'all')), '; Q_c: ', num2str(mean(Q_c,'all')), '; Q_s: ', num2str(mean(Q_s,'all'))]);

x = region(1):d0:(region(1)+region(3));
y = region(2):d0:(region(2)+region(4));
if mod(numel(x),2) % odd number of elements
    x=x(1:end-1)+d0/2; % be sure to get an even number
end
if mod(numel(y),2) % odd number of elements
    y=y(1:end-1)+d0/2; % be sure to get an even number
end

%
% Plot single frame of Q
RedBlue = redblue_v2(384,0.5);
% RedBlue = RedBlue(64:320,:);

fig = figure(); hold on

imagesc(x,y,Q_s);
% quiver(px,py,pvx*5,pvy*5,'off','r-');
% quiver(px(good),py(good),pvx(good)*5,pvy(good)*5,'off','k-','LineWidth',0.2);
quiver(x_plot(Down_Samp),y_plot(Down_Samp),mean_u_cut(Down_Samp,Down_Samp),mean_v_cut(Down_Samp,Down_Samp),0.8,'LineWidth',0.2,'Color',color1,'ShowArrowHead','on');

colormap(RedBlue);
% clb = colorbar;
clim([-4e-3 4e-3]);
% clb.LineWidth = 0.2;
% clb.Ticks = -6e-4:2e-4:8e-4;

axis([50 950 50 950]);
% axis([45 955 45 955]);

ax = gca;
set(gcf,'MenuBar','figure',...
    'Units','centimeters',...
    'Position',[10,10,4,4],...
    'Resize',0,...
    'PaperUnits','centimeters',...
    'PaperSize',[4 4]);

set(gca,'FontName','Arial',...
    'FontSize',7,... 
    'FontUnits','Points',...
    'XTick',[],...
    'YTick',[],...
    'XTicklabel',[],...
    'YTicklabel',[],...
    'Visible','on',...
    'YDir','normal',...
    'Box','on',...
    'linewidth',0.3,...
    'Units','centimeters',...
    'looseInset',[0,0,0,0]...
    );
ax.PlotBoxAspectRatio = [1,1,1];