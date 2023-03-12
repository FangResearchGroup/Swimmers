clear;
load('set1_mean.mat');
load('set1_eigenVV.mat')
load('c_string.mat');
color1 = c_string{50}{1};
color2 = c_string{24}{4};
color3 = c_string{3}{6};

d0 = 10;
ppcm = 1024/8.4;
region_window = [50 940 50 940]; %[left,right,bottom,top]


x_select = x_grid>=region_window(1) & x_grid<=region_window(2);
y_select = y_grid>=region_window(3) & y_grid<=region_window(4);
mean_v_cut = mean_v(y_select',x_select);
mean_u_cut = mean_u(y_select',x_select);
% v_cut = v(y_select',x_select,:);
% u_cut = u(y_select',x_select,:);
% v_cut = num2cell(v_cut,[1,2]);
% v_cut = vertcat(v_cut{:});

eigenvx = eigvec_x(y_select',x_select);
eigenvy = eigvec_y(y_select',x_select);

x_plot = x(x_select);
x_plot = x_plot - x_plot(1);
y_plot = y(y_select);
y_plot = y_plot - y_plot(1);

Down_Samp = 5:8:85;


fig = figure(); hold on

% need down sample
yyaxis right
quiver(x_plot(Down_Samp),y_plot(Down_Samp),mean_u_cut(Down_Samp,Down_Samp),mean_v_cut(Down_Samp,Down_Samp),0.7,'-','LineWidth',1,'Color',color1);
quiver(x_plot(Down_Samp),y_plot(Down_Samp),eigenvx(Down_Samp,Down_Samp),eigenvy(Down_Samp,Down_Samp),0.5,'-','LineWidth',1,'Color',color3);
axis([0 7 0 7])
xlb = xlabel('x (cm)');
ylb = ylabel('y (cm)');
set(gca,'Ytick',0:1:7,'ycolor','k');


fillpatch_x = [x_plot,flip(x_plot)];
fillpatch_y = [mean(mean_v_cut)+v_cut_std,flip(mean(mean_v_cut)-v_cut_std)];
yyaxis left
fill(fillpatch_x,fillpatch_y,'b','FaceColor',color2,'FaceAlpha',0.2,'EdgeColor','none');
plot(x_plot,mean(mean_v_cut),'-','linewidth',1,'Color',color2);
% errorbar(x_plot,mean(mean_v_cut),std(v_cut),'linewidth',0.5,'Color',color2);

axis([0 7 -1.5 1.5]);

p = polyfit(x(x_select),mean(mean_v_cut),1);
P1 = p(1);
mean_v_y = mean(mean_v_cut);
P_mean = mean(diff(mean_v_y)/(d0/ppcm));
ylb2 = ylabel('v (cm/s)');
set(gca,'Ytick',-1.5:0.5:1.5,'ycolor',color2);

set(gcf,'MenuBar','figure',...
    'Units','centimeters',...
    'Position',[3,3,6.6,6],...
    'Resize',0,...
    'PaperUnits','centimeters',...
    'PaperSize',[6.6,6]);

ax = gca;
set(gca,'FontName','Arial',...
    'FontSize',7,... 
    'FontUnits','Points',...       
    'Box','on',...
    'LineWidth',0.5,...
    'Units','centimeters',...
    'looseInset',[0,0,0,0]...
    );
ax.PlotBoxAspectRatio = [1,1,1];