%% load toolbox path
clear;
load('data.mat');
load('c_string.mat');
color1 = c_string{50}{1};

iqr_cut = 0; % if 1, apply iqr process to remove outliers
vel_outlier_ratio=3; % drop particles w/ velocity > 3*iqr

d0 = 5;
region = [0 0 1024 1024];

x_grid = region(1):d0:(region(1) + region(3));
y_grid = region(2):d0:(region(2) + region(4));

if mod(numel(x_grid),2) % odd number of elements
    x_grid=x_grid(1:end-1)+d0/2; % be sure to get an even number
    x = x_grid;
end

if mod(numel(y_grid),2) % odd number of elements
    y_grid=y_grid(1:end-1)+d0/2; % be sure to get an even number
    y = y_grid;
end

[x0,y0] = meshgrid(x,y);

ppcm = 1024/8.4;

if iqr_cut
    pvv = sqrt(pvx.^2 + pvx.^2); % velocity magnitude of this frame
    pviqr = iqr(pvv); % interquarter range of pv
    good = pvv < vel_outlier_ratio * pviqr;
    if sum(~good)>0
        disp(['Dropping ' num2str(sum(~good)) ...
            ' outlier particle(s) from frame ' num2str(Frame) '.'])
    end
    pvx_grid = griddata(px(good),py(good),pvx(good),x0,y0,'cubic');
    pvy_grid = griddata(px(good),py(good),pvy(good),x0,y0,'cubic');
else
    pvx_grid = griddata(px,py,pvx,x0,y0,'cubic');
    pvy_grid = griddata(px,py,pvy,x0,y0,'cubic');
end
pvx_grid(isnan(pvx_grid)) = 0;
pvy_grid(isnan(pvy_grid)) = 0;

% pvx_grid = pvx_grid/ppcm*fps; % change to cm/s
% pvy_grid = pvy_grid/ppcm*fps; % change to cm/s

% save pvx and pvy
u = pvx_grid; % in cm/s
v = pvy_grid; % in cm/s
%
L = ppcm*4;
Wp=2*d0/L;
sigma = 0.8;

u = [u(1,:).*ones(size(u)); u; u(end,:).*ones(size(u))];
v = [v(1,:).*ones(size(v)); v; v(end,:).*ones(size(v))];

u = [u(:,1).*ones(size(u)), u, u(:,end).*ones(size(u))];
v = [v(:,1)+d0*ones(size(v,1),1).*(-size(v,2):1:-1)*P1, v, v(:,end)+d0 *ones(size(v,1),1).*(1:1:size(v,2))*P1];

padding = 0;

u2=sharpfilt2(u,Wp,sigma,1,padding);
v2=sharpfilt2(v,Wp,sigma,1,padding);

u = u(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
v = v(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
u2=u2(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
v2=v2(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
%
fps = 60;
region_window = [50 940 50 940]; %[left,right,bottom,top]
x_select = x_grid>=region_window(1) & x_grid<=region_window(2);
y_select = y_grid>=region_window(3) & y_grid<=region_window(4);
v_cut = v(y_select',x_select);
u_cut = u(y_select',x_select);
v_cut2 = v2(y_select',x_select);
u_cut2 = u2(y_select',x_select);

x_plot = x(x_select);
x_plot = x_plot - x_plot(1);
y_plot = y(y_select);
y_plot = y_plot - y_plot(1);

Down_Samp = 1:2:178;


fig = figure(); hold on

% plot (a)
% quiver(x_plot(Down_Samp)/ppcm,y_plot(Down_Samp)/ppcm,...
%     u_cut(Down_Samp,Down_Samp)/ppcm*fps,v_cut(Down_Samp,Down_Samp)/ppcm*fps,...
%     2.5,'-','LineWidth',1,'Color',color1);

% plot (b) and (c)
% quiver(x_plot(Down_Samp)/ppcm,y_plot(Down_Samp)/ppcm,...
%     u_cut2(Down_Samp,Down_Samp)/ppcm*fps,v_cut2(Down_Samp,Down_Samp)/ppcm*fps,...
%     2.5,'-','LineWidth',1,'Color',color1);

% plot (d)
quiver(x_plot(Down_Samp)/ppcm,y_plot(Down_Samp)/ppcm,...
    u_cut(Down_Samp,Down_Samp)/ppcm*fps - u_cut2(Down_Samp,Down_Samp)/ppcm*fps,...
    v_cut(Down_Samp,Down_Samp)/ppcm*fps - v_cut2(Down_Samp,Down_Samp)/ppcm*fps,...
    2.5,'-','LineWidth',1,'Color',color1);

axis([0 7 2.5 4.5]);
xlb = xlabel('x (cm)');
ylb = ylabel('y (cm)');
set(gca,'Ytick',3:4,'ycolor','k');

set(gcf,'MenuBar','figure',...
    'Units','centimeters',...
    'Position',[3,3,16,9],...
    'Resize',0,...
    'PaperUnits','centimeters',...
    'PaperSize',[8,8]);

ax = gca;
set(gca,'FontName','Arial',...
    'FontSize',7,... 
    'FontUnits','Points',...       
    'Box','on',...
    'LineWidth',0.5,...
    'Units','centimeters',...
    'looseInset',[0,0,0,0]...
    );
ax.PlotBoxAspectRatio = [7,2,1];