function [Q,Q_L,Q_c,Q_s,x,y] = fst_extrapolation_demo(pu,pv,px,py,L,region,d0,RO,ppcm,fps,nv,sigma,vel_outlier_ratio,P1,filtertype)

% Calculate the FST of a flow field as well as the three terms of Q.
% For the energy transfer (Q), positive values represent down-scale
% transfer, i.e., transfer to smaller length scales (as in a forward cascade).

%********************************************************************
% This version is used for shear flow with major velocity gradient in
% dV/dx. The flow field is extrapolated in order to avoid boundary
% artifacts due to zero padding.
% Since the horizontal velocity is small, extrapolated velocity field
% has horizontal velocty 0.
%********************************************************************

% Input:
% filename: file location where velocity field information saves
% framerange: [min max]
% L: length scale of filtering in unit of pixel
% region: [minx miny rangex rangey]
% d0: grid space for 2D interpolation
% RO: remove velocity outlier or not, if 1, removes outliers with 3 iqr
% ppcm: pixel per cm
% fps: frame per second
% nv: kinematic viscosity
% sigma: window width for gauss (or whatever) window
% vel_outlier_ratio: velocity outlier ratio
% P1: extrapolation coefficient in x direction: y = y0 + P1 * x
% filtertype: 'sharp' or 'gaussian'

% Output [all in units of cm and s]:
% Q: energy flux through length scale L
% Q_L: Q corresponding to the Leonard term
% Q_c: Q correponding to the cross term
% Q_s: Q corresponding to the subgrid term
% x: x location
% y: y location

%% Pre-determination
sigma_default = 0.5; % window width for sharpfilt2.m
d0default = 20; % pixels; gives roughly 10k grid pts
vel_outlier_ratio_default = 3; % drop particles w/ velocity > 3*iqr
nv_default = 1.25e-2; %cm^2/s
padding = 0;

if ~exist('d0','var') || isempty(d0)
    d0=d0default;
elseif numel(d0)>1
    error('Sorry, d0 must be a single scalar.');
end
if ~exist('nv','var') || isempty(nv)
    nv = nv_default;
end
if ~exist('sigma','var') || isempty(sigma)
    sigma = sigma_default;
end
if ~exist('vel_outlier_ratio','var') || isempty(vel_outlier_ratio)
    vel_outlier_ratio = vel_outlier_ratio_default;
end

%% main function
Wp=2*d0/L; % cutoff freq, norm by Nyquist

region_x = [region(1), region(1)+region(3), region(1)+region(3), region(1), region(1)];
region_y = [region(2), region(2), region(2)+region(4), region(2)+region(4), region(2)];

% setup mesh grid
x=region(1):d0:(region(1)+region(3)); % set up x coords of grid
if mod(numel(x),2) % odd number of elements
    x=x(1:end-1)+d0/2; % be sure to get an even number
end

y=region(2):d0:(region(2)+region(4)); % set up y coords of grid
if mod(numel(y),2) % odd number of elements
    y=y(1:end-1)+d0/2; % be sure to get an even number
end

% Need to use matlab version 2012b or higher
begins = 1;
ends = length(pu);
Nt=1; % number of total frames

% pre-allocation
QQ = nan*zeros(numel(y),numel(x),Nt);
QQ_L = nan*zeros(numel(y),numel(x),Nt);
QQ_c = nan*zeros(numel(y),numel(x),Nt);
QQ_s = nan*zeros(numel(y),numel(x),Nt);

%% main loop

% loop over each single frame
for ii=1:Nt
    
    % -=- extract a single frame
    ind=begins(ii):ends(ii);
    if numel(ind)<3 % need at least 3 pts for griddata
        continue
    end

    % -=- interpolate the velocity data into gridded form -=-
    if RO % remove outlier velocity vectors
        pv_mag = sqrt(pu(ind).^2 + pv(ind).^2); % velocity magnitude of this frame
        pviqr = iqr(pv_mag); % interquarter range of pv
        good = pv_mag < vel_outlier_ratio * pviqr;
        if any(~good) && sum(good)>=3 % need at least 3 pts for griddata
            good_inpolygon = inpolygon(px(ind(~good)),py(ind(~good)),region_x,region_y); % count number of outliers in DOI
            disp(['Dropping ' num2str(sum(~good)) ...
                ' outlier particle(s),' num2str(sum(good_inpolygon)) ' in selected region.']);
            u=griddata(px(ind(good)),py(ind(good)),pu(ind(good)),x,y'); % grid velocity
            v=griddata(px(ind(good)),py(ind(good)),pv(ind(good)),x,y');
        else
            u=griddata(px(ind),py(ind),pu(ind),x,y'); % grid velocity
            v=griddata(px(ind),py(ind),pv(ind),x,y');
        end
    else % directly perform the interpolation
        u=griddata(px(ind),py(ind),pu(ind),x,y'); % grid velocity
        v=griddata(px(ind),py(ind),pv(ind),x,y');
    end
    u(isnan(u))=0; % kill off the NaNs
    v(isnan(v))=0;

    % -=- extrapolation
    x_temp = d0:d0:3*length(x)*d0; % extend each side by one domain width
    y_temp = d0:d0:3*length(y)*d0; % extend each side by one domain width

    u = [u(1,:).*ones(size(u)); u; u(end,:).*ones(size(u))];
    v = [v(1,:).*ones(size(v)); v; v(end,:).*ones(size(v))];

    u = [u(:,1).*ones(size(u)), u, u(:,end).*ones(size(u))];
    v = [v(:,1)+d0*ones(size(v,1),1).*(-size(v,2):1:-1)*P1, v, v(:,end)+d0 *ones(size(v,1),1).*(1:1:size(v,2))*P1];

    
    % -=- Filter everything -=-
    switch filtertype
        case 'sharp'
            u2=sharpfilt2(u,Wp,sigma,1,padding);
            v2=sharpfilt2(v,Wp,sigma,1,padding);

            uu2=sharpfilt2(u.^2,Wp,sigma,1,padding);
            vv2=sharpfilt2(v.^2,Wp,sigma,1,padding);
            uv2=sharpfilt2(u.*v,Wp,sigma,1,padding);

            % for Leonard term
            uu3 = sharpfilt2(u2.*u2,Wp,sigma,1,padding);
            vv3 = sharpfilt2(v2.*v2,Wp,sigma,1,padding);
            uv3 = sharpfilt2(u2.*v2,Wp,sigma,1,padding);

            % for cross term
            uu4 = sharpfilt2(u2.*(u - u2),Wp,sigma,1,padding);
            vv4 = sharpfilt2(v2.*(v - v2),Wp,sigma,1,padding);
            uv4 = sharpfilt2(u2.*(v - v2),Wp,sigma,1,padding);
            vu4 = sharpfilt2(v2.*(u - u2),Wp,sigma,1,padding);

            % for subgrid term
            uu5 = sharpfilt2((u - u2).*(u - u2),Wp,sigma,1,padding);
            vv5 = sharpfilt2((v - v2).*(v - v2),Wp,sigma,1,padding);
            uv5 = sharpfilt2((u - u2).*(v - v2),Wp,sigma,1,padding);
        case 'gaussian'
            u2=gaussfilt2(u,Wp,padding);
            v2=gaussfilt2(v,Wp,padding);

            uu2=gaussfilt2(u.^2,Wp,padding);
            vv2=gaussfilt2(v.^2,Wp,padding);
            uv2=gaussfilt2(u.*v,Wp,padding);

            % for Leonard term
            uu3 = gaussfilt2(u2.*u2,Wp,padding);
            vv3 = gaussfilt2(v2.*v2,Wp,padding);
            uv3 = gaussfilt2(u2.*v2,Wp,padding);

            % for cross term
            uu4 = gaussfilt2(u2.*(u - u2),Wp,padding);
            vv4 = gaussfilt2(v2.*(v - v2),Wp,padding);
            uv4 = gaussfilt2(u2.*(v - v2),Wp,padding);
            vu4 = gaussfilt2(v2.*(u - u2),Wp,padding);

            % for subgrid term
            uu5 = gaussfilt2((u - u2).*(u - u2),Wp,padding);
            vv5 = gaussfilt2((v - v2).*(v - v2),Wp,padding);
            uv5 = gaussfilt2((u - u2).*(v - v2),Wp,padding);
    end

    % rate of strain term
    [u2x,u2y]=gradient(u2,x_temp,y_temp);
    [v2x,v2y]=gradient(v2,x_temp,y_temp);

    % -=- remove the extrapolated frame
    u = u(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    v = v(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    u2=u2(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    v2=v2(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    uu2=uu2(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    vv2=vv2(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    uv2=uv2(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    u2x = u2x(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    u2y = u2y(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    v2x = v2x(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    v2y = v2y(length(x)+1:length(x)*2,length(y)+1:length(y)*2);

    uu3 = uu3(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    vv3 = vv3(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    uv3 = uv3(length(x)+1:length(x)*2,length(y)+1:length(y)*2);

    uu4 = uu4(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    vv4 = vv4(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    uv4 = uv4(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    vu4 = vu4(length(x)+1:length(x)*2,length(y)+1:length(y)*2);

    uu5 = uu5(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    vv5 = vv5(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    uv5 = uv5(length(x)+1:length(x)*2,length(y)+1:length(y)*2);
    
    % -=- Calc energy transfer between scales -=-
    QQ(:,:,ii) = (u2.^2-uu2).*u2x + (v2.^2-vv2).*v2y + ...
        (u2.*v2-uv2).*(u2y+v2x); % - tau_ij S_ij in frames and pixels
    
    % Leonard term
    QQ_L(:,:,ii) = (u2.*u2 - uu3).*u2x + (v2.*v2 - vv3).*v2y + ...
        (u2.*v2 - uv3).*(u2y+v2x);

    % Cross term
    QQ_c(:,:,ii) = -1 * (2*uu4.*u2x + 2*vv4.*v2y + (uv4 + vu4).*(u2y + v2x));

    % Subgrid term
    QQ_s(:,:,ii) = -1 * (uu5.*u2x + vv5.*v2y + uv5.*(u2y+v2x));

%     % -=- Calc energy dissipation in larger-than-filter length scale -=-
%     % Epsilon(:,:,ii) = 2*nv*(u2x.*u2x + 0.25*(u2y+v2x).*(u2y+v2x) + v2y.*v2y); % 2 nv S_ij S_ij in frames and pixels [pope]
%     Epsilon(:,:,ii) = nv*(u2x.*u2x + u2y.*u2y + v2x.*v2x + v2y.*v2y); % nv du_i/dx_j du_i/dx_j [Eyink 1995]
%     
%     % -=- Calculate the stress strain eigenvalue and angle -=-
%     trace_of_stress = 0.5*(uu2-u2.*u2)+0.5*(vv2-v2.*v2); % mean of trace of stress tensor (tau_ij)
%     
% %     eigval_stress=zeros(size(QQ));
%     for i = 1:size(QQ,1)
%         for j=1:size(QQ,2)
%             stress=[uu2(i,j)-u2(i,j)*u2(i,j)-trace_of_stress(i,j), uv2(i,j)-u2(i,j)*v2(i,j);...
%                     vu2(i,j)-v2(i,j)*u2(i,j), vv2(i,j)-v2(i,j)*v2(i,j)-trace_of_stress(i,j)]; % traceless
%             eigval_stress(i,j,ii)=max(abs(eig(stress)));
%             
%             [V,D] = eig(stress);
%             eig_vec = V(:,diag(D)==max(diag(D)));
%             eigvec_stress(i,j,ii) = atan(eig_vec(2)/eig_vec(1));
%         end
%     end
%        
% %     eigval_strain=zeros(size(QQ));
%     for k=1:size(QQ,1)
%         for l=1:size(QQ,2)
%             strain=[u2x(k,l), 0.5*(u2y(k,l)+v2x(k,l));...
%                     0.5*(u2y(k,l)+v2x(k,l)), v2y(k,l)];
%             eigval_strain(k,l,ii)=max(abs(eig(strain)));
%             % eigval_strain(k,l,ii)=sqrt(0.25*(strain(1,1)-strain(2,2))^2 + strain(1,2)^2);
% 
%             [V,D] = eig(strain);
%             eig_vec = V(:,diag(D)==max(diag(D)));
%             eigvec_strain(k,l,ii) = atan(eig_vec(2)/eig_vec(1));
%         end
%     end
%     
%     Cos_2theta(:,:,ii) = -(1/2)*QQ(:,:,ii)./eigval_stress./eigval_strain;
%     thetatheta(:,:,ii) = 0.5*acos(Cos_2theta);        
end

% Unit conversion
Q = QQ/(ppcm^2)*(fps^3);
Q_L = QQ_L/(ppcm^2)*(fps^3);
Q_c = QQ_c/(ppcm^2)*(fps^3);
Q_s = QQ_s/(ppcm^2)*(fps^3);
x = x/ppcm;
y = y/ppcm;

if numel(Q)==0
    warning('MATLAB:energyxfer:tooFewParticles', ...
        'Found too few particles to calculate energy transfer.');
end

end