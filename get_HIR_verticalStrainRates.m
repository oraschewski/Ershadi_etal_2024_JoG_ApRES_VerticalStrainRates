% Deriving vertical strain rates from repeat Autonomous phase-sensitive
% Radio Echo Sounder (ApRES) measurements on Hammarryggen ice rise
%
% Supplement to:
% Ershadi et al. (2024), Investigating the Dynamic History of a
% Promontory Ice Rise using Radar Data, Journal of Glaciology
%
% Falk Oraschewski, 2024
close all
clear
addpath('fmcw_scripts');

%%% Settings
% Number of data points
Nsites = 15;
% Plot range
range_plot = 0:600;
% Define fit intervals
fitIntervals = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7] * 549;

% Path to data
dirDataHIR = '.\data';
% Subdirectory for output files
subdirOutput = 'processed';
% Filename for output file
fileOutput = 'StrainRatesHIR.mat';
% Subdirectory with ApRES data
subdirApRES = 'renamed';
% Note: To make the ApRES files identifiable, they were collected in one
% folder and renamed by adding "pN_" at the beginning of the filenames with
% the site number N.

% Name of file with study site information
fileStudySite = 'ApRESposReordered.csv';
% Note:T The sites in 'ApRESposReordered.csv' were reordered to match the
% following order of ApRES sites on the HIR profile
siteOrder = [7:-1:0, 8:1:14];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make directory for output data
mkdir(fullfile(dirDataHIR, 'processed'))

% Load study site information
dataHIR = readtable(fullfile(dirDataHIR, fileStudySite));
lon = dataHIR.X;
lat = dataHIR.Y;
% Convert coordinates to polar stereographic
proj = projcrs(3031, 'Authority', 'EPSG');
[psX, psY] = projfwd(proj, lat, lon);
% Compute distance along profile
dist = sqrt((psX-psX(1)).^2 + (psY-psY(1)).^2)/1000;

% Create empty matrizes
vsr_pts = zeros(length(range_plot),Nsites);
vsr_pts_mm = zeros(length(range_plot),Nsites);
vsr_fit = zeros(6,Nsites);
vsre_fit = zeros(6,Nsites);
vsr_mean = zeros(6,Nsites);
vsre_mean = zeros(6,Nsites);
vsre_mean_e = zeros(6,Nsites);
vsr_mean_mm = zeros(6,Nsites);
vsre_mean_mm = zeros(6,Nsites);
vsre_mean_e_mm = zeros(6,Nsites);
vsr_mean_interp = zeros(6,Nsites);
vsre_mean_interp = zeros(6,Nsites);

% Create counter
nCount = 0;
% Iterate over ApRES sites
for n = siteOrder
    nCount = nCount+1;
    % Get files
    files = dir(fullfile(dirDataHIR, subdirApRES, ['p' num2str(n) '_*.dat']));
    file1 = fullfile(files(1).folder, files(1).name);
    file2 = fullfile(files(2).folder, files(2).name);

    % Compute deformation
    site = fmcw_deformation(file1, file2);
    
    % Load strain rates with moving mean filter
    mm_range = site.sr_pt.range_gn_mm;
    mm_vsr = site.sr_pt.vsr_mm;
    mm_vsre = site.sr_pt.vsr_mm_std;
    vsr_pt_mm = interp1(mm_range,mm_vsr,range_plot,'nearest',NaN);
    vsr_pts_mm(:,nCount) = vsr_pt_mm;
    
    % Load strain rates between strong reflections
    gn_range = site.sr_pt.range_gn;
    gn_vsr = site.sr_pt.vsr;
    gn_vsre = site.sr_pt.vsre;
    vsr_pt = interp1(gn_range,gn_vsr,range_plot,'nearest',NaN);
    vsr_pts(:,nCount) = vsr_pt;

    % Load strain rate fits
    vsr_fit(:,nCount) = site.sr_fit.vsr;
    vsre_fit(:,nCount) = site.sr_fit.vsre;

    % Split data into depth intervals
    for m = 1:6
        vsr_mean(m,nCount) = mean(gn_vsr(gn_range > fitIntervals(m) & gn_range <= fitIntervals(m+1)));
        vsre_mean(m,nCount) = std(gn_vsr(gn_range > fitIntervals(m) & gn_range <= fitIntervals(m+1)));
        vsre_mean_e(m,nCount) = sqrt(sum(gn_vsre(gn_range > fitIntervals(m) & gn_range <= fitIntervals(m+1)).^2));
        vsr_mean_mm(m,nCount) = mean(mm_vsr(mm_range > fitIntervals(m) & mm_range <= fitIntervals(m+1)));
        vsre_mean_mm(m,nCount) = std(mm_vsr(mm_range > fitIntervals(m) & mm_range <= fitIntervals(m+1)));
        vsre_mean_e_mm(m,nCount) = sqrt(sum(mm_vsre(mm_range > fitIntervals(m) & mm_range <= fitIntervals(m+1)).^2));
        vsr_mean_interp(m,nCount) = mean(vsr_pt(range_plot > fitIntervals(m) & range_plot <= fitIntervals(m+1)));
        vsre_mean_interp(m,nCount) = std(vsr_pt(range_plot > fitIntervals(m) & range_plot <= fitIntervals(m+1)));
    end

    % Collect data into data structure
    data.dist = dist(nCount);
    data.psX = psX(nCount);
    data.psY = psY(nCount);
    data.depth = gn_range;
    data.vsr = gn_vsr;
    data.vsre = gn_vsre;
    data.depth_movemean = mm_range;
    data.vsr_movemean = mm_vsr;
    data.vsre_movemean = mm_vsre;
    data.vsr_fit = vsr_fit;
    data.vsre_fit = vsre_fit;
    data.vsr_mean = vsr_mean;
    data.vsre_mean = vsre_mean;
    data.vsre_mean_e = vsre_mean_e;
    data.vsr_mean_interp = vsr_mean_interp;
    data.vsre_mean_interp = vsre_mean_interp;
    
    % Save data for ApRES site
    %file_out = fullfile(dirDataHIR, subdirOutput, ['p' num2str(n) '_vsr.mat']);
    %save(file_out, 'data', '-mat')
end

% Save output data for all ApRES sites
data.fitIntervals = fitIntervals;
% save(fullfile(dirDataHIR, subdirOutput, fileOutput),'data','-mat')

%% Prepare plots

% Process distance information
dist_plot = dist-dist(8); % Distance relative to ice-divide
dist_plot = dist_plot/0.456; % Normalize with thickness
dist_mid = [dist_plot(1)-diff(dist_plot(1:2))/2; dist_plot(2:end)-diff(dist_plot)/2; dist_plot(end)+diff(dist_plot(end-1:end))/2];
[X,Y] = meshgrid(dist_mid,range_plot);

% Adapt strain rate data to meshgrid
vsr_pts_plot = vsr_pts;
vsr_pts_plot(:,end+1) = vsr_pts_plot(:,end);
vsr_pts_plot_mm = vsr_pts_mm;
vsr_pts_plot_mm(:,end+1) = vsr_pts_plot_mm(:,end);

%% Plot strain rates (pointwise between strong reflections)
f = figure;
f.Position = [100 100 900 700];
clear ax
ax(1) = subplot(3,1,1);
errorbar(dist_plot,vsr_fit',vsre_fit',-vsre_fit')
hold on
xticks([-4.5 -3 -2 -1 0 1 2 3 4.5])
ylabel('$\dot\varepsilon_{zz}\ \left(\mathrm{yr}^{-1}\right)$','Interpreter','latex')
box on
title('Vertical strain rates along ApRES line','Interpreter','latex')
leg = legend('1-0.8','0.8-0.7','0.7-0.6','0.6-0.5','0.5-0.4','0.4-0.3','Location','southwest','NumColumns',2);
leg.ItemTokenSize = [15,9];

ax(2) = subplot(3,1,2:3);
h = pcolor(X,Y,vsr_pts_plot);
hold on
colormap(parula);
set(h, 'EdgeColor', 'none')
clim([-0.007,0])
ch = colorbar('Southoutside');
ch.Label.String = '$\dot\varepsilon_{zz}\ \left(\mathrm{yr}^{-1}\right)$';
ch.Label.Interpreter = 'latex';
axis ij
xticks([-4.5 -3 -2 -1 0 1 2 3 4.5])
xlabel('Distance along line $(\mathrm{km})$','Interpreter','latex')
ylabel('Depth $(\mathrm{m})$','Interpreter','latex')
box on
linkaxes(ax,'x')

%% Plot smoothed strain rates (after moving mean filtering)

f = figure;
f.Position = [100 100 900 700];
clear ax
ax(1) = subplot(3,1,1);
errorbar(dist_plot,vsr_fit',vsre_fit',-vsre_fit')
hold on
xticks([-4.5 -3 -2 -1 0 1 2 3 4.5])
ylabel('$\dot\varepsilon_{zz}\ (\mathrm{yr}^{-1})$','Interpreter','latex')
box on
title('Vertical strain rates along ApRES line','Interpreter','latex')
leg = legend('1-0.8','0.8-0.7','0.7-0.6','0.6-0.5','0.5-0.4','0.4-0.3','Location','southwest','NumColumns',2);
leg.ItemTokenSize = [15,9];

ax(2) = subplot(3,1,2:3);
h = pcolor(X,Y,vsr_pts_plot_mm);
hold on
colormap(parula);
set(h, 'EdgeColor', 'none')
clim([-0.007,0])
ch = colorbar('Southoutside');
ch.Label.String = '$\dot\varepsilon_{zz}\ \left(\mathrm{yr}^{-1}\right)$';
ch.Label.Interpreter = 'latex';
axis ij
xticks([-4.5 -3 -2 -1 0 1 2 3 4.5])
xlabel('Distance along line $(\mathrm{km})$','Interpreter','latex')
ylabel('Depth $(\mathrm{m})$','Interpreter','latex')
box on
linkaxes(ax,'x')

