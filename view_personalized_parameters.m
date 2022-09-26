
clear all;
clc; 

figure(234);

% Load inferred parameter information
load('RESULTS/CVI_RESULTS');

% Get parameter sizes
mat_size = size(RESULTS.TTA_MU);
N_subjects = mat_size(1);
N_params = mat_size(2);
sidxs = 1:N_subjects;

% Get uncertainty
e_PHI_SG = exp(RESULTS.PHI_SG);
e_TTA_SG = exp(RESULTS.TTA_SG);

% Get number of plots required
params_to_plot = 1:N_params;
n_plot = length(params_to_plot);

% Get covariance projections for pairwise parameter confidence intervals
ind_offdiag = find(tril(ones(N_params,N_params),-1));
L = zeros(N_params,N_params);
L(ind_offdiag) = RESULTS.PHI_COV;
ind_diag = find(eye(N_params));
L(ind_diag) = e_PHI_SG;
L_cut = L(params_to_plot,params_to_plot);

% Colors
color_ssp_hex = '#4E79A7' ;
color_gen_hex = '#E15759' ;
color_ssp = sscanf(color_ssp_hex(2:end),'%2x%2x%2x',[1 3])/255;
color_gen = sscanf(color_gen_hex(2:end),'%2x%2x%2x',[1 3])/255;

% Plot ALL Subject-Specific Params ---------------------------------------

for i = 1:n_plot
    for j = 1:n_plot
        
        subplot(n_plot, n_plot, (i-1)*n_plot+j);
        ppt = [params_to_plot(i) params_to_plot(j)];
        
        if ppt(1)==ppt(2)
            
            xlim([0,1]); ylim([0,10]);
            set(gca,'XTick',[]); set(gca,'YTick',[]); box on;
            
        else

            for sidx = sidxs 
                
                % Get confidence ellipse characteristics for subject
                S_all = diag(e_TTA_SG(sidx,:));
                this_center = RESULTS.TTA_MU(sidx,:)';
                this_COV = (S_all^2);
                
                % Plot confidence ellipse
                p = 0.95;
                s = -2 * log(1 - p);
                [V, D] = eig(this_COV(ppt,ppt) * s);
                t = linspace(0, 2 * pi);
                a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
                plot(a(1, :) + this_center(ppt(1)), a(2, :) + this_center(ppt(2)), 'Color', color_ssp);
                hold on;
                
            end
            
            % Get confidence ellipse characteristics for generator
            this_center = RESULTS.PHI_MU';
            this_COV = L*L';

            % Plot confidence ellipse
            p = 0.95;
            s = -2 * log(1 - p);
            [V, D] = eig(this_COV(ppt,ppt) * s);
            t = linspace(0, 2 * pi);
            a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
            plot(a(1, :) + this_center(ppt(1)), a(2, :) + this_center(ppt(2)), 'LineWidth', 1, 'Color', color_gen);

            set(gca,'XTick',[]); set(gca,'YTick',[]);
            xlim([0,1]); ylim([0,1]);
            
        end
        
    end
end

% Plot TWO selected parameters -------------------------------------------

i = 11;
j = 13;

param_names = {'alpha_Gain','alpha_Loss','Kp','Ka','Kv','a','b','Ktp','beta_v','Kco','BV0','HCT0','CO0','BP0'};

rlb = [ 0    0    0.01   10    1    0     5      0     0    0.001  1   0.1  3   70];
rub = [ 6    6    0.5    100   10   2     100    500   3    0.05   5   0.5  6   130];

fprintf('i: ');
fprintf(param_names{i});
fprintf(', j: ');
fprintf(param_names{j});
fprintf('\n');

figure;
set(gcf, 'Position', [100, 100, 400, 300]);
        
ppt = [params_to_plot(i) params_to_plot(j)];

for sidx = sidxs
    
    % Get confidence ellipse characteristics for subject
    S_all = diag(e_TTA_SG(sidx,:));
    this_center = RESULTS.TTA_MU(sidx,:)';
    this_COV = (S_all^2);
    
    % Plot confidence ellipse
    p = 0.95;
    s = -2 * log(1 - p);
    [V, D] = eig(this_COV(ppt,ppt) * s);
    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    plot_x = a(1, :) + this_center(ppt(1));
    plot_y = a(2, :) + this_center(ppt(2));
    plot_x_scaled = rlb(i)+(rub(i)-rlb(i))*plot_x;
    plot_y_scaled = rlb(j)+(rub(j)-rlb(j))*plot_y;
    plot(plot_x_scaled, plot_y_scaled, 'Color', color_ssp);
    hold on;

end

% Get confidence ellipse characteristics for generator
this_center = RESULTS.PHI_MU';
this_COV = L*L';

% Plot confidence ellipse
p = 0.95;
s = -2 * log(1 - p);
[V, D] = eig(this_COV(ppt,ppt) * s);
t = linspace(0, 2 * pi);
a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
plot_x = a(1, :) + this_center(ppt(1));
plot_y = a(2, :) + this_center(ppt(2));
plot_x_scaled = rlb(i)+(rub(i)-rlb(i))*plot_x;
plot_y_scaled = rlb(j)+(rub(j)-rlb(j))*plot_y;
plot(plot_x_scaled, plot_y_scaled, 'Color', color_gen);

xlabel(param_names{i}); ylabel(param_names{j});
