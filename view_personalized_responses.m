
clear all;
clc;

sidx = 1; % Select Subject Number to View
ci_res = 100; % confidence interval resolution: No. of samples from the personalized approximate posterior
n_params = 14; % Number of model parameters

% Load inference results
load('RESULTS/CVI_RESULTS');
VARS = RESULTS.VARS;
var_sizes = RESULTS.var_sizes;

% Load dataset
DATASET = prepare_data;
Data = DATASET{sidx};
n_subjects = numel(DATASET);

% Get all most-likely personalized model parameters P from raw results stored in VARS
% Alternatively we can read these values from RESULTS.TTA_MU.
P = model_variational( VARS, var_sizes, zeros(n_subjects,n_params) );

% Simulate the selected subject with most-likely personalized parameters
Outputs = HR_run_model(Data.Inputs, P(sidx,:));

% Start parallel cluster for higher speed
if isempty(gcp)
    this_cluster = parcluster('local');
    this_cluster.NumWorkers = 6;
    this_pool = parpool(this_cluster, this_cluster.NumWorkers);
end

% Calculate Output Confidence Intervals ----------------------------------

% Get simulated signal length for different measured variables
hc_length = length(Outputs.HCT.Times);
co_length = length(Outputs.CO.Times);
bp_length = length(Outputs.MAP.Times);

% Initialize variables to store sampled responses
S_HC = zeros(ci_res, hc_length);
S_CO = zeros(ci_res, co_length);
S_BP = zeros(ci_res, bp_length);

% Get subject specific signal quality (noise) variables
e_NOISE_all = exp(RESULTS.NOISE);
e_NOISE = e_NOISE_all(sidx,:);

parfor i = 1:ci_res
    
    % Sample from the personalized approximate posterior
    P_sample = model_variational( VARS, var_sizes, randn(n_subjects,n_params) );
    SampleOutputs = HR_run_model(Data.Inputs, P_sample(sidx,:));
    
    S_HC(i,:) = SampleOutputs.HCT.Values;
    S_CO(i,:) = SampleOutputs.CO.Values;
    S_BP(i,:) = SampleOutputs.MAP.Values;
    
    S_HC_n(i,:) = SampleOutputs.HCT.Values + e_NOISE(1) * randn(hc_length,1);
    S_CO_n(i,:) = SampleOutputs.CO.Values  + e_NOISE(2) * randn(co_length,1);
    S_BP_n(i,:) = SampleOutputs.MAP.Values + e_NOISE(3) * randn(bp_length,1);
    
    fprintf('Sample: %d\n', i);
    
end

mu_HC = mean(S_HC);
mu_CO = mean(S_CO);
mu_BP = mean(S_BP);

sd_HC = std(S_HC);
sd_CO = std(S_CO);
sd_BP = std(S_BP);

sd_HC_n = std(S_HC_n);
sd_CO_n = std(S_CO_n);
sd_BP_n = std(S_BP_n);

%% ------------------------------------------------------------------------

disp_colors = { '-k', '-b' };
shade_colors = { [194 224 255]/255 };
outer_colors = { [113 140 227]/255 };

figure; set(gcf, 'Position', [100, 100, 250, 700]);

% Plot Inputs ------------------------------------------------------------

subplot(4,1,1);
plot(Data.Inputs.Infusion.Times, Data.Inputs.Infusion.Values,'-k','LineWidth',1); hold on;
plot(Data.Inputs.Hemorrhage.Times, -Data.Inputs.Hemorrhage.Values,'-.k','LineWidth',1); hold on;
plot(Data.Inputs.UO.Times,-Data.Inputs.UO.Values,'--k','LineWidth',1);
set(gca,'XTick',[0:30:180]); xlim([0 180]);
xlabel('Time (min)'); ylabel('Fluid I/O (mL/min)');
grid on;

% Plot 2-Sigma Confidence Intervals --------------------------------------

subplot(4,1,2);
hold on;
patch([Outputs.HCT.Times; flipud(Outputs.HCT.Times)], [(mu_HC-2*sd_HC)'*100; flipud((mu_HC+2*sd_HC)'*100)], ...
    shade_colors{1}, 'EdgeColor', outer_colors{1} );
plot(Outputs.HCT.Times, (mu_HC-2*sd_HC_n)'*100, 'Color', outer_colors{1});
plot(Outputs.HCT.Times, (mu_HC+2*sd_HC_n)'*100, 'Color', outer_colors{1});

subplot(4,1,3);
hold on;
patch([Outputs.CO.Times; flipud(Outputs.CO.Times)], [(mu_CO-2*sd_CO)'; flipud((mu_CO+2*sd_CO)')], ...
    shade_colors{1}, 'EdgeColor', outer_colors{1} );
plot(Outputs.CO.Times, (mu_CO-2*sd_CO_n)', 'Color', outer_colors{1});
plot(Outputs.CO.Times, (mu_CO+2*sd_CO_n)', 'Color', outer_colors{1});

subplot(4,1,4);
hold on;
patch([Outputs.MAP.Times; flipud(Outputs.MAP.Times)], [(mu_BP-2*sd_BP)'; flipud((mu_BP+2*sd_BP)')], ...
    shade_colors{1}, 'EdgeColor', outer_colors{1} );
plot(Outputs.MAP.Times, (mu_BP-2*sd_BP_n)', 'Color', outer_colors{1});
plot(Outputs.MAP.Times, (mu_BP+2*sd_BP_n)', 'Color', outer_colors{1});

% Plot Most-Likely Personalized Model ------------------------------------

subplot(4,1,2);
hold on;
plot(Outputs.HCT.Times,Outputs.HCT.Values*100, disp_colors{1}, 'LineWidth',1);
plot(Data.Measurements.HCT.Times, Data.Measurements.HCT.Values*100, 'o','MarkerSize',4,'MarkerEdgeColor','k');
set(gca,'XTick',[0:30:180]);
ylabel('Hematocrit (%)'); xlabel('Time (min)');
set(gca,'fontsize',10); set(gca,'XLim',[0 180],'XTick',[0:30:180])
axis on; grid on; box on;

subplot(4,1,3);
hold on;
plot(Outputs.CO.Times, Outputs.CO.Values, disp_colors{1}, 'LineWidth',1);
plot(Data.Measurements.CO.Times, Data.Measurements.CO.Values, 'o','MarkerSize',4,'MarkerEdgeColor','k');
set(gca,'XTick',[0:30:180]); ylim([0 7]);
ylabel('Cardiac Output (L/min)'); xlabel('Time (min)');
set(gca,'fontsize',10); set(gca,'XLim',[0 180],'XTick',[0:30:180])
axis on; grid on; box on;

subplot(4,1,4);
hold on;
plot(Outputs.MAP.Times, Outputs.MAP.Values, disp_colors{1}, 'LineWidth',1);
plot(Data.Measurements.MAP.Times, Data.Measurements.MAP.Values, 'o','MarkerSize',4,'MarkerEdgeColor','k');
set(gca,'XTick',[0:30:180]); ylim([0 130]);
ylabel('MAP (mmHg)'); xlabel('Time (min)');
set(gca,'fontsize',10); set(gca,'XLim',[0 180],'XTick',[0:30:180])
axis on; grid on; box on;
