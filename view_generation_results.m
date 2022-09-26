
clear all;
clc;

% Settings
sidx = 4; % Select input to give the virtual subjects for visualization
gen_range_hi = 5; % How many standard deviations to randomize for subject generation
ci_res = 1000; % How many subjects to generate
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
P = model_variational(VARS, var_sizes, zeros(n_subjects,n_params));

% Simulate the selected subject with most-likely personalized parameters
Outputs = HR_run_model(Data.Inputs, P(sidx,:));

% Start parallel cluster for higher speed
if isempty(gcp)
    this_cluster = parcluster('local');
    this_cluster.NumWorkers = 6;
    this_pool = parpool(this_cluster, this_cluster.NumWorkers);
end

% Plot Output Confidence Intervals ---------------------------------------

% Get simulated signal length for different measured variables
hc_length = length(Outputs.HCT.Times);
co_length = length(Outputs.CO.Times);
bp_length = length(Outputs.MAP.Times);

% Initialize variables to store sampled responses
S_HC = zeros(ci_res, hc_length);
S_CO = zeros(ci_res, co_length);
S_BP = zeros(ci_res, bp_length);

% Generate 2-sigma random numbers to be later used for subject generation
RANDS = zeros(ci_res, n_params);
i = 1;
while i <= ci_res
    random_no = randn(1, n_params);
    if (all(abs(random_no) < gen_range_hi))
        RANDS(i,:) = random_no;
        i = i+1;
    end
end

% Simulate the generated subjects
parfor i = 1:ci_res
    
    [P_sample] = model_generative( VARS, var_sizes, RANDS(i,:));
    SampleOutputs = HR_run_model(Data.Inputs, P_sample);
    
    S_HC(i,:) = SampleOutputs.HCT.Values;
    S_CO(i,:) = SampleOutputs.CO.Values;
    S_BP(i,:) = SampleOutputs.MAP.Values;
    
    fprintf('Sample: %d\n', i);
    
end

P_HC = zeros(size(P,1), hc_length);
P_CO = zeros(size(P,1), co_length);
P_BP = zeros(size(P,1), bp_length);

% Simulate the personalized subject
parfor i = 1:size(P,1)

    SampleOutputs = HR_run_model(Data.Inputs, P(i,:));
    
    P_HC(i,:) = SampleOutputs.HCT.Values;
    P_CO(i,:) = SampleOutputs.CO.Values;
    P_BP(i,:) = SampleOutputs.MAP.Values;
    
    fprintf('Subject: %d\n', i);
    
end

%% -----------------------------------------------------------------------

figure; set(gcf, 'Position', [100, 100, 250, 700]);

% Plot Inputs ------------------------------------------------------------

subplot(4,1,1);
plot(Data.Inputs.Infusion.Times, Data.Inputs.Infusion.Values,'-k','LineWidth',1); hold on;
plot(Data.Inputs.Hemorrhage.Times, -Data.Inputs.Hemorrhage.Values,'-.k','LineWidth',1); hold on;
plot(Data.Inputs.UO.Times,-Data.Inputs.UO.Values,'--k','LineWidth',1);
set(gca,'XTick',[0:30:180]); xlim([0 180]);
xlabel('Time (min)'); ylabel('Fluid I/O (mL/min)');
grid on;

% Plot Virtual Subject Responses -----------------------------------------

subplot(4,1,2);
hold on;
for i = 1:ci_res
    p = plot(Outputs.HCT.Times, S_HC(i,:)*100, '-b'); p.Color(4) = 0.2;
end
for i = 1:size(P,1)
    plot(Outputs.HCT.Times, P_HC(i,:)*100, '-r', 'LineWidth', 1);
end
set(gca,'XTick',[0:30:180]);
ylabel('Hematocrit (%)'); xlabel('Time (min)');
set(gca,'fontsize',10); set(gca,'XLim',[0 180],'XTick',[0:30:180])
axis on; grid on; box on;

subplot(4,1,3);
hold on;
for i = 1:ci_res
    p = plot(Outputs.CO.Times, S_CO(i,:), '-b');
    p.Color(4) = 0.2;
end
for i = 1:size(P,1)
    plot(Outputs.HCT.Times, P_CO(i,:), '-r', 'LineWidth', 1);
end
set(gca,'XTick',[0:30:180]); ylim([0 9]);
ylabel('Cardiac Output (L/min)'); xlabel('Time (min)');
set(gca,'fontsize',10); set(gca,'XLim',[0 180],'XTick',[0:30:180])
axis on; grid on; box on;

subplot(4,1,4);
hold on;
for i = 1:ci_res
    p = plot(Outputs.MAP.Times, S_BP(i,:), '-b');
    p.Color(4) = 0.2;
end
for i = 1:size(P,1)
    plot(Outputs.HCT.Times, P_BP(i,:), '-r', 'LineWidth', 1);
end
set(gca,'XTick',[0:30:180]); ylim([0 130]);
ylabel('MAP (mmHg)'); xlabel('Time (min)');
set(gca,'fontsize',10); set(gca,'XLim',[0 180],'XTick',[0:30:180]);
axis on; grid on; box on;
