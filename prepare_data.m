
function [ DATASET ] = prepare_data()

    j = 0;           % Experiment index
    DATASET = {};    % Container for the dataset

    % 0. Load data 
    raw_data = load('SDATA/synthetic_hemorrhage_resuscitation_dataset');
    synthetic_data = raw_data.GEN_DATA;
    
    for SN = 1:numel(synthetic_data)

        % Increment experiment index
        j = j+1;

        % 1. Name the experiment type. For example, these data belong to
        % hemorrhage resuscitation experiments, so we choose 'HR' as the name.
        Experiment = 'HR';

        % 2. Load input data into an appropriate data structure.
        Inputs.Infusion   = struct('Values', synthetic_data{SN}.Inputs.Infusion.Values, 'Times', synthetic_data{SN}.Inputs.Infusion.Times);
        Inputs.Hemorrhage = struct('Values', synthetic_data{SN}.Inputs.Hemorrhage.Values, 'Times', synthetic_data{SN}.Inputs.Hemorrhage.Times);
        Inputs.UO         = struct('Values', synthetic_data{SN}.Inputs.UO.Values, 'Times', synthetic_data{SN}.Inputs.UO.Times);

        % 3. Load output data into an appropriate data structure.
        Measurements.HCT = struct('Values', synthetic_data{SN}.Measurements.HCT.Values, 'Times', synthetic_data{SN}.Measurements.HCT.Times);
        Measurements.CO  = struct('Values', synthetic_data{SN}.Measurements.CO.Values, 'Times', synthetic_data{SN}.Measurements.CO.Times);
        Measurements.MAP = struct('Values', synthetic_data{SN}.Measurements.MAP.Values, 'Times', synthetic_data{SN}.Measurements.MAP.Times);

        % 4. Save the data structure into a cell array.
        DATASET{j} = struct('Experiment', Experiment, 'Inputs', Inputs, 'Measurements', Measurements);

    end

end
