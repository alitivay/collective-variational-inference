function [ lpy ] = log_py( z, data, lognoise )
            
    % Run the appropriate likelihood function (according to experiment type)
    this_log_py = str2func(strcat(data.Experiment, '_log_likelihood'));
    lpy = this_log_py( z, data, lognoise );

end

