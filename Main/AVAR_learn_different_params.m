function [Results, title_cell] = AVAR_learn_different_params(  params_struct, params_type, params_vec)
% Running VAR learning experiments at different scenarios

if nargin<2
    params_vec = linspace(0,1,10);
end

if isfield(params_struct, 'sr')               % sampling rate
    sr = params_struct.sr;
else
    sr = 0.8;
end

if isfield(params_struct, 'std_n')            % noise level
    std_n = params_struct.std_n;
else
    std_n = 0.1;
end


if isfield(params_struct, 'Ntrials')          % Number of trials
    Ntrials = params_struct.Ntrials;
else
    Ntrials = 20;
end



noise_fixed = 1;
sampling_fixed = 1;

switch(params_type)
    case 'noise'
        noise_fixed = 0;
        param_name = 'std_n';
    case 'sr'
        sampling_fixed = 0;
        param_name = 'sr';
    case 'alpha_0'
        param_name = 'alpha_0';
    case 'alpha_1'
        param_name = 'alpha_1';
    case 'tau'
        param_name = 'tau';
    otherwise
        warning('no valid param type')
end



%% Running the VAR learning algorithms for Ntrials

N_params = length(params_vec);


RelativeEr_cell_all = cell(N_params, Ntrials);
Fscore_cell_all = cell(N_params, Ntrials);
times_all = cell(N_params, Ntrials);
AUC_cell_all = cell(N_params, Ntrials);
NumEdge_cell_all = cell(N_params, Ntrials);

%--------
N = params_struct.N;
T = params_struct.T;



data_struct = synthetic_data( params_struct );


if sampling_fixed
    SampleNum = floor(N*sr);     % the number of sampled data at each time index
    SampleMatrix = zeros(N,T);
    
    for i = 1:T
        SampleMatrix(randperm(N, SampleNum),i) = 1;
    end
end

if noise_fixed
    noise = std_n * randn(N,T);  % measurement noise
end


%%
for n = 1:N_params
    
    param = params_vec(n);

%     if strcmp(param_name, "alpha_0")
%         param = T*param;
%     end

    if strcmp(param_name, "noise")
        std_n = param;
    end

    if strcmp(param_name, "sr")
        sr = param;
    end
    
    for trial = 1:Ntrials
            
        
        X = data_struct.X;
        
        
        if ~sampling_fixed
            SampleNum = floor(N*sr); 
            SampleMatrix = zeros(N,T);
            
            for i = 1:T
                SampleMatrix(randperm(N, SampleNum),i) = 1;
            end
        end
        
        if ~noise_fixed
            noise = std_n * randn(size(X));   
        end
            
        Y = SampleMatrix.*(X + noise);
        data_struct.SampleMatrix = SampleMatrix;
        params_struct.(param_name) = param;

        [Result, title_cell] = AVAR_learning_algorithms( Y, data_struct, params_struct);


    
        % Save the results
        RelativeEr_cell_all{n, trial} = Result.RelativeEr;
        Fscore_cell_all{n, trial} = Result.Fscore;
        times_all{n, trial} = Result.times;
        AUC_cell_all{n, trial} = Result.AUC;
        NumEdge_cell_all{n, trial} = Result.NumEdge_cell;
        
    end
    
end

Results.RelativeEr = RelativeEr_cell_all;
Results.Fscore = Fscore_cell_all;
Results.times = times_all;
Results.AUC = AUC_cell_all;
Results.NumEdge_cell = NumEdge_cell_all;

end

