function [Result, title_cell] = X_recovery_algorithms( Y, data_struct, params_struct, methods)
% Evaluation of different signal recovery algorithms

if isfield(params_struct, 'show_graph_results')              % show_graph_results
    show_graph_results = params_struct.show_graph_results;
else
    show_graph_results = 0;
end

if isfield(params_struct, 'Normalization')   % Normalization = 'max','trace';
    Normalization = params_struct.Normalization;
else
    Normalization = 'trace';
end

if isfield(params_struct, 'W_thr')   % W threshold value
    W_thr = params_struct.W_thr;
else
    W_thr = 5e-2;
end

Xtrue = data_struct.X;

N = params_struct.N;


if isfield(data_struct, 'SampleMatrix')
    Mask = data_struct.SampleMatrix;
else
    Mask = ones(N);
end



if isfield(params_struct, 'std_n')
    std_n = params_struct.std_n;
else
    std_n = 0.1;
end

if isfield(params_struct, 'alpha_0')
    alpha_0 = params_struct.alpha_0;
else
    alpha_0 = 0;
end


if isfield(params_struct, 'alpha_1')
    alpha_1 = params_struct.alpha_1;
else
    alpha_1 = 20;
end

if isfield(params_struct, 'tau')
    tau = params_struct.tau;
else
    tau = 50;
end

X_cell = cell(1);
title_cell = cell(1);
times = cell(1);



%% Signal recovery algorithms
%==============================

if nargin<5
    methods = {
                'STSRGL', ...
                'STSRGL X-sub'                
              };
end

ind = 0;

%-------------------------------------------------------
if any(strcmp(methods,'STSRGL'))
    ind = ind + 1;    
    t0 = tic;
    params = struct;
    params.std_n = std_n;
    params.tau = tau;
    params.alpha_0 = alpha_0;
    params.alpha_1 = alpha_1;
    params.Normalization = Normalization;
    output = learn_STSRGL( Y, Mask, params );
    Xhat = output.X;
    times{ind} = toc(t0);
    X_cell{ind} = Xhat;
    title_cell{ind} = 'STSRGL (Proposed)';
end

%-------------------------------------------------------
if any(strcmp(methods,'STSRGL X-sub'))
    ind = ind + 1;    
    t0 = tic;
    params = struct;
    params.std_n = std_n;
    params.tau = tau;
    params.alpha_0 = alpha_0;
    params.alpha_1 = alpha_1;
    params.Optimize_A = 0;
    params.Optimize_L = 0;
    params.Normalization = Normalization;
    output = learn_STSRGL( Y, Mask, params );
    Xhat = output.X;
    times{ind} = toc(t0);
    X_cell{ind} = Xhat;
    title_cell{ind} = 'STSRGL X-sub (Proposed)';
end


%-------------------------------------------------------
%% Error calculation

N_alg = length(X_cell);
RelativeErs = zeros(1,N_alg);
SNRs = zeros(1,N_alg);


for i = 1:N_alg
    
    Xhat = X_cell{i};
        
    RelativeErs(i) = norm(Xhat-Xtrue,'fro')/norm(Xtrue,'fro');
    SNRs(i) = 20*log10(1/RelativeErs(i));
    RelativeErs(i) = NMSE(Xtrue,Xhat);
       
end

Result.RelativeEr = RelativeErs;
Result.SNRs = SNRs;
Result.times = cell2mat(times);

