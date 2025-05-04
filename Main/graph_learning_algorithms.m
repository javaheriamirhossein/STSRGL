function [Result, title_cell] = graph_learning_algorithms( Y, data_struct, params_struct, methods)
% Evaluation of different (undirected) graph learning algorithms

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

Ltrue = data_struct.Ltrue;
Wtrue = data_struct.Wtrue;
W_connectivity_true = data_struct.Wconn_true;

N = params_struct.N;
S_cov = cov(Y',1);

if strcmp(Normalization, 'trace')==1
    Ltrue = Ltrue/trace(Ltrue)*N;
    Wtrue = Wtrue/trace(Ltrue)*N;
elseif strcmp(Normalization, 'max')==1
    Wtrue = Wtrue/abs(max(max(Wtrue)));
    Ltrue = diag(sum(Wtrue,2))-Wtrue;
end


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


L_cell = cell(1);
title_cell = cell(1);
times = cell(1);




%% (Undirected) Graph learning algorithms
%==============================

if nargin<5
    methods = {
                'STSRGL', ...
                'STSRGL L-sub', ...
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
    Laplacian = output.L;
    times{ind} = toc(t0);
    L_cell{ind} = Laplacian;
    title_cell{ind} = 'STSRGL (Proposed)';
end

%-------------------------------------------------------
if any(strcmp(methods,'STSRGL L-sub'))
    ind = ind + 1;    
    t0 = tic;
    params = struct;
    params.std_n = std_n;
    params.tau = tau;
    params.alpha_0 = alpha_0;
    params.alpha_1 = alpha_1;
    params.Optimize_A = 0;
    params.Optimize_X = 0;
    params.Normalization = Normalization;
    output = learn_STSRGL( Y, Mask, params );
    Laplacian = output.L;
    times{ind} = toc(t0);
    L_cell{ind} = Laplacian;
    title_cell{ind} = 'STSRGL L-sub (Proposed)';
end




%% Error calculation
%-------------------------------------------------------

if show_graph_results
    close all
    Laplacian = Ltrue;
    Laplacian(1:N+1:end) = 0;
    Laplacian = Laplacian/max(max(abs(Laplacian)));
    imagesc(abs(Laplacian))
    colorbar
    colormap hot
    title('$L$ true', 'interpreter', 'latex')
end

N_alg = length(L_cell);
RelativeEr = zeros(1,N_alg);
Fscore = zeros(1,N_alg);
AUC = zeros(1,N_alg);
NumEdge = nan(1,N_alg);


for i = 1:N_alg
    Laplacian = L_cell{i};
    
    W = abs(Laplacian);
    W(1:N+1:end) = 0;
    W_norm = W/max(max(W)); % Map W to [0,1] for ROC
    [~,~,~,AUC(i)] = perfcurve(W_connectivity_true(:),W_norm(:),1);

    

    % Scale the Laplacian
    if strcmp(Normalization, 'trace')==1
        Laplacian = Laplacian/trace(Laplacian)*N;
    elseif strcmp(Normalization, 'max')==1
        Laplacian = Laplacian/max(max(W));
    end
     
    % Threshold the weights
    W = abs(Laplacian);
    W(1:N+1:end) = 0; 
    W(W<W_thr) = 0;
    Laplacian = diag(sum(W,2))-W;

    Adjacency = Laplacian < 0;
    NumEdge(i) = sum(sum(W>0))/2 ;
    RelativeEr(i) = norm(Laplacian-Ltrue,'fro')/norm(Ltrue,'fro');
    tp = sum(sum(Adjacency.*W_connectivity_true));
    fp = sum(sum(Adjacency.*(~W_connectivity_true)));
    fn = sum(sum((~Adjacency).*W_connectivity_true));
    Fscore(i) = 2*tp/(2*tp+fn+fp);
    
    if show_graph_results
        figure
        Laplacian(1:N+1:end) = 0;
        Laplacian = Laplacian/max(max(abs(Laplacian)));
        imagesc(abs(Laplacian))
        colormap hot
        colorbar
        title(title_cell{i}, 'interpreter', 'latex')
    end
    
end

Result.AUC = AUC;
Result.RelativeEr = RelativeEr;
Result.Fscore = Fscore;
Result.times = cell2mat(times);
Result.NumEdge = NumEdge;

