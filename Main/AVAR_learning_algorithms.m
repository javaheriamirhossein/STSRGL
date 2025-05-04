function [Result, title_cell] = AVAR_learning_algorithms( Y, data_struct, params_struct )
% Evaluation of different (directed graph) VAR model learning algorithms


if isfield(params_struct, 'show_graph_results')              % show_graph_results
    show_graph_results = params_struct.show_graph_results;
else
    show_graph_results = 0;
end

if isfield(params_struct, 'VAR_Normalization')   % VAR matrix Normalization = 'norm', 'max','trace';
    VAR_Normalization = params_struct.VAR_Normalization;
else
    VAR_Normalization = 'norm';
end

if isfield(params_struct, 'Normalization')   % Laplacian Normalization = 'max','trace';
    Normalization = params_struct.Normalization;
else
    Normalization = 'trace';
end

if isfield(params_struct, 'W_thr')   % matrix threshold value
    W_thr = params_struct.W_thr;
else
    W_thr = 5e-2;
end

AVAR_true_cell = data_struct.AVAR_true_cell;
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


AVAR_cell = cell(1);
title_cell = cell(1);
times = cell(1);



%% VAR (directed graph) learning algorithms
%==============================

if nargin<5
    methods = {
                'STSRGL', ...
                'STSRGL A-sub', ...
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
    AVAR = output.A;
    times{ind} = toc(t0);
    AVAR_cell{ind} = AVAR;
    title_cell{ind} = 'STSRGL (Proposed)';
end

%-------------------------------------------------------
if any(strcmp(methods,'STSRGL A-sub'))
    ind = ind + 1;    
    t0 = tic;
    params = struct;
    params.std_n = std_n;
    params.tau = tau;
    params.alpha_0 = alpha_0;
    params.alpha_1 = alpha_1;
    params.Optimize_L = 0;
    params.Optimize_X = 0;
    params.Normalization = Normalization;
    output = learn_STSRGL( Y, Mask, params );
    AVAR = output.A;
    times{ind} = toc(t0);
    AVAR_cell{ind} = AVAR;
    title_cell{ind} = 'STSRGL A-sub (Proposed)';
end



%% Error calculation
%-------------------------------------------------------

N_alg = length(AVAR_cell);
order = length(AVAR_true_cell);

RelativeEr_cell = cell(1, order);
Fscore_cell = cell(1, order);
AUC_cell = cell(1, order);
NumEdge_cell = cell(1, order);

for p = 1:order
    AVAR_true = AVAR_true_cell{p};
    if strcmp(VAR_Normalization, 'norm')==1
        AVAR_true = AVAR_true/norm(AVAR_true);
    elseif strcmp(VAR_Normalization, 'trace')==1
        AVAR_true = AVAR_true/trace(AVAR_true)*N;
    elseif strcmp(VAR_Normalization, 'max')==1
        AVAR_true = AVAR_true/(max(max(AVAR_true)));
    end


    if show_graph_results
        close all
        AVAR = AVAR_true;
        AVAR = AVAR/max(max(abs(AVAR)));
        imagesc(abs(AVAR))
        colorbar
        colormap hot
        title('$A$ true', 'interpreter', 'latex')
    end


    RelativeEr = nan(1,N_alg);
    Fscore = nan(1,N_alg);
    AUC = nan(1,N_alg);
    NumEdge = nan(1,N_alg);

    AVAR_connectivity_true = AVAR_true~=0;

    for i = 1:N_alg
        AVAR = AVAR_cell{i};
        
        if iscell(AVAR) 
            if p <= length(AVAR)
                AVAR = AVAR{p};
            end
        end
        
    
        AVAR = AVAR/max(max(abs(AVAR)));        
        [~,~,~,AUC(i)] = perfcurve(AVAR_connectivity_true(:), abs(AVAR(:)),1);


        AVAR(abs(AVAR)<W_thr) = 0;        
        
        if strcmp(VAR_Normalization, 'norm')==1
            AVAR = AVAR/norm(AVAR);
        elseif strcmp(VAR_Normalization, 'trace')==1
            AVAR = AVAR/trace(AVAR)*N;
        elseif strcmp(VAR_Normalization, 'max')==1
            AVAR = AVAR/(max(max(AVAR)));
        end        
    
    
        Adjacency = AVAR~=0; 
        NumEdge(i) = sum(sum(Adjacency)) ;
        RelativeEr(i) = norm(AVAR-AVAR_true,'fro')/norm(AVAR_true,'fro');       
        tp = sum(sum(Adjacency.*AVAR_connectivity_true));
        fp = sum(sum(Adjacency.*(~AVAR_connectivity_true)));
        fn = sum(sum((~Adjacency).*AVAR_connectivity_true));
        Fscore(i) = 2*tp/(2*tp+fn+fp);
        
        if show_graph_results
            figure
            AVAR = AVAR/max(max(abs(AVAR)));
            imagesc(abs(AVAR))
            colormap hot
            colorbar 
            title(title_cell{i}, 'interpreter', 'latex')
        end


    end
    
    RelativeEr_cell{p} = RelativeEr;
    Fscore_cell{p} = Fscore ;
    AUC_cell{p} = AUC ;
    NumEdge_cell{p} = NumEdge ;



end


Result.times = cell2mat(times);
Result.RelativeEr_cell = RelativeEr_cell;
Result.Fscore_cell = Fscore_cell;
Result.AUC_cell = AUC_cell;
Result.NumEdge_cell = NumEdge_cell;

if order>1
    Result.RelativeEr = RelativeEr_cell;
    Result.Fscore = Fscore_cell;
    Result.AUC = AUC_cell;
    Result.NumEdge = NumEdge_cell;
else
    Result.RelativeEr = RelativeEr;
    Result.Fscore = Fscore;
    Result.AUC = AUC;
    Result.NumEdge = NumEdge;
end

end

