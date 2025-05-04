% Demo for visualization of the learned weight matrix for undirected graph
% modeling
%---------------------------------------------------

savedir = [cd, '\Results\Graph_Learning\Visual\', date, '\'];
addpath(genpath(cd))
rng(1);

type = 1;       % graph type 
N = 100;        % N nodes
T = 10*N;       % N measurments
sr = 0.8;       % sampling rate
std_n = 0.1;    % noise level      
C_n_inv = 100*eye(N);    % noise inverse cov matrix
AVAR_type = 'randl';
Normalization = 'trace';
W_thr = 0.05;    % threshold on the weights 
A_thr = 0.7;     % VAR matrix threshold
epsilon = 1;
Ntrials = 1;  
alpha_0 = 0.1*T;
alpha_1 = 20;
tau = 50;
order = 1; 
noise_type = 'randn';

params_struct.show_graph = 0;
params_struct.N = N;
params_struct.T = T;
params_struct.type = type;
params_struct.AVAR_type = AVAR_type;
params_struct.mask_AVAR = 1;
params_struct.Normalization = Normalization;
params_struct.W_thr = W_thr;
params_struct.A_thr = A_thr;

data_struct = synthetic_data( params_struct );

Wtrue = data_struct.Wtrue;
Ltrue = data_struct.Ltrue;
X = data_struct.X;

[N,T] = size(X);
SampleNum = floor(N*sr); % the number of sampled points at each time
SampleMatrix = zeros(N,T);
for i = 1:T
    SampleMatrix(randperm(N, SampleNum),i) = 1;
end 


noise = std_n * randn(size(X)); % measurement noise
Y = SampleMatrix.*(X + noise);
 
Mask = SampleMatrix;



[N, T] = size(Y);
S_data = cov(Y',1);


%% Variables  


L_cell = cell(0);
title_cell = cell(0);
times = cell(0);


if strcmp(Normalization, 'trace')==1
    Wtrue = Wtrue/sum(sum(Wtrue))*N;
else
    Wtrue = Wtrue/abs(max(max(Wtrue)));
end
W_connectivity_true = Wtrue > 0;

%% Graph learning algorithms
%-------------------------------------------------------

params = struct;
params.std_n = std_n;
params.Normalization = Normalization;
t0 = tic;
output = learn_STSRGL( Y, Mask, params );
times{end+1} = toc(t0);
L_cell{end+1} = output.L;
title_cell{end+1} = 'STSRGL (Proposed)';


%-------------------------------------------------------
params = struct;
params.std_n = std_n;
params.Optimize_A = 0;
params.Optimize_X = 0;
params.Normalization = Normalization;
t0 = tic;
output = learn_STSRGL( Y, Mask, params );
times{end+1} = toc(t0);
L_cell{end+1} = output.L;
title_cell{end+1} = 'STSRGL L-sub (Proposed)';



%% Plot the learned graphs

close all
if ~isfolder(savedir)
    mkdir(savedir)
end

figure
Wtrue = Wtrue/max(max(abs(Wtrue)));
imagesc(abs(Wtrue))
colorbar 
colormap hot
title('$W$ true', 'interpreter', 'latex')
savefig(gcf,[savedir, 'Wtrue ',num2str(N),'.fig'])
% saveas(gcf,[savedir, 'Wtrue ',num2str(N),'.pdf'])


N_alg = length(L_cell);
RelativeEr = zeros(1,N_alg);
Fscore = zeros(1,N_alg);
AUC = zeros(1,N_alg);
X_roc_cell = cell(1,N_alg);
Y_roc_cell{i} = cell(1,N_alg);
T_roc_cell{i} = cell(1,N_alg);


for i = 1:N_alg
    Laplacian = L_cell{i};
     
   
    W = abs(Laplacian);
    W(1:N+1:end) = 0;
    W_norm = W/max(max(W));        
    [X_roc_cell{i},Y_roc_cell{i},T_roc_cell{i},AUC(i)] = perfcurve_graph(W_connectivity_true(:),W_norm(:));
    
    
    if strcmp(Normalization, 'trace')==1
        Laplacian = Laplacian/trace(Laplacian)*N;
    elseif strcmp(Normalization, 'max')==1
        Laplacian = Laplacian/max(max(Laplacian));
    end
    
    W = -Laplacian;
    W(1:N+1:end) = 0;
    W(W<W_thr) = 0;        
    Laplacian = diag(sum(W,2))-W;   
      
    Adjacency = W > 0; 
    RelativeEr(i) = norm(Laplacian-Ltrue,'fro')/norm(Ltrue,'fro');       
    Fscore(i) = Fscore_metric(W_connectivity_true, Adjacency);
    
    figure        
    W = W/max(max(abs(W)));
    imagesc(abs(W))
    colormap hot
    colorbar 
    title(title_cell{i}, 'interpreter', 'latex')
    savefig(gcf,[savedir, title_cell{i},' ',num2str(N),'.fig'])
    % saveas(gcf,[savedir, title_cell{i},' ',num2str(N),'.pdf'])
end



%% =====================================

Markers = {'h', '+','o','*' ,'x','s','p','d','^','v','>','<','.'};
% Colors = {'y', 'r', 'g', 'b', 'm', 'k', 'c'};
LineTypes = { '-', '-.', '--', '--', '--', ':', '-.', '-'};
Colors = {[1 1 0], [1 0 0], [0 0 1], [0.466666666666667 0.674509803921569 0.188235294117647], [0 0 0], [0.870588235294118 0.490196078431373 0], ...
          [0.494117647058824 0.184313725490196 0.556862745098039], [0.294117647058824 0.184313725490196 0.556862745098039], ...
          [0.301960784313725 0.745098039215686 0.933333333333333], [0.929411764705882 0.694117647058824 0.125490196078431]};
      
      
figure
grid on
hold on
for i=1:N_alg
    marker = Markers{mod(i,numel(Markers))+1};
    color = Colors{mod(i,numel(Colors))+1};
    line = LineTypes{mod(i,numel(LineTypes))+1};
    y_roc = Y_roc_cell{i};
    plot(X_roc_cell{i},y_roc,'Color',color,'LineStyle',line, 'LineWidth',2)
end

hold off
legend(title_cell,'location','best','interpreter','latex','FontSize',13)
xlabel('FPR','interpreter','latex')
ylabel('TPR','interpreter','latex')
title('ROC curve')
set(gcf, 'Position', [100 200 400 400])
axis equal
axis([0 1 0 1]) 
savefig(gcf,[savedir, 'AUC_performance.fig'])
% saveas(gcf,[savedir, 'AUC_performance.pdf'])
