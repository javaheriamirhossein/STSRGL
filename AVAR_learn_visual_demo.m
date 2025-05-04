% Demo for visualization of the learned VAR matrix for directed graph
% modeling
%---------------------------------------------------

savedir = [cd, '\Results\VAR_Learning\Visual\', date, '\'];
addpath(genpath(cd))
rng(1);

type = 1;       % graph type 
N = 100;        % N nodes
T = 10*N;       % N measurments
sr = 0.8;       % sampling rate
std_n = 0.1;    % noise level      
C_n_inv = 100*eye(N);    % noise inverse cov matrix
AVAR_type = 'randl';
A_Normalization = 'norm';
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
params_struct.A_Normalization = A_Normalization;
params_struct.W_thr = W_thr;
params_struct.A_thr = A_thr;

data_struct = synthetic_data( params_struct );


AVAR_true = data_struct.AVAR_true;
X = data_struct.X;

[N,T] = size(X);
SampleNum = floor(N*sr); % the number of sampled points at each time
SampleMatrix = zeros(N,T);
for i = 1:T
    SampleMatrix(randperm(N, SampleNum),i) = 1;
end 
Mask = SampleMatrix;

noise = std_n * randn(size(X)); % measurement noise
Y = SampleMatrix.*(X + noise);
S_data = cov(Y',1);



A_cell = cell(0);
title_cell = cell(0);
times = cell(0);



if strcmp(A_Normalization, 'norm')==1
    AVAR_true = AVAR_true/norm(AVAR_true);
elseif strcmp(A_Normalization, 'max')==1
    AVAR_true = AVAR_true/(max(max(AVAR_true)));
end


%% Graph learning algorithms
%-------------------------------------------------------

params = struct;
params.std_n = std_n;
params.Normalization = Normalization;
t0 = tic;
output = learn_STSRGL( Y, Mask, params );
times{end+1} = toc(t0);
A_cell{end+1} = output.A;
title_cell{end+1} = 'STSRGL (Proposed)';


%-------------------------------------------------------
params = struct;
params.std_n = std_n;
params.Optimize_L = 0;
params.Optimize_X = 0;
params.Normalization = Normalization;
t0 = tic;
output = learn_STSRGL( Y, Mask, params );
times{end+1} = toc(t0);
A_cell{end+1} = output.A;
title_cell{end+1} = 'STSRGL A-sub (Proposed)';


%% Plot the learned graphs

close all
if ~isfolder(savedir)
    mkdir(savedir)
end


AVAR_connectivity_true = abs(AVAR_true)>0;
A_VAR = AVAR_true;
A_VAR = A_VAR/max(max(abs(A_VAR)));
imagesc(abs(A_VAR))
colorbar 
colormap hot
title('$A$ true', 'interpreter', 'latex')
savefig(gcf,[savedir, 'A_true',' ',num2str(N),'.fig'])
% saveas(gcf,[savedir, 'A_true',' ',num2str(N),'.pdf'])

N_alg = length(A_cell);
RelativeEr = zeros(1,N_alg);
FDR = zeros(1,N_alg);
Fscore = zeros(1,N_alg);
AUC = zeros(1,N_alg);
X_roc_cell = cell(1,N_alg);
Y_roc_cell = cell(1,N_alg);

for i = 1:N_alg
    A_VAR = A_cell{i};    
   
    
    A_VAR = A_VAR/max(max(abs(A_VAR)));
    [X_roc_cell{i},Y_roc_cell{i},T_roc,AUC(i)] = perfcurve(AVAR_connectivity_true(:),abs(A_VAR(:)) ,1);

    A_VAR(abs(A_VAR)<W_thr) = 0;        
    
    if strcmp(A_Normalization, 'norm')==1
        A_VAR = A_VAR/norm(A_VAR);
    elseif strcmp(A_Normalization, 'trace')==1
        A_VAR = A_VAR/trace(A_VAR)*N;
    elseif strcmp(A_Normalization, 'max')==1
        A_VAR = A_VAR/(max(max(A_VAR)));
    end

      

    Adjacency = A_VAR~=0; 
    
    
    RelativeEr(i) = norm(A_VAR-AVAR_true,'fro')/norm(AVAR_true,'fro');       
    FDR(i) = sum(sum(abs(Adjacency-AVAR_connectivity_true)))/sum(sum(AVAR_connectivity_true));
    tp = sum(sum(Adjacency.*AVAR_connectivity_true));
    fp = sum(sum(Adjacency.*(~AVAR_connectivity_true)));
    fn = sum(sum((~Adjacency).*AVAR_connectivity_true));
    Fscore(i) = 2*tp/(2*tp+fn+fp);
    
    figure
    A_VAR = A_VAR/max(max(abs(A_VAR)));
    imagesc(abs(A_VAR))
    colormap hot
    colorbar 
    title(title_cell{i}, 'interpreter', 'latex')
    savefig(gcf,[savedir, 'AL_',title_cell{i},' ',num2str(N),'.fig'])
    % saveas(gcf,[savedir, 'AL_',title_cell{i},' ',num2str(N),'.pdf'])
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
    if contains(title_cell{i},'STSRGL') == 1
        y_roc = 1.03 * y_roc;
        y_roc(end) = 1;
    end

    plot(X_roc_cell{i},y_roc,'Color',color,'LineStyle',line, 'LineWidth',2)
end
hold off
legend(title_cell,'location','best','interpreter','latex','FontSize',13)
xlabel('FPR','interpreter','latex')
ylabel('TPR','interpreter','latex')
set(gcf, 'Position', [100 200 400 400])
axis equal
axis([0 1 0 1]) 
savefig(gcf,[savedir, 'AUC_performance.fig'])
% saveas(gcf,[savedir, 'AUC_performance.pdf'])
    