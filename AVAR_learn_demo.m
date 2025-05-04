% Demo for examining the VAR learning performance
%---------------------------------------------------
maindir = [cd, '\Results\VAR_Learning\', date, '\'];
addpath(genpath(cd))

experiments = {'sr'};


%% Synthetic data generation
%-------------------------------------------------------
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

params_struct_orig = struct;
params_struct_orig.noise_type = noise_type;
params_struct_orig.order = order;
params_struct_orig.type = type;
params_struct_orig.N = N;
params_struct_orig.T = T;
params_struct_orig.sr = sr;
params_struct_orig.std_n = std_n;    
params_struct_orig.C_n_inv = C_n_inv;    
params_struct_orig.alpha_1 = alpha_1; 
params_struct_orig.alpha_0 = alpha_0; 
params_struct_orig.tau = tau; 
params_struct_orig.AVAR_type = AVAR_type;
params_struct_orig.Normalization = Normalization;
params_struct_orig.W_thr = W_thr;
params_struct_orig.A_thr = A_thr;
params_struct_orig.epsilon = epsilon;
params_struct_orig.Ntrials = Ntrials;
params_struct_orig.sampling_fixed = 1;
params_struct_orig.noise_fixed = 1;
params_struct_orig.mask_AVAR = 1;
params_struct_orig.AVAR_symmetric = 0;

%%
t0 = tic;
params_struct = params_struct_orig;
AVAR_learn_different_experiments( params_struct, experiments, maindir);



time_taken = toc(t0);
fprintf('Time taken = %f\n', time_taken);
