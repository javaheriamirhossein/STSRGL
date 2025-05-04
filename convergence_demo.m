% Demo for examining the convergene of the graph learning methods
%---------------------------------------------------

savedir = [cd, '\Results\Demo_Convergence\', date, '\'];
addpath(genpath(cd))


%% Synthetic data generation
%-------------------------------------------------------

rng(1);
type = 1;
N = 100;
T = 10*N;
sr = 1;         % sampling rate
std_n = 0.1;    % noise std
Normalization = 'trace';
AVAR_type = 'randl';
W_thr = 0.05;   % weight threshold
A_thr = 0.7;    % VAR matrix threshold


params_struct.N = N;
params_struct.T = T;
params_struct.type = type;
params_struct.AVAR_type = AVAR_type;
params_struct.mask_AVAR = 1;
params_struct.Normalization = Normalization;
params_struct.W_thr = W_thr;
params_struct.A_thr = A_thr;

data_struct = synthetic_data( params_struct );

Xtrue = data_struct.X;
Ltrue = data_struct.Ltrue;
Wtrue = -Ltrue;
Wtrue(1:N+1:end) = 0;
AVAR_true = data_struct.AVAR_true;
AVAR_true_cell = {AVAR_true, zeros(size(AVAR_true))};


% Missing samples and adding noise 
%-------------------------------------------------------

X = Xtrue;
SampleNum = floor(N*sr); % the number of sampled points at each time
SampleMatrix = zeros(N,T);
for i = 1:T
    SampleMatrix(randperm(N, SampleNum),i) = 1;
end 


noise = std_n * randn(size(X)); % measurement noise
Y = SampleMatrix.*(X + noise);
Mask = SampleMatrix;




%% Parameters of the algorithms
%-------------------------------------------------------

tau = 50;
std_n = 0.1;
W_thr = 0.01;
alpha_1 = 20;
alpha_0 = 0.1*T;

maxIter = 1000;
error_cell = cell(0);
title_cell = cell(0);
times = cell(0);


%% Running algorithms
%-------------------------------------------------------

params = struct;
params.std_n = std_n;
params.W_thr = W_thr;
params.maxIter = maxIter;
params.Normalization = Normalization;
params.Xtrue = Xtrue;
params.Atrue = AVAR_true;
params.Ltrue = Ltrue;
params.alpha_1 = alpha_1 ;
params.tau = tau ;
params.alpha_0 = alpha_0 ;
t0 = tic;
params.W_thr = W_thr;
output = learn_STSRGL( Y, Mask, params );
times{end+1} = toc(t0);
error_cell{end+1} = output;
title_cell{end+1} = 'STSRGL (Proposed)';


%% Error calculation
%-----------------------------
close all;
N_alg = length(error_cell);


description = ['type ', num2str(type), ' ', num2str(N), ' ', num2str(T), ' ', num2str(sr), ' ', num2str(std_n), ' ', AVAR_type, ...
                ' ', Normalization, ' ', num2str(W_thr), ' '];

            

if ~isfolder(savedir)
    mkdir(savedir)
end



Markers = {'h', '+','o','*' ,'x','s','p','d','^','v','>','<','.'};
LineTypes = { '-', '-.', '--', '--', '--', ':', '-.', '-'};
Colors = {[1 1 0], [1 0 0], [0 0 1], [0.466666666666667 0.674509803921569 0.188235294117647], [0 0 0], [0.870588235294118 0.490196078431373 0], ...
          [0.494117647058824 0.184313725490196 0.556862745098039], [0.294117647058824 0.184313725490196 0.556862745098039], ...
          [0.301960784313725 0.745098039215686 0.933333333333333], [0.929411764705882 0.694117647058824 0.125490196078431]};
 

label_x = 'time';

figure
grid on
hold on
legend_cell = cell(1);
cnt = 0;

for i = 1:N_alg
    if isfield( error_cell{i}, 'error_X') && isfield( error_cell{i}, 'times_X')
       cnt = cnt+1;
       y_data = error_cell{i}.error_X/ norm(Xtrue, 'fro');
       x_data = error_cell{i}.times_X; 
       marker = Markers{mod(cnt,numel(Markers))+1};
       color = Colors{mod(cnt,numel(Colors))+1};
       line = LineTypes{mod(cnt,numel(LineTypes))+1};
%        plot( x_data, y_data, 'Marker',marker,'Color',color,'LineStyle',line)
       plot( x_data, y_data,'--s','Color',color)
       legend_cell{cnt} = title_cell{i};
    end
end
title('X Error vs time', 'interpreter', 'latex','FontSize',15);
xlabel(label_x, 'interpreter', 'latex','FontSize',13);
ylabel('Relative Error', 'interpreter', 'latex','FontSize',13);
legend( legend_cell, 'interpreter', 'latex','location','northeast' );
figsavename = [savedir, 'Error X ', description, '.fig'];                 
savefig(gcf,figsavename);
       
 
figure
grid on
hold on
legend_cell = cell(1);
cnt = 0;

for i = 1:N_alg
    if isfield( error_cell{i}, 'error_A') && isfield( error_cell{i}, 'times_A')
       cnt = cnt+1;
       y_data = error_cell{i}.error_A/ norm(AVAR_true, 'fro');
       x_data = error_cell{i}.times_A;       
       marker = Markers{mod(cnt,numel(Markers))+1};
       color = Colors{mod(cnt,numel(Colors))+1};
       line = LineTypes{mod(cnt,numel(LineTypes))+1};
%        plot( x_data, y_data, 'Marker',marker,'Color',color,'LineStyle',line)
       plot( x_data, y_data,'--s','Color',color)
       legend_cell{cnt} = title_cell{i};
    end
end
title('A Error vs time', 'interpreter', 'latex','FontSize',15);
xlabel(label_x, 'interpreter', 'latex','FontSize',13);
ylabel('Relative Error', 'interpreter', 'latex','FontSize',13);
legend( legend_cell, 'interpreter', 'latex','location','northeast' );
figsavename = [savedir, 'Error A ', description, '.fig'];                 
savefig(gcf,figsavename);


figure
grid on
hold on
legend_cell = cell(1);
cnt = 0;

for i = 1:N_alg
    if isfield( error_cell{i}, 'Fscore_A') && isfield( error_cell{i}, 'times_A')
       cnt = cnt+1;
       y_data = error_cell{i}.Fscore_A;
       x_data = error_cell{i}.times_A;       
       marker = Markers{mod(cnt,numel(Markers))+1};
       color = Colors{mod(cnt,numel(Colors))+1};
       line = LineTypes{mod(cnt,numel(LineTypes))+1};
%        plot( x_data, y_data, 'Marker',marker,'Color',color,'LineStyle',line)
       plot( x_data, y_data,'--s','Color',color)
       legend_cell{cnt} = title_cell{i};
    end
end
title('A F-score vs time', 'interpreter', 'latex','FontSize',15);
xlabel(label_x, 'interpreter', 'latex','FontSize',13);
ylabel('F-score', 'interpreter', 'latex','FontSize',13);
legend( legend_cell, 'interpreter', 'latex','location','northeast' );
figsavename = [savedir, 'F-scroe A ', description, '.fig'];                 
savefig(gcf,figsavename);


figure
grid on
hold on
legend_cell = cell(1);
cnt = 0;

for i = 1:N_alg
    if isfield( error_cell{i}, 'error_L') && isfield( error_cell{i}, 'times_L')
       cnt = cnt+1;
       y_data = error_cell{i}.error_L/ norm(Ltrue, 'fro');
       x_data = error_cell{i}.times_L;       
       marker = Markers{mod(cnt,numel(Markers))+1};
       color = Colors{mod(cnt,numel(Colors))+1};
       line = LineTypes{mod(cnt,numel(LineTypes))+1};
%        plot( x_data, y_data, 'Marker',marker,'Color',color,'LineStyle',line)
       plot( x_data, y_data,'--s','Color',color)
       legend_cell{cnt} = title_cell{i};
    end
end
title('L Error vs time', 'interpreter', 'latex','FontSize',15);
xlabel(label_x, 'interpreter', 'latex','FontSize',13);
ylabel('Relative Error', 'interpreter', 'latex','FontSize',13);
legend( legend_cell, 'interpreter', 'latex','location','northeast' );
figsavename = [savedir, 'Error L ', description, '.fig'];                 
savefig(gcf,figsavename);


figure
grid on
hold on
legend_cell = cell(1);
cnt = 0;

for i = 1:N_alg
    if isfield( error_cell{i}, 'Fscore_L') && isfield( error_cell{i}, 'times_L')
       cnt = cnt+1;
       y_data = error_cell{i}.Fscore_L;
       x_data = error_cell{i}.times_L;       
       marker = Markers{mod(cnt,numel(Markers))+1};
       color = Colors{mod(cnt,numel(Colors))+1};
       line = LineTypes{mod(cnt,numel(LineTypes))+1};
%        plot( x_data, y_data, 'Marker',marker,'Color',color,'LineStyle',line)
       plot( x_data, y_data,'--s','Color',color)
       legend_cell{cnt} = title_cell{i};
    end
end
title('L F-score vs time', 'interpreter', 'latex','FontSize',15);
xlabel(label_x, 'interpreter', 'latex','FontSize',13);
ylabel('F-score', 'interpreter', 'latex','FontSize',13);
legend( legend_cell, 'interpreter', 'latex','location','northeast' );
figsavename = [savedir, 'F-scroe L ', description, '.fig'];                 
savefig(gcf,figsavename);

