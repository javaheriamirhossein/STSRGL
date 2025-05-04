function [ params_struct, all_results ] = AVAR_learn_different_experiments( params_struct, experiments, savedir)
clc
close all

if nargin<2
    experiments = {'noise', 'sr'};
end


if nargin<3
    savedir = [cd, '\Results\VAR_Learning\', date, '\'];
end


if ~isfolder(savedir)
    mkdir(savedir);
end


if isfield(params_struct,'default') &&   params_struct.default ==1
    type = 1;       % graph type 
    N = 100;        % N nodes
    T = 10*N;       % N measurments
    sr = 0.8;       % sampling rate
    std_n = 0.1;    % noise level        
    AVAR_type = 'randl';
    Normalization = 'trace';
    W_thr = 5e-2;
    Ntrials = 50;        
    
    params_struct.type = type;
    params_struct.N = N;
    params_struct.T = T;
    params_struct.sr = sr;
    params_struct.std_n = std_n;    
    params_struct.AVAR_type = AVAR_type;
    params_struct.Normalization = Normalization;
    params_struct.W_thr = W_thr;
    params_struct.epsilon = epsilon;
    params_struct.Ntrials = Ntrials;
    params_struct.sampling_fixed = 1;
    params_struct.noise_fixed = 1;
end


type = params_struct.type;
N = params_struct.N;
T = params_struct.T;
sr = params_struct.sr;
std_n = params_struct.std_n;
AVAR_type = params_struct.AVAR_type;
Normalization = params_struct.Normalization; 
W_thr = params_struct.W_thr;
noise_type = params_struct.noise_type;
order = params_struct.order;


task = 'VAR learning';

Nexp = length(experiments);
all_results = cell(1,Nexp);


code = [num2str(params_struct.sampling_fixed), num2str(params_struct.noise_fixed)];
 
W_thr_str = num2str(W_thr);
W_thr_str = W_thr_str(3:end);

description = ['type ', num2str(type), ' ', num2str(N), ' ', num2str(T), ' ', num2str(10*sr), ' ', num2str(10*std_n), ' ', AVAR_type, ...
                ' ', Normalization, ' ', W_thr_str, ' ', noise_type, ' ', num2str(order), ' ', code];
 


x_axis_log = 0; 
for n=1:Nexp
    
    experiment = experiments{n};

    switch experiment
        case 'noise'
            input_vec = linspace(0,1,11);
            labelx = 'Noise Variance';
            param_type = 'noise';
            [Results, title_cell] = AVAR_learn_different_params(  params_struct, param_type, input_vec );            x_axis_log = 0;
            
        case 'sr'
            input_vec = linspace(0.2,1,10);
            labelx = 'Sampling Rate';        
            param_type = 'sr';
            [Results, title_cell] = AVAR_learn_different_params(  params_struct, param_type, input_vec );

        case 'alpha_0'
            input_vec = logspace(-1,3,10);
            labelx = '$\alpha_0$';        
            param_type = 'alpha_0';
            [Results, title_cell] = AVAR_learn_different_params(  params_struct, param_type, input_vec );
            x_axis_log = 1;


        case 'tau'
            input_vec = logspace(0,3,10);
            labelx = '$\tau$';        
            param_type = 'tau';
            [Results, title_cell] = AVAR_learn_different_params(  params_struct, param_type, input_vec );
            x_axis_log = 1;


        case 'alpha_1'
            input_vec = logspace(-1,2,10);
            labelx = '$\alpha_1$';        
            param_type = 'alpha_1';
            [Results, title_cell] = AVAR_learn_different_params(  params_struct, param_type, input_vec );            x_axis_log = 1;
            x_axis_log = 1;


    end


    [Results_mat_cell, field_names] = Results_nestedcell2tensor( Results );

    titleStr = [task, ' different ', experiment];
    legend_cell = title_cell;

    N_figs = length(Results_mat_cell);
    ylabls = {'Relative Error', 'F-score', 'time (s)', 'AUC', 'Number of Edges'};
    
    dim = params_struct.Ntrials>1;
    for k = 1:N_figs
        if length(size(Results_mat_cell{k}))>(2+dim) 
            Result_tensor = mean(Results_mat_cell{k},4);
            order = size(Results_mat_cell{k},3);
            for p = 1:order

                Result = Result_tensor(:,:,p);
                figure
                figname = field_names{k};
                if x_axis_log
                    semilogx(input_vec, Result, '--s', 'LineWidth',1.2)
                else
                    plot(input_vec, Result, '--s', 'LineWidth',1.2)
                end
                
                xlabel(labelx, 'interpreter', 'latex','FontSize',13)
                ylabel(ylabls{k},'interpreter','latex','FontSize',13)
                legend(legend_cell, 'interpreter', 'latex','location','northeast')  % 'northeastoutside'
        
                xlim([input_vec(1), input_vec(end)*1.1])
                % axis([-inf inf -inf inf])
        
                set(gcf, 'Position', [100 100 500 400])
                titleStr_p = [titleStr, 'order ', num2str(p), ''];
                title(titleStr_p , 'interpreter', 'latex','FontSize',15)
                grid on
                
                
        
                figsavename = [savedir, 'Figure AL order ', num2str(p),' ', description, ' diff ', experiment, ' ', figname, '.fig'];                 
                savefig(gcf,figsavename);
                figsavename = [savedir, 'Figure AL order ', num2str(p),' ', description, ' diff ', experiment, ' ', figname, '.png'];
                saveas(gcf,figsavename);
                close(gcf)

            end

        else
            Result = mean(Results_mat_cell{k},3);
            figure
            figname = field_names{k};
            if x_axis_log
                semilogx(input_vec, Result, '--s', 'LineWidth',1.2)
            else
                plot(input_vec, Result, '--s', 'LineWidth',1.2)
            end
            
            xlabel(labelx, 'interpreter', 'latex','FontSize',13)
            ylabel(ylabls{k},'interpreter','latex','FontSize',13)
            legend(legend_cell, 'interpreter', 'latex','location','northeast') 
    
            xlim([input_vec(1), input_vec(end)*1.1])
            % axis([-inf inf -inf inf])
    
            set(gcf, 'Position', [100 100 500 400])    
            title(titleStr, 'interpreter', 'latex','FontSize',15)
            grid on
            
            
    
            figsavename = [savedir, 'Figure AL ', description, ' diff ', experiment, ' ', figname, '.fig'];                 
            savefig(gcf,figsavename);
            figsavename = [savedir, 'Figure AL ', description, ' diff ', experiment, ' ', figname, '.png'];
            saveas(gcf,figsavename);
            close(gcf)

        end

    end
    
    all_results{n} = Results_mat_cell;
    
end

matsavename = [savedir, 'Results all AL ', description, '.mat']; 
save(matsavename, 'all_results', 'experiments') 


end

