function [Results_mat_struct, field_names] = Results_nestedcell2tensor( Results )

if ~isstruct(Results)
    error('the input should be a struct')
end

field_names = fieldnames(Results);

N_fields = length(field_names);
Results_mat_struct = cell(1,N_fields);

for i=1:N_fields
    result_cell = Results.(field_names{i});
    [L, N_trial] = size(result_cell);
    if iscell(result_cell{1})
        N_alg = length(result_cell{1}{1});
        order = length(result_cell{1});
        result_mat = zeros(L, N_alg, order, N_trial);
        is_4d = 1;
    else
        N_alg = length(result_cell{1});
        result_mat = zeros(L, N_alg, N_trial);
        is_4d = 0;
    end
    
    if is_4d
        for p=1:order
            for r=1:L
                for c=1:N_trial
                    result_mat(r,:,p,c) = result_cell{r,c}{p};
                end
            end
        end
            
    else

        for r=1:L
            for c=1:N_trial
                result_mat(r,:,c) = result_cell{r,c};
            end
        end
    end

    Results_mat_struct{i} = result_mat;    

end

