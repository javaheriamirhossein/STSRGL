function [Results_mat_struct, field_names] = Results_cell2tensor( Results )

if ~isstruct(Results)
    error('the input should be a struct')
end

field_names = fieldnames(Results);

N_fields = length(field_names);
Results_mat_struct = cell(1,N_fields);

for i=1:N_fields
    result_cell = Results.(field_names{i});
    [L, N_trial] = size(result_cell);
    N_alg = numel(result_cell{1});
    
    result_mat = zeros(L, N_alg, N_trial);
    
    for r=1:L
        for c=1:N_trial
            result_mat(r,:,c) = result_cell{r,c};
        end
    end
            
    Results_mat_struct{i} = result_mat;    

end

