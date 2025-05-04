function [ diag_ind, lowerdiag_ind] = diag_lowerdiag_index( N )
% Gives the diagonal and the lower-diagonal indices of a symmetric matrix
% of size NxN

diag_ind = zeros(N,1);
N_lower = 0.5*N*(N-1);
lowerdiag_ind = zeros(N_lower,2);

for i=1:N
    diag_ind(i) = row_col_2_id( i,i,N );
end

cnt = 0;

for j=1:N
    for i=1:N    
        if i>j
            cnt = cnt+1;
            lowerdiag_ind(cnt,:) = [row_col_2_id( i,j,N ), row_col_2_id( j,i,N )];
        end
    end
end 


end


