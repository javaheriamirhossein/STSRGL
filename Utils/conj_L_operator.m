function [ v ] = conj_L_operator( A, row_col_id )
% Conjugate of the Laplacian operator

N = size(A,1);


if nargin<2
    [~, lowerdiag_ind] = diag_lowerdiag_index(N);
    Nw = size(lowerdiag_ind,1);    
    v = zeros(Nw,1);
    for k=1:Nw
        [i,j] = id_2_row_col( lowerdiag_ind(k,1), N );
        v(k) = A(i,i) + A(j,j) - A(i,j) - A(j,i);
    end
    
else
    
    Nw = size(row_col_id,1);    
    v = zeros(Nw,1);
    for k=1:Nw
        i = row_col_id(k,1);
        j = row_col_id(k,2);
        v(k) = A(i,i) + A(j,j) - A(i,j) - A(j,i);
    end
    
end

end
