function [ id ] = row_col_2_id( i,j,N )
% Convert row and col to vector index
id = (j-1)*N+i;


end