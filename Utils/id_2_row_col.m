function [i,j] = id_2_row_col ( id,N )
% convert the vector index to matrix row col 
i = mod(id,N);

if i==0
    i = N;
    j = floor(id/N);    
else
    j = floor(id/N) + 1;

end