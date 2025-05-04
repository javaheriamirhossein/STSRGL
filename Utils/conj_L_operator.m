function [ v ] = conj_L_operator_simple( A )
% the conjugate of the Laplacian operator

N = size(A,1);


Nw = fix(0.5*N*(N-1));    
v = zeros(Nw,1);

k = 0;
for j=1:N
    for i=1:N    
        if i>j
            k = k+1;
            v(k) = A(i,i) + A(j,j) - A(i,j) - A(j,i);
        end
    end
end 
    


end
