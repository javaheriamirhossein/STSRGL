function [ v ] = conj_d_operator( a )
% the conjugate of the degree operator

N = numel(a);
Nw = fix(0.5*N*(N-1));    
v = zeros(Nw,1);

k = 0;
for j=1:N
    for i=1:N    
        if i>j
            k = k+1;
            v(k) = a(i) + a(j);
        end
    end
end 



end
