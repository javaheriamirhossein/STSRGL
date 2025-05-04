function [ L ] = L_operator( w, N )
% The Laplacian operator: maps the vector of edge weights w to Laplacian matrix L

Nw = length(w);
if nargin<2
    N = fix( (1+ sqrt(1+8*Nw))/2 );
end

W = zeros(N);
cnt = 0;
for i=1:N
    for j=(i+1):N
        cnt = cnt+1;
        W(i,j) = w(cnt);
    end
end
        
W = (W + W');

L = diag(sum(W,2)) - W;

end

