function [ W ] = W_from_w( w, N )
% Convert w (weight vector) to W (weight matrix)

n = length(w);
if nargin<2
    N = fix(0.5 * (1+ sqrt(1+8*n) ) );
end

W = zeros(N);
l = 1;
for i=1:N
    for j=(i+1):N
        W(i,j) = w(l);
        l = l+1;
    end
end
W = (W + W');

end

