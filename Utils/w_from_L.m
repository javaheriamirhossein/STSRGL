function [ w ] = w_from_L( L )
% Convert the Laplacian matrix to w (weight vector)

N = size(L,1);

W = -L;
W(1:N+1:end) = 0;
[ ~, lowerdiag_ind] = diag_lowerdiag_index( N );

w = W(lowerdiag_ind(:,1));
w = max( w,0 );
w = w(:);
end

