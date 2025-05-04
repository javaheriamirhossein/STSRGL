function A_arma = generate_AVAR( N, A_arma_type )
% Generate random VAR evolution (transition) matrix
if nargin<2
    A_arma_type = 'randl';
end

switch(A_arma_type)
    case 'rand'
        A_arma = rand(N);
    case 'exp'
        A_arma = random('exp',1,N,N);
    case 'logn'
        A_arma = random('logn',1,1,N,N);
    case 'randl'
        A_arma = randl(N, N);
    case 'randn'
        A_arma = randn(N);
    case 'eye'
        A_arma = eye(N);
    otherwise
        error('unidentified A_arma tyoe')
end



end

