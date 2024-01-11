function X = generate_random_matrix( Nr, Nc, params )

if nargin<3
    params = struct;
end

if isfield(params, 'type')
    type = params.type;
else
    type = 'randn';
end

if isfield(params, 'mu')
    mu = params.mu;
else
    mu = 1;
end

if isfield(params, 'sigma')
    sigma = params.sigma;
else
    sigma = 1;
end

if isfield(params, 'nu')
    nu = params.nu;
else
    nu = 2.5;
end


switch(type)
    case 'rand'
        X = rand(Nr,Nc);
    case 'exp'
        X = random('exp',mu,Nr,Nc);
    case 'logn'
        X = random('logn',mu,sigma,Nr,Nc);
    case 'randl'
        X = randl(Nr, Nc);
    case 'randn'
        X = randn(Nr,Nc);
    case 't'
        X = random('t',nu,Nr,Nc);
    otherwise
        error('unidentified  type')
end

end

