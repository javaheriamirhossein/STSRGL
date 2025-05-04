function  output  = learn_STSRGL( Y, Mask, params )

% learn_STSRGL: Learn graph from noisy and incomplete
%               observations 
%   Usage:  output = learn_STSRGL(Y, SampleMask)
%           output = learn_STSRGL(Y, SampleMask, params)
%
%   Inputs:
%         Y         : The noisy and missing data matrix 
%         Mask      : The sampling mask matrix
%         params    : Optional parameters
%
%   Outputs:
%         output    : Output struct
%         output.X  : The reconstructed data matrix
%         output.L  : The inferred Laplacian mtrix
%         output.A  : The inferred state transition matrix
%
%
%
%   params
%   ---------------------
%
%    Normalization   : Specify how the Laplacian or the state transition
%                      matrix should be normalized
%    alpha_0         : The regularization parameter for the Laplacian
%                      matrix (weight) norm
%    alpha_1         : The regularization parameter for the state
%                      transition matrix norm
%    maxit           : Maximum number of iterations. Default: 1000
%    W_thr           : Weight pruning threshold
%    tol             : Tolerance for stopping criterion. Defaul: 1e-5
%    tau             : parameter of the algorithm (controlling the speed of
%                      convergence)
%    C_n_inv         : Noise inverse covariance matrix
%    X0              : Initial value of X
%    L0              : Initial value of L
%    A0              : Initial value of A
%    Optimize_X      : Whether to run the X-subproblem (optimize over X) 
%    Optimize_L      : Whether to run the L-subproblem (optimize over L) 
%    Optimize_A      : Whether to run the A-subproblem (optimize over A) 
%    Xtrue           : Ground-truth X
%    Ltrue           : Ground-truth L
%    Atrue           : Ground-truth A 
%    Compute_f       : Whether to compute the value of the objective
%                      function


% Author: Amirhossein Javaheri
% Date: Dec 2023


%% ============================

[N, T] = size(Y);
S_cov = cov(Y',1);
S_inv = pinv(S_cov);



% Useful matrices
J = ones(N, N) / N;
Hoff = eye(N) - ones(N);



%============================
% Algorithm parameters

if nargin<3
    params = struct;
end


if isfield(params, 'Normalization')    % Normalization should be 'max' or 'trace' ;
    Normalization = params.Normalization;
    Normalize = 1;
else
    Normalization = 'none';
    Normalize = 0;
end


if isfield(params,'alpha_0')
    alpha_0 = params.alpha_0;
else
    alpha_0 = 0.1*T;
end



if isfield(params, 'W_thr')   
    W_thr = params.W_thr;
else
    W_thr = 5e-2;
end


if isfield(params,'tau')
    tau = params.tau;
else
    tau = 50;
end


if isfield(params,'C_n_inv')
    C_n_inv = params.C_n_inv;
elseif isfield(params,'std_n')
    std_n =  params.std_n;
    if std_n < 0.1
        std_n = 0.1;
    end
    C_n_inv = 1/(std_n^2) *eye(N);
else
    C_n_inv = 100*eye(N); 
end


if isfield(params,'alpha_1')
    alpha_1 = params.alpha_1;
else
    alpha_1 = 20;  
end


if isfield(params, 'maxIter')   
    maxIter = params.maxIter;
else
    maxIter = 200;
end


compute_error_A = 0;
compute_error_L = 0;
compute_error_X = 0;

if isfield(params,'Atrue')
    Atrue = params.Atrue;
    errors_A = nan(1,maxIter);
    Fscore_A = nan(1,maxIter);
    times_A = nan(1,maxIter);
    compute_error_A = 1;
end

if isfield(params,'Ltrue')
    Ltrue = params.Ltrue;
    errors_L = nan(1,maxIter);
    Fscore_L = nan(1,maxIter);
    times_L = nan(1,maxIter);
    compute_error_L = 1;
end

if isfield(params,'Xtrue')
    Xtrue = params.Xtrue;
    errors_X = nan(1,maxIter);
    times_X = nan(1,maxIter);
    compute_error_X = 1;
end

Optimize_X = 1;
Optimize_L = 1;
Optimize_A = 1;
if isfield(params,'Optimize_X')
    Optimize_X = params.Optimize_X; end
if isfield(params,'Optimize_L')
    Optimize_L = params.Optimize_L; end
if isfield(params,'Optimize_A')
    Optimize_A = params.Optimize_A; end

if isfield(params,'Compute_f')
    Compute_f = params.Compute_f;
else
    Compute_f = 0;
end

if Compute_f
    f_X = 0;
    f_A = 0;
    f_w = 0;
    f_arr = nan(1,maxIter+1);
    f_arr(1) = Inf;
end

rel_errors_X = nan(1,maxIter);
rel_errors_A = nan(1,maxIter);
rel_errors_w = nan(1,maxIter);

has_converged_w = 0;
has_converged_A = 0;
has_converged_X = 0;



% ==============================
% Initialization

if isfield(params,'X0')
    X = params.X0;
else
    X = Y;
end



if isfield(params,'L0')
    L = params.L0;
    w = w_from_L(L);
else

    w = w_from_L(S_inv-J);
    w = w + 0.01*mean(w);
    
    if strcmp(Normalization, 'trace')==1
        w = w/sum(w)*N/2;
    elseif strcmp(Normalization, 'max')==1
        w = w/max(w);
    end
    L = L_operator(w);

end

if isfield(params,'A0')
    A = params.A0;
else
    A = zeros(N);
end

C_n_inv_norm = norm(C_n_inv,2);
L_norm = norm(L,2);
A_norm = norm(A,2);
X_1 = [zeros(N,1), X(:,1:end-1)];
X1_norm = norm(X_1);
B = Mask.*(C_n_inv*Y);




reltol = 1e-5;
Nloop = 1;
t0 = tic;

%=========================
% Iterative updates

for i=1:maxIter

    if Compute_f
        Cx = A*X_1;
        Sx = (X*X');
        Bx = X*Cx';
        
        Dx = (Cx*Cx');
        K = ( Sx - Bx - Bx' + Dx + 2*alpha_0 * Hoff)/T;
    
        r = conj_L_operator( K );
        f_joint = T*sum(w.*r);
    end


    if Optimize_X
        for loop = 1:Nloop
            if Compute_f
                f_X = C_n_inv_norm * norm(Y-Mask.*X,"fro")^2;
            end
            mu = C_n_inv_norm + 2*L_norm*(1+ A_norm^2);
            
            X_1 = [zeros(N,1), X(:,1:end-1)];
            
            U = X - A*X_1;
            G1 = L*U;
            G2 = A'*[G1(:,2:end), zeros(N,1)];
            G3 = -B + Mask.*(C_n_inv*X);
            X_new = X - 1/mu*(G1 -G2 + G3);

            
            rel_err = norm(X_new - X,'fro')/norm(X,'fro');
            rel_errors_X(i) = rel_err;
            has_converged_X = rel_err< reltol;
            if (has_converged_X) &&  (i>1)
                Optimize_X = 0;
                break;
            end
            X = X_new;
            
            if compute_error_X
                errors_X(i) = norm(X - Xtrue, 'fro');
                times_X(i) = toc(t0);
            end
            
            
            X_1 = [zeros(N,1), X(:,1:end-1)];
            X1_norm = norm(X_1);
            

        end
    end
    


    if Optimize_A
        for loop = 1:Nloop            
            

            beta = X1_norm^2* L_norm;
            
            if Compute_f
                f_A = 2*alpha_1*sum(sum(abs(A)));
            end
            
            Diff = A*X_1-X;
            Diff = Diff*X_1';
            temp = L*Diff;
            
            gk = A - (1/beta)*temp ;
            thr = alpha_1/beta;
            A_new = Threshold_Soft(gk,thr) ;
            
            
            rel_err = norm(A_new-A,'fro') / norm(A, 'fro');
            rel_errors_A(i) = rel_err;
            has_converged_A = rel_err< reltol;
            if (has_converged_A) &&  (i>1)
                Optimize_A = 0;
                break;
            end
            
            
            A = A_new;
            
            if compute_error_A
                if Normalize
                    Ahat = A/max(max(abs(A)));
                    Ahat(abs(Ahat)<W_thr) = 0;
                    Ahat = Ahat/norm(Ahat);
                else
                    Ahat = A; 
                end 
                errors_A(i) = norm(Ahat- Atrue,'fro');
                Fscore_A(i) = Fscore_metric(Atrue, Ahat);
                times_A(i) = toc(t0);
            end
            
            
            A_norm = norm(A,2);
            
        end
    end



    if Optimize_L
        for loop = 1:Nloop

            U = X - A*X_1;
            K = ( (U*U')  + 2*alpha_0 * Hoff)/T;    
            r = conj_L_operator( K );
            
             
            [U,Lamb] = eig(L + J);
            lambs = diag(Lamb);

            if Compute_f
                log_det = log(prod(lambs));
            end

            lambs = 1./lambs;
            L_inv = (U*diag(lambs)*U');        
    
            q =  conj_L_operator( L_inv );
            zeta = tau.*w.^2;
            ratio = (zeta + 1).* q ./ (zeta.*q+r);
            w_new = w.* sqrt(ratio);
            
            if Compute_f
                f_w = -T* log_det + 2*alpha_0*sum(sum(L.*Hoff));
            end

            rel_err = norm(w_new - w) / norm(w);
            rel_errors_w(i) = rel_err;
            has_converged_w = rel_err< reltol;
            if (has_converged_w) &&  (i>1)
                Optimize_L = 0;
                break;
            end
           
            
            w = w_new;
            L = L_operator(w);
            

            if compute_error_L
                what = w;
                what( what< (W_thr*max(w)) ) = 0;
                Lhat = L_operator(what);
                if strcmp(Normalization, 'trace')==1
                    Lhat = Lhat/trace(Lhat)*N;
                elseif strcmp(Normalization, 'max')==1
                    Lhat = Lhat/max(what);
                end
                errors_L(i) = norm(Lhat - Ltrue, 'fro');
                Fscore_L(i) = Fscore_metric(Ltrue, Lhat);
                times_L(i) = toc(t0);
            end
            
            L_norm = norm(L,2);
            
        end
    end
    
    if Compute_f
        f_arr(i+1) = f_w + f_A + f_X + f_joint;
    end

    if (has_converged_w) && (has_converged_A) && (has_converged_X)
         break;
    end

       
end

fprintf('STSRGL finished at iteration %d\n',i);

L = L_operator(w);

if compute_error_X
    output.error_X = errors_X;
    output.times_X = times_X;
end

if compute_error_L
    output.error_L = errors_L;
    output.Fscore_L = Fscore_L;
    output.times_L = times_L;
end

if compute_error_A
    output.error_A = errors_A;
    output.Fscore_A = Fscore_A;
    output.times_A = times_A;
end



% Output matrices
output.X = X;
output.L = L;
output.A = A;

output.rel_error_A = rel_errors_A;
output.rel_error_X = rel_errors_X;
output.rel_error_w = rel_errors_w;

if Compute_f
    output.f = f_arr;
end


end

% Soft thresholding function
function Ahat = Threshold_Soft(A,thr)
    N = size(A,1);
    a = A(:);
    temp = abs(a)-thr;
    supp = temp>=0;
    ahat = zeros(size(a));
    ahat(supp) = sign(a(supp)).*temp(supp);
    Ahat = reshape(ahat, N, N);
end
