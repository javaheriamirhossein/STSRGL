function [FPR,TPR,Thr,AUC] = perfcurve_graph(A_true,A,Npt)

% Plots the ROC 
if nargin<3
    Npt = 100;
end

Thr = linspace(1,0,Npt);
TPR = zeros(1,Npt);
FPR = zeros(1,Npt);

P = sum(A_true);
N = sum(~A_true);
AUC = 0;

for i = 1:Npt
    Adjacency = A>=Thr(i);
    TPR(i) = sum(sum(Adjacency.*A_true))/P;
    FPR(i) = sum(sum(Adjacency.*(~A_true)))/N;
    if i>1
        AUC = AUC + (FPR(i)-FPR(i-1))*TPR(i-1);
    end
end

end