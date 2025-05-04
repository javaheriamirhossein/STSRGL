function [Fscore_val,FDR_val] = Fscore_metric(Mtrue,M)

A_connectivity_true = Mtrue~= 0;
Adjacency = M~= 0;
    

FDR_val = 1-sum(sum(abs(Adjacency-A_connectivity_true)))/sum(sum(A_connectivity_true));
tp = sum(sum(Adjacency.*A_connectivity_true));
fp = sum(sum(Adjacency.*(~A_connectivity_true)));
fn = sum(sum((~Adjacency).*A_connectivity_true));
Fscore_val = 2*tp/(2*tp+fn+fp);
end