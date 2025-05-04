function  Error  = NMSE( Xtrue, Xhat )
[N, T] = size(Xtrue);
Error = 0;
for i=1:T
    Error = Error + norm( Xtrue(:,i) - Xhat(:,i) )^2 / norm( Xtrue(:,i) )^2;
end
Error = Error/T;
end

