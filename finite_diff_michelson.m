function [DF]=finite_diff_michelson(x,lambda)

h=1e-8;
M=length(x);
E=eye(M);
DF=zeros(M);
for j=1:M
    xh=x+h*E(:,j);
    DF(:,j)=(F_michelson(xh,lambda)-F_michelson(x,lambda))/h;
end
end