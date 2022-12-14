function [rho,C] = exp_decay_b_least_square(b)

abs_b=abs(b);

m=length(abs_b);

Y=log(abs_b);

k=(1:m)';
sum_k=sum(k);
sum_k2=sum(k.^2);

sum_Y=sum(Y);
sum_kY=sum(k.*Y);
denom=m*sum_k2-sum_k^2;
A=(m*sum_kY-sum_k*sum_Y)/denom;
B=(sum_Y*sum_k2-sum_k*sum_kY)/denom;

C=exp(B);
rho=exp(-A);

display(['C = ', num2str(C),', rho = ', num2str(rho)])

plot(k,A*k+B,'color',[0 0 0],'linewidth',1)
hold on
plot(Y,'*','color',[1 0 0])

end






