function [s,C] = alg_decay_b_least_square(b)

m=length(b);

abs_b=abs(b);

Y=log(abs_b+realmin); %% We add the value realmin to avoid getting -Inf

k=(1:m)';
log_k=log(k);
sum_log_k=sum(log_k);
sum_log_k2=sum(log_k.^2);

sum_Y=sum(Y);
sum_log_kY=sum(log_k.*Y);
denom=m*sum_log_k2-sum_log_k^2;
A=(m*sum_log_kY-sum_log_k*sum_Y)/denom;
B=(sum_Y*sum_log_k2-sum_log_k*sum_log_kY)/denom;

C=exp(B);
s=-A;

k(1)=[];
b(1)=[];

plot(k,C*k.^(-s),'color',[1 0 0],'linewidth',2)
hold on
plot(k,b,'color',[0 0 0])

end






