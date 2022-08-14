close all

load S2

steps=40;
delta=.01;

N=length(x);

x_tot=zeros(N,steps);
lambda_tot=zeros(1,steps);
norm_tot=zeros(1,steps);

lambda_tot(1)=lambda;
x_tot(:,1)=x;
norm_tot(1)=norm(x);

%%% Oth order predictor-corrector
for k=1:steps
    lambda_tot(k+1)=lambda_tot(k)+delta;
    [x_tot(:,k+1)]=newton_michelson(x_tot(:,k),lambda_tot(k+1));
    norm_tot(k+1)=norm(x_tot(:,k+1));
end

figure
plot(lambda_tot,norm_tot)

