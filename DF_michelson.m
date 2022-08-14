function [DF]=DF_michelson(x)

% The Galerkin projection
% Inputs
%  1) x = (L,a_0,a_1,b_1,a_2,b_2,...,a_{m-1},b_{m-1}).
%  2) lambda : parameter

m=length(x)/2;
DF=zeros(2*m);
L=x(1);
a=[x(2);x(2*(1:m-1)+1)];
b=[0;x(2*(2:m))];

%%% Derivative of the phase condition %%%

% tempJ1=[0 1 2*ones(1,2*m-2)]; tempJ1(1,(2*(2:m)))=zeros(m-1,1);
% DF(1,:)=tempJ1;

DF(1,3)=1; %% a_1=0

%%% Second row %%%
DF(2,2)=a(1); 
DF(2,2*(1:m-1)+1)=2*a(2:end)'; 
DF(2,2*(1:m-1)+2)=2*b(2:end)';

%%% First column %%%
k=(1:m-1)';
DF(2*(1:m-1)+1,1)=(3*k.^3*L^2-k).*b(2:end);
DF(2*(2:m),1)=-(3*k.^3*L^2-k).*a(2:end);

%%% Second column %%%
DF(2*(1:m-1)+1,2)=a(2:end);
DF(2*(2:m),2)=b(2:end);

%%% Remaining rows %%% 

for k=1:m-1  
    R=[[0 k^3*L^3-k*L];[-(k^3*L^3-k*L) 0]];
    for l=1:m-1
        tempJ_kl=isequal(k,l)*R+[[coeff_a(a,k-l) -coeff_b(b,k-l)];[coeff_b(b,k-l) coeff_a(a,k-l)]];
        if k+l>=m
            DF([2*k+1;2*k+2],[2*l+1;2*l+2])=tempJ_kl;
        else
            DF([2*k+1;2*k+2],[2*l+1;2*l+2])=tempJ_kl+[[coeff_a(a,k+l) coeff_b(b,k+l)];[coeff_b(b,k+l) -coeff_a(a,k+l)]];
        end
    end
    
end

end

function [coeff]=coeff_a(a,n1)
    m=length(a);
    if abs(n1)<m
        coeff=a(abs(n1)+1,1);
    else
        coeff=0;
    end
end

function [coeff]=coeff_b(b,n1)
    m=length(b);
    if abs(n1)<m
        if n1<0
            coeff=-(b(abs(n1)+1,1));
        else
            coeff=b(n1+1,1);
        end
    else
        coeff=0;
    end
end
