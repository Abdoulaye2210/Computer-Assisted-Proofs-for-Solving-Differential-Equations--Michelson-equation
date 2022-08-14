function [x]=newton_michelson(x,lambda)

% The Galerkin projection
% Inputs
%  1) x = (L,a_0,a_1,b_1,a_2,b_2,...,a_{m-1},b_{m-1}).
%  2) lambda : parameter

tol=4e-16;
F=F_michelson(x,lambda);
DF_inv=inv(DF_michelson(x));

k=0;
while (k<=20) && (norm(F)> tol),
    x = x - DF_inv*F_michelson(x,lambda);
    DF_inv = inv(DF_michelson(x));
    F=F_michelson(x,lambda);
    %display(['La norme de F = ',num2str(norm(F))])
    k=k+1;
end

if k>20
    display('Newton did not converge !')
else
    if (norm(F)< tol)
        display(['Newton converged after ',num2str(k),' iterations'])
    end
end

if (norm(F)< tol) && (abs(x(1))>.001)
    figure
    plot_periodic(x)
end

return
