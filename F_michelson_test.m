function [F]=F_michelson_test(m,lambda)

% The Galerkin projection
% Inputs
%  1) x = (L,a_0,a_1,b_1,a_2,b_2,...,a_{m-1},b_{m-1}).
%  2) lambda : parameter

x=rand(6*m-2,1);
F=zeros(6*m-2,1);
L=x(1);
a=[x(2);x(2*(1:3*m-2)+1)];
b=[0;x(2*(2:3*m-1))];

aa=quadratic_sumFFTee(a,a);
ab=quadratic_sumFFTeo(a,b);
bb=quadratic_sumFFToo(b,b);

%%% Phase Condition %%%
%u0=0; 
%F(1)=a(1)+2*sum(a(2:m))-u0; %%% The periodic solution satisfies u(0)=u0

F(1)=a(2); %%% a_1=0
%%% Second term
F(2)=(1/2)*(aa(1)-bb(1))-lambda^2;

%%% Remaining terms
k=(1:3*m-2)';
F(2*(1:3*m-2)+1)=(k.^3*L^3-k*L).*b(2:end)+(1/2)*(aa(2:end)-bb(2:end));
F(2*(2:3*m-1))=(-k.^3*L^3+k*L).*a(2:end)+ab(2:end);

end

function [s]=quadratic_sumFFTee(a1,a2)

%%% Inputs: a1 and a2 that are thought of as being even, that is 
%%%         a1_{-k} = a1_k and a2_{-k} = a2_k

m=length(a1);

b1=flipdim(a1(2:m),1);
ta1=[zeros(m,1);b1;a1;zeros(m,1)];
tu1=ifft(ifftshift(ta1));

b2=flipdim(a2(2:m),1);
ta2=[zeros(m,1);b2;a2;zeros(m,1)];
tu2=ifft(ifftshift(ta2));

F=fftshift(fft(tu1.*tu2));

s=real((4*m-1)*F(2*m:3*m-1));

end

function [s]=quadratic_sumFFToo(a1,a2)

%%% Inputs: a1 and a2 that are thought of as being odd, that is 
%%%         a1_{-k} = -a1_k and a2_{-k} = -a2_k

m=length(a1);

b1=flipdim(a1(2:m),1);
ta1=[zeros(m,1);-b1;a1;zeros(m,1)];
tu1=ifft(ifftshift(ta1));

b2=flipdim(a2(2:m),1);
ta2=[zeros(m,1);-b2;a2;zeros(m,1)];
tu2=ifft(ifftshift(ta2));

F=fftshift(fft(tu1.*tu2));

s=real((4*m-1)*F(2*m:3*m-1));

end

function [s]=quadratic_sumFFTeo(a1,a2)

%%% Inputs: a1 is thought of as being even, that is a1_{-k} = a1_k
%%%          and a2 is thought of as being odd, that is a2_{-k} = -a2_k

m=length(a1);

b1=flipdim(a1(2:m),1);
ta1=[zeros(m,1);b1;a1;zeros(m,1)];
tu1=ifft(ifftshift(ta1));

b2=flipdim(a2(2:m),1);
ta2=[zeros(m,1);-b2;a2;zeros(m,1)];
tu2=ifft(ifftshift(ta2));

F=fftshift(fft(tu1.*tu2));

s=real((4*m-1)*F(2*m:3*m-1));

end