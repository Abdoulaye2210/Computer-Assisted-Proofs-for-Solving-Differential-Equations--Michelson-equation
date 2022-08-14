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