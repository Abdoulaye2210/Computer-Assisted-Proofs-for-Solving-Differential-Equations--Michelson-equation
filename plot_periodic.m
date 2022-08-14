function []=plot_periodic(x)

period=2*pi/abs(x(1));
plot((0:.001:period),periodic(x,(0:.001:period)),'Color',[0 0 1],'LineWidth',2)
box off
hold on
plot([0 period],[0 0],'color',[0 0 0])

end

function [periodic]=periodic(x,t)

m1=length(x);
m=m1/2;
L=x(1,1);
a=[x(2,1);x((2*(1:m-1)'+1),1)];
b=[0;x(2*(2:m)',1)];

sum=0;
for k=1:m-1
    sum=sum+a(k+1,1)*cos(k*L*t)-b(k+1,1)*sin(k*L*t);
end
periodic=a(1,1)+2*sum;
end
