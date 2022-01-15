function dxdt = sveqhird(~,x,r)
dxdt = zeros(8,1);

dxdt(1)=r(1)-(x(4)-r(4))*x(4)-r(5)*x(1)-r(3)*(x(1)-r(5))*x(1);
dxdt(2)=r(5)*x(1)-r(2)*x(2);
dxdt(3)=r(3)*(x(1)-r(5))*x(1)-x(3)^2;
dxdt(4)=r(2)*x(2)+(x(3)-r(6))*x(3)-x(4)^2;
dxdt(5)=r(6)*x(3)+r(4)*x(4)-x(5)^2;
dxdt(6)=r(7)*x(5)-x(6)^2;
dxdt(7)=(x(5)-r(7)-r(8)*(x(5)-r(7)))*x(5)+(x(6)-r(9))*x(6);
dxdt(8)=(r(8)*(x(5)-r(7)))*x(5)+r(9)*x(6);
end