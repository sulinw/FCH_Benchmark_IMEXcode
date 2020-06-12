function W = Wz(u,qtype)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Keith's Well - with a foot at u = bm.
% Standard W(u)=(u-a)^2/2 * [ (u-b)^2/2 + gam/3*(u-(3b-a)/2) ];
% New Wq(u)=[ (u-a)^2/2 *(1-2*Coef*(u-(a+b)/2)^2) + qtype*eps*(1-sech((u-bm)/eps))]
%            * [ (u-b)^2/2 + gam/3*(u-(3b-a)/2) ];
% {Wz, Wzp, Wzpp} can be sort together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global bm bp gam eps CoefE
Q = eps*(1-sech((u-bm)/eps));

W1 = (u-bm).^2.*(1/2-CoefE*(u-(bm+bp)/2).^2)+qtype.*Q;
W2 = (u-bp).^2/2+(gam/3)*(u-(3*bp-bm)/2);
W = W1.*W2;
end