function Wp = Wzp(u,qtype)
global bm bp gam eps CoefE
term = sech((u-bm)/eps);
%
Q = eps*(1-term);
Qp = term.*tanh((u-bm)/eps);
term1 = (u-bm).^2;        
term1p = 2*(u-bm);
term2 = 1/2-CoefE*(u-(bm+bp)/2).^2;
term2p = -2*CoefE*(u-(bm+bp)/2);

W1 = term1.*term2+qtype.*Q;
W1p = (term1.*term2p+term1p.*term2)+qtype.*Qp;
W2 = (u-bp).^2/2+(gam/3)*(u-(3*bp-bm)/2);
W2p = (u-bp)+(gam/3);

Wp = W1p.*W2+W1.*W2p;
end