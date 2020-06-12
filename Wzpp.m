function Wpp = Wzpp(u,qtype)
global bm bp gam eps CoefE
term = sech((u-bm)/eps);
Q = eps*(1-term);
Qp = term.*tanh((u-bm)/eps);
Qpp = 1/eps*term.*(2*term.^2-1);

term1 = (u-bm).^2;        
term1p = 2*(u-bm);
term1pp = 2;

term2 = 1/2-CoefE*(u-(bm+bp)/2).^2;
term2p = -2*CoefE*(u-(bm+bp)/2);
term2pp = -2*CoefE;

W1 = term1.*term2+qtype.*Q;
W1p = (term1.*term2p+term1p.*term2)+qtype.*Qp;
W1pp = (term1.*term2pp+2*term1p.*term2p+term1pp.*term2)+qtype.*Qpp;
W2 = (u-bp).^2/2+(gam/3)*(u-(3*bp-bm)/2);
W2p = (u-bp)+(gam/3);
W2pp = 1;

Wpp = W1pp.*W2+2*W1p.*W2p+W1.*W2pp;
end