function EgyIntgd = EnergyIntegrand(phi,lap,derivx,derivy,eps,eta1,eta2,qtype)
%% the Energy functional of FCH equatio, note that eta's have eps already
Uh = fft2(phi);
LapU = ifft2(lap.*Uh);

Ux = ifft2(derivx.*Uh);   % derivx:= 1i*alpha2/ratio and derivy:= 1i*beta2/ratio
Uy = ifft2(derivy.*Uh);   % look at init_operators.m and FCH_Imex_MCP_New.m

H2f = 1/2*(eps^2*LapU-Wzp(phi,qtype)).^2; % the first part of integrand
Eta = 1/2*eps^2*eta1.*(Ux.^2+Uy.^2)+eta2.*Wz(phi,qtype);  % second part of integrand
EgyIntgd = H2f-Eta;       % the whole integrand
end