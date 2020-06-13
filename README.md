# FCH_Benchmark_IMEXcode
% This code is for Benchmark problems by using a BDF2-IMEX scheme to solve the FCH equation 
%                            u_t = G dE/dt,
% described in the paper 
%         "BENCHMARK COMPUTATION OF MORPHOLOGICAL COMPLEXITY IN THE FCH GRADIENT FLOW"
% by Andrew Christlieb, Keith Promislowm Zengqiang Tan, Sulin Wang, Brian Wetton, and Steven Wise.
% This paper has been submitted to Elsevier on June 6, 2020
% You could download the paper: https://arxiv.org/abs/2006.04784

% where the operator G := Lap/(1-rho*Lap), dE/dt is the variational derivative of the Energy
%    E = Energy(u) := INTEGRAL(1/2*(eps^2 Lap u - Wq'(u))^2 - (eps^2/2*eta1*|Grad u|^2+eta2*Wq(u))) 
% where eta1 and eta2 are constants about eps,
%    Wq(u)=[(u-a)^2 *(1/2-Coef*(u-(a+b)/2)^2) + qtype*eps*(1-sech((u-a)/eps))]
%          * [(u-b)^2/2+gam/3*(u-(3b-a)/2)];
% where Coef is some constant to make sure u=b is a local minimum. Here a:=bm=-1, b:=bp=1, gam=0.3,
% To get the Wq(u) in the FCH-Benchmark paper (Promislow, etc), please let CoefE = 0.

% The scheme is as follows:
%       a*u^{n+1}+b*u^{n}+c*u^{n-1} = G VarDerivE(u^{n+1},u^{*,n+1})
% where in the VarDerivE(u^{n+1},u^{*,n+1}),
% the linear operator about u^{n+1} is (eps2^2*lap-alpham).^2, (beta1,beta2)=(2,1) in the paper. 
% for other nonlinear and eta2 terms, u^{*,n+1}=interpolation(u^{n-1},u^n).

% To start, please try examples like follows:
% [phi,T0,FFT2] = BDF2IMEX;
% or 
% [phi,T0,FFT2] = BDF2IMEX(256,0,0,20,0.2,3,1e-5,0.2);

% Verision 2. Updated on June 12, 2020 by Sulin Wang.
%% For any questions, please contact Sulin Wang at Email: wangsuli@msu.edu
