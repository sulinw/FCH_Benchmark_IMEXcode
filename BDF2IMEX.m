function [phi,T0,FFT2] = BDF2IMEX(N,phi,T0,T,dcoef,delta,tol,qtype)
% This code uses an BDF2-IMEX scheme to solve the FCH equation 
%                            u_t = G dE/dt,
% described in the paper 
%         "BENCHMARK COMPUTATION OF MORPHOLOGICAL COMPLEXITY IN THE FCH GRADIENT FLOW"
% by Andrew Christlieb, Keith Promislowm Zengqiang Tan, Sulin Wang, Brian Wetton, and Steven Wise.

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

% Selected Important Variables
% alpham := Wq''(bm,qtype) controls the foot at u=bm;
% dcoef is used in constructing the initial data;
% FFT2 is to calculate the FFT2 calls (EnergyIntegrand(...) is excluded);
% {bm,bp,gam,qtype,CoefE} are used in Wq(u) and its derivatives;
% {derivx,derivy,X,getFull} are used in EnergyIntegrand to compute Energy(time);
% {dt_min,dt_max,safty,VSUbd} are for updating the temporal step-sizes;
% tMCLdtE := [time; Max(u); u(center); u(left-corner); dt; Energy];
% Main.{ttime(k),SolutionT(:,:,k)} is a strcture data.

% To start, please try examples like follows:
% [phi,T0,FFT2] = BDF2IMEX;
% or 
% [phi,T0,FFT2] = BDF2IMEX(256,0,0,20,0.2,3,1e-5,0.2);

% Verision 2. Updated on June 12, 2020 by Sulin Wang.
%% For any questions, please contact Sulin Wang at Email: wangsuli@msu.edu

%% ====================== Basic settings start ==============================
global bm bp gam eps alpham CoefE
L = 4*pi;
bm = -1;
bp = 1;
gam = 0.3;
eps = 0.1;  
FFT2 = 0;  % FFT2 counts (EnergyIntegrand(...) is excluded)

tic
if nargin == 0    % the default parameters {N,phi,T0,T,dcoef,delta,tol,qtype}
    N = 2^8;
    T0 = 0;
    T = 10;
    dcoef = 0.2;
	delta = 3;    % the default eta2 = 3*eps
    tol = 1e-5;
    qtype = 0;
end
eta1 = 1.45*eps;  eta2 = delta*eps; 

% The constant CoefE is used to make sure u=bp=1 is a local min for Wq(u),
% i,e, Wq'(bm) = 0. So CoefE = 0 for qtype=0; however, CoefE is nonzero for other qtype
CoefE = qtype/(bm-bp)^4*(-2*eps+sech((-bm+bp)/eps)*(2*eps+(-bm+bp)*tanh((-bm+bp)/eps)));
alpham = Wzpp(bm,qtype);
fprintf('(Wq''(%d,q))^2+(Wq''(%d,q))^2=%.2e\n',bm,bp,(Wzp(-1,qtype))^2+(Wzp(1,qtype))^2);

% for the temporal step-size
dt_max = 100;
dt_min = 1e-9;
safty = 0.9;
VSUbd = 1+sqrt(2)-1e-10;  % to ensure zero-stability, the upper bound for variable BDF2

%
[derivx,derivy,lap,~] = init_operators(N,L/(2*pi),1); % derivx,derivy is used for Energy(time)
eps2lap = eps^2*lap;
rho = 0;
G = lap./(1-rho*lap);   % this is for G=lap/(1-rho*lap), for rho>=0
Sigma = G.*(eps2lap-alpham).^2;
OpeF = @(oper,f)real(ifft2(oper.*fft2(f)));   % gives operator(func), like Lap(u)

%% set up the initial data
if T0 == 0
    load('phi256_aSym2.mat','phi256','sym_Sym');
    h0 = L/256;       [xx0,yy0] = meshgrid(h0:h0:L);
    h = L/N;            [xx,yy] = meshgrid(h:h:L);
    phi = interp2(xx0,yy0,phi256,xx,yy,'spline');
    clear h0 xx0 yy0;
    phi = phi+dcoef*eps/(1.7^2);   % add mass to the background, here 1.7 = Wq(bm,q=0)
end
endT = T0+T;
FileName1 = ['IMEX2_',sym_Sym,'_{q=',num2str(qtype),',d=',num2str(dcoef),',tol=',num2str(tol)];
FileName2 = [',eta2=',num2str(delta),'eps,N=',num2str(N),',T=',num2str(T0),'-',num2str(endT),'}'];
FileName = [FileName1,FileName2];

%% Main code for the adaptive BDF2IMEX-AM3 scheme
% getFull and X are used in computing Energy(time) since the data is periodic
getFull = @(A)[A,A(:,1);A(1,:),A(1,1)];  
X = linspace(0,L,N+1);    % starts on Feb 15th
Nm = round(N/2);
%
k = 1;
time = T0;
dt = dt_min;
tMCLdtE = [];
tMCLdtE(1,k) = time;
tMCLdtE(2,k) = max(phi(:));
tMCLdtE(3,k) = phi(Nm,Nm)-bm;
tMCLdtE(4,k) = phi(1,1)-bm;
tMCLdtE(5,k) = dt;
EgyIntgd = EnergyIntegrand(phi,lap,derivx,derivy,eps,eta1,eta2,qtype);
tMCLdtE(6,k) = trapz(X,trapz(X,getFull(EgyIntgd)));
%
Main.ttime(k) = time;
Main.SolutionT(:,:,k) = phi;
phioo = phi;
RHSoo = RightHandSide(phi);

%% 
phi = IMEX1(phi,dt);   % IMEX1 solution for the first step
k = k+1;        % now k is 2
dt = dt_min;
time = time+dt; % now the time = T0+dt_min
tMCLdtE(1,k) = time;
tMCLdtE(2,k) = max(phi(:));
tMCLdtE(3,k) = phi(Nm,Nm)-bm;
tMCLdtE(4,k) = phi(1,1)-bm;
tMCLdtE(5,k) = dt;
EgyIntgd = EnergyIntegrand(phi,lap,derivx,derivy,eps,eta1,eta2,qtype);
tMCLdtE(6,k) = trapz(X,trapz(X,getFull(EgyIntgd)));
%
phio = phi;
RHSo = RightHandSide(phi);

%%
% below: parameters for save(Main).
p = [0.1,     1, max(1,1/20*T), max(1,1/10*T);
     0.002, 0.1,             1,            2];
k1 = (p(2,2)-p(2,1))/(p(1,2)-p(1,1));
k2 = (p(2,3)-p(2,2))/(p(1,3)-p(1,2));
k3 = (p(2,4)-p(2,3))/(p(1,4)-p(1,3));
w1 = (k2-k1)/2;  w2 = (k3-k2)/2;  w3 = k1+w1+w2;
gw = @(t)w1*abs(t-p(1,2))+w2*abs(t-p(1,3))+w3*t+p(2,2)-((w3-w2)*p(1,2)+w2*p(1,3));
dt_Save = @(t)max(p(2,1),min(gw(t),p(2,4)));
%
Ns = diff([p(1,:),endT])./p(2,:);  
fprintf('It needs at most %.2f MB\n',sum(Ns)*0.5243*(N/256)^2);
% tt = T0:0.01:endT; figure, 
% plot(tt,dt_Save(tt),'-r',tt,gw(tt),'--b',p(1,:),p(2,:),'k*','LineWidth',1.5);
flagT = T0;  idSave = 1;  flagMax = bm; dto = dt;
% above: parameters for save(Main).
while time < endT
    if time+dt > endT     % Make sure we can get the solution precisely at some time T
        dt = endT-time;   % the equivalent one is dt = min(dt,T0+T-time);
    end
    phi2 = IMEX2(phio,phioo,dt,dto);   % IMEX2 solutions
    RHS = RightHandSide(phi2);
    %
    rto = dt/dto;  % the ratio of temporal stepsizes (dto,dt)
    phi3 = phio+dt/6*( (3+2*rto)/(1+rto)*RHS+(3+rto)*RHSo-rto^2/(1+rto)*RHSoo ); % AM3
    
    rErr = norm(phi2(:)-phi3(:),2)/norm(phi3(:),2);   % the relative error
    % disp([dt,rto,rErr]);
    if rErr <= tol || dt == dt_min
        k = k+1;   % this k starts from 3
        time = time+dt;
        phi = phi2;
        
        tMCLdtE(1,k) = time;
        tMCLdtE(2,k) = max(phi(:));
        tMCLdtE(3,k) = phi(Nm,Nm)-bm;
        tMCLdtE(4,k) = phi(1,1)-bm;
        tMCLdtE(5,k) = dt;
        EgyIntgd = EnergyIntegrand(phi,lap,derivx,derivy,eps,eta1,eta2,qtype);
        tMCLdtE(6,k) = trapz(X,trapz(X,getFull(EgyIntgd)));
        
        % save the structure: Main.(ttime,Solution)
        if time-flagT > dt_Save(time) && abs(tMCLdtE(2,k)-flagMax)>=5e-5
            idSave = idSave+1;   flagT = time;
            flagMax = tMCLdtE(2,k);
            Main.ttime(idSave) = time;
            Main.SolutionT(:,:,idSave) = phi;
            fprintf('IMEX2(q=%.1f, d=%.4f, eta2=%.2f eps, tol=%.0e)=%.2f%%: max(T=%.3f) = %.4f, dt = %.3e\n',...
                     qtype,dcoef,delta,tol,100*(time-T0)/T,time,flagMax,dt);
        end
        
        % 
        phioo = phio;   RHSoo = RHSo;
        phio = phi;     RHSo = RHS;
        dto = dt;
    end
    
    % To ensure zero-stability, bound "safty*...(1/3)" by "1+sqrt(2)" from above
    dt = dt*min(VSUbd,safty*(tol/rErr)^(1/3));
    dt = min(max(dt,dt_min),dt_max);
end

%%
T0 = time;
runtime = toc;
fprintf('When {q=%.1f,d=%.4f,tol=%.0e}: max(phi) = %.4f, FFT2 = %d, runtime = %.2f min\n',qtype,dcoef,tol,max(phi(:)),FFT2,runtime/60);
%
clear FileName1 FileName2 phi256 phi2 phi3 phio RHS RHSo RHSoo k1 k2 k3 w1 w2 w3 flagMax flagT;

%% save all data in a mat file (format v7.3) 
save([FileName,'.mat'],'-v7.3'); % 

%% Plot the solution at T
load('myjet.mat','mycolor');
figure, pcolor(xx,yy,phi); shading interp; colormap(mycolor), caxis([-1,1]); colorbar, 
title(['q=',num2str(qtype),', d=',num2str(dcoef),', \eta_2=',num2str(delta),'\epsilon, \sigma=',num2str(tol),', T = ',num2str(time),', max(u) = ',num2str(max(phi(:)))]); 
set(gcf,'unit','centimeters','position',[2,0.5,15.4,13]);
print('-dpng','-r400',[FileName,'.png']);

%% Plot the Energy verse time=0~T
figure,
subplot(1,2,1), 
semilogx(tMCLdtE(1,:),tMCLdtE(6,:),'-b','LineWidth',1.25); 
title('Energy v.s. time');
subplot(1,2,2), 
loglog(tMCLdtE(1,:),tMCLdtE(6,:)-tMCLdtE(6,end),'-.b','LineWidth',1.25); 
title('Energy-Energy(end) v.s. time');
set(gcf,'position',[20 140 1260 470]);
%
disp('program done!!!')

%% ========================== Other functions ==================================
    function [u] = IMEX1(phio,dt)
		WpU0 = Wzp(phio,qtype);
        WppU0 = Wzpp(phio,qtype);
        G0 = (2*alpham+eta1-WppU0).*OpeF(eps2lap,phio)-OpeF(eps2lap,WpU0);
		G1 = (WppU0-eta2).*WpU0-alpham^2*phio;
		GG = OpeF(G,G0+G1);
        %
        Coefinv = 1./(1-dt*Sigma);
        u = OpeF(Coefinv,phio+dt*GG);  % this is IMEX1 solution
        FFT2 = FFT2+4;
    end
    function [u] = IMEX2(phio,phioo,dt,dto)
        phi_p2 = (phio-phioo)/dto*dt+phio;
		WpU0 = Wzp(phi_p2,qtype);
        WppU0 = Wzpp(phi_p2,qtype);
        G0 = (2*alpham+eta1-WppU0).*OpeF(eps2lap,phi_p2)-OpeF(eps2lap,WpU0);
		G1 = (WppU0-eta2).*WpU0-alpham^2*phi_p2;
		GG = OpeF(G,G0+G1);
        %
        Coefinv = 1./(1/dt+1/(dt+dto)-Sigma);
        RHS0 = (1/dt+1/dto)*phio+(1/(dt+dto)-1/dto)*phioo + GG;
        u = OpeF(Coefinv,RHS0);
        FFT2 = FFT2+4;
    end

	%% compute the Variational Derivative: dE/du
    function RHS = RightHandSide(phi)
        Wpu = Wzp(phi,qtype);    % this is a matrix!!!
        Wppu = Wzpp(phi,qtype);  % this is a matrix!!!
        omega = OpeF(eps2lap,phi)-Wpu;
        Right = OpeF(eps2lap,omega)+(eta1-Wppu).*omega+(eta1-eta2)*Wpu;
        RHS = OpeF(G,Right);
        FFT2 = FFT2+3;
    end
end