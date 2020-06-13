function [derivx,derivy,lap,lap0] = init_operators(N,ratio,index)
if (mod(N,2)==0)
    alpha_alias  = [0:N/2-1, N/2, -N/2+1:-1];
    alpha_alias2 = [0:N/2-1,   0, -N/2+1:-1];
else
    alpha_alias  = [0:(N-1)/2, -(N-1)/2:-1];
    alpha_alias2 = [0:(N-1)/2, -(N-1)/2:-1];
end
%
[alpha,beta] = meshgrid(alpha_alias);
[alpha2,beta2] = meshgrid(alpha_alias2);
%
derivx = 1i*alpha2/ratio;
derivy = 1i*beta2/ratio;
%
if index == 1
    lap = - (alpha .^2+beta .^2)/(ratio^2);
else
    lap = - (alpha2.^2+beta2.^2)/(ratio^2);
end
lap0 = lap; lap0(1,1) = 1.0;   % this is necessary for PSD scheme for lap.^(-1)