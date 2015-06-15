function [ thecovariance ] = Cov_theorical(K, L, M)
% Author: Wacha Bounliphone - wacha.bounliphone@centralesupelec.fr
% Copyright (c) 2014-2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
%
% If you use this software in your research, please cite:
%
%Bounliphone, W., Gretton, A.,  Tenenhaus, A., Blaschko, M.B. (2015).  
%A low variance consistent test of relative dependency.  
%International Conference on Machine Learning, Jul 2015, Lille, France. 



%numSamples
m = size(K,1);

%R_xyxz
constante  = (1/(4*m)) * (1/((m-1)*(m-2)*(m-3))^2);

K = centerKern(K);
L = centerKern(L);
M = centerKern(M);
Ks = K-diag(diag(K));
Ls = L-diag(diag(L));
Ms = M-diag(diag(M));

h_xy = computehvec(Ks,Ls);
h_xz = computehvec(Ks,Ms);

R_xyxz = constante * h_xy' * h_xz;



%HSIC_xy
HSIC_xy = HSICunbiased(K,L);

%HSIC_xz
HSIC_xz = HSICunbiased(K,M);


%thecovariance
thecovariance = (16/m) *  (R_xyxz - (HSIC_xy*HSIC_xz));

end


function h_xy = computehvec(Ks,Ls);
m = size(Ks,1);
e = ones(m,1);
t1 = ((m-2)*(m-2)).*(Ks.*Ls)*e;
t2 = (m-2)*( ((sum(sum(Ks.*Ls)))*e) - (Ks*(Ls*e)) - (Ls*(Ks * e)) );
t3 = m * (Ks *e) .* (Ls * e);
t4 = (sum(sum(Ls))) * (Ks*e);
t5 = (sum(sum(Ks))) * (Ls*e);
t6 = ( (e'*Ks) * (Ls * e) )* e;

h_xy = t1 +t2 - t3 + t4 + t5 -t6;
end 



