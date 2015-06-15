function pvalue=relative_dependency_test(K,L,M);
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

% K : target
% L : source 1
% M : source 2
Cov(1,1) = Cov_theorical(K, L, L);
Cov(2,2) = Cov_theorical(K, M, M);
Cov(2,1) = Cov_theorical(K, L, M);
Cov(1,2) = Cov(2,1);


H = [HSICunbiased(K,L);HSICunbiased(K,M)];


theta=pi/4;
R = [cos(theta) -sin(theta);sin(theta) cos(theta)];

H = R*H;
Cov = R*Cov*R';

t = H(1);
sig = sqrt(Cov(1,1));
t = t/sig;
pvalue=1-normcdf(t);


end
