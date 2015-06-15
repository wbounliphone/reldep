function [HSIC_xy,HSIC_xy_b] = HSICunbiased(K,L);
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

% Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., Smola, A. J. (2007). 
%A kernel statistical test of independence. 
%In Advances in Neural Information Processing Systems, 585-592.

m = size(K,1);
e = ones(m,1);

Ks = centerKern(K);
Ls = centerKern(L);
clear K L
Ks = Ks-diag(diag(Ks));
Ls = Ls-diag(diag(Ls));

t1 = sum(sum(Ks.*Ls));
t2 = sum(sum(Ks)).*sum(sum(Ls))./((m-1)*(m-2));
t3 = (e'*Ks)*(Ls*e)*2/(m-2);

HSIC_xy = (t1+t2-t3)/(m*(m-3));

if(nargout>1)
  HSIC_xy_b = trace(Ks*Ls)/(m*(m-1));
end

end

