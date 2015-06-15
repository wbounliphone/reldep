function K = centerKern(K);
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


m = size(K,1);
e = ones(m,1);
K = K - e*(e'*K)./m; %premultiply by H
K = K - (K*e)*e'./m; %postmultiply by H
end
