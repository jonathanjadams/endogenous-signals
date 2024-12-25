function SLk = lagmulti(S,k) %this function left-multiplies an operator by the lag operator of order k
if k>=0
SLk = cat(3,zeros(size(S,1),size(S,2),k),S(:,:,1:(end-k)));
else
SLk = cat(3,S(:,:,(1-k):end),zeros(size(S,1),size(S,2),-k));
end
end

%   SIGNAL_OP Package: Version 1.0
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu