function Z_inv=inv_signal(S,N)
Z_inv=t2s(s2t(S,N)^-1,N); %this inverts the toeplitz operator if possible
end

%   SIGNAL_OP Package: Version 1.0
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu