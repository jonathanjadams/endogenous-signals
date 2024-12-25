function V = s2v(S,N,column)
if nargin<3
    column = 1; %default: make a column vector
end
if column==1
    V=reshape(permute(flip(S,3),[2 1 3]),size(S,2),[])'; 
    %flip: because causal signals correspond to upper triangular operators
else
    V=reshape(S,size(S,1),[]); %create a block row vector
end
end


%   SIGNAL_OP Package: Version 1.1
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu