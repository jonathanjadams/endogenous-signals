function S = v2s(V,N,column)
if nargin<3
    column=0; 
end
    if column ==0 %i.e. its a row: V= [S0, S1, S2...]
    S=reshape(V,size(V,1),size(V,2)/N,[]);
    else %otherwise  V = [... S2; S1; S0] 
    S=flip(permute(reshape(V',size(V,2),size(V,1)/N,[]),[2 1 3]),3);
    end
        
end


%   SIGNAL_OP Package: Version 1.1
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu