function S3=smulti(S1,S2,N) %this function "multiplies" 2 operators

if size(S1,1)==0 %if S1 has an empty row dimension, return an empty operator instead of an error
    S3 = zeros(0,size(S2,2),max(size(S1,3),size(S2,3))); 
else

    if size(S1,3)<=N && size(S2,3)<=N
        noncausal=0;
    else
        noncausal=1;  
    end

    if noncausal>0 %need to pad both ends with zeros for accurate convolution
        if size(S1,3)==1
            S1c = cat(3,S1,zeros(size(S1,1),size(S1,2),N-1));
        else
            S1c = S1(:,:,1:N);
        end
        if size(S2,3)==1
            S2c = cat(3,S2,zeros(size(S2,1),size(S2,2),N-1));
        else
            S2c = S2(:,:,1:N);
        end
        if size(S1,3)<=N
            S1nc = cat(3,zeros(size(S1,1),size(S1,2),N-1),S1c);
        else
            S1nc = S1;
        end
        if size(S2,3)<=N
            S2nc = cat(3,zeros(size(S2,1),size(S2,2),N-1),S2c);
        else
            S2nc = S2;
        end
        S1extended = cat(3,cat(3,zeros(size(S1,1),size(S1,2),N-1),S1nc),zeros(size(S1,1),size(S1,2),N-1));

        V3 = s2t(S1extended,2*N-1)*s2v(S2nc);
        S3 = v2s(V3,2*N-1,1);
    else
        T3=s2t(S1,N)*s2t(S2,N);
        S3 = t2s(T3,N);
    end

end
end


%   SIGNAL_OP Package: Version 1.2
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu