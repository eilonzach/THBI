function [ allpair_invCe ] = pairwise_Ce( nstas,r,sigma )
%[ allpair_invCe ] = pairwise_Ce( nstas,r,sigma )
%   Detailed explanation goes here

npair = handshake(nstas);
allpair_invCe = struct();

kk = 0;
for is1 = 1:nstas
    for is2 = is1+1:nstas
        kk = kk+1;
        
        
    ruse = mean(r([is1,is2]));
    [ allpair_invCe(kk).invCe ] = calc_Ce( sigma,ruse,length(dat1) );


    end
end

end

