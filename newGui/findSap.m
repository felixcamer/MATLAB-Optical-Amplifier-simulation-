function [ Sap ] = findSap(lbd_p,lambda,S_ap)


for ii=1:size(lambda,2)
    if lambda(1,ii)<=lbd_p
        %inc=ii
        Sap=S_ap(1,ii+1); % (m2)	Section efficace d’absorption de la pompe 
        
    end
    
end


end

