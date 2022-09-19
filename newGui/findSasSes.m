function [ Sas,Ses ] = findSasSes(lbd_s,lambda,S_as,S_es)


for ii=1:size(lambda,2)
    if lambda(1,ii)<=lbd_s
        %inc=ii
        Sas=S_as(1,ii+1); % (m2)	Section efficace d’absorption de la sonde 
        Ses=S_es(1,ii+1); % (m2)	Section efficace d’émission de la sonde
    end
    
end


end

