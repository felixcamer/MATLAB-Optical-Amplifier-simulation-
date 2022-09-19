function [ind]=findech(tt,t)
%--------------------------------------------------------------------------
% Cette fonction recherche l'indice (ind) du vecteur temps (tt) 
% dont la valeur est la plus proche de la valeur (t)
% tt (ps)
% t (ps)
% ind (sd)
%--------------------------------------------------------------------------

for ii=1:size(tt,2)
   if tt(1,ii)<t
       ind=ii;
   end       
end    

end

