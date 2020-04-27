function result = fermi(E,Ef)
 % Fermi dirac distribition at T=0 K
if E>Ef
    result=0;
else if E==Ef;
        result=1/2;
    else 
        result=1;
    end 
end 
end