// **** Purpose ****
// calculate the chemcial potential of a given energy spectrum
// **** Variables ****
// [E_level]: n x 1, real 
// <= the full energy spectrum of discrete DOS
// [e_tot]: 1x1, integer
// <= the total electrons
// [temp]: 1x1, real > 0
// <= the temperature
// [accuracy]: 1x1, real >0, default=10^-6
// <= the accuracy of the chemcial potential
// [Ef]: 1x1, real
// => the chemical potential, in unit of eV
// **** Version ****
// 05/01/2014
// **** Comment ****

function Ef=PIL_Ef(E_level,e_tot,temp,accuracy)
    [lhs,rhs]=argn();
    if rhs==3
        accuracy=1e-4;
    end
    E_level=gsort(real(E_level(:)),'g','i')
    
    // round e_tot
    if abs(round(e_tot)-e_tot) >= 0.01 then
        disp('Warning! e_tot has been rounded to integer!');
    end
    e_tot=round(e_tot);
    
    // check e_tot consistency
    if e_tot >=length(E_level)
        disp('Error! e_tot >= E_level, no Fermi level!');
        abort        
    end
    
    k_Bolt=(25.6/298)/100; //298*k_bolt=25.6meV (Wiki "kT(energy)")
    E_upper=E_level(e_tot+1);
    E_lower=E_level(e_tot);
    E_error=1;
    count=0
    while E_error > accuracy
        count=count+1;
        Ef=(E_upper+E_lower)/2
        n_E=sum((exp((E_level-Ef)/(k_Bolt*temp))+1).^(-1));
        if n_E > e_tot
            E_error=Ef-E_lower;
            E_upper=Ef;
            Ef=(E_upper+E_lower)/2;
        elseif n_E < e_tot
            E_error=E_upper-Ef;
            E_lower=Ef;
            Ef=(E_upper+E_lower)/2;
        else 
            break;
        end
        if count >= 15
            disp('Error: PIL_Ef, Ef cannot converge!');
        end
    end
endfunction

//examples:
//E_level=zeros(1,1000);
//E_level(1:500)=linspace(3,0.5,500);
//E_level(501:1000)=linspace(-0.9,-3,500);
//e_tot=500;
//temp=300;
//Ef=mylib_Ef(E_level,500,300);
//Result:
//Ef=-0.100055
