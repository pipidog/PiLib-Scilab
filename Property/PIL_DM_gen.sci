// **** Purpose ****
// generates the density matrix  with T specitied.
// **** Variables ****
// [E_level]: n x 1, real
// <= the energy levels
// [temperature]: 1x1, real > 0
// <= the temperature in Kelvin
// [Ef]: 1x1, real
// <= the Fermi level 
// [tot_e]: 1x1, integer >0
// <= the total charge
// [V_trans]: n x n, complex, defualt=eye(n,n)
// <= the transformation matrix of diagonalizing H(k)
// [DM]: n x n, complex
// => the density matrix
// **** Version ****
// 05/01/2014 first built
// 06/02/2014 bug fix: error when tot_e=0
// **** Comment ****
// This function put your input E_levels in diagonal form with Fermi-Dirac
// weighting. Then transform the DM to your original basis. Therefore, the
// input E_level can be at single k-poing ( |D(q)><D(q)| ) or in real space
// with all space sites (|D><D|) 

function DM=PIL_DM_gen(E_level,temperature,Ef,tot_e,V_trans)
    if abs(tot_e) >= 10^-3 then
        k_Bolt=(25.6/298)/100;
        E_level=gsort(real(E_level),'g','i')
        [lhs,rhs]=argn();
        if rhs==4
            V_trans=eye(length(E_level),length(E_level));
        end
        DM=diag((exp((E_level-Ef)/(k_Bolt*temperature))+1).^(-1));
        DM=DM*(tot_e/sum(diag(DM))); // useful for degenerate case!
        DM=V_trans*DM*V_trans';
    else 
        DM=diag(zeros(E_level));
    end
endfunction
