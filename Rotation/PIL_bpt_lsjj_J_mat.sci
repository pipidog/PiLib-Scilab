// **** Purpose ****
// generates the two-particle Jx, Jy, Jz matrix with LS or JJ coupling 
// **** Variables ****
// [L1],[S1],[L2],[S2]: 1x1, integer or half-integer
// <= the quantum number
// [coup]: 1x1, string, 'jj' or 'ls'
// <= the coupling type
// [Jx_j],[Jy_j],[Jz_j]: n x n, real
// => the angular momentum matrix in spin-orbital coupled basis
// [state_list]: n x 4, integer
// => the state list
// **** Version ****
// 05/01/2014
// **** Comment ****
// For a given [L1,S1,L2,S2] in LS coupling or JJ coupling, they forms different
// possible J. This code generates the Jx(J), Jy(J), Jz(J) with the Jz(J) diagonal
// basis and direct sum them. To get Jx, Jy, Jz in direct product basis, use:
// [U_jm,state_list]=PIL_ls_jj_trans(L1,S1,L2,S2,coup)
// Jx_m=U_jm*Jx_j*U_jm'; and so on. 

function[Jx_j,Jy_j,Jz_j,state_list]=PIL_bpt_lsjj_J_mat(L1,S1,L2,S2,coup)
select coup
case 'jj'
    J1_list=linspace(abs(L1-S1),L1+S1,L1+S1-abs(L1-S1)+1); // J1
    J2_list=linspace(abs(L2-S2),L2+S2,L2+S2-abs(L2-S2)+1); // J2
case 'ls'
    J1_list=linspace(abs(L1-L2),L1+L2,L1+L2-abs(L1-L2)+1); // L
    J2_list=linspace(abs(S1-S2),S1+S2,S1+S2-abs(S1-S2)+1); // S
end

dim=[(2*L1+1)*(2*S1+1)*(2*L2+1)*(2*S2+1),(2*L1+1)*(2*S1+1)*(2*L2+1)*(2*S2+1)];
U_jm=zeros(dim(1),dim(2));
state_list=zeros(dim(1),4);
count_ini=0;
count_end=0;
for n1=1:length(J1_list) //J1 or L
    for n2=1:length(J2_list) //J2 or S
        J1=J1_list(n1);
        J2=J2_list(n2);
        J_list=linspace(abs(J1-J2),J1+J2,J1+J2-abs(J1-J2)+1);
        for n3=1:length(J_list) // J
            J=J_list(n3)
            count_ini=count_end+1;
            count_end=count_end+2*J+1;
            [jx,jy,jz]=PIL_J_mat(J);
            Jx_j(count_ini:count_end,count_ini:count_end)=jx;
            Jy_j(count_ini:count_end,count_ini:count_end)=jy;
            Jz_j(count_ini:count_end,count_ini:count_end)=jz;
            state_list(count_ini:count_end,1:3)=repmat([J1,J2,J],(count_end-count_ini+1),1);
            state_list(count_ini:count_end,4)=[-J:+J]';
        end
    end
end
endfunction
