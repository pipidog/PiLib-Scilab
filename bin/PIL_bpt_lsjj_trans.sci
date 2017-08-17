// **** Purpose ****
// calculate the two-particle CG coefficients of LS and JJ coupling
// **** Variables ****
// [L1,S1,L2,S2]: integer or half-integer
// <= the quantum numbers
// [coup]: char, 'ls', 'jj'
// <= specify the type of coupling
// [U_jm]: nxn, complex
// => unitary transformation, U_jm*|J1,J2,J,MJ> -> |L1,S1,ML1,MS1;L2,S2,ML2,MS2>
// [state_list]: nx8, integer or half-integer
// => recrods the order of the state of the unitary transformation [L,S,J,MJ,ML1,MS1,ML2,MS2]
// **** Version ****
// 05/01/2014
// **** Comment ****
// This code generates the basis transformation matrix between LS or JJ coupling
// basis and the direct product two particle basis. The unitary matrix satisifies
// the relation: U_jm*|J1,J2,J,MJ> -> |L1,S1,ML1,MS1;L2,S2,ML2,MS2>
// state_list records [L,S,J,MJ,ML1,MS1,ML2,MS2]

function [U_jm,state_list]=PIL_bpt_lsjj_trans(L1,S1,L2,S2,coup)
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
    state_list=zeros(dim(1),8);
    count_j=0;
    for n1=1:length(J1_list) //J1 or L
        for n2=1:length(J2_list) //J2 or S
            J1=J1_list(n1);
            J2=J2_list(n2);
            J_list=linspace(abs(J1-J2),J1+J2,J1+J2-abs(J1-J2)+1);
            for n3=1:length(J_list) // J
                J=J_list(n3)
                for n4=1:2*J+1 //MJ
                    MJ=n4-1-J;
                    count_j=count_j+1;
                    count_m=0; 
                    state_list(count_j,1:4)=[J1,J2,J,MJ];
                    // start expansion loop ======================================== 
                    for n5=1:2*L1+1 // ML1 
                        for n6=1:2*S1+1 // MS1
                            for n7=1:2*L2+1 // ML2
                                for n8=1:2*S2+1 // MS2
                                    ML1=n5-1-L1;
                                    MS1=n6-1-S1;
                                    ML2=n7-1-L2;
                                    MS2=n8-1-S2;
                                    count_m=count_m+1;
                                    state_list(count_m,5:8)=[ML1,MS1,ML2,MS2];
                                    U_jm(count_m,count_j)=PIL_bpt_lsjj_cg_coff(L1,S1,ML1,MS1,L2,S2,ML2,MS2,J1,J2,J,MJ,coup);
                                end
                            end
                        end
                    end
                    // end expansion loop ======================================
                end
            end
        end
    end
endfunction
