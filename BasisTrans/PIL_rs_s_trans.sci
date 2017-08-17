// **** Purpose ****
// generate the unitary transformaton form relasitivsitc spherical to spherical
// **** Variables ****
// [L]: 1x1, integer
// <= orbital quantum number
// [S]: 1x1, integer or half-integer
// <= spin quantum number
// [U_rs_s]: n x n, real or complex
// => unitary transformation from 'rs' basis to 's' 
// **** Version ****
// 05/01/2014
// **** Comment ****
// This program transforms |J,Mj> into |ML,MS>==> U_jm*|L,S,J,MJ>=|L,S,MS,ML>
// Note: this is not biparticle problem. This is just single particle transfomation.
// The state order: 
// |L,S,J,MJ> => J: |L-S| ~ L+S, ML: -J ~ +J
// |L,S,MS,ML> => MS:-S ~ +S, ML: -L ~ +L

function [U_rs_s]=PIL_rs_s_trans(L,S)
    if L<0 | S<0 then
        disp('L, S must positive!');
    end
    J_list=linspace(abs(L-S),L+S,L+S-abs(L-S)+1);
    // construct U_rs_s
    U_rs_s=zeros((2*L+1)*(2*S+1),(2*L+1)*(2*S+1));
    state_list=zeros((2*L+1)*(2*S+1));
    count_c=0;
    for n1=1:length(J_list) // run all J=L+S
        J=J_list(n1);
        for n2=1:2*J+1 // run all MJ
            MJ=n2-1-J;
            count_c=count_c+1;
            count_r=0;
            state_list(count_c,1:2)=[J,MJ];
            for n3=1:2*S+1 // run all MS
                for n4=1:2*L+1 // run all ML
                    count_r=count_r+1;
                    MS=n3-1-S;
                    ML=n4-1-L;
                    if count_c==1 then
                        state_list(count_r,3:4)=[ML,MS];
                    end
                    if ML+MS==MJ then
                        U_rs_s(count_r,count_c)=PIL_cg_coff(L,S,ML,MS,J,MJ);
                    end
                end
            end

        end
    end
endfunction

