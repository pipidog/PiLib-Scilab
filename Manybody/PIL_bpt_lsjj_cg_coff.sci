// **** Purpose ****
// calculate the two-particle CG coefficients of LS and JJ coupling
// **** Variables ****
// [L1,S1,ML1,MS1,L2,S2,ML2,MS2,J1,J2,J,MJ]: integer or half-integer
// <= the quantum numbers
// [coup]: char, 'ls', 'jj'
// <= specify the type of coupling
// [CG]: 1x1, real or complex
// => cg coefficient
// **** Version ****
// 05/01/2014
// **** Comment ****
// This program gives you the cg cofficients of LS or JJ coulping:
//
// LS: L=L1+L2,S=S1+S2 => J=L+S
// <L1,S1,ML1,MS1;L2,S2,ML2,MS2|L,S,J,MJ>
// =C(L,S,ML,MS|J,MJ)*C(L1,L2,ML1,ML2|L,ML)*C(S1,S2,MS1,MS2|S,MS)
//
// JJ: J1=L1+S1, J2=L2+S2 => J=J1+J2
// <L1,S1,ML1,MS1;L2,S2,ML2,MS2|J1,J2,J,MJ>
// =C(J1,J2,MJ1,MJ2|J,MJ)*C(L1,S1,ML1,MS1|J1,MJ1)*C(L2,S2,ML2,MS2|J2,MJ2)

function [CG]=PIL_bpt_lsjj_cg_coff(L1,S1,ML1,MS1,L2,S2,ML2,MS2,J1,J2,J,MJ,coup)
    if sign(L1)<0 | sign(L2)<0 | sign(S1)<0 | sign(S2)<0 | sign(J1)<0 | sign(J2) <0 | sign(J) <0 then
        disp('Error! L1,L2,S1,S2,J1(L),J2(S),J must be non-negative!');
        pause;
    end
    CG=0;
    for MJ1=-J1:+J1
        for MJ2=-J2:+J2
            select coup
            case 'ls'
                if MJ1+MJ2==MJ & ML1+ML2==MJ1 & MS1+MS2==MJ2 then
                    CG=CG+PIL_cg_coff(J1,J2,MJ1,MJ2,J,MJ)*PIL_cg_coff(L1,L2,ML1,ML2,J1,MJ1)*PIL_cg_coff(S1,S2,MS1,MS2,J2,MJ2);
                end
            case 'jj'
                if MJ1+MJ2==MJ & ML1+MS1==MJ1 & ML2+MS2==MJ2 then
                    CG=CG+PIL_cg_coff(J1,J2,MJ1,MJ2,J,MJ)*PIL_cg_coff(L1,S1,ML1,MS1,J1,MJ1)*PIL_cg_coff(L2,S2,ML2,MS2,J2,MJ2);
                end 
            end
        end
    end
endfunction
