// **** Purpose ****
// calculates the Fermion creation or annihilation operator
// **** Variables ****
// [F_op]: 2 x total operation, integer
// <= the fermion operatros (see comment)
// [Fock_state_in]: 1 x total single particle states, integer
// <= the inpur Fock state, elements should be only 0 or 1
// [Fock_phase_out]: 1x1, -1 or +1
// => the phase of the Fock_state_out
// [Fock_state_out]: 1 x total single particle states, integer
// => the output Fock state  
// **** Version ****
// 05/01/2014
// **** Comment ****
// the input operation is a Nx2 matrix: [2,+1;1,-1;3,+1;5,-1] which represents
// C_{5}C_{3}^{+}C_{1}C_{2}^{+}  (Note: order reversed!!!)

function [Fock_phase_out,Fock_state_out]=PIL_Fermion_op(F_op,Fock_state_in)
    // check operation has correct form
    if find(F_op(:,2)~=1 & F_op(:,2)~=-1)~=[] then
        disp('Error! Fermion operator can be +1 or -1 on 2nd column!');
        break;
    end
    if max(F_op(:,1))>length(Fock_state(1,:)) | min(F_op(:,1))<0  then
        disp('Error! State label not in range!');
        break;
    end

    tot_op=length(F_op(:,1));
    Fock_state_out=Fock_state_in;
    Fock_phase_out=1;
    for n=1:tot_op
        select F_op(n,2)
        case +1
            if Fock_state_out(F_op(n,1))==0 then
                Fock_state_out(F_op(n,1))=1;
                if F_op(n,1)~=1 then
                    Fock_phase_out=(-1)^sum(Fock_state_out(1:F_op(n,1)-1))*Fock_phase_out;
                end 
            else
                Fock_state_out=zeros(Fock_state_out);
                Fock_phase_out=0;
                break;
            end
        case -1
            if Fock_state_out(F_op(n,1))==1 then
                Fock_state_out(F_op(n,1))=0;
                if F_op(n,1)~=1 then
                    Fock_phase_out=(-1)^sum(Fock_state_out(1:F_op(n,1)-1))*Fock_phase_out;
                end 
            else
                Fock_state_out=zeros(Fock_state_out);
                Fock_phase_out=0;
                break;
            end
        end
    end
endfunction
