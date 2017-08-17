// **** Purpose ****
// Slaster-Koster Table
// **** Variables ****
// [orb1]: 1x1, integer, 1~32
// <= orbital index 1 in cubic harmonics basis
// [orb2]: 1x1, integer, 1~32
// <= orbital index 2 in cubic harmonics basis
// [r]: 1x3, real
// <= position vector between two ions
// [SK_paremeter]: 1x4, real
// <= the Slaster-Koster paremeter sigma, pi, delta, phi
// [SK_int]: 1x1, real
// => Slaster-Koster value
// **** Version ****
// 05/11/2014 first built
// 05/20/2014 use consistent notations (with spin) for orb1 and orb2
// **** Comment ****
// 1. orb1 and orb2 are defined using the PiLib_basis_trans
//      check the manual file. it contains spin, range 1~32
// 2. in-code spinless state label:
//      <s> s:1 
//      <p> x:2,y:3,z:4 
//      <d> xy:5, yz:6, xz:7, x2-y2:8, 3z2-r2:9 
//      <f> xyz:10, x(5x2-3r2):11, y(5y2-3r2):12, z(5z2-3r2):13, x(y2-z2):14, y(z2-x2):15, z(x2-y2):16 
//      Note: this order is different from PIL_basis_trans! 
// 3. SK_parameter
//      there are four SK parameters, sigma, pi, delta, phi
//      altohugh not all of them are needed for all calculation, this code still follows
// a four-paremeter formate. So you can leave the uncessary paremeters zero. They
//      will not be used in the calculation. 
// 4. See User manual "Slaster-Koster Table"

function SK_int=PIL_SK_int(orb1,orb2,r,SK_parameter)
    // spin to no-spin index    
    function val=orb_spin2nspin(orb)
        if orb==0 then
            val=0
        elseif orb==1 | orb==2 
            val=1
        elseif orb==3 | orb==6
            val=2
        elseif orb==4 | orb==7
            val=3
        elseif orb==5 | orb==8
            val=4
        elseif orb==9 | orb==14
            val=5
        elseif orb==10 | orb==15
            val=6
        elseif orb==11 | orb==16
            val=7
        elseif orb==12 | orb==17
            val=8
        elseif orb==13 | orb==18
            val=9
        elseif orb==19 | orb==26
            val=10
        elseif orb==20 | orb==27
            val=11
        elseif orb==21 | orb==28
            val=12
        elseif orb==22 | orb==29
            val=13
        elseif orb==23 | orb==30
            val=14
        elseif orb==24 | orb==31
            val=15
        elseif orb==25 | orb==32
            val=16
        end
    endfunction
    
    // Define Calculateed Slaster-Koster Table
    function SK_int=SK_table(SK_index,method)
        select method
        case 'n'
            L=r(1)/norm(r); M=r(2)/norm(r); N=r(3)/norm(r);
        case 'r'
            L=r(2)/norm(r); M=r(3)/norm(r); N=r(1)/norm(r);
        case 'l'
            L=r(3)/norm(r); M=r(1)/norm(r); N=r(2)/norm(r);
        end

        select SK_index
        case [1,1]    
            SK_int=ts;
        case [1,2]
            SK_int=L*ts
        case [1,5]    
            SK_int=sqrt(3)*L*M*ts
        case [1,8]    
            SK_int=sqrt(3)/2*(L^2-M^2)*ts
        case [1,9]    
            SK_int=(N^2-(1/2)*(L^2+M^2))*ts  
        case [1,10]    
            SK_int=sqrt(15)*L*M*N*ts
        case [1,11]    
            SK_int=(1/2)*L*(5*L^2-3)*ts 
        case [1,14]    
            SK_int=(1/2)*sqrt(15)*L*(M^2-N^2)*ts 
        case [2,2]    
            SK_int=L^2*ts+(1-L^2)*tp
        case [2,3]
            SK_int=L*M*ts-L*M*tp
        case [2,4]
            SK_int=L*N*ts-L*N*tp
        case [2,5]
            SK_int=sqrt(3)*L^2*M*ts+M*(1-2*L^2)*tp
        case [2,6]
            SK_int=sqrt(3)*L*M*N*ts-2*L*M*N*tp
        case [2,7]
            SK_int=sqrt(3)*L^2*N*ts+N*(1-2*L^2)*tp
        case [2,8]
            SK_int=sqrt(3)/2*L*(L^2-M^2)*ts+L*(1-L^2+M^2)*tp
        case [2,9]
            SK_int=L*(N^2-(1/2)*(L^2+M^2))*ts-sqrt(3)*L*N^2*tp
        case [2,10]
            SK_int=sqrt(15)*L^2*M*N*ts-sqrt(5/2)*(3*L^2-1)*M*N*tp
        case [2,11]
            SK_int=(1/2)*L^2*(5*L^2-3)*ts-sqrt(3/8)*(5*L^2-1)*(L^2-1)*tp
        case [2,12]
            SK_int=(1/2)*L*M*(5*M^2-3)*ts-sqrt(3/8)*L*M*(5*M^2-1)*tp
        case [2,13]
            SK_int=(1/2)*L*N*(5*N^2-3)*ts-sqrt(3/8)*L*N*(5*N^2-1)*tp
        case [2,14]
            SK_int=(1/2)*sqrt(15)*L^2*(M^2-N^2)*ts-sqrt(5/8)*(3*L^2-1)*(M^2-N^2)*tp
        case [2,15]
            SK_int=(1/2)*sqrt(15)*L*M*(N^2-L^2)*ts-sqrt(5/8)*L*M*(3*(N^2-L^2)+2)*tp
        case [2,16]
            SK_int=(1/2)*sqrt(15)*L*N*(L^2-M^2)*ts-sqrt(5/8)*L*N*(3*(L^2-M^2)-2)*tp        
        case [3,8]    
            SK_int=sqrt(3)/2*M*(L^2-M^2)*ts-M*(1+L^2-M^2)*tp
        case [3,9]    
            SK_int=M*(N^2-(1/2)*(L^2+M^2))*ts-sqrt(3)*M*N^2*tp        
        case [4,8]
            SK_int=sqrt(3)/2*N*(L^2-M^2)*ts-N*(L^2-M^2)*tp
        case [4,9]
            SK_int=N*(N^2-(1/2)*(L^2+M^2))*ts+sqrt(3)*N*(L^2+M^2)*tp        
        case [5,5]
            SK_int=3*L^2*M^2*ts+(L^2+M^2-4*L^2*M^2)*tp+(N^2+L^2*M^2)*td
        case [5,6]
            SK_int=3*L*M^2*N*ts+L*N*(1-4*M^2)*tp+L*N*(M^2-1)*td
        case [5,7]
            SK_int=3*L^2*M*N*ts+M*N*(1-4*L^2)*tp+M*N*(L^2-1)*td
        case [5,8]
            SK_int=(3/2)*L*M*(L^2-M^2)*ts+2*L*M*(M^2-L^2)*tp+(1/2)*L*M*(L^2-M^2)*td
        case [5,9]
            SK_int=sqrt(3)*L*M*(N^2-(1/2)*(L^2+M^2))*ts-2*sqrt(3)*L*M*N^2*tp+(1/2)*sqrt(3)*L*M*(1+N^2)*td
        case [5,10]
            SK_int=sqrt(45)*L^2*M^2*N*ts-sqrt(5/2)*N*(6*L^2*M^2+N^2-1)*tp+N*(3*L^2*M^2+2*N^2-1)*td
        case [5,11]
            SK_int=(1/2)*sqrt(3)*L^2*M*(5*L^2-3)*ts-sqrt(3/8)*M*(5*L^2-1)*(2*L^2-1)*tp+(1/2)*sqrt(15)*L^2*M*(L^2-1)*td
        case [5,12]
            SK_int=(1/2)*sqrt(3)*L*M^2*(5*M^2-3)*ts-sqrt(3/8)*L*(5*M^2-1)*(2*M^2-1)*tp+(1/2)*sqrt(15)*L*M^2*(M^2-1)*td
        case [5,13]
            SK_int=(1/2)*sqrt(3)*L*M*N*(5*N^2-3)*ts-sqrt(3/2)*L*M*N*(5*N^2-1)*tp+(1/2)*sqrt(15)*L*M*N*(N^2-1)*td
        case [5,14]
            SK_int=(3/2)*sqrt(5)*L^2*M*(M^2-N^2)*ts-sqrt(5/8)*M*((6*L^2-1)*(M^2-N^2)-2*L^2)*tp+(1/2)*M*(3*L^2*(M^2-N^2)-4*N^2-2*L^2)*td
        case [5,15]
            SK_int=(3/2)*sqrt(5)*L*M^2*(N^2-L^2)*ts-sqrt(5/8)*L*((6*M^2-1)*(N^2-L^2)+2*M^2)*tp+(1/2)*L*(3*M^2*(N^2-L^2)-4*N^2+2*M^2)*td
        case [5,16]    
            SK_int=(3/2)*sqrt(5)*L*M*N*(L^2-M^2)*ts-3*sqrt(5/2)*L*M*N*(L^2-M^2)*tp+(3/2)*L*M*N*(L^2-M^2)*td        
        case [6,8]
            SK_int=(3/2)*M*N*(L^2-M^2)*ts-M*N*(1+2*(L^2-M^2))*tp+M*N*(1+(1/2)*(L^2-M^2))*td
        case [6,9]
            SK_int=sqrt(3)*M*N*(N^2-(1/2)*(L^2+M^2))*ts+sqrt(3)*M*N*(L^2+M^2-N^2)*tp-(1/2)*sqrt(3)*M*N*(L^2+M^2)*td        
        case [7,8]
            SK_int=(3/2)*N*L*(L^2-M^2)*ts+N*L*(1-2*(L^2-M^2))*tp-N*L*(1-(1/2)*(L^2-M^2))*td
        case [7,9]
            SK_int=sqrt(3)*L*N*(N^2-(1/2)*(L^2+M^2))*ts+sqrt(3)*L*N*(L^2+M^2-N^2)*tp-(1/2)*sqrt(3)*L*N*(L^2+M^2)*td        
        case [8,8]
            SK_int=(3/4)*(L^2-M^2)^2*ts+(L^2+M^2-(L^2-M^2)^2)*tp+(N^2+(1/4)*(L^2-M^2)^2)*td
        case [8,9]
            SK_int=(1/2)*sqrt(3)*(L^2-M^2)*(N^2-(1/2)*(L^2+M^2))*ts+sqrt(3)*N^2*(M^2-L^2)*tp+(1/4)*sqrt(3)*(1+N^2)*(L^2-M^2)*td
        case [8,10]
            SK_int=(3/2)*sqrt(5)*L*M*N*(L^2-M^2)*ts-3*sqrt(5/2)*L*M*N*(L^2-M^2)*tp+(3/2)*L*M*N*(L^2-M^2)*td
        case [8,11]
            SK_int=(1/4)*sqrt(3)*L*(L^2-M^2)*(5*L^2-3)*ts-sqrt(3/8)*L*(L^2-M^2-1)*(5*L^2-1)*tp-(1/4)*sqrt(15)*L*((L^2-M^2)*(1-L^2)-2*N^2)*td
        case [8,12]
            SK_int=(1/4)*sqrt(3)*M*(L^2-M^2)*(5*M^2-3)*ts-sqrt(3/8)*M*(L^2-M^2+1)*(5*M^2-1)*tp-(1/4)*sqrt(15)*M*((L^2-M^2)*(1-M^2)+2*N^2)*td
        case [8,13]
            SK_int=(1/4)*sqrt(3)*N*(L^2-M^2)*(5*N^2-3)*ts-sqrt(3/8)*N*(L^2-M^2)*(5*N^2-1)*tp+(1/4)*sqrt(15)*N*(N^2+1)*(L^2-M^2)*td
        case [8,14] 
            SK_int=(3/4)*sqrt(5)*L*(L^2-M^2)*(M^2-N^2)*ts-sqrt(5/8)*L*(3*(L^2-M^2)*(M^2-N^2)-L^2+1)*tp-(1/4)*L*(3*(L^2-M^2)*(M^2-N^2)-4*L^2+2)*td
        case [8,15]
            SK_int=(3/4)*sqrt(5)*M*(L^2-M^2)*(N^2-L^2)*ts-sqrt(5/8)*M*(3*(L^2-M^2)*(N^2-L^2)-M^2+1)*tp+(1/4)*M*(3*(L^2-M^2)*(N^2-L^2)-4*M^2+2)*td
        case [8,16]
            SK_int=(3/4)*sqrt(5)*N*(L^2-M^2)^2*ts-sqrt(5/8)*N*(3*(L^2-M^2)^2+2*N^2-2)*tp+(1/4)*N*(3*(L^2-M^2)^2+8*N^2-4)*td        
        case [9,9]
            SK_int=(N^2-(1/2)*(L^2+M^2))^2*ts+3*N^2*(L^2+M^2)*tp+(3/4)*(L^2+M^2)*td
        case [9,10]
            SK_int=(1/2)*sqrt(15)*L*M*N*(3*N^2-1)*ts-sqrt(15/2)*L*M*N*(3*N^2-1)*tp+(1/2)*sqrt(3)*L*M*N*(3*N^2-1)*td
        case [9,11]
            SK_int=(1/4)*L*(3*N^2-1)*(5*L^2-3)*ts-(3/4)*sqrt(2)*L*N^2*(5*L^2-1)*tp+(3/4)*sqrt(5)*L*(L^2*N^2-M^2)*td
        case [9,12]
            SK_int=(1/4)*M*(3*N^2-1)*(5*M^2-3)*ts-(3/4)*sqrt(2)*M*N^2*(5*M^2-1)*tp+(3/4)*sqrt(5)*M*(M^2*N^2-L^2)*td
        case [9,13]
            SK_int=(1/4)*N*(3*N^2-1)*(5*N^2-3)*ts-(3/4)*sqrt(2)*N*(5*N^2-1)*(N^2-1)*tp+(3/4)*sqrt(5)*N*(N^2-1)^2*td
        case [9,14]
            SK_int=(1/4)*sqrt(15)*L*(M^2-N^2)*(3*N^2-1)*ts-sqrt(15/8)*L*N^2*(3*(M^2-N^2)+2)*tp+(1/4)*sqrt(3)*L*((3*N^2-1)*(M^2-N^2)-4*L^2+2)*td
        case [9,15]
            SK_int=(1/4)*sqrt(15)*M*(N^2-L^2)*(3*N^2-1)*ts-sqrt(15/8)*M*N^2*(3*(N^2-L^2)-2)*tp+(1/4)*sqrt(3)*M*((3*N^2-1)*(N^2-L^2)+4*M^2-2)*td
        case [9,16]
            SK_int=(1/4)*sqrt(15)*N*(3*N^2-1)*(L^2-M^2)*ts-sqrt(15/8)*N*(3*N^2-1)*(L^2-M^2)*tp+(1/4)*sqrt(3)*N*(3*N^2-1)*(L^2-M^2)*td        
        case [10,10]
            SK_int=15*L^2*M^2*N^2*ts+(5/2)*(L^2*M^2+M^2*N^2+N^2*L^2-9*L^2*M^2*N^2)*tp...
            +(1-4*(L^2*M^2+M^2*N^2+N^2*L^2)+9*L^2*M^2*N^2)*td+(3/2)*(1-L^2)*(1-M^2)*(1-N^2)*tf
        case [10,13]
            SK_int=(1/2)*sqrt(15)*L*M*N^2*(5*N^2-3)*ts-(1/4)*sqrt(15)*L*M*(3*N^2-1)*(5*N^2-1)*tp...
            +(1/2)*sqrt(15)*L*M*N^2*(3*N^2-1)*td+(1/4)*sqrt(15)*L*M*(1-N^4)*tf
        case [10,16]
            SK_int=(15/2)*L*M*N^2*(L^2-M^2)*ts-(5/4)*L*M*(L^2-M^2)*(9*N^2-1)*tp...
            +(1/2)*L*M*(L^2-M^2)*(9*N^2-4)*td+(3/4)*L*M*(L^2-M^2)*(1-N^2)*tf        
        case [11,13]
            SK_int=(1/4)*L*N*(5*L^2-3)*(5*N^2-3)*ts-(3/8)*L*N*(5*L^2-1)*(5*N^2-1)*tp...
            +(15/4)*L*N*(L^2*N^2-M^2)*td+(5/8)*L*N*(3*M^2-L^2*N^2)*tf
        case [11,16]
            SK_int=(1/4)*sqrt(15)*L*N*(L^2-M^2)*(5*L^2-3)*ts+(1/8)*sqrt(15)*L*N*(2-3*(L^2-M^2))*(5*L^2-1)*tp...
            +(1/4)*sqrt(15)*L*N*(3*(1+L^2)*(L^2-M^2)-8*L^2+2)*td+(1/8)*sqrt(15)*L*N*(-(L^2+3)*(L^2-M^2)+6*L^2-2)*tf        
        case [12,13]
            SK_int=(1/4)*M*N*(5*M^2-3)*(5*N^2-3)*ts-(3/8)*M*N*(5*M^2-1)*(5*N^2-1)*tp...
            +(15/4)*M*N*(M^2*N^2-L^2)*td+(5/8)*M*N*(3*L^2-M^2*N^2)*tf
        case [12,16]    
            SK_int=(1/4)*sqrt(15)*M*N*(L^2-M^2)*(5*M^2-3)*ts-(1/8)*sqrt(15)*M*N*(2+3*(L^2-M^2))*(5*M^2-1)*tp...
            +(1/4)*sqrt(15)*M*N*(3*(1+M^2)*(L^2-M^2)+8*M^2-2)*td-(1/8)*sqrt(15)*M*N*((M^2+3)*(L^2-M^2)+6*M^2-2)*tf        
        case [13,13]    
            SK_int=(1/4)*N^2*(5*N^2-3)^2*ts+(3/8)*(5*N^2-1)^2*(1-N^2)*tp+(15/4)*N^2*(1-N^2)^2*td+(5/8)*(1-N^2)^3*tf
        case [13,16]    
            SK_int=(1/4)*sqrt(15)*(L^2-M^2)*N^2*(5*N^2-3)*ts+(1/8)*sqrt(15)*(L^2-M^2)*(5*N^2-1)*(3*N^2-1)*tp...
            +(1/4)*sqrt(15)*(L^2-M^2)*N^2*(3*N^2-1)*td-(1/8)*sqrt(15)*(L^2-M^2)*(1-N^4)*tf        
        case [14,16]    
            SK_int=(15/4)*L*N*(L^2-M^2)*(M^2-N^2)*ts-(5/8)*L*N*(9*(L^2-M^2)*(M^2-N^2)-2*M^2+2)*tp...
            +(1/4)*L*N*(9*(L^2-M^2)*(M^2-N^2)-8*M^2+2)*td+(3/8)*L*N*(-(L^2-M^2)*(M^2-N^2)+2*M^2+2)*tf        
        case [15,16]    
            SK_int=(15/4)*M*N*(L^2-M^2)*(N^2-L^2)*ts-(5/8)*M*N*(9*(L^2-M^2)*(N^2-L^2)-2*L^2+2)*tp...
            +(1/4)*M*N*(9*(L^2-M^2)*(N^2-L^2)-8*L^2+2)*td+(3/8)*M*N*(-(L^2-M^2)*(N^2-L^2)+2*L^2+2)*tf     
        case [16,16]    
            SK_int=(15/4)*N^2*(L^2-M^2)^2*ts+(5/8)*(4*N^2*(1-N^2)+(L^2-M^2)^2*(1-9*N^2))*tp...
            +(1/4)*((L^2-M^2)^2*(9*N^2-4)+4*(1-2*N^2)^2)*td+(3/8)*(1-N^2)*((1+N^2)-4*L^2*M^2)*tf
        end
    endfunction
    // check input
    [lhs,rhs]=argn();
    if (rhs~=4) | (length(r)~=3) | (length(SK_parameter)~=4) then
        disp('Error: PIL_SK_int, input variable dimension wrong!');
    end
    // determin spin
    dn_state=[1,3:5,9:13,19:25];
    up_state=[2,6:8,14:18,26:32];
    if (length(find(dn_state==orb1))+length(find(up_state==orb2)))==2 then
        orb1=0; orb2=0;
    elseif (length(find(up_state==orb1))+length(find(dn_state==orb2)))==2
        orb1=0; orb2=0;
    else
        orb1=orb_spin2nspin(orb1);
        orb2=orb_spin2nspin(orb2)
    end
    
    // Define SK integrals by calling SK table
    if orb1 > orb2
        r=-r;
        orb=gsort([orb1,orb2],'g','i');
    else
        orb=[orb1,orb2];
    end    
    ts=SK_parameter(1); tp=SK_parameter(2); td=SK_parameter(3); tf=SK_parameter(4);
    
    select orb
    case [1,1]
        SK_int=SK_table([1,1],'n');
    case [1,2]
        SK_int=SK_table([1,2],'n');
    case [1,3]
        SK_int=SK_table([1,2],'r');
    case [1,4]
        SK_int=SK_table([1,2],'l');
    case [1,5]
        SK_int=SK_table([1,5],'n');
    case [1,6]
        SK_int=SK_table([1,5],'r');
    case [1,7]
        SK_int=SK_table([1,5],'l');
    case [1,8]
        SK_int=SK_table([1,8],'n');
    case [1,9]
        SK_int=SK_table([1,9],'n');
    case [1,10]
        SK_int=SK_table([1,10],'n');
    case [1,11]
        SK_int=SK_table([1,11],'n');
    case [1,12]
        SK_int=SK_table([1,11],'r');
    case [1,13]
        SK_int=SK_table([1,11],'l');
    case [1,14]
        SK_int=SK_table([1,14],'n');
    case [1,15]
        SK_int=SK_table([1,14],'r');
    case [1,16]
        SK_int=SK_table([1,14],'l');
    case [2,2]
        SK_int=SK_table([2,2],'n');
    case [2,3]
        SK_int=SK_table([2,3],'n');
    case [2,4]
        SK_int=SK_table([2,4],'n');
    case [2,5]
        SK_int=SK_table([2,5],'n');
    case [2,6]
        SK_int=SK_table([2,6],'n');
    case [2,7]
        SK_int=SK_table([2,7],'n');
    case [2,8]
        SK_int=SK_table([2,8],'n');
    case [2,9]
        SK_int=SK_table([2,9],'n');
    case [2,10]
        SK_int=SK_table([2,10],'n');
    case [2,11]
        SK_int=SK_table([2,11],'n');
    case [2,12]
        SK_int=SK_table([2,12],'n');
    case [2,13]
        SK_int=SK_table([2,13],'n');
    case [2,14]
        SK_int=SK_table([2,14],'n');
    case [2,15]
        SK_int=SK_table([2,15],'n');
    case [2,16]
        SK_int=SK_table([2,16],'n');
    case [3,3]
        SK_int=SK_table([2,2],'r');
    case [3,4]
        SK_int=SK_table([2,3],'r');
    case [3,5]
        SK_int=SK_table([2,4],'r');
    case [3,6]
        SK_int=SK_table([2,5],'r');
    case [3,7]
        SK_int=SK_table([2,6],'r');
    case [3,8]
        SK_int=SK_table([3,8],'n');
    case [3,9]
        SK_int=SK_table([3,9],'n');
    case [3,10]
        SK_int=SK_table([2,10],'r');
    case [3,11]
        SK_int=SK_table([2,13],'r');
    case [3,12]
        SK_int=SK_table([2,11],'r');
    case [3,13]
        SK_int=SK_table([2,12],'r');
    case [3,14]
        SK_int=SK_table([2,16],'r');
    case [3,15]
        SK_int=SK_table([2,14],'r');
    case [3,16]
        SK_int=SK_table([2,15],'r');
    case [4,4]
        SK_int=SK_table([2,2],'l');
    case [4,5]
        SK_int=SK_table([2,6],'l');
    case [4,6]
        SK_int=SK_table([2,7],'l');
    case [4,7]
        SK_int=SK_table([2,5],'l');
    case [4,8]
        SK_int=SK_table([4,8],'n');
    case [4,9]
        SK_int=SK_table([4,9],'n');
    case [4,10]
        SK_int=SK_table([2,10],'l');
    case [4,11]
        SK_int=SK_table([2,12],'l');
    case [4,12]
        SK_int=SK_table([2,13],'l');
    case [4,13]
        SK_int=SK_table([2,11],'l');
    case [4,14]
        SK_int=SK_table([2,15],'l');
    case [4,15]
        SK_int=SK_table([2,16],'l');
    case [4,16]
        SK_int=SK_table([2,14],'l');
    case [5,5]
        SK_int=SK_table([5,5],'n');
    case [5,6]
        SK_int=SK_table([5,6],'n');
    case [5,7]
        SK_int=SK_table([5,7],'n');
    case [5,8]
        SK_int=SK_table([5,8],'n');
    case [5,9]
        SK_int=SK_table([5,9],'n');
    case [5,10]
        SK_int=SK_table([5,10],'n');
    case [5,11]
        SK_int=SK_table([5,11],'n');
    case [5,12]
        SK_int=SK_table([5,12],'n');
    case [5,13]
        SK_int=SK_table([5,13],'n');
    case [5,14]
        SK_int=SK_table([5,14],'n');
    case [5,15]
        SK_int=SK_table([5,15],'n');
    case [5,16]
        SK_int=SK_table([5,16],'n');
    case [6,6]
        SK_int=SK_table([5,5],'r');
    case [6,7]
        SK_int=SK_table([5,6],'r');
    case [6,8]
        SK_int=SK_table([6,8],'n');
    case [6,9]
        SK_int=SK_table([6,9],'n');
    case [6,10]
        SK_int=SK_table([5,10],'r');
    case [6,11]
        SK_int=SK_table([5,13],'r');
    case [6,12]
        SK_int=SK_table([5,11],'r');
    case [6,13]
        SK_int=SK_table([5,12],'r');
    case [6,14]
        SK_int=SK_table([5,16],'r');
    case [6,15]
        SK_int=SK_table([5,14],'r');
    case [6,16]
        SK_int=SK_table([5,15],'r');
    case [7,7]
        SK_int=SK_table([5,5],'l');
    case [7,8]
        SK_int=SK_table([7,8],'n');
    case [7,9]
        SK_int=SK_table([7,9],'n');
    case [7,10]
        SK_int=SK_table([5,10],'l');
    case [7,11]
        SK_int=SK_table([5,12],'l');
    case [7,12]
        SK_int=SK_table([5,13],'l');
    case [7,13]
        SK_int=SK_table([5,11],'l');
    case [7,14]
        SK_int=SK_table([5,15],'l');
    case [7,15]
        SK_int=SK_table([5,16],'l');
    case [7,16]
        SK_int=SK_table([5,14],'l');
    case [8,8]
        SK_int=SK_table([8,8],'n');
    case [8,9]
        SK_int=SK_table([8,9],'n');
    case [8,10]
        SK_int=SK_table([8,10],'n');
    case [8,11]
        SK_int=SK_table([8,11],'n');
    case [8,12]
        SK_int=SK_table([8,12],'n');
    case [8,13]
        SK_int=SK_table([8,13],'n');
    case [8,14]
        SK_int=SK_table([8,14],'n');
    case [8,15]
        SK_int=SK_table([8,15],'n');
    case [8,16]
        SK_int=SK_table([8,16],'n');
    case [9,9]
        SK_int=SK_table([9,9],'n');
    case [9,10]
        SK_int=SK_table([9,10],'n');
    case [9,11]
        SK_int=SK_table([9,11],'n');
    case [9,12]
        SK_int=SK_table([9,12],'n'); 
    case [9,13]
        SK_int=SK_table([9,13],'n');
    case [9,14]
        SK_int=SK_table([9,14],'n');
    case [9,15]
        SK_int=SK_table([9,15],'n');
    case [9,16]
        SK_int=SK_table([9,16],'n');
    case [10,10]
        SK_int=SK_table([10,10],'n');
    case [10,11]
        SK_int=SK_table([10,13],'r');
    case [10,12]
        SK_int=SK_table([10,13],'l');
    case [10,13]
        SK_int=SK_table([10,13],'n');
    case [10,14]
        SK_int=SK_table([10,16],'r');
    case [10,15]
        SK_int=SK_table([10,16],'l');
    case [10,16]
        SK_int=SK_table([10,16],'n');
    case [11,11]
        SK_int=SK_table([13,13],'r');
    case [11,12]
        SK_int=SK_table([11,13],'r');
    case [11,13]
        SK_int=SK_table([11,13],'n'); 
    case [11,14]
        SK_int=SK_table([13,16],'r');
    case [11,15]
        SK_int=SK_table([12,16],'l');
    case [11,16]
        SK_int=SK_table([11,16],'n');
    case [12,12]
        SK_int=SK_table([13,13],'l');
    case [12,13]
        SK_int=SK_table([12,13],'n'); 
    case [12,14]
        SK_int=SK_table([11,16],'r');
    case [12,15]
        SK_int=SK_table([13,16],'l');
    case [12,16]
        SK_int=SK_table([12,16],'n');
    case [13,13]
        SK_int=SK_table([13,13],'n');
    case [13,14]
        SK_int=SK_table([12,16],'r');
    case [13,15]
        SK_int=SK_table([11,16],'l');
    case [13,16]
        SK_int=SK_table([13,16],'n');
    case [14,14]
        SK_int=SK_table([16,16],'r');
    case [14,15]
        SK_int=SK_table([15,16],'l');
    case [14,16]
        SK_int=SK_table([14,16],'n');
    case [15,15]
        SK_int=SK_table([16,16],'l');
    case [15,16]
        SK_int=SK_table([15,16],'n'); 
    case [16,16]
        SK_int=SK_table([16,16],'n');
    else
        SK_int=0;
    end
        
endfunction
