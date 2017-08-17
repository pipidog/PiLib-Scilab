// **** Purpose ****
// This function generates the lattice vectors to construct a conventional
// cell. 
// **** Variables ****
// [a_vec]: 3x3, real
// <= lattice row vectors 
// [new_b_red]: 3x3, real
// <= new reciprocal row vectors in reduced coordinate
// [new_a_red]: 3x3, real
// => new lattice row vectors in terms of a_vec 
// [orth_chk]: 1x3, real
// => the inner product of new lattice vectors between (1,2), (1,3), (2,3)
//    if a3 is supposed to be the perpendicular axis, then the last two 
//    should be close to zero. 
// **** Version ****
// 5/16/2016: first built
// **** Comment ****
// 1. How to use:
// 1). open xsf file using xcrysden, then open k-path selection. If in 
//     real space, there is an axis perpendicular to the other two, then
//     in k-space, the reciprocal lattice should have the same behavior.
//     Therefore, you can define a set of new axis in BZ where b3 is 
//     perperndicular to the other two. Then you can construct a slab
//     structure that has the same BZ as b1 and b2 but b3 is compressed.
// 2). The choosen new reciprocal lattice vector should be labeled using
//     reduced coordinate (usually high-symmetry points). Then convert 
//     these to new lattice vector in real space, you will be able to get
//     a conventional lattice. 
// 3). Once the conventional lattice vector is obtained, use PIL_conv_cell_vec
//     to obtain its corresponding sublattice

function [new_a_red,orth_chk]=PIL_conv_cell_gen(a_vec,new_b_red)
    b_cart=PIL_recip_vec(a_vec);
    new_b_cart=new_b_red*b_cart

    new_a_cart=PIL_recip_vec(new_b_cart);
    new_a_red=zeros(3,3);
    for n=1:3
        // construct new reduced lattice vector
        new_a_red(n,:)=clean(PIL_linexpan(new_a_cart(n,:),a_vec')',1e-4);
        gcd_factor=min(abs(new_a_red(n,abs(new_a_red(n,:))>1e-4)));
        new_a_red(n,:)=new_a_red(n,:)/gcd_factor;
        if sum(abs(new_a_red(n,:)-fix(new_a_red(n,:)))) > 1e-4 then
            disp('Error: PIL_conv_cell_gen, new_a_red are not integers!');
            abort
        end
    end

    new_a_cart=new_a_red*a_vec;
    orth_chk=zeros(1,3);
    orth_chk(1)=new_a_cart(1,:)*new_a_cart(2,:)';
    orth_chk(2)=new_a_cart(1,:)*new_a_cart(3,:)';
    orth_chk(3)=new_a_cart(2,:)*new_a_cart(3,:)';
endfunction


// example:
//a_vec=..
//[6.141900    0.000000    0.000000
//3.070950    5.319041    0.000000
//3.070950    1.773014    5.014840]
//
//new_b_red=[0 0 1/2; 1/2 0 0; 3/8 3/4 3/8];
//[new_a_red,orth_chk]=PIL_conv_cell_gen(a_vec,new_b_red)

