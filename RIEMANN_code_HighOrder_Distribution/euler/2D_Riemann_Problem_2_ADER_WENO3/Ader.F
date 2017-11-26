! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#include "directives"

#define SKIP_LINEAR_ALGEBRA 1
! Setting SKIP_LINEAR_ALGEBRA to 1 will cause us to circumvent the
! matrix inversion. (Might be useful for testing purposes.)

! Zone-centered densities go as:

! rho (x,y,z,t) = rho_1(t) + rho_2(t) x + rho_3(t) y + rho_4(t) z
!                      (second order terms till here)

!          + rho_5(t) (x^2-1/12) + rho_6(t) (y^2-1/12) + rho_7(t) (z^2-1/12)
!          + rho_8(t) x y + rho_9(t) y z + rho_10(t) z x
!                      (third order terms till here)

!          + rho_11(t) (x^3-3x/20) + rho_12(t) (y^3-3y/20)
!          + rho_13(t) (z^3-3z/20)
!          + rho_14(t) (x^2-1/12) y + rho_15(t) (x^2-1/12) z
!          + rho_16(t) (y^2-1/12) x + rho_17(t) (y^2-1/12) z
!          + rho_18(t) (z^2-1/12) x + rho_19(t) (z^2-1/12) y
!          + rho_20(t) x y z
!                      (fourth order terms till here)

!          + rho_21(t) (x^4-3x^2/14+3/560) + rho_22(t) (y^4-3y^2/14+3/560)
!          + rho_23(t) (z^4-3z^2/14+3/560)
!          + rho_24(t) (x^3-3x/20) y + rho_25(t) (x^3-3x/20) z
!          + rho_26(t) x (y^3-3y/20) + rho_27(t) (y^3-3y/20) z
!          + rho_28(t) x (z^3-3z/20) + rho_29(t) y (z^3-3z/20)
!          + rho_30(t) (x^2-1/12) (y^2-1/12) + rho_31(t) (y^2-1/12) (z^2-1/12)
!          + rho_32(t) (x^2-1/12) (z^2-1/12)
!          + rho_33(t) (x^2-1/12) y z + rho_34(t) x (y^2-1/12) z
!          + rho_35(t) x y (z^2-1/12)
!                      (fifth order terms till here)

!          + rho_36(t) (x^3-3x/20) y z + rho_37(t) x (y^3-3y/20) z
!          + rho_38(t) x y (z^3-3z/20)
!          + rho_39(t) (x^4-3x^2/14+3/560) y + rho_40(t) (x^4-3x^2/14+3/560) z
!          + rho_41(t) x (y^4-3y^2/14+3/560) + rho_42(t) (y^4-3y^2/14+3/560) z
!          + rho_43(t) x (z^4-3z^2/14+3/560) + rho_44(t) y (z^4-3z^2/14+3/560)
!                      (terms only needed for AMR-MHD)

! The above sequence is also in the same order as the basis functions
! themselves. i.e. the coefficients have the same ordinal number as the
! basis functions as programmed here.

! ------------------------------

! Face-centered Bx, By & Bz fields go as:
! ( See eqns. (2.2), (2.3) & (2.4) of Balsara 2006)

! bfx (y,z,t) = bfx_1(t) + bfx_2(t) y + bfx_3(t) z
!                      (second order terms till here)

!          + bfx_4(t) (y^2-1/12) + bfx_5(t) (z^2-1/12) + bfx_6(t) y z
!                      (third order terms till here)

!          + bfx_7(t) (y^3-3y/20) + bfx_8(t) (z^3-3z/20)
!          + bfx_9(t) (y^2-1/12) z + bfx_10(t) (z^2-1/12) y
!                      (fourth order terms till here)

!          + bfx_11(t) (y^3-3y/20) z + bfx_12(t) y (z^3-3z/20)
!                      (terms only needed for AMR)

! and in terms of basis functions it becomes

! bfx (y,z,t) = bfx_1(t) phi_1 + bfx_2(t) phi_3 + bfx_3(t) phi_4
!                      (second order terms till here)

!          + bfx_4(t) phi_6 + bfx_5(t) phi_7 + bfx_6(t) phi_9
!                      (third order terms till here)

!          + bfx_7(t) phi_12 + bfx_8(t) phi_13
!          + bfx_9(t) phi_17 + bfx_10(t) phi_19
!                      (fourth order terms till here)

!          + bfx_11(t) phi_27 + bfx_12(t) phi_29
!                      (terms only needed for AMR)

! ----------

! bfy (x,z,t) = bfy_1(t) + bfy_2(t) x + bfy_3(t) z
!                      (second order terms till here)

!          + bfy_4(t) (x^2-1/12) + bfy_5(t) (z^2-1/12) + bfy_6(t) x z
!                      (third order terms till here)

!          + bfy_7(t) (x^3-3x/20) + bfy_8(t) (z^3-3z/20)
!          + bfy_9(t) (x^2-1/12) z + bfy_10(t) (z^2-1/12) x
!                      (fourth order terms till here)

!          + bfy_11(t) (x^3-3x/20) z + bfy_12(t) x (z^3-3z/20)
!                      (terms only needed for AMR)

! and in terms of basis functions it becomes

! bfy (x,z,t) = bfy_1(t) phi_1 + bfy_2(t) phi_2 + bfy_3(t) phi_4
!                      (second order terms till here)

!          + bfy_4(t) phi_5 + bfy_5(t) phi_7 + bfy_6(t) phi_10
!                      (third order terms till here)

!          + bfy_7(t) phi_11 + bfy_8(t) phi_13 
!          + bfy_9(t) phi_15 + bfy_10(t) phi_18
!                      (fourth order terms till here)

!          + bfy_11(t) phi_25 + bfy_12(t) phi_28
!                      (terms only needed for AMR)

! ----------

! bfz (x,y,t) = bfz_1(t) + bfz_2(t) x + bfz_3(t) y
!                      (second order terms till here)

!          + bfz_4(t) (x^2-1/12) + bfz_5(t) (y^2-1/12) + bfz_6(t) x y
!                      (third order terms till here)

!          + bfz_7(t) (x^3-3x/20) + bfz_8(t) (y^3-3y/20)
!          + bfz_9(t) (x^2-1/12) y + bfz_10(t) (y^2-1/12) x
!                      (fourth order terms till here)

!          + bfz_11(t) (x^3-3x/20) y + bfz_12(t) x (y^3-3y/20)
!                      (terms only needed for AMR)

! and in terms of basis functions it becomes

! bfz (x,y,t) = bfz_1(t) phi_1 + bfz_2(t) phi_2 + bfz_3(t) phi_3
!                      (second order terms till here)

!          + bfz_4(t) phi_5 + bfz_5(t) phi_6 + bfz_6(t) phi_8
!                      (third order terms till here)

!          + bfz_7(t) phi_11 + bfz_8(t) phi_12
!          + bfz_9(t) phi_14 + bfz_10(t) phi_16
!                      (fourth order terms till here)

!          + bfz_11(t) phi_24 + bfz_12(t) phi_26
!                      (terms only needed for AMR)

! ------------------------------

! Zone-centered reconstruction of face-centered Bx, By & Bz fields go as:
! ( See eqns. (2.5), (2.6) & (2.7) of Balsara 2009)

! bfx (x,y,z,t) = bfx_1(t) + bfx_2(t) x + bfx_3(t) y + bfx_4(t) z
!          + bfx_5(t) (x^2-1/12) + bfx_6(t) x y + bfx_7(t) x z
!                      (second order terms till here)

!          + bfx_8(t) (y^2-1/12) + bfx_9(t) (z^2-1/12) + bfx_10(t) y z
!          + bfx_11(t) (x^3-3x/20) + bfx_12(t) (x^2-1/12) y
!          + bfx_13(t) (x^2-1/12) z + bfx_14(t) (y^2-1/12) x
!          + bfx_15(t) (z^2-1/12) x + bfx_16(t) x y z
!                      (third order terms till here)

!          + bfx_17(t) (y^3-3y/20) + bfx_18(t) (z^3-3z/20)
!          + bfx_19(t) (y^2-1/12) z + bfx_20(t) y (z^2-1/12)
!          + bfx_21(t) (x^4-3x^2/14+3/560) + bfx_22(t) (x^3-3x/20) y
!          + bfx_23(t) (x^3-3x/20) z + bfx_24(t) x (y^3-3y/20)
!          + bfx_25(t) x (z^3-3z/20) + bfx_26(t) (x^2-1/12) (y^2-1/12)
!          + bfx_27(t) (x^2-1/12) (z^2-1/12) + bfx_28(t) x (y^2-1/12) z
!          + bfx_29(t) x y (z^2-1/12)
!                      (fourth order terms till here)

!          + bfx_30(t) (y^3-3y/20) z + bfx_32(t) x (y^3-3y/20) z
!          + bfx_31(t) y (z^3-3z/20) + bfx_33(t) x y (z^3-3z/20)
!          + bfx_34(t) (x^4-3x^2/14+3/560) y + bfx_35(t) (x^4-3x^2/14+3/560) z
!                      (terms only needed for AMR)

! and in terms of basis functions it becomes

! bfx (x,y,z,t) = bfx_1 phi_1 + bfx_2 phi_2 + bfx_3 phi_3 + bfx_4 phi_4
!          + bfx_5 phi_5 + bfx_6 phi_8 + bfx_7 phi_10
!                      (second order terms till here)

!          + bfx_8 phi_6 + bfx_9 phi_7 + bfx_10 phi_9
!          + bfx_11 phi_11 + bfx_12 phi_14
!          + bfx_13 phi_15 + bfx_14 phi_16
!          + bfx_15 phi_18 + bfx_16 phi_20
!                      (third order terms till here)

!          + bfx_17 phi_12 + bfx_18 phi_13
!          + bfx_19 phi_17 + bfx_20 phi_19
!          + bfx_21 phi_21 + bfx_22 phi_24
!          + bfx_23 phi_25 + bfx_24 phi_26
!          + bfx_25 phi_28 + bfx_26 phi_30
!          + bfx_27 phi_32 + bfx_28 phi_34
!          + bfx_29 phi_35
!                      (fourth order terms till here)

!          + bfx_30(t) phi_27 + bfx_32(t) phi_37
!          + bfx_31(t) phi_29 + bfx_33(t) phi_38
!          + bfx_34(t) phi_39 + bfx_35(t) phi_40
!                      (terms only needed for AMR)

! ----------

! bfy (x,y,z,t) = bfy_1(t) + bfy_2(t) x + bfy_3(t) y + bfy_4(t) z
!          + bfy_5(t) (y^2-1/12) + bfy_6(t) x y + bfy_7(t) y z
!                      (second order terms till here)

!          + bfy_8(t) (x^2-1/12) + bfy_9(t) (z^2-1/12) + bfy_10(t) x z
!          + bfy_11(t) (y^3-3y/20) + bfy_12(t) (x^2-1/12) y
!          + bfy_13(t) (y^2-1/12) x + bfy_14(t) (y^2-1/12) z
!          + bfy_15(t) (z^2-1/12) y + bfy_16(t) x y z
!                      (third order terms till here)

!          + bfy_17(t) (x^3-3x/20) + bfy_18(t) (z^3-3z/20)
!          + bfy_19(t) (x^2-1/12) z + bfy_20(t) x (z^2-1/12)
!          + bfy_21(t) (y^4-3y^2/14+3/560) + bfy_22(t) (x^3-3x/20) y
!          + bfy_23(t) x (y^3-3y/20) + bfy_24(t) (y^3-3y/20) z
!          + bfy_25(t) y (z^3-3z/20) + bfy_26(t) (x^2-1/12) (y^2-1/12)
!          + bfy_27(t) (y^2-1/12) (z^2-1/12) + bfy_28 (x^2-1/12) y z
!          + bfy_29(t) x y (z^2-1/12)
!                      (fourth order terms till here)

!          + bfy_30(t) (x^3-3x/20) z + bfy_32(t) (x^3-3x/20) y z
!          + bfy_31(t) x (z^3-3z/20) + bfy_33(t) x y (z^3-3z/20)
!          + bfy_34(t) x (y^4-3y^2/14+3/560) + bfy_35(t) (y^4-3y^2/14+3/560) z
!                      (terms only needed for AMR)

! and in terms of basis functions it becomes

! bfy (x,y,z,t) = bfy_1 phi_1 + bfy_2 phi_2 + bfy_3 phi_3 + bfy_4 phi_4
!          + bfx_5 phi_6 + bfy_6 phi_8 + bfy_7 phi_9
!                      (second order terms till here)

!          + bfy_8 phi_5 + bfy_9 phi_7 + bfy_10 phi_10
!          + bfy_11 phi_12 + bfy_12 phi_14
!          + bfy_13 phi_16 + bfy_14 phi_17
!          + bfy_15 phi_19 + bfy_16 phi_20
!                      (third order terms till here)

!          + bfy_17 phi_11 + bfy_18 phi_13
!          + bfy_19 phi_15 + bfy_20 phi_18
!          + bfy_21 phi_22 + bfy_22 phi_24
!          + bfy_23 phi_26 + bfy_24 phi_27
!          + bfy_25 phi_29 + bfy_26 phi_30
!          + bfy_27 phi_31 + bfy_28 phi_33
!          + bfy_29 phi_35
!                      (fourth order terms till here)

!          + bfy_30(t) phi_25 + bfy_32(t) phi_36
!          + bfy_31(t) phi_28 + bfy_33(t) phi_38
!          + bfy_34(t) phi_41 + bfy_35(t) phi_42
!                      (terms only needed for AMR)

! ----------

! bfz (x,y,z,t) = bfz_1(t) + bfz_2(t) x + bfz_3(t) y + bfz_4(t) z
!          + bfz_5(t) (z^2-1/12) + bfz_6(t) y z + bfz_7(t) x z
!                      (second order terms till here)

!          + bfz_8(t) (x^2-1/12) + bfz_9(t) (y^2-1/12) + bfz_10(t) x y
!          + bfz_11(t) (z^3-3z/20) + bfz_12(t) (x^2-1/12) z
!          + bfz_13(t) (y^2-1/12) z + bfz_14(t) (z^2-1/12) x
!          + bfz_15(t) (z^2-1/12) y + bfz_16(t) x y z
!                      (third order terms till here)

!          + bfz_17(t) (x^3-3x/20) + bfz_18(t) (y^3-3y/20)
!          + bfz_19(t) (x^2-1/12) y + bfz_20(t) x (y^2-1/12)
!          + bfz_21(t) (z^4-3z^2/14+3/560) + bfz_22(t) (x^3-3x/20) z
!          + bfz_23(t) (y^3-3y/20) z + bfz_24(t) x (z^3-3z/20)
!          + bfz_25(t) y (z^3-3z/20) + bfz_26(t) (y^2-1/12) (z^2-1/12)
!          + bfz_27(t) (x^2-1/12) (z^2-1/12) + bfz_28(t) (x^2-1/12) y z
!          + bfz_29(t) x (y^2-1/12) z
!                      (fourth order terms till here)

!          + bfz_30(t) (x^3-3x/20) y + bfz_32(t) (x^3-3x/20) y z
!          + bfz_31(t) x (y^3-3y/20) + bfz_33(t) x (y^3-3y/20) z
!          + bfz_34(t) x (z^4-3z^2/14+3/560) + bfz_35(t) y (z^4-3z^2/14+3/560)
!                      (terms only needed for AMR)

! and in terms of basis functions it becomes

! bfz (x,y,z,t) = bfz_1 phi_1 + bfz_2 phi_2 + bfz_3 phi_3 + bfz_4 phi_4
!          + bfz_5 phi_7 + bfz_6 phi_9 + bfz_7 phi_10
!                      (second order terms till here)

!          + bfz_8 phi_5 + bfz_9 phi_6 + bfz_10 phi_8
!          + bfz_11 phi_13 + bfz_12 phi_15
!          + bfz_13 phi_17 + bfz_14 phi_18
!          + bfz_15 phi_19 + bfz_16 phi_20
!                      (third order terms till here)

!          + bfz_17 phi_11 + bfz_18 phi_12
!          + bfz_19 phi_14 + bfz_20 phi_16
!          + bfz_21 phi_23 + bfz_22 phi_25
!          + bfz_23 phi_27 + bfz_24 phi_28
!          + bfz_25 phi_29 + bfz_26 phi_31
!          + bfz_27 phi_32 + bfz_28 phi_33
!          + bfz_29 phi_34
!                      (fourth order terms till here)

!          + bfz_30(t) phi_24 + bfz_32(t) phi_36
!          + bfz_31(t) phi_26 + bfz_33(t) phi_37
!          + bfz_34(t) phi_43 + bfz_35(t) phi_44
!                      (terms only needed for AMR)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADER_CG_MODAL_3D_MHD_PTWISE

     1           ( dtmove, igeom, pi, rhofloor, prsfloor, smallnum,
     1             gamma, gasconst,
     1             ader_tolerance, max_iter_ader_cg,

     1             corrupted_zone,

     1             xc, yc, zc, dx, dy, dz,

     1             w, entropy_1d, gamma_eos_quadpt, molewt_quadpt,

     1             u, f, g, h, s, s_int)

! ----------------------------------------------------------------------

! Written by Dinshaw Balsara

! Incorporates ADER CG with source terms as contained
! in the paper by Balsara et al 2009, JCP. The present version
! accomodated Galerkin elements that are continuous in time and is,
! therefore, unsuitable for use with stiff source terms.
! It is substantially faster than ADER DG. It is useful for
! situations where the hyperbolic system is large and/or the Jacobian of
! the flux matrix is difficult to evaluate. Current version is for
! single-fluid MHD with multiple species of 2-fluid MHD. It scales linearly
! as the number of species and is, therefore, very useful for reactive flow.

! One MHD fluid and "N_EULER" number of Euler fluids (which may also be
! repurposed as radiation fields in various frequency bins) are treated here.

! Takes the conserved variables and their spatial moments "w" as input
! and returns the space-time moments of the conserved variables as "u".
! This is all done on a reference space-time element in 3d that spans
! (x,y,z,t) in the range [-.5,.5]X[-.5,.5]X[-.5,.5]X[0,1].

! "dtmove, dx, dy, dz" are timestep and x,y,z zone sizes. The update has
! to be done in physical variables so various combinations of these
! variables are multiplied into the flux and source terms.

! "gamma, pi, rhofloor, prsfloor, smallnum" are the usual inputs.

! gamma law gas is assumed, but allowance is made for more general eos if
! the spatial moments are provided in "gamma_eos_spacetime" with the time terms zeroed.

! "max_iter_ader_cg" the maximum number of permitted iterations.

! "f, g, h, s" are fluxes in x,y,z-directions and the source term respectively.
! They are not saved here. 

! "s_int" is the space-time integral of the source term which is also returned.

! ----------------------------------------------------------------------

! It helps to know the space-time modal basis set into which the solution is
! being projected. Here it is for orders 2, 3, and 4.

! P0[x_] := 1 ; P1[x_] := x ; P2[x_] := ( x^2 - 1/12) ;
! P3[x_] := ( x^3 - 3 * x / 20)

! Q0[t_] := 1 ; Q1[t_] := 2 * t - 1 ; Q2[t_] := 6 * t^2 - 6 * t + 1 ;
! Q3[t_] := 20 * t^3 - 30 * t^2 + 12 * t - 1 ; DQ0[t_] := D[Q0[t],t] ;
! DQ1[t_] := D[Q1[t],t] ; DQ2[t_] := D[Q2[t],t] ; DQ3[t_] := D[Q3[t],t]

! ------------------------------

! Second order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P0[x] * P0[y] * P0[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
! [1/2,0,0,0], [-1/2,0,0,0], [0,1/2,0,0], [0,-1/2,0,0],
! [0,0,1/2,0], [0,0,-1/2,0],
!   [0,0,0,1]

! ------------------------------

! Third order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P2[x] * P0[y] * P0[z] * Q0[t] +
!                   u6 * P0[x] * P2[y] * P0[z] * Q0[t] +
!                   u7 * P0[x] * P0[y] * P2[z] * Q0[t] +
!                   u8 * P1[x] * P1[y] * P0[z] * Q0[t] +
!                   u9 * P0[x] * P1[y] * P1[z] * Q0[t] +
!                   u10 * P1[x] * P0[y] * P1[z] * Q0[t] +
!                         u11 * P0[x] * P0[y] * P0[z] * Q1[t] +
!                         u12 * P0[x] * P0[y] * P0[z] * Q2[t] +
!                   u13 * P1[x] * P0[y] * P0[z] * Q1[t] +
!                   u14 * P0[x] * P1[y] * P0[z] * Q1[t] +
!                   u15 * P0[x] * P0[y] * P1[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
!   [0,0,0,0], [1/2,0,0,0], [-1/2,0,0,0], [0,1/2,0,0], [0,-1/2,0,0],
!   [0,0,1/2,0], [0,0,-1/2,0],
! [1/2,1/2,1/2,0], [-1/2,1/2,1/2,0], [1/2,-1/2,1/2,0], [-1/2,-1/2,1/2,0],
! [1/2,1/2,-1/2,0], [-1/2,1/2,-1/2,0], [1/2,-1/2,-1/2,0], [-1/2,-1/2,-1/2,0],
!   [1/2,0,0,1/2], [-1/2,0,0,1/2], [0,1/2,0,1/2], [0,-1/2,0,1/2],
!   [0,0,1/2,1/2], [0,0,-1/2,1/2],
! [0,0,0,1]

! ------------------------------

! Fourth order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P2[x] * P0[y] * P0[z] * Q0[t] +
!                   u6 * P0[x] * P2[y] * P0[z] * Q0[t] +
!                   u7 * P0[x] * P0[y] * P2[z] * Q0[t] +
!                   u8 * P1[x] * P1[y] * P0[z] * Q0[t] +
!                   u9 * P0[x] * P1[y] * P1[z] * Q0[t] +
!                   u10 * P1[x] * P0[y] * P1[z] * Q0[t] +
!                         u11 * P3[x] * P0[y] * P0[z] * Q0[t] +
!                         u12 * P0[x] * P3[y] * P0[z] * Q0[t] +
!                         u13 * P0[x] * P0[y] * P3[z] * Q0[t] +
!                         u14 * P2[x] * P1[y] * P0[z] * Q0[t] +
!                         u15 * P2[x] * P0[y] * P1[z] * Q0[t] +
!                         u16 * P1[x] * P2[y] * P0[z] * Q0[t] +
!                         u17 * P0[x] * P2[y] * P1[z] * Q0[t] +
!                         u18 * P1[x] * P0[y] * P2[z] * Q0[t] +
!                         u19 * P0[x] * P1[y] * P2[z] * Q0[t] +
!                         u20 * P1[x] * P1[y] * P1[z] * Q0[t] +
!                   u21 * P0[x] * P0[y] * P0[z] * Q1[t] +
!                   u22 * P0[x] * P0[y] * P0[z] * Q2[t] +
!                   u23 * P0[x] * P0[y] * P0[z] * Q3[t] +
!                         u24 * P1[x] * P0[y] * P0[z] * Q1[t] +
!                         u25 * P0[x] * P1[y] * P0[z] * Q1[t] +
!                         u26 * P0[x] * P0[y] * P1[z] * Q1[t] +
!                         u27 * P1[x] * P0[y] * P0[z] * Q2[t] +
!                         u28 * P0[x] * P1[y] * P0[z] * Q2[t] +
!                         u29 * P0[x] * P0[y] * P1[z] * Q2[t] +
!                   u30 * P2[x] * P0[y] * P0[z] * Q1[t] +
!                   u31 * P0[x] * P2[y] * P0[z] * Q1[t] +
!                   u32 * P0[x] * P0[y] * P2[z] * Q1[t] +
!                   u33 * P1[x] * P1[y] * P0[z] * Q1[t] +
!                   u34 * P0[x] * P1[y] * P1[z] * Q1[t] +
!                   u35 * P1[x] * P0[y] * P1[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
! [0,0,0,0], [1/2,0,0,0], [1/4,0,0,0], [-1/4,0,0,0], [-1/2,0,0,0],
! [0,1/2,0,0], [0,1/4,0,0], [0,-1/4,0,0], [0,-1/2,0,0],
! [0,0,1/2,0], [0,0,1/4,0], [0,0,-1/4,0], [0,0,-1/2,0],
!   [1/2,1/2,1/2,0], [0,1/2,1/2,0], [-1/2,1/2,1/2,0], [1/2,-1/2,1/2,0],
!   [0,-1/2,1/2,0], [-1/2,-1/2,1/2,0], [1/2,1/2,-1/2,0], [0,1/2,-1/2,0],
!   [-1/2,1/2,-1/2,0], [1/2,-1/2,-1/2,0], [0,-1/2,-1/2,0], [-1/2,-1/2,-1/2,0],
! [1/2,0,1/2,0], [-1/2,0,1/2,0], [1/2,0,-1/2,0], [-1/2,0,-1/2,0],
! [1/2,1/2,0,0], [-1/2,1/2,0,0], [1/2,-1/2,0,0], [-1/2,-1/2,0,0],
!   [0,0,0,1/3], [1/2,0,0,1/3], [-1/2,0,0,1/3], [0,1/2,0,1/3], [0,-1/2,0,1/3],
!   [0,0,1/2,1/3], [0,0,-1/2,1/3],
! [1/2,1/2,1/2,1/3], [-1/2,1/2,1/2,1/3], [1/2,-1/2,1/2,1/3], [-1/2,-1/2,1/2,1/3],
! [1/2,1/2,-1/2,1/3], [-1/2,1/2,-1/2,1/3], [1/2,-1/2,-1/2,1/3], [-1/2,-1/2,-1/2,1/3],
!   [1/2,0,0,2/3], [-1/2,0,0,2/3], [0,1/2,0,2/3], [0,-1/2,0,2/3],
!   [0,0,1/2,2/3], [0,0,-1/2,2/3],
! [0,0,0,1]

! ----------------------------------------------------------------------

      IMPLICIT NONE

! Variables declared here go through the subroutine interface.

      INTEGER max_iter_ader_cg, igeom,
     1        corrupted_zone

      REAL dtmove, gamma, gamma_euler, pi, gasconst, ader_tolerance, 
     1     rhofloor, prsfloor, smallnum,
     1     xc, yc, zc, dx, dy, dz

      REAL, DIMENSION ( NVAR, NUM_ZONECEN_BASES) :: w

      REAL, DIMENSION ( 5) :: entropy_1d

      REAL, DIMENSION ( NUM_SPACETIME_QUADPTS) ::
     1                gamma_eos_quadpt, molewt_quadpt

      REAL, DIMENSION ( NVAR, NUM_SPACETIME_BASES) :: u, f, g, h, s

      REAL, DIMENSION ( NVAR) :: s_int

! ----------------------------------------------------------------------

! Variables declared past this point don't go through the subroutine interface.

      INTEGER iter_ader, iquadpt, ibases, ivar

      REAL fourpi, ob4pi, ob8pi, ob3, tempaa, tempab, tempac,
     1     x1, x2, x3, obx1, cos_x2, sin_x2, cot_x2

! These are the flow variables (and their spatial derivatives, if needed)
! at the quadrature points.

      REAL, DIMENSION ( NUM_SPACETIME_QUADPTS) ::

     1      rho, prs, eng_t, vlx, vly, vlz, bfx, bfy, bfz,
     1      entropy

      REAL, DIMENSION ( 0: NFLUID, NUM_SPACETIME_QUADPTS) :: specie

! These are the conserved variables, fluxes and source terms at the
! quadrature points.

      REAL, DIMENSION ( NVAR, NUM_SPACETIME_QUADPTS) ::

     1      u_quadpt, f_quadpt, g_quadpt, h_quadpt, s_quadpt,
     1      tu, tf, tg, th, ts,

     1      du_quadpt_dx, du_quadpt_dy, du_quadpt_dz

! List of nodal points in each zone.

      REAL, DIMENSION ( 1:4, NUM_SPACETIME_QUADPTS) :: nodal_points

! ----------------------------------------------------------------------

! Make some sanity checks here.

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_MODAL_3D_MHD_PTWISE, NFI","ELD .NE. 7",
     1              NFIELD
      STOP
      END IF

! ----------

      ivar = NFIELD + 1 + NFLUID

      IF ( NVAR .NE. ivar) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1  "NVAR .NE. ivar","NVA","R = ", NVAR," ivar = ", ivar
      STOP
      END IF

! ----------

      IF ( max_iter_ader_cg .LT. NUMORDER - 1 ) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1  "max_iter_ader_cg .LT. NUMORDER ; max_iter_ader_cg = ",
     1   max_iter_ader_cg, NUMORDER - 1
      STOP
      END IF

! ----------

#if ( NUMORDER == 2)
      IF ( NUM_ZONECEN_BASES .NE. 4) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_ZONECEN","_BASES .NE. 4", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 5) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_SPACETIME","_BASES .NE. 5", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 6) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 6", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 7", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

! ----------

#if ( NUMORDER == 3)
      IF ( NUM_ZONECEN_BASES .NE. 10) THEN
      WRITE ( 6, *)
     1  "Stopping in third order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_ZONECEN","_BASES .NE. 10", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 15) THEN
      WRITE ( 6, *)
     1  "Stopping in third order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_SPACETIME","_BASES .NE. 15", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 13) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 15", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 20) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 22", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

! ----------

#if ( NUMORDER == 4)
      IF ( NUM_ZONECEN_BASES .NE. 20) THEN
      WRITE ( 6, *)
     1  "Stopping in fourth order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_ZONECEN","_BASES .NE. 20", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 35) THEN
      WRITE ( 6, *)
     1  "Stopping in fourth order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_SPACETIME","_BASES .NE. 35", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 21) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 33", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 42) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_MODAL_3D_MHD_PTWISE, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 55", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

! ----------------------------------------------------------------------

! Set up some constants that are used all over in this subroutine.

      fourpi = 4.0 * pi
      ob4pi = 1.0 / ( 4.0 * pi)
      ob8pi = 1.0 / ( 8.0 * pi)
      ob3 = 1.0 / 3.0

! ----------------------------------------------------------------------

! STEP 0) Set up the nodal points that are to be used at each of the orders.

! ------------------------------
#if ( NUMORDER == 2)
! cg, 2nd order nodal_points

      nodal_points ( 1:4, 01) = (/  0.5 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( 1:4, 02) = (/ -0.5 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( 1:4, 03) = (/  0.0 ,  0.5 ,  0.0 , 0.0  /) 
      nodal_points ( 1:4, 04) = (/  0.0 , -0.5 ,  0.0 , 0.0  /) 
      nodal_points ( 1:4, 05) = (/  0.0 ,  0.0 ,  0.5 , 0.0  /) 
      nodal_points ( 1:4, 06) = (/  0.0 ,  0.0 , -0.5 , 0.0  /) 

      nodal_points ( 1:4, 07) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /) 

#endif
! End "NUMORDER == 2".
! ------------------------------
#if ( NUMORDER == 3)
! cg, 3rd order nodal_points

      nodal_points ( 1:4, 01) = (/  0.0 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( 1:4, 02) = (/  0.5 ,  0.5 ,  0.0 , 0.0  /) 
      nodal_points ( 1:4, 03) = (/  0.5 , -0.5 ,  0.0 , 0.0  /) 
      nodal_points ( 1:4, 04) = (/  0.5 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( 1:4, 05) = (/  0.5 ,  0.0 , -0.5 , 0.0  /)
      nodal_points ( 1:4, 06) = (/ -0.5 ,  0.5 ,  0.0 , 0.0  /)
      nodal_points ( 1:4, 07) = (/ -0.5 , -0.5 ,  0.0 , 0.0  /)
      nodal_points ( 1:4, 08) = (/ -0.5 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( 1:4, 09) = (/ -0.5 ,  0.0 , -0.5 , 0.0  /)
      nodal_points ( 1:4, 10) = (/  0.0 ,  0.5 ,  0.5 , 0.0  /)
      nodal_points ( 1:4, 11) = (/  0.0 , -0.5 ,  0.5 , 0.0  /)
      nodal_points ( 1:4, 12) = (/  0.0 ,  0.5 , -0.5 , 0.0  /)
      nodal_points ( 1:4, 13) = (/  0.0 , -0.5 , -0.5 , 0.0  /)

      nodal_points ( 1:4, 14) = (/  0.5 ,  0.0 ,  0.0 , 0.5  /)
      nodal_points ( 1:4, 15) = (/ -0.5 ,  0.0 ,  0.0 , 0.5  /)
      nodal_points ( 1:4, 16) = (/  0.0 ,  0.5 ,  0.0 , 0.5  /)
      nodal_points ( 1:4, 17) = (/  0.0 , -0.5 ,  0.0 , 0.5  /)
      nodal_points ( 1:4, 18) = (/  0.0 ,  0.0 ,  0.5 , 0.5  /)
      nodal_points ( 1:4, 19) = (/  0.0 ,  0.0 , -0.5 , 0.5  /)

      nodal_points ( 1:4, 20) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /)

#endif
! End "NUMORDER == 3".
! ------------------------------
#if ( NUMORDER == 4)
! cg, 4th order nodal_points

      nodal_points ( 1:4, 01) = (/  0.0 ,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( 1:4, 02) = (/  0.5 ,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( 1:4, 03) = (/  0.25,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( 1:4, 04) = (/ -0.25,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( 1:4, 05) = (/ -0.5 ,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( 1:4, 06) = (/  0.0 ,  0.5 ,  0.0 , 0.0  /)
      nodal_points ( 1:4, 07) = (/  0.0 ,  0.25,  0.0 , 0.0  /)
      nodal_points ( 1:4, 08) = (/  0.0 , -0.25,  0.0 , 0.0  /)
      nodal_points ( 1:4, 09) = (/  0.0 , -0.5 ,  0.0 , 0.0  /)
      nodal_points ( 1:4, 10) = (/  0.0 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( 1:4, 11) = (/  0.0 ,  0.0 ,  0.25, 0.0  /)
      nodal_points ( 1:4, 12) = (/  0.0 ,  0.0 , -0.25, 0.0  /)
      nodal_points ( 1:4, 13) = (/  0.0 ,  0.0 , -0.5 , 0.0  /)
      nodal_points ( 1:4, 14) = (/  0.5 ,  0.5 ,  0.5 , 0.0  /)
      nodal_points ( 1:4, 15) = (/  0.5 , -0.5 ,  0.5 , 0.0  /)
      nodal_points ( 1:4, 16) = (/  0.5 ,  0.5 , -0.5 , 0.0  /)
      nodal_points ( 1:4, 17) = (/  0.5 , -0.5 , -0.5 , 0.0  /)
      nodal_points ( 1:4, 18) = (/ -0.5 ,  0.5 ,  0.5 , 0.0  /)
      nodal_points ( 1:4, 19) = (/ -0.5 , -0.5 ,  0.5 , 0.0  /)
      nodal_points ( 1:4, 20) = (/ -0.5 ,  0.5 , -0.5 , 0.0  /)
      nodal_points ( 1:4, 21) = (/ -0.5 , -0.5 , -0.5 , 0.0  /)

      nodal_points ( 1:4, 22) = (/  0.0 ,  0.0 ,  0.0 , 1.0/3.0 /)
      nodal_points ( 1:4, 23) = (/  0.5 ,  0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( 1:4, 24) = (/  0.5 , -0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( 1:4, 25) = (/  0.5 ,  0.0 ,  0.5 , 1.0/3.0 /)
      nodal_points ( 1:4, 26) = (/  0.5 ,  0.0 , -0.5 , 1.0/3.0 /)
      nodal_points ( 1:4, 27) = (/ -0.5 ,  0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( 1:4, 28) = (/ -0.5 , -0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( 1:4, 29) = (/ -0.5 ,  0.0 ,  0.5 , 1.0/3.0 /)
      nodal_points ( 1:4, 30) = (/ -0.5 ,  0.0 , -0.5 , 1.0/3.0 /)
      nodal_points ( 1:4, 31) = (/  0.0 ,  0.5 ,  0.5 , 1.0/3.0 /)
      nodal_points ( 1:4, 32) = (/  0.0 , -0.5 ,  0.5 , 1.0/3.0 /)
      nodal_points ( 1:4, 33) = (/  0.0 ,  0.5 , -0.5 , 1.0/3.0 /)
      nodal_points ( 1:4, 34) = (/  0.0 , -0.5 , -0.5 , 1.0/3.0 /)

      nodal_points ( 1:4, 35) = (/  0.0 ,  0.0 ,  0.0 , 2.0/3.0 /)
      nodal_points ( 1:4, 36) = (/  0.5 ,  0.0 ,  0.0 , 2.0/3.0 /)
      nodal_points ( 1:4, 37) = (/ -0.5 ,  0.0 ,  0.0 , 2.0/3.0 /)
      nodal_points ( 1:4, 38) = (/  0.0 ,  0.5 ,  0.0 , 2.0/3.0 /)
      nodal_points ( 1:4, 39) = (/  0.0 , -0.5 ,  0.0 , 2.0/3.0 /)
      nodal_points ( 1:4, 40) = (/  0.0 ,  0.0 ,  0.5 , 2.0/3.0 /)
      nodal_points ( 1:4, 41) = (/  0.0 ,  0.0 , -0.5 , 2.0/3.0 /)

      nodal_points ( 1:4, 42) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /)

#endif
! End "NUMORDER == 4".
! ------------------------------

#if ( ENFORCE_POSITIVITY == 1)

! While "nodal_points" still represents positions in the reference space-time
! hypercube, we use it to build the "entropy" at all the nodes.

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      entropy ( iquadpt) = entropy_1d ( 01)
     1    + nodal_points ( 1, iquadpt) * entropy_1d ( 02)
     1    + nodal_points ( 2, iquadpt) * entropy_1d ( 03)
     1    + nodal_points ( 3, iquadpt) * entropy_1d ( 04)
     1    + nodal_points ( 4, iquadpt) * entropy_1d ( 05)

      END DO ! End of "DO iquadpt =" loop.

#endif
! End of "ENFORCE_POSITIVITY == 1".

! ------------------------------

! Set up the nodal points in a global space-time. These are coordinate
! values.

      nodal_points ( 1, :) = nodal_points ( 1, :) * dx + xc 
      nodal_points ( 2, :) = nodal_points ( 2, :) * dy + yc
      nodal_points ( 3, :) = nodal_points ( 3, :) * dz + zc
      nodal_points ( 4, :) = nodal_points ( 4, :) * dtmove 

! ----------------------------------------------------------------------

! STEP 1) Having come into the subroutine with variables loaded into "w",
! load the appropriate spatial variation of "w" into "u" while setting
! the temporal variation of "u" as well as "f,g,h & s" to zero initially.

      u ( :, 1: NUM_ZONECEN_BASES) = w ( :, 1: NUM_ZONECEN_BASES)
      u ( :, NUM_ZONECEN_BASES + 1: NUM_SPACETIME_BASES) = 0.0

      f ( :, :) = 0.0
      g ( :, :) = 0.0
      h ( :, :) = 0.0
      s ( :, :) = 0.0

! ----------------------------------------------------------------------

! STEP 2) Use the spatial modes to make "u_quadpt ( :, 1: NUM_ZONECEN_QUADPTS)"
! from "u ( :, 1: NUM_ZONECEN_BASES)" at time t = 0.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : modal -> nodal

      u_quadpt ( :, 01) = u ( :, 01) + 0.5 * u ( :, 02)

      u_quadpt ( :, 02) = u ( :, 01) - 0.5 * u ( :, 02)

      u_quadpt ( :, 03) = u ( :, 01) + 0.5 * u ( :, 03)

      u_quadpt ( :, 04) = u ( :, 01) - 0.5 * u ( :, 03)

      u_quadpt ( :, 05) = u ( :, 01) + 0.5 * u ( :, 04)

      u_quadpt ( :, 06) = u ( :, 01) - 0.5 * u ( :, 04)

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : modal -> nodal

      u_quadpt ( :, 01) = u ( :, 01) - ( u ( :, 05)
     1 + u ( :, 06) + u ( :, 07) ) / 12.0

      u_quadpt ( :, 02) = u_quadpt ( :, 01)
     1 + 0.5 * ( u ( :, 02) + u ( :, 03) )
     1 + 0.25 * ( u ( :, 05) + u ( :, 06) + u ( :, 08) )

      u_quadpt ( :, 03) = u_quadpt ( :, 02) - u ( :, 03)
     1 - 0.5 * u ( :, 08)

      u_quadpt ( :, 04) = u_quadpt ( :, 01)
     1 + 0.5 * ( u ( :, 02) + u ( :, 04) )
     1 + 0.25 * ( u ( :, 05) + u ( :, 07) + u ( :, 10) )

      u_quadpt ( :, 05) = u_quadpt ( :, 04) - u ( :, 04)
     1 - 0.5 * u ( :, 10)

      u_quadpt ( :, 06) = u_quadpt ( :, 03) - u ( :, 02)
     1 + u ( :, 03)

      u_quadpt ( :, 07) = u_quadpt ( :, 02) - u ( :, 02)
     1 - u ( :, 03)

      u_quadpt ( :, 08) = u_quadpt ( :, 05) - u ( :, 02)
     1 + u ( :, 04)

      u_quadpt ( :, 09) = u_quadpt ( :, 04) - u ( :, 02)
     1 - u ( :, 04)

      u_quadpt ( :, 10) = u_quadpt ( :, 01)
     1 + 0.5 * ( u ( :, 03) + u ( :, 04) )
     1 + 0.25 * ( u ( :, 06) + u ( :, 07) + u ( :, 09) )

      u_quadpt ( :, 11) = u_quadpt ( :, 10) - u ( :, 03)
     1 - 0.5 * u ( :, 09)

      u_quadpt ( :, 12) = u_quadpt ( :, 11) + u ( :, 03)
     1 - u ( :, 04)

      u_quadpt ( :, 13) = u_quadpt ( :, 10) - u ( :, 03)
     1 - u ( :, 04)

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : modal -> nodal

      u_quadpt ( :, 01) = u ( :, 01) - ( u ( :, 05)
     1 + u ( :, 06) + u ( :, 07) ) / 12.0

      u_quadpt ( :, 02) = u_quadpt ( :, 01)
     1 + ( 60.0 * u ( :, 02) + 30.0 * u ( :, 05)
     1 + 6.0 * u ( :, 11) - 5.0 * u ( :, 16)
     1 - 5.0 * u ( :, 18) ) / 120.0

      u_quadpt ( :, 03) = 0.5 * ( u_quadpt ( :, 01)
     1 + u_quadpt ( :, 02) ) - 0.015625 * ( 4 * u ( :, 05)
     1 + 3.0 * u ( :, 11) )

      u_quadpt ( :, 04) = 3.0 * u_quadpt ( :, 01)
     1 + u_quadpt ( :, 02) - 3.0 * u_quadpt ( :, 03)
     1 - 0.09375 * u ( :, 11)

      u_quadpt ( :, 05) = 6.0 * u_quadpt ( :, 01)
     1 + 3.0 * u_quadpt ( :, 02) - 8.0 * u_quadpt ( :, 03)
     1 - 0.375 * u ( :, 11)

      u_quadpt ( :, 06) = u_quadpt ( :, 01)
     1 + ( 60.0 * u ( :, 03) + 30.0 * u ( :, 06)
     1 + 6.0 * u ( :, 12) - 5.0 * u ( :, 14)
     1 - 5.0 * u ( :, 19) ) / 120.0

      u_quadpt ( :, 07) = 0.5 * ( u_quadpt ( :, 01)
     1 + u_quadpt ( :, 06) ) - 0.015625 * ( 4.0 * u ( :, 06)
     1 + 3.0 * u ( :, 12) )

      u_quadpt ( :, 08) = 3.0 * u_quadpt ( :, 01)
     1 + u_quadpt ( :, 06) - 3.0 * u_quadpt ( :, 07)
     1 - 0.09375 * u ( :, 12)

      u_quadpt ( :, 09) = 6.0 * u_quadpt ( :, 01)
     1 + 3.0 * u_quadpt ( :, 06) - 8.0 * u_quadpt ( :, 07)
     1 - 0.375 * u ( :, 12)

      u_quadpt ( :, 10) = u_quadpt ( :, 01)
     1 + ( 60.0 * u ( :, 04) + 30.0 * u ( :, 07)
     1 + 6.0 * u ( :, 13) - 5.0 * u ( :, 15)
     1 - 5.0 * u ( :, 17) ) / 120.0

      u_quadpt ( :, 11) = 0.5 * ( u_quadpt ( :, 01)
     1 + u_quadpt ( :, 10) ) - 0.015625 * ( 4.0 * u ( :, 07)
     1 + 3.0 * u ( :, 13) )

      u_quadpt ( :, 12) = 3.0 * u_quadpt ( :, 01)
     1 + u_quadpt ( :, 10) - 3.0 * u_quadpt ( :, 11)
     1 - 0.09375 * u ( :, 13)

      u_quadpt ( :, 13) = 6.0 * u_quadpt ( :, 01)
     1 + 3.0 * u_quadpt ( :, 10) - 8.0 * u_quadpt ( :, 11)
     1 - 0.375 * u ( :, 13)

      u_quadpt ( :, 14) = u_quadpt ( :, 02)
     1 + u_quadpt ( :, 07) + u_quadpt ( :, 10)
     1 - 2.0 * u_quadpt ( :, 01) + 0.25 * ( u ( :, 03)
     1 + u ( :, 08) + u ( :, 09) + u ( :, 10) )
     1 + 0.125 * ( u ( :, 15) + u ( :, 16) + u ( :, 17)
     1 + u ( :, 18) + u ( :, 20) ) + 0.1875 * u ( :, 06)
     1 + 5.0 * ( u ( :, 14) + u ( :, 19) ) / 48.0
     1 + 0.071875 * u ( :, 12)

      u_quadpt ( :, 15) = u_quadpt ( :, 14)
     1 + u_quadpt ( :, 09) - u_quadpt ( :, 06)
     1 - 0.5 * ( u ( :, 08) + u ( :, 09) )
     1 - 0.25 * ( u ( :, 14) + u ( :, 19) + u ( :, 20) )

      u_quadpt ( :, 16) = u_quadpt ( :, 14)
     1 + u_quadpt ( :, 13) - u_quadpt ( :, 10)
     1 - 0.5 * ( u ( :, 09) + u ( :, 10) )
     1 - 0.25 * ( u ( :, 15) + u ( :, 17) + u ( :, 20) )

      u_quadpt ( :, 17) = u_quadpt ( :, 15)
     1 + u_quadpt ( :, 16) - u_quadpt ( :, 14) + u ( :, 09)
     1 + 0.5 * u ( :, 20)

      u_quadpt ( :, 18) = u_quadpt ( :, 14)
     1 + u_quadpt ( :, 05) - u_quadpt ( :, 02)
     1 - 0.5 * ( u ( :, 08) + u ( :, 10) )
     1 - 0.25 * ( u ( :, 16) + u ( :, 18) + u ( :, 20) )

      u_quadpt ( :, 19) = u_quadpt ( :, 18)
     1 + u_quadpt ( :, 09) - u_quadpt ( :, 06)
     1 + 0.5 * ( u ( :, 08) - u ( :, 09) )
     1 - 0.25 * ( u ( :, 14) + u ( :, 19) - u ( :, 20) )

      u_quadpt ( :, 20) = u_quadpt ( :, 18)
     1 + u_quadpt ( :, 13) - u_quadpt ( :, 10)
     1 - 0.5 * ( u ( :, 09) - u ( :, 10) )
     1 - 0.25 * ( u ( :, 15) + u ( :, 17) - u ( :, 20) )

      u_quadpt ( :, 21) = u_quadpt ( :, 19)
     1 + u_quadpt ( :, 20) - u_quadpt ( :, 18) + u ( :, 09)
     1 - 0.5 * u ( :, 20)

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------

! STEP 3) Get the primitive variables from the conserved variables at
! the spatial nodes at t=0. The density and pressure variables are
! given positive floor values to ensure that they can be used for physically
! meaningful source terms.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

! ----------

      u_quadpt ( 1, iquadpt) =
     1                   AMAX1 ( u_quadpt ( 1, iquadpt), rhofloor)

      rho ( iquadpt) = u_quadpt ( 1, iquadpt)

      tempaa = 1.0 / rho ( iquadpt)

      vlx ( iquadpt) = u_quadpt ( 2, iquadpt) * tempaa
      vly ( iquadpt) = u_quadpt ( 3, iquadpt) * tempaa
      vlz ( iquadpt) = u_quadpt ( 4, iquadpt) * tempaa

      bfx ( iquadpt) = u_quadpt ( 6, iquadpt)
      bfy ( iquadpt) = u_quadpt ( 7, iquadpt)
      bfz ( iquadpt) = u_quadpt ( 8, iquadpt)

      tempab = 0.5 * rho ( iquadpt) * ( vlx ( iquadpt)**2
     1         + vly ( iquadpt)**2 + vlz ( iquadpt)**2)

      tempac = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      eng_t ( iquadpt) = u_quadpt ( 5, iquadpt)
     1                     - tempab - tempac

      prs ( iquadpt) = eng_t ( iquadpt)
     1                   * ( gamma_eos_quadpt ( iquadpt) - 1.0)

#if ( ENFORCE_POSITIVITY == 1)
      IF ( prs ( iquadpt) .LE. prsfloor) THEN

      prs ( iquadpt) = entropy ( iquadpt)
     1* rho ( iquadpt)**( gamma_eos_quadpt ( iquadpt) - 1.0)
      prs ( iquadpt) = AMAX1 ( prs ( iquadpt), prsfloor)

      eng_t ( iquadpt) = prs ( iquadpt)
     1                 / ( gamma_eos_quadpt ( iquadpt) - 1.0)

      u_quadpt ( 5, iquadpt) = eng_t ( iquadpt)
     1                           + tempab + tempac

      END IF
#else
      IF ( prs ( iquadpt) .LE. prsfloor) THEN

      prs ( iquadpt) = prsfloor

      eng_t ( iquadpt) = prs ( iquadpt)
     1                 / ( gamma_eos_quadpt ( iquadpt) - 1.0)

      u_quadpt ( 5, iquadpt) = eng_t ( iquadpt)
     1                           + tempab + tempac

      END IF
#endif

#if ( NFLUID >= 1)
      specie ( 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( 1: NFLUID, iquadpt) * rho ( iquadpt)
#endif

! ----------

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 4.a) Build fluxes ( f_quadpt, g_quadpt, h_quadpt) from the primitive
! variables that have now become available at spatial nodes at t=0.

! Multiply "f_quadpt" by "dtmove/dx" and so on.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

! -----

      tempaa = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      tempab = ( vlx ( iquadpt) * bfx ( iquadpt)
     1         + vly ( iquadpt) * bfy ( iquadpt)
     1         + vlz ( iquadpt) * bfz ( iquadpt)) * ob4pi

! -----

      f_quadpt ( 1, iquadpt) = u_quadpt ( 2, iquadpt)

      f_quadpt ( 2, iquadpt) = u_quadpt ( 2, iquadpt)
     1         * vlx ( iquadpt)
     1         + prs ( iquadpt)
     1         + tempaa - bfx ( iquadpt)**2 * ob4pi

      f_quadpt ( 3, iquadpt) = u_quadpt ( 2, iquadpt)
     1         * vly ( iquadpt)
     1         - bfx ( iquadpt) * bfy ( iquadpt) * ob4pi

      f_quadpt ( 4, iquadpt) = u_quadpt ( 2, iquadpt)
     1         * vlz ( iquadpt)
     1         - bfx ( iquadpt) * bfz ( iquadpt) * ob4pi

      f_quadpt ( 5, iquadpt) = ( u_quadpt ( 5, iquadpt)
     1         + prs ( iquadpt) + tempaa) * vlx ( iquadpt)
     1         - bfx ( iquadpt) * tempab

      f_quadpt ( 6, iquadpt) = 0.0

      f_quadpt ( 7, iquadpt) =
     1           vlx ( iquadpt) * bfy ( iquadpt)
     1         - vly ( iquadpt) * bfx ( iquadpt)

      f_quadpt ( 8, iquadpt) =
     1           vlx ( iquadpt) * bfz ( iquadpt)
     1         - vlz ( iquadpt) * bfx ( iquadpt)

#if ( NFLUID >= 1)
      f_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlx ( iquadpt)
#endif

! -----

      g_quadpt ( 1, iquadpt) = u_quadpt ( 3, iquadpt)

      g_quadpt ( 2, iquadpt) = u_quadpt ( 3, iquadpt)
     1         * vlx ( iquadpt)
     1         - bfy ( iquadpt) * bfx ( iquadpt) * ob4pi

      g_quadpt ( 3, iquadpt) = u_quadpt ( 3, iquadpt)
     1         * vly ( iquadpt)
     1         + prs ( iquadpt)
     1         + tempaa - bfy ( iquadpt)**2 * ob4pi

      g_quadpt ( 4, iquadpt) = u_quadpt ( 3, iquadpt)
     1         * vlz ( iquadpt)
     1         - bfy ( iquadpt) * bfz ( iquadpt) * ob4pi

      g_quadpt ( 5, iquadpt) = ( u_quadpt ( 5, iquadpt)
     1         + prs ( iquadpt) + tempaa) * vly ( iquadpt)
     1                       - bfy ( iquadpt) * tempab

      g_quadpt ( 6, iquadpt) =
     1           vly ( iquadpt) * bfx ( iquadpt)
     1         - vlx ( iquadpt) * bfy ( iquadpt)

      g_quadpt ( 7, iquadpt) = 0.0

      g_quadpt ( 8, iquadpt) =
     1           vly ( iquadpt) * bfz ( iquadpt)
     1         - vlz ( iquadpt) * bfy ( iquadpt)

#if ( NFLUID >= 1)
      g_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vly ( iquadpt)
#endif

! -----

      h_quadpt ( 1, iquadpt) = u_quadpt ( 4, iquadpt)

      h_quadpt ( 2, iquadpt) = u_quadpt ( 4, iquadpt)
     1         * vlx ( iquadpt)
     1         - bfz ( iquadpt) * bfx ( iquadpt) * ob4pi

      h_quadpt ( 3, iquadpt) = u_quadpt ( 4, iquadpt)
     1         * vly ( iquadpt)
     1         - bfz ( iquadpt) * bfy ( iquadpt) * ob4pi

      h_quadpt ( 4, iquadpt) = u_quadpt ( 4, iquadpt)
     1         * vlz ( iquadpt)
     1         + prs ( iquadpt)
     1         + tempaa - bfz ( iquadpt)**2 * ob4pi

      h_quadpt ( 5, iquadpt) = ( u_quadpt ( 5, iquadpt)
     1         + prs ( iquadpt) + tempaa) * vlz ( iquadpt)
     1                       - bfz ( iquadpt) * tempab

      h_quadpt ( 6, iquadpt) =
     1           vlz ( iquadpt) * bfx ( iquadpt)
     1         - vlx ( iquadpt) * bfz ( iquadpt)

      h_quadpt ( 7, iquadpt) =
     1           vlz ( iquadpt) * bfy ( iquadpt)
     1         - vly ( iquadpt) * bfz ( iquadpt)

      h_quadpt ( 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      h_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlz ( iquadpt)
#endif

! -----

      END DO
! Finished "DO iquadpt ="

! ------------------------------

! To make cylindrical and spherical geometries accessible to non-rel.
! flow, we put this part after either flux construction.

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      x1 = nodal_points ( 1, iquadpt)
      x2 = nodal_points ( 2, iquadpt)
      x3 = nodal_points ( 3, iquadpt)

      obx1 = 1.0 / x1

      IF ( igeom == 2 .or. igeom == 3) THEN

         g_quadpt ( :, iquadpt) = g_quadpt ( :, iquadpt) * obx1

      END IF

      IF ( igeom == 3) THEN

         sin_x2 = 1.0 / ( x1 * sin ( x2))

         h_quadpt ( :, iquadpt) = h_quadpt ( :, iquadpt) * sin_x2

      END IF

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

      tempaa = dtmove / dx

      f_quadpt ( :, 1: NUM_ZONECEN_QUADPTS) = tempaa
     1         * f_quadpt ( :, 1: NUM_ZONECEN_QUADPTS)

      tempab = dtmove / dy

      g_quadpt ( :, 1: NUM_ZONECEN_QUADPTS) = tempab
     1         * g_quadpt ( :, 1: NUM_ZONECEN_QUADPTS)

      tempac = dtmove / dz

      h_quadpt ( :, 1: NUM_ZONECEN_QUADPTS) = tempac
     1         * h_quadpt ( :, 1: NUM_ZONECEN_QUADPTS)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 4.b) Now project the fluxes "( f_quadpt, g_quadpt, h_quadpt)"
! from the spatial nodes at t=0 back to the t=0 part of the modal
! space time variables for those variables, i.e. to "( f, g, h)" .

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal

! ----------

!      f ( :, 01) = ( f_quadpt ( :, 01) + f_quadpt ( :, 02)
!     1 + f_quadpt ( :, 03) + f_quadpt ( :, 04)
!     1 + f_quadpt ( :, 05) + f_quadpt ( :, 06)) / 6.0

      f ( :, 02) = f_quadpt ( :, 01) - f_quadpt ( :, 02)

!      f ( :, 03) = f_quadpt ( :, 03) - f_quadpt ( :, 04)

!      f ( :, 04) = f_quadpt ( :, 05) - f_quadpt ( :, 06)

! ----------

!      g ( :, 01) = ( g_quadpt ( :, 01) + g_quadpt ( :, 02)
!     1 + g_quadpt ( :, 03) + g_quadpt ( :, 04)
!     1 + g_quadpt ( :, 05) + g_quadpt ( :, 06)) / 6.0

!      g ( :, 02) = g_quadpt ( :, 01) - g_quadpt ( :, 02)

      g ( :, 03) = g_quadpt ( :, 03) - g_quadpt ( :, 04)

!      g ( :, 04) = g_quadpt ( :, 05) - g_quadpt ( :, 06)

! ----------

!      h ( :, 01) = ( h_quadpt ( :, 01) + h_quadpt ( :, 02)
!     1 + h_quadpt ( :, 03) + h_quadpt ( :, 04)
!     1 + h_quadpt ( :, 05) + h_quadpt ( :, 06)) / 6.0

!      h ( :, 02) = h_quadpt ( :, 01) - h_quadpt ( :, 02)

!      h ( :, 03) = h_quadpt ( :, 03) - h_quadpt ( :, 04)

      h ( :, 04) = h_quadpt ( :, 05) - h_quadpt ( :, 06)

! ----------

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal

! ----------

      tf ( :, 01) = f_quadpt ( :, 02) + f_quadpt ( :, 03)
     1           + f_quadpt ( :, 04) + f_quadpt ( :, 05)

      tf ( :, 02) = f_quadpt ( :, 06) + f_quadpt ( :, 07)
     1           + f_quadpt ( :, 08) + f_quadpt ( :, 09)

      tf ( :, 03) = f_quadpt ( :, 10) + f_quadpt ( :, 11)
     1           + f_quadpt ( :, 12) + f_quadpt ( :, 13)

!      tf ( :, 04) = f_quadpt ( :, 02) + f_quadpt ( :, 06)
!     1           + f_quadpt ( :, 10) + f_quadpt ( :, 12)

!      tf ( :, 05) = f_quadpt ( :, 03) + f_quadpt ( :, 07)
!     1           + f_quadpt ( :, 11) + f_quadpt ( :, 13)

!      tf ( :, 06) = f_quadpt ( :, 04) + f_quadpt ( :, 05)
!     1           + f_quadpt ( :, 08) + f_quadpt ( :, 09)

!      tf ( :, 07) = f_quadpt ( :, 04) + f_quadpt ( :, 08)
!     1           + f_quadpt ( :, 10) + f_quadpt ( :, 11)

!      tf ( :, 08) = f_quadpt ( :, 05) + f_quadpt ( :, 09)
!     1           + f_quadpt ( :, 12) + f_quadpt ( :, 13)

!      tf ( :, 09) = f_quadpt ( :, 02) + f_quadpt ( :, 03)
!     1           + f_quadpt ( :, 06) + f_quadpt ( :, 07)

      f ( :, 02) = 0.25 * ( tf ( :, 01) - tf ( :, 02) )

!      f ( :, 03) = 0.25 * ( tf ( :, 04) - tf ( :, 05) )

!      f ( :, 04) = 0.25 * ( tf ( :, 07) - tf ( :, 08) )

      f ( :, 05) = 0.5 * ( tf ( :, 01) + tf ( :, 02)
     1 - tf ( :, 03) ) - 2.0 * f_quadpt ( :, 01)

!      f ( :, 06) = 0.5 * ( tf ( :, 04) + tf ( :, 05)
!     1 - tf ( :, 06) ) - 2.0 * f_quadpt ( :, 01)

!      f ( :, 07) = 0.5 * ( tf ( :, 07) + tf ( :, 08)
!     1 - tf ( :, 09) ) - 2.0 * f_quadpt ( :, 01)

      f ( :, 08) = f_quadpt ( :, 02) - f_quadpt ( :, 03)
     1 - f_quadpt ( :, 06) + f_quadpt ( :, 07)

!      f ( :, 09) = f_quadpt ( :, 10) - f_quadpt ( :, 11)
!     1 - f_quadpt ( :, 12) + f_quadpt ( :, 13)

      f ( :, 10) = f_quadpt ( :, 04) - f_quadpt ( :, 05)
     1 - f_quadpt ( :, 08) + f_quadpt ( :, 09)

!      f ( :, 01) = f_quadpt ( :, 01) + ( f ( :, 05)
!     1 + f ( :, 06) + f ( :, 07) ) / 12.0

! ----------

!      tg ( :, 01) = g_quadpt ( :, 02) + g_quadpt ( :, 03)
!     1           + g_quadpt ( :, 04) + g_quadpt ( :, 05)

!      tg ( :, 02) = g_quadpt ( :, 06) + g_quadpt ( :, 07)
!     1           + g_quadpt ( :, 08) + g_quadpt ( :, 09)

!      tg ( :, 03) = g_quadpt ( :, 10) + g_quadpt ( :, 11)
!     1           + g_quadpt ( :, 12) + g_quadpt ( :, 13)

      tg ( :, 04) = g_quadpt ( :, 02) + g_quadpt ( :, 06)
     1           + g_quadpt ( :, 10) + g_quadpt ( :, 12)

      tg ( :, 05) = g_quadpt ( :, 03) + g_quadpt ( :, 07)
     1           + g_quadpt ( :, 11) + g_quadpt ( :, 13)

      tg ( :, 06) = g_quadpt ( :, 04) + g_quadpt ( :, 05)
     1           + g_quadpt ( :, 08) + g_quadpt ( :, 09)

!      tg ( :, 07) = g_quadpt ( :, 04) + g_quadpt ( :, 08)
!     1           + g_quadpt ( :, 10) + g_quadpt ( :, 11)

!      tg ( :, 08) = g_quadpt ( :, 05) + g_quadpt ( :, 09)
!     1           + g_quadpt ( :, 12) + g_quadpt ( :, 13)

!      tg ( :, 09) = g_quadpt ( :, 02) + g_quadpt ( :, 03)
!     1           + g_quadpt ( :, 06) + g_quadpt ( :, 07)

!      g ( :, 02) = 0.25 * ( tg ( :, 01) - tg ( :, 02) )

      g ( :, 03) = 0.25 * ( tg ( :, 04) - tg ( :, 05) )

!      g ( :, 04) = 0.25 * ( tg ( :, 07) - tg ( :, 08) )

!      g ( :, 05) = 0.5 * ( tg ( :, 01) + tg ( :, 02)
!     1 - tg ( :, 03) ) - 2.0 * g_quadpt ( :, 01)

      g ( :, 06) = 0.5 * ( tg ( :, 04) + tg ( :, 05)
     1 - tg ( :, 06) ) - 2.0 * g_quadpt ( :, 01)

!      g ( :, 07) = 0.5 * ( tg ( :, 07) + tg ( :, 08)
!     1 - tg ( :, 09) ) - 2.0 * g_quadpt ( :, 01)

      g ( :, 08) = g_quadpt ( :, 02) - g_quadpt ( :, 03)
     1 - g_quadpt ( :, 06) + g_quadpt ( :, 07)

      g ( :, 09) = g_quadpt ( :, 10) - g_quadpt ( :, 11)
     1 - g_quadpt ( :, 12) + g_quadpt ( :, 13)

!      g ( :, 10) = g_quadpt ( :, 04) - g_quadpt ( :, 05)
!     1 - g_quadpt ( :, 08) + g_quadpt ( :, 09)

!      g ( :, 01) = g_quadpt ( :, 01) + ( g ( :, 05)
!     1 + g ( :, 06) + g ( :, 07) ) / 12.0

! ----------

!      th ( :, 01) = h_quadpt ( :, 02) + h_quadpt ( :, 03)
!     1           + h_quadpt ( :, 04) + h_quadpt ( :, 05)

!      th ( :, 02) = h_quadpt ( :, 06) + h_quadpt ( :, 07)
!     1           + h_quadpt ( :, 08) + h_quadpt ( :, 09)

!      th ( :, 03) = h_quadpt ( :, 10) + h_quadpt ( :, 11)
!     1           + h_quadpt ( :, 12) + h_quadpt ( :, 13)

!      th ( :, 04) = h_quadpt ( :, 02) + h_quadpt ( :, 06)
!     1           + h_quadpt ( :, 10) + h_quadpt ( :, 12)

!      th ( :, 05) = h_quadpt ( :, 03) + h_quadpt ( :, 07)
!     1           + h_quadpt ( :, 11) + h_quadpt ( :, 13)

!      th ( :, 06) = h_quadpt ( :, 04) + h_quadpt ( :, 05)
!     1           + h_quadpt ( :, 08) + h_quadpt ( :, 09)

      th ( :, 07) = h_quadpt ( :, 04) + h_quadpt ( :, 08)
     1           + h_quadpt ( :, 10) + h_quadpt ( :, 11)

      th ( :, 08) = h_quadpt ( :, 05) + h_quadpt ( :, 09)
     1           + h_quadpt ( :, 12) + h_quadpt ( :, 13)

      th ( :, 09) = h_quadpt ( :, 02) + h_quadpt ( :, 03)
     1           + h_quadpt ( :, 06) + h_quadpt ( :, 07)

!      h ( :, 02) = 0.25 * ( th ( :, 01) - th ( :, 02) )

!      h ( :, 03) = 0.25 * ( th ( :, 04) - th ( :, 05) )

      h ( :, 04) = 0.25 * ( th ( :, 07) - th ( :, 08) )

!      h ( :, 05) = 0.5 * ( th ( :, 01) + th ( :, 02)
!     1 - th ( :, 03) ) - 2.0 * h_quadpt ( :, 01)

!      h ( :, 06) = 0.5 * ( th ( :, 04) + th ( :, 05)
!     1 - th ( :, 06) ) - 2.0 * h_quadpt ( :, 01)

      h ( :, 07) = 0.5 * ( th ( :, 07) + th ( :, 08)
     1 - th ( :, 09) ) - 2.0 * h_quadpt ( :, 01)

!      h ( :, 08) = h_quadpt ( :, 02) - h_quadpt ( :, 03)
!     1 - h_quadpt ( :, 06) + h_quadpt ( :, 07)

      h ( :, 09) = h_quadpt ( :, 10) - h_quadpt ( :, 11)
     1 - h_quadpt ( :, 12) + h_quadpt ( :, 13)

      h ( :, 10) = h_quadpt ( :, 04) - h_quadpt ( :, 05)
     1 - h_quadpt ( :, 08) + h_quadpt ( :, 09)

!      h ( :, 01) = h_quadpt ( :, 01) + ( h ( :, 05)
!     1 + h ( :, 06) + h ( :, 07) ) / 12.0

! ----------

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal

! ----------

      f ( :, 05) = 2.0 * ( f_quadpt ( :, 02)
     1 - 2.0 * f_quadpt ( :, 01) + f_quadpt ( :, 05) )

!      f ( :, 06) = 2.0 * ( f_quadpt ( :, 06)
!     1 - 2.0 * f_quadpt ( :, 01) + f_quadpt ( :, 09) )

!      f ( :, 07) = 2.0 * ( f_quadpt ( :, 10)
!     1 - 2.0 * f_quadpt ( :, 01) + f_quadpt ( :, 13) )

      f ( :, 08) = 0.5 * ( f_quadpt ( :, 14)
     1 - f_quadpt ( :, 15) + f_quadpt ( :, 16)
     1 - f_quadpt ( :, 17) - f_quadpt ( :, 18)
     1 + f_quadpt ( :, 19) - f_quadpt ( :, 20)
     1 + f_quadpt ( :, 21) )

!      f ( :, 09) = 0.5 * ( f_quadpt ( :, 14)
!     1 - f_quadpt ( :, 16) + f_quadpt ( :, 18)
!     1 - f_quadpt ( :, 20) - f_quadpt ( :, 15)
!     1 + f_quadpt ( :, 17) - f_quadpt ( :, 19)
!     1 + f_quadpt ( :, 21) )

      f ( :, 10) = 0.5 * ( f_quadpt ( :, 14)
     1 - f_quadpt ( :, 18) + f_quadpt ( :, 15)
     1 - f_quadpt ( :, 19) - f_quadpt ( :, 16)
     1 + f_quadpt ( :, 20) - f_quadpt ( :, 17)
     1 + f_quadpt ( :, 21) )

      f ( :, 11) = 16 * ( f_quadpt ( :, 02)
     1 - 2.0 * f_quadpt ( :, 03) + 2.0 * f_quadpt ( :, 04)
     1 - f_quadpt ( :, 05) ) / 3.0

!      f ( :, 12) = 16.0 * ( f_quadpt ( :, 06)
!     1 - 2.0 * f_quadpt ( :, 07) + 2.0 * f_quadpt ( :, 08)
!     1 - f_quadpt ( :, 09) ) / 3.0

!      f ( :, 13) = 16.0 * ( f_quadpt ( :, 10)
!     1 - 2.0 * f_quadpt ( :, 11) + 2.0 * f_quadpt ( :, 12)
!     1 - f_quadpt ( :, 13) ) / 3.0

      tf ( :, 01) = f_quadpt ( :, 14) + f_quadpt ( :, 15)
     1 + f_quadpt ( :, 16) + f_quadpt ( :, 17)

      tf ( :, 02) = f_quadpt ( :, 18) + f_quadpt ( :, 19)
     1 + f_quadpt ( :, 20) + f_quadpt ( :, 21)

      tf ( :, 03) = f_quadpt ( :, 14) + f_quadpt ( :, 16)
     1 + f_quadpt ( :, 18) + f_quadpt ( :, 20)

      tf ( :, 04) = f_quadpt ( :, 15) + f_quadpt ( :, 17)
     1 + f_quadpt ( :, 19) + f_quadpt ( :, 21)

      tf ( :, 05) = f_quadpt ( :, 14) + f_quadpt ( :, 15)
     1 + f_quadpt ( :, 18) + f_quadpt ( :, 19)

      tf ( :, 06) = f_quadpt ( :, 16) + f_quadpt ( :, 17)
     1 + f_quadpt ( :, 20) + f_quadpt ( :, 21)

      f ( :, 02) = 2.0 * ( f_quadpt ( :, 02)
     1 - f_quadpt ( :, 05) ) / 3.0 - 0.1 * f ( :, 11)
     1 + 0.25 * ( tf ( :, 01) - tf ( :, 02) ) / 3.0

!      f ( :, 03) = 2.0 * ( f_quadpt ( :, 06)
!     1 - f_quadpt ( :, 09) ) / 3.0 - 0.1 * f ( :, 12)
!     1 + 0.25 * ( tf ( :, 03) - tf ( :, 04) ) / 3.0

!      f ( :, 04) = 2.0 * ( f_quadpt ( :, 10)
!     1 - f_quadpt ( :, 13) ) / 3.0 - 0.1 * f ( :, 13)
!     1 + 0.25 * ( tf ( :, 05) - tf ( :, 06) ) / 3.0

      f ( :, 20) = f_quadpt ( :, 14) - f_quadpt ( :, 18)
     1 - f_quadpt ( :, 16) + f_quadpt ( :, 20)
     1 - f_quadpt ( :, 15) + f_quadpt ( :, 19)
     1 + f_quadpt ( :, 17) - f_quadpt ( :, 21)

!      f ( :, 01) = 0.125 * ( f_quadpt ( :, 14)
!     1 + f_quadpt ( :, 15) + f_quadpt ( :, 16)
!     1 + f_quadpt ( :, 17) + f_quadpt ( :, 18)
!     1 + f_quadpt ( :, 19) + f_quadpt ( :, 20)
!     1 + f_quadpt ( :, 21) ) - 0.5 * ( f ( :, 05)
!     1 + f ( :, 06) + f ( :, 07) ) / 3.0

! ----------

!      g ( :, 05) = 2.0 * ( g_quadpt ( :, 02)
!     1 - 2.0 * g_quadpt ( :, 01) + g_quadpt ( :, 05) )

      g ( :, 06) = 2.0 * ( g_quadpt ( :, 06)
     1 - 2.0 * g_quadpt ( :, 01) + g_quadpt ( :, 09) )

!      g ( :, 07) = 2.0 * ( g_quadpt ( :, 10)
!     1 - 2.0 * g_quadpt ( :, 01) + g_quadpt ( :, 13) )

      g ( :, 08) = 0.5 * ( g_quadpt ( :, 14)
     1 - g_quadpt ( :, 15) + g_quadpt ( :, 16)
     1 - g_quadpt ( :, 17) - g_quadpt ( :, 18)
     1 + g_quadpt ( :, 19) - g_quadpt ( :, 20)
     1 + g_quadpt ( :, 21) )

      g ( :, 09) = 0.5 * ( g_quadpt ( :, 14)
     1 - g_quadpt ( :, 16) + g_quadpt ( :, 18)
     1 - g_quadpt ( :, 20) - g_quadpt ( :, 15)
     1 + g_quadpt ( :, 17) - g_quadpt ( :, 19)
     1 + g_quadpt ( :, 21) )

!      g ( :, 10) = 0.5 * ( g_quadpt ( :, 14)
!     1 - g_quadpt ( :, 18) + g_quadpt ( :, 15)
!     1 - g_quadpt ( :, 19) - g_quadpt ( :, 16)
!     1 + g_quadpt ( :, 20) - g_quadpt ( :, 17)
!     1 + g_quadpt ( :, 21) )

!      g ( :, 11) = 16 * ( g_quadpt ( :, 02)
!     1 - 2.0 * g_quadpt ( :, 03) + 2.0 * g_quadpt ( :, 04)
!     1 - g_quadpt ( :, 05) ) / 3.0

      g ( :, 12) = 16.0 * ( g_quadpt ( :, 06)
     1 - 2.0 * g_quadpt ( :, 07) + 2.0 * g_quadpt ( :, 08)
     1 - g_quadpt ( :, 09) ) / 3.0

!      g ( :, 13) = 16.0 * ( g_quadpt ( :, 10)
!     1 - 2.0 * g_quadpt ( :, 11) + 2.0 * g_quadpt ( :, 12)
!     1 - g_quadpt ( :, 13) ) / 3.0

      tg ( :, 01) = g_quadpt ( :, 14) + g_quadpt ( :, 15)
     1 + g_quadpt ( :, 16) + g_quadpt ( :, 17)

      tg ( :, 02) = g_quadpt ( :, 18) + g_quadpt ( :, 19)
     1 + g_quadpt ( :, 20) + g_quadpt ( :, 21)

      tg ( :, 03) = g_quadpt ( :, 14) + g_quadpt ( :, 16)
     1 + g_quadpt ( :, 18) + g_quadpt ( :, 20)

      tg ( :, 04) = g_quadpt ( :, 15) + g_quadpt ( :, 17)
     1 + g_quadpt ( :, 19) + g_quadpt ( :, 21)

      tg ( :, 05) = g_quadpt ( :, 14) + g_quadpt ( :, 15)
     1 + g_quadpt ( :, 18) + g_quadpt ( :, 19)

      tg ( :, 06) = g_quadpt ( :, 16) + g_quadpt ( :, 17)
     1 + g_quadpt ( :, 20) + g_quadpt ( :, 21)

!      g ( :, 02) = 2.0 * ( g_quadpt ( :, 02)
!     1 - g_quadpt ( :, 05) ) / 3.0 - 0.1 * g ( :, 11)
!     1 + 0.25 * ( tg ( :, 01) - tg ( :, 02) ) / 3.0

      g ( :, 03) = 2.0 * ( g_quadpt ( :, 06)
     1 - g_quadpt ( :, 09) ) / 3.0 - 0.1 * g ( :, 12)
     1 + 0.25 * ( tg ( :, 03) - tg ( :, 04) ) / 3.0

!      g ( :, 04) = 2.0 * ( g_quadpt ( :, 10)
!     1 - g_quadpt ( :, 13) ) / 3.0 - 0.1 * g ( :, 13)
!     1 + 0.25 * ( tg ( :, 05) - tg ( :, 06) ) / 3.0

      g ( :, 20) = g_quadpt ( :, 14) - g_quadpt ( :, 18)
     1 - g_quadpt ( :, 16) + g_quadpt ( :, 20)
     1 - g_quadpt ( :, 15) + g_quadpt ( :, 19)
     1 + g_quadpt ( :, 17) - g_quadpt ( :, 21)

!      g ( :, 01) = 0.125 * ( g_quadpt ( :, 14)
!     1 + g_quadpt ( :, 15) + g_quadpt ( :, 16)
!     1 + g_quadpt ( :, 17) + g_quadpt ( :, 18)
!     1 + g_quadpt ( :, 19) + g_quadpt ( :, 20)
!     1 + g_quadpt ( :, 21) ) - 0.5 * ( g ( :, 05)
!     1 + g ( :, 06) + g ( :, 07) ) / 3.0

! ----------

!      h ( :, 05) = 2.0 * ( h_quadpt ( :, 02)
!     1 - 2.0 * h_quadpt ( :, 01) + h_quadpt ( :, 05) )

!      h ( :, 06) = 2.0 * ( h_quadpt ( :, 06)
!     1 - 2.0 * h_quadpt ( :, 01) + h_quadpt ( :, 09) )

      h ( :, 07) = 2.0 * ( h_quadpt ( :, 10)
     1 - 2.0 * h_quadpt ( :, 01) + h_quadpt ( :, 13) )

!      h ( :, 08) = 0.5 * ( h_quadpt ( :, 14)
!     1 - h_quadpt ( :, 15) + h_quadpt ( :, 16)
!     1 - h_quadpt ( :, 17) - h_quadpt ( :, 18)
!     1 + h_quadpt ( :, 19) - h_quadpt ( :, 20)
!     1 + h_quadpt ( :, 21) )

      h ( :, 09) = 0.5 * ( h_quadpt ( :, 14)
     1 - h_quadpt ( :, 16) + h_quadpt ( :, 18)
     1 - h_quadpt ( :, 20) - h_quadpt ( :, 15)
     1 + h_quadpt ( :, 17) - h_quadpt ( :, 19)
     1 + h_quadpt ( :, 21) )

      h ( :, 10) = 0.5 * ( h_quadpt ( :, 14)
     1 - h_quadpt ( :, 18) + h_quadpt ( :, 15)
     1 - h_quadpt ( :, 19) - h_quadpt ( :, 16)
     1 + h_quadpt ( :, 20) - h_quadpt ( :, 17)
     1 + h_quadpt ( :, 21) )

!      h ( :, 11) = 16 * ( h_quadpt ( :, 02)
!     1 - 2.0 * h_quadpt ( :, 03) + 2.0 * h_quadpt ( :, 04)
!     1 - h_quadpt ( :, 05) ) / 3.0

!      h ( :, 12) = 16.0 * ( h_quadpt ( :, 06)
!     1 - 2.0 * h_quadpt ( :, 07) + 2.0 * h_quadpt ( :, 08)
!     1 - h_quadpt ( :, 09) ) / 3.0

      h ( :, 13) = 16.0 * ( h_quadpt ( :, 10)
     1 - 2.0 * h_quadpt ( :, 11) + 2.0 * h_quadpt ( :, 12)
     1 - h_quadpt ( :, 13) ) / 3.0

      th ( :, 01) = h_quadpt ( :, 14) + h_quadpt ( :, 15)
     1 + h_quadpt ( :, 16) + h_quadpt ( :, 17)

      th ( :, 02) = h_quadpt ( :, 18) + h_quadpt ( :, 19)
     1 + h_quadpt ( :, 20) + h_quadpt ( :, 21)

      th ( :, 03) = h_quadpt ( :, 14) + h_quadpt ( :, 16)
     1 + h_quadpt ( :, 18) + h_quadpt ( :, 20)

      th ( :, 04) = h_quadpt ( :, 15) + h_quadpt ( :, 17)
     1 + h_quadpt ( :, 19) + h_quadpt ( :, 21)

      th ( :, 05) = h_quadpt ( :, 14) + h_quadpt ( :, 15)
     1 + h_quadpt ( :, 18) + h_quadpt ( :, 19)

      th ( :, 06) = h_quadpt ( :, 16) + h_quadpt ( :, 17)
     1 + h_quadpt ( :, 20) + h_quadpt ( :, 21)

!      h ( :, 02) = 2.0 * ( h_quadpt ( :, 02)
!     1 - h_quadpt ( :, 05) ) / 3.0 - 0.1 * h ( :, 11)
!     1 + 0.25 * ( th ( :, 01) - th ( :, 02) ) / 3.0

!      h ( :, 03) = 2.0 * ( h_quadpt ( :, 06)
!     1 - h_quadpt ( :, 09) ) / 3.0 - 0.1 * h ( :, 12)
!     1 + 0.25 * ( th ( :, 03) - th ( :, 04) ) / 3.0

      h ( :, 04) = 2.0 * ( h_quadpt ( :, 10)
     1 - h_quadpt ( :, 13) ) / 3.0 - 0.1 * h ( :, 13)
     1 + 0.25 * ( th ( :, 05) - th ( :, 06) ) / 3.0

      h ( :, 20) = h_quadpt ( :, 14) - h_quadpt ( :, 18)
     1 - h_quadpt ( :, 16) + h_quadpt ( :, 20)
     1 - h_quadpt ( :, 15) + h_quadpt ( :, 19)
     1 + h_quadpt ( :, 17) - h_quadpt ( :, 21)

!      h ( :, 01) = 0.125 * ( h_quadpt ( :, 14)
!     1 + h_quadpt ( :, 15) + h_quadpt ( :, 16)
!     1 + h_quadpt ( :, 17) + h_quadpt ( :, 18)
!     1 + h_quadpt ( :, 19) + h_quadpt ( :, 20)
!     1 + h_quadpt ( :, 21) ) - 0.5 * ( h ( :, 05)
!     1 + h ( :, 06) + h ( :, 07) ) / 3.0

! ----------

#endif
! End of "NUMORDER == 4"

! ----------------------------------------------------------------------
! Default is to build source terms at t=0 . In some ADER-CG updates
! (on Cartesian meshes), we may not require any source terms.
#if ( BUILD_SOURCE_TERMS == 1)
! ----------------------------------------------------------------------

! STEP 5.a) Build source terms "s_quadpt" at the spatial quadrature points at t=0.
! Remember that the primitive variables are still available from the previous
! loop at each of the t=0 quadrature points for constructing the source terms.

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      s_quadpt ( :, iquadpt) = 0.0

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

#if ( 0 == 1)

! For single-fluid adiabatic MHD turbulence, we may want to provide a gradual
! (explicit) heating/cooling term which ensures that pressure always
! tends to a mean pressure "prs_bar" in a time "tau_bar" (comparable
! to a turnover time). Such a source term can be implemented explicitly
! and its implementation is given below.
! Remember to declare "prs_bar, tau_bar".

      prs_bar = 1.0
      tau_bar = 1.0

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

! ----------

      s_quadpt ( 5, iquadpt) = - ( prs ( iquadpt) - prs_bar)
     1     / ( tau_bar * ( gamma_eos_quadpt ( iquadpt) - 1.0) )

! ----------

      END DO
! Finished "DO iquadpt ="

      END IF

#endif

! ------------------------------

! Cylindrical Geometry

      IF ( igeom .EQ. 2) THEN

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

! ----------

      x1 = nodal_points ( 1, iquadpt)
      x2 = nodal_points ( 2, iquadpt)
      x3 = nodal_points ( 3, iquadpt)

      obx1 = 1.0 / x1

      tempaa = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      tempab = ( vlx ( iquadpt) * bfx ( iquadpt)
     1         + vly ( iquadpt) * bfy ( iquadpt)
     1         + vlz ( iquadpt) * bfz ( iquadpt)) * ob4pi

      tempac = u_quadpt ( 5, iquadpt) + prs ( iquadpt) + tempaa

! ----------

      s_quadpt ( 1, iquadpt) = - rho ( iquadpt)
     1                           * vlx ( iquadpt) * obx1

      s_quadpt ( 2, iquadpt) =
     1 ( rho ( iquadpt) * vly ( iquadpt)**2
     1 - rho ( iquadpt) * vlx ( iquadpt)**2
     1 + bfx ( iquadpt)**2 * ob4pi
     1 - bfy ( iquadpt)**2 * ob4pi) * obx1

      s_quadpt ( 3, iquadpt) =
     1 - 2.0 * ( rho ( iquadpt) * vlx ( iquadpt) 
     1                              * vly ( iquadpt)
     1 - bfx ( iquadpt) * bfy ( iquadpt) * ob4pi) * obx1

      s_quadpt ( 4, iquadpt) =
     1 - ( rho ( iquadpt) * vlx ( iquadpt) 
     1                        * vlz ( iquadpt)
     1 - bfx ( iquadpt) * bfz ( iquadpt) * ob4pi) * obx1

      s_quadpt ( 5, iquadpt) = - ( tempac * vlx ( iquadpt)
     1 - bfx ( iquadpt) * tempab) * obx1

      s_quadpt ( 6, iquadpt) = 0.0

      s_quadpt ( 7, iquadpt) = 0.0

      s_quadpt ( 8, iquadpt) = - f_quadpt ( 8, iquadpt) * obx1

! ----------

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Spherical Geometry

      IF ( igeom .EQ. 3) THEN

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

! ----------

      x1 = nodal_points ( 1, iquadpt)
      x2 = nodal_points ( 2, iquadpt)
      x3 = nodal_points ( 3, iquadpt)

      obx1 = 1.0 / x1

      cos_x2 = cos ( x2)
      sin_x2 = sin ( x2)
      cot_x2 = cos_x2 / sin_x2

      tempaa = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      tempab = ( vlx ( iquadpt) * bfx ( iquadpt)
     1         + vly ( iquadpt) * bfy ( iquadpt)
     1         + vlz ( iquadpt) * bfz ( iquadpt)) * ob4pi

      tempac = u_quadpt ( 5, iquadpt) + prs ( iquadpt) + tempaa

! ----------

      s_quadpt ( 1, iquadpt) = - rho ( iquadpt)
     1 * ( 2.0 * vlx ( iquadpt)
     1 + cot_x2 * vly ( iquadpt) ) * obx1

      s_quadpt ( 2, iquadpt) = rho ( iquadpt)
     1 * ( - 2.0 * vlx ( iquadpt)**2 + vly ( iquadpt)**2
     1   + vlz ( iquadpt)**2) * obx1
     1 - cot_x2
     1 * ( rho ( iquadpt) * vlx ( iquadpt) * vly ( iquadpt)
     1   - bfx ( iquadpt) * bfy ( iquadpt) * ob4pi) * obx1
     1 - ( - 2.0 * bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1   + bfz ( iquadpt)**2) * ob4pi * obx1

      s_quadpt ( 3, iquadpt) = - rho ( iquadpt)
     1 * ( 3.0 * vlx ( iquadpt) * vly ( iquadpt)
     1 + cot_x2
     1 * ( vly ( iquadpt)**2 - vlz ( iquadpt)**2)) * obx1
     1 + ( 3.0 * bfx ( iquadpt) * bfy ( iquadpt)
     1   + cot_x2
     1   * ( bfy ( iquadpt)**2 - bfz ( iquadpt)**2))
     1 * ob4pi * obx1

      s_quadpt ( 4, iquadpt) =
     1 - rho ( iquadpt) * vlz ( iquadpt)
     1 * ( 3.0 * vlx ( iquadpt) + 2.0 * cot_x2
     1   * vly ( iquadpt) ) * obx1
     1 + bfz ( iquadpt)
     1 * ( 3.0 * bfx ( iquadpt) + 2.0 * cot_x2
     1   * bfy ( iquadpt) ) * ob4pi * obx1

      s_quadpt ( 5, iquadpt) =
     1 - 2.0 * ( tempac * vlx ( iquadpt)
     1         - bfx ( iquadpt) * tempab) * obx1
     1 - cot_x2 * ( tempac * vly ( iquadpt)
     1 - bfy ( iquadpt) * tempab) * obx1

      s_quadpt ( 6, iquadpt) = - cot_x2
     1                           * g_quadpt ( 6, iquadpt) * obx1

      s_quadpt ( 7, iquadpt) = g_quadpt ( 6, iquadpt) * obx1

      s_quadpt ( 8, iquadpt) = - f_quadpt ( 8, iquadpt) * obx1

! ----------

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Multiply "s_quadpt" by "dtmove".

      s_quadpt ( :, 1: NUM_ZONECEN_QUADPTS) = dtmove
     1         * s_quadpt ( :, 1: NUM_ZONECEN_QUADPTS)

! ----------------------------------------------------------------------

! STEP 5.b) Now project the sources "s_quadpt" from the spatial nodes at t=0
! back to the t=0 part of the modal space time variables for those variables,
! i.e. to "s".

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal

! ----------

      s ( :, 01) = ( s_quadpt ( :, 01) + s_quadpt ( :, 02)
     1 + s_quadpt ( :, 03) + s_quadpt ( :, 04)
     1 + s_quadpt ( :, 05) + s_quadpt ( :, 06)) / 6.0

      s ( :, 02) = s_quadpt ( :, 01) - s_quadpt ( :, 02)

      s ( :, 03) = s_quadpt ( :, 03) - s_quadpt ( :, 04)

      s ( :, 04) = s_quadpt ( :, 05) - s_quadpt ( :, 06)

! ----------

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal

! ----------

      ts ( :, 01) = s_quadpt ( :, 02) + s_quadpt ( :, 03)
     1           + s_quadpt ( :, 04) + s_quadpt ( :, 05)

      ts ( :, 02) = s_quadpt ( :, 06) + s_quadpt ( :, 07)
     1           + s_quadpt ( :, 08) + s_quadpt ( :, 09)

      ts ( :, 03) = s_quadpt ( :, 10) + s_quadpt ( :, 11)
     1           + s_quadpt ( :, 12) + s_quadpt ( :, 13)

      ts ( :, 04) = s_quadpt ( :, 02) + s_quadpt ( :, 06)
     1           + s_quadpt ( :, 10) + s_quadpt ( :, 12)

      ts ( :, 05) = s_quadpt ( :, 03) + s_quadpt ( :, 07)
     1           + s_quadpt ( :, 11) + s_quadpt ( :, 13)

      ts ( :, 06) = s_quadpt ( :, 04) + s_quadpt ( :, 05)
     1           + s_quadpt ( :, 08) + s_quadpt ( :, 09)

      ts ( :, 07) = s_quadpt ( :, 04) + s_quadpt ( :, 08)
     1           + s_quadpt ( :, 10) + s_quadpt ( :, 11)

      ts ( :, 08) = s_quadpt ( :, 05) + s_quadpt ( :, 09)
     1           + s_quadpt ( :, 12) + s_quadpt ( :, 13)

      ts ( :, 09) = s_quadpt ( :, 02) + s_quadpt ( :, 03)
     1           + s_quadpt ( :, 06) + s_quadpt ( :, 07)

      s ( :, 02) = 0.25 * ( ts ( :, 01) - ts ( :, 02) )

      s ( :, 03) = 0.25 * ( ts ( :, 04) - ts ( :, 05) )

      s ( :, 04) = 0.25 * ( ts ( :, 07) - ts ( :, 08) )

      s ( :, 05) = 0.5 * ( ts ( :, 01) + ts ( :, 02)
     1 - ts ( :, 03) ) - 2.0 * s_quadpt ( :, 01)

      s ( :, 06) = 0.5 * ( ts ( :, 04) + ts ( :, 05)
     1 - ts ( :, 06) ) - 2.0 * s_quadpt ( :, 01)

      s ( :, 07) = 0.5 * ( ts ( :, 07) + ts ( :, 08)
     1 - ts ( :, 09) ) - 2.0 * s_quadpt ( :, 01)

      s ( :, 08) = s_quadpt ( :, 02) - s_quadpt ( :, 03)
     1 - s_quadpt ( :, 06) + s_quadpt ( :, 07)

      s ( :, 09) = s_quadpt ( :, 10) - s_quadpt ( :, 11)
     1 - s_quadpt ( :, 12) + s_quadpt ( :, 13)

      s ( :, 10) = s_quadpt ( :, 04) - s_quadpt ( :, 05)
     1 - s_quadpt ( :, 08) + s_quadpt ( :, 09)

      s ( :, 01) = s_quadpt ( :, 01) + ( s ( :, 05)
     1 + s ( :, 06) + s ( :, 07) ) / 12.0

! ----------

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal

! ----------

      s ( :, 05) = 2.0 * ( s_quadpt ( :, 02)
     1 - 2.0 * s_quadpt ( :, 01) + s_quadpt ( :, 05) )

      s ( :, 06) = 2.0 * ( s_quadpt ( :, 06)
     1 - 2.0 * s_quadpt ( :, 01) + s_quadpt ( :, 09) )

      s ( :, 07) = 2.0 * ( s_quadpt ( :, 10)
     1 - 2.0 * s_quadpt ( :, 01) + s_quadpt ( :, 13) )

      s ( :, 08) = 0.5 * ( s_quadpt ( :, 14)
     1 - s_quadpt ( :, 15) + s_quadpt ( :, 16)
     1 - s_quadpt ( :, 17) - s_quadpt ( :, 18)
     1 + s_quadpt ( :, 19) - s_quadpt ( :, 20)
     1 + s_quadpt ( :, 21) )

      s ( :, 09) = 0.5 * ( s_quadpt ( :, 14)
     1 - s_quadpt ( :, 16) + s_quadpt ( :, 18)
     1 - s_quadpt ( :, 20) - s_quadpt ( :, 15)
     1 + s_quadpt ( :, 17) - s_quadpt ( :, 19)
     1 + s_quadpt ( :, 21) )

      s ( :, 10) = 0.5 * ( s_quadpt ( :, 14)
     1 - s_quadpt ( :, 18) + s_quadpt ( :, 15)
     1 - s_quadpt ( :, 19) - s_quadpt ( :, 16)
     1 + s_quadpt ( :, 20) - s_quadpt ( :, 17)
     1 + s_quadpt ( :, 21) )

      s ( :, 11) = 16 * ( s_quadpt ( :, 02)
     1 - 2.0 * s_quadpt ( :, 03) + 2.0 * s_quadpt ( :, 04)
     1 - s_quadpt ( :, 05) ) / 3.0

      s ( :, 12) = 16.0 * ( s_quadpt ( :, 06)
     1 - 2.0 * s_quadpt ( :, 07) + 2.0 * s_quadpt ( :, 08)
     1 - s_quadpt ( :, 09) ) / 3.0

      s ( :, 13) = 16.0 * ( s_quadpt ( :, 10)
     1 - 2.0 * s_quadpt ( :, 11) + 2.0 * s_quadpt ( :, 12)
     1 - s_quadpt ( :, 13) ) / 3.0

      ts ( :, 01) = s_quadpt ( :, 14) + s_quadpt ( :, 15)
     1 + s_quadpt ( :, 16) + s_quadpt ( :, 17)

      ts ( :, 02) = s_quadpt ( :, 18) + s_quadpt ( :, 19)
     1 + s_quadpt ( :, 20) + s_quadpt ( :, 21)

      ts ( :, 03) = s_quadpt ( :, 14) + s_quadpt ( :, 16)
     1 + s_quadpt ( :, 18) + s_quadpt ( :, 20)

      ts ( :, 04) = s_quadpt ( :, 15) + s_quadpt ( :, 17)
     1 + s_quadpt ( :, 19) + s_quadpt ( :, 21)

      ts ( :, 05) = s_quadpt ( :, 14) + s_quadpt ( :, 15)
     1 + s_quadpt ( :, 18) + s_quadpt ( :, 19)

      ts ( :, 06) = s_quadpt ( :, 16) + s_quadpt ( :, 17)
     1 + s_quadpt ( :, 20) + s_quadpt ( :, 21)

      s ( :, 02) = 2.0 * ( s_quadpt ( :, 02)
     1 - s_quadpt ( :, 05) ) / 3.0 - 0.1 * s ( :, 11)
     1 + 0.25 * ( ts ( :, 01) - ts ( :, 02) ) / 3.0

      s ( :, 03) = 2.0 * ( s_quadpt ( :, 06)
     1 - s_quadpt ( :, 09) ) / 3.0 - 0.1 * s ( :, 12)
     1 + 0.25 * ( ts ( :, 03) - ts ( :, 04) ) / 3.0

      s ( :, 04) = 2.0 * ( s_quadpt ( :, 10)
     1 - s_quadpt ( :, 13) ) / 3.0 - 0.1 * s ( :, 13)
     1 + 0.25 * ( ts ( :, 05) - ts ( :, 06) ) / 3.0

      s ( :, 20) = s_quadpt ( :, 14) - s_quadpt ( :, 18)
     1 - s_quadpt ( :, 16) + s_quadpt ( :, 20)
     1 - s_quadpt ( :, 15) + s_quadpt ( :, 19)
     1 + s_quadpt ( :, 17) - s_quadpt ( :, 21)

      s ( :, 01) = 0.125 * ( s_quadpt ( :, 14)
     1 + s_quadpt ( :, 15) + s_quadpt ( :, 16)
     1 + s_quadpt ( :, 17) + s_quadpt ( :, 18)
     1 + s_quadpt ( :, 19) + s_quadpt ( :, 20)
     1 + s_quadpt ( :, 21) ) - 0.5 * ( s ( :, 05)
     1 + s ( :, 06) + s ( :, 07) ) / 3.0

! ----------

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
#endif
! End of "BUILD_SOURCE_TERMS == 1"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Start the outer ADER iteration here.

      DO iter_ader = 1, max_iter_ader_cg

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 6) Make the conserved variables at the t>0 space-time nodes.
! We put in the flexibility of having different "gamma_eos_spacetime" at each
! space-time node, but that is not used.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : modal -> nodal

! ----------

      u_quadpt ( :, 07) = u ( :, 01) + u ( :, 05)

! ----------

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : modal -> nodal

! ----------

      u_quadpt ( :, 14) = 0.25 * ( u_quadpt ( :, 02)
     1 + u_quadpt ( :, 03) + u_quadpt ( :, 04)
     1 + u_quadpt ( :, 05) + u ( :, 12) + u ( :, 13) )
     1 + 0.5 * u ( :, 11) - 0.125 * ( u ( :, 06)
     1 + u ( :, 07) )

      u_quadpt ( :, 15) = u_quadpt ( :, 14)
     1 - u ( :, 02) - 0.5 * u ( :, 13)

      u_quadpt ( :, 16) = 0.25 * ( u_quadpt ( :, 02)
     1 + u_quadpt ( :, 06) + u_quadpt ( :, 10)
     1 + u_quadpt ( :, 12) + u ( :, 12) + u ( :, 14) )
     1 + 0.5 * u ( :, 11) - 0.125 * ( u ( :, 05)
     1 + u ( :, 07) )

      u_quadpt ( :, 17) = u_quadpt ( :, 16)
     1 - u ( :, 03) - 0.5 * u ( :, 14)

      u_quadpt ( :, 18) = 0.25 * ( u_quadpt ( :, 04)
     1 + u_quadpt ( :, 08) + u_quadpt ( :, 10)
     1 + u_quadpt ( :, 11) + u ( :, 12) + u ( :, 15) )
     1 + 0.5 * u ( :, 11) - 0.125 * ( u ( :, 05)
     1 + u ( :, 06) )

      u_quadpt ( :, 19) = u_quadpt ( :, 18) - u ( :, 04)
     1 - 0.5 * u ( :, 15)

      u_quadpt ( :, 20) = u_quadpt ( :, 01) + u ( :, 11)
     1 + u ( :, 12)

! ----------

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : modal -> nodal

! ----------

      u_quadpt ( :, 22) = u_quadpt ( :, 01)
     1 + ( 36.0 * u ( :, 21) + 12.0 * u ( :, 22)
     1 + 4.0 * u ( :, 23) - 3.0 * u ( :, 30)
     1 - 3.0 * u ( :, 31) - 3.0 * u ( :, 32) ) / 108.0

      u_quadpt ( :, 23) = u_quadpt ( :, 22)
     1 + u_quadpt ( :, 02) + u_quadpt ( :, 06)
     1 - 2.0 * u_quadpt ( :, 01) + ( 18.0 * u ( :, 08)
     1 + 9.0 * u ( :, 14) + 9.0 * u ( :, 16)
     1 + 12.0 * u ( :, 24) + 12.0 * u ( :, 25)
     1 + 4.0 * u ( :, 27) + 4.0 * u ( :, 28)
     1 + 6.0 * u ( :, 30) + 6.0 * u ( :, 31)
     1 + 6.0 * u ( :, 33) ) / 72.0

      u_quadpt ( :, 24) = u_quadpt ( :, 23)
     1 - u_quadpt ( :, 06) + u_quadpt ( :, 09)
     1 - 0.5 * u ( :, 08) - 0.25 * u ( :, 14)
     1 - ( 6.0 * u ( :, 25) + 2.0 * u ( :, 28)
     1 + 3.0 * u ( :, 33) ) / 18.0

      u_quadpt ( :, 25) = 0.5 * ( u_quadpt ( :, 23)
     1 + u_quadpt ( :, 24) + u_quadpt ( :, 14)
     1 + u_quadpt ( :, 15) ) - u_quadpt ( :, 02)
     1 - 0.5 * u ( :, 06) - 0.25 * u ( :, 16)
     1 - 0.125 * u ( :, 17) + ( 6.0 * u ( :, 26)
     1 + 2.0 * u ( :, 29) - 3.0 * u ( :, 31)
     1 + 3.0 * u ( :, 32) + 3.0 * u ( :, 35) ) / 36.0

      u_quadpt ( :, 26) = u_quadpt ( :, 23)
     1 + u_quadpt ( :, 24) - u_quadpt ( :, 25)
     1 - 0.5 * ( u ( :, 06) - u ( :, 07) )
     1 - 0.25 * ( u ( :, 16) - u ( :, 18) )
     1 - ( u ( :, 31) - u ( :, 32) ) / 6.0

      u_quadpt ( :, 27) = u_quadpt ( :, 23)
     1 + u_quadpt ( :, 05) - u_quadpt ( :, 02)
     1 - 0.5 * u ( :, 08) - 0.25 * u ( :, 16)
     1 - ( 6.0 * u ( :, 24) + 2.0 * u ( :, 27)
     1 + 3.0 * u ( :, 33) ) / 18.0

      u_quadpt ( :, 28) = u_quadpt ( :, 27)
     1 - u_quadpt ( :, 06) + u_quadpt ( :, 09)
     1 + 0.5 * u ( :, 08) - 0.25 * u ( :, 14)
     1 - ( 6.0 * u ( :, 25) + 2.0 * u ( :, 28)
     1 - 3.0 * u ( :, 33) ) / 18.0

      u_quadpt ( :, 29) = 0.5 * ( u_quadpt ( :, 27)
     1 + u_quadpt ( :, 28) + u_quadpt ( :, 18)
     1 + u_quadpt ( :, 19) ) - u_quadpt ( :, 05)
     1 - 0.5 * u ( :, 06) + 0.25 * u ( :, 16)
     1 - 0.125 * u ( :, 17) + ( 6.0 * u ( :, 26)
     1 + 2.0 * u ( :, 29) - 3.0 * u ( :, 31)
     1 + 3.0 * u ( :, 32)  - 3.0 * u ( :, 35) ) / 36.0

      u_quadpt ( :, 30) = u_quadpt ( :, 27)
     1 + u_quadpt ( :, 28) - u_quadpt ( :, 29)
     1 - 0.5 * ( u ( :, 06) - u ( :, 07) )
     1 + 0.25 * ( u ( :, 16) - u ( :, 18) )
     1 - ( u ( :, 31) - u ( :, 32) ) / 6.0

      u_quadpt ( :, 31) = 0.5 * ( u_quadpt ( :, 25)
     1 + u_quadpt ( :, 29) + u_quadpt ( :, 14)
     1 + u_quadpt ( :, 18) ) - u_quadpt ( :, 10)
     1 - 0.5 * u ( :, 05) - 0.25 * u ( :, 15)
     1 - 0.125 * u ( :, 14) + ( 6.0 * u ( :, 25)
     1 + 2.0 * u ( :, 28) - 3.0 * u ( :, 30)
     1 + 3.0 * u ( :, 31) + 3.0 * u ( :, 34) ) / 36.0

      u_quadpt ( :, 32) = u_quadpt ( :, 31)
     1 - u_quadpt ( :, 06) + u_quadpt ( :, 09)
     1 - 0.5 * u ( :, 09) - 0.25 * u ( :, 19)
     1 - ( 6.0 * u ( :, 25) + 2.0 * u ( :, 28)
     1 + 3.0 * u ( :, 34) ) / 18.0

      u_quadpt ( :, 33) = u_quadpt ( :, 23)
     1 + u_quadpt ( :, 27) - u_quadpt ( :, 31)
     1 - 0.5 * ( u ( :, 05) - u ( :, 07) )
     1 - 0.25 * ( u ( :, 14) - u ( :, 19) )
     1 - ( u ( :, 30) - u ( :, 32) ) / 6.0

      u_quadpt ( :, 34) = u_quadpt ( :, 24)
     1 + u_quadpt ( :, 28) - u_quadpt ( :, 32)
     1 - 0.5 * ( u ( :, 05) - u ( :, 07) )
     1 + 0.25 * ( u ( :, 14) - u ( :, 19) )
     1 - ( u ( :, 30) - u ( :, 32) ) / 6.0


      tu ( :, 01) = ( 36.0 * u ( :, 21) + 24.0 * u ( :, 22)
     1 + 16.0 * u ( :, 23) ) / 54.0

      tu ( :, 02) = ( 18.0 * u ( :, 24)
     1 + 12.0 * u ( :, 27) ) / 54.0

      tu ( :, 03) = ( 6.0 * u ( :, 30) - 3.0 * u ( :, 31)
     1 - 3.0 * u ( :, 32) ) / 54.0

      tu ( :, 04) = ( 18.0 * u ( :, 25)
     1 + 12.0 * u ( :, 28) ) / 54.0

      tu ( :, 05) = ( - 3.0 * u ( :, 30) + 6.0 * u ( :, 31)
     1 - 3.0 * u ( :, 32) ) / 54.0

      tu ( :, 06) = ( 18.0 * u ( :, 26)
     1 + 12.0 * u ( :, 29) ) / 54.0

      tu ( :, 07) = ( - 3.0 * u ( :, 30) - 3.0 * u ( :, 31)
     1 + 6.0 * u ( :, 32) ) / 54.0

      u_quadpt ( :, 35) = ( 2.0 * u_quadpt ( :, 22)
     1 - u_quadpt ( :, 01) ) + 2.0 * ( u ( :, 22)
     1 + u ( :, 23) ) / 9.0

      u_quadpt ( :, 36) = u_quadpt ( :, 02) + tu ( :, 01)
     1 + tu ( :, 02) + tu ( :, 03)

      u_quadpt ( :, 37) = u_quadpt ( :, 05) + tu ( :, 01)
     1 - tu ( :, 02) + tu ( :, 03)

      u_quadpt ( :, 38) = u_quadpt ( :, 06) + tu ( :, 01)
     1 + tu ( :, 04) + tu ( :, 05)

      u_quadpt ( :, 39) = u_quadpt ( :, 09) + tu ( :, 01)
     1 - tu ( :, 04) + tu ( :, 05)

      u_quadpt ( :, 40) = u_quadpt ( :, 10) + tu ( :, 01)
     1 + tu ( :, 06) + tu ( :, 07)

      u_quadpt ( :, 41) = u_quadpt ( :, 13) + tu ( :, 01)
     1 - tu ( :, 06) + tu ( :, 07)


      u_quadpt ( :, 42) = u_quadpt ( :, 01) + u ( :, 21)
     1 + u ( :, 22) + u ( :, 23) - ( u ( :, 30)
     1 + u ( :, 31) + u ( :, 32) ) / 12.0

! ----------

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 7) Get the primitive variables from the conserved variables at
! the spatial nodes at t>0. The density and pressure variables are
! given positive floor values to ensure that they can be used for physically
! meaningful source terms.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

! ----------

      u_quadpt ( 1, iquadpt) =
     1                   AMAX1 ( u_quadpt ( 1, iquadpt), rhofloor)

      rho ( iquadpt) = u_quadpt ( 1, iquadpt)

      tempaa = 1.0 / rho ( iquadpt)

      vlx ( iquadpt) = u_quadpt ( 2, iquadpt) * tempaa
      vly ( iquadpt) = u_quadpt ( 3, iquadpt) * tempaa
      vlz ( iquadpt) = u_quadpt ( 4, iquadpt) * tempaa

      bfx ( iquadpt) = u_quadpt ( 6, iquadpt)
      bfy ( iquadpt) = u_quadpt ( 7, iquadpt)
      bfz ( iquadpt) = u_quadpt ( 8, iquadpt)

      tempab = 0.5 * rho ( iquadpt) * ( vlx ( iquadpt)**2
     1         + vly ( iquadpt)**2 + vlz ( iquadpt)**2)

      tempac = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      eng_t ( iquadpt) = u_quadpt ( 5, iquadpt)
     1                     - tempab - tempac

      prs ( iquadpt) = eng_t ( iquadpt)
     1                   * ( gamma_eos_quadpt ( iquadpt) - 1.0)

#if ( ENFORCE_POSITIVITY == 1)
      IF ( prs ( iquadpt) .LE. prsfloor) THEN

      prs ( iquadpt) = entropy ( iquadpt)
     1* rho ( iquadpt)**( gamma_eos_quadpt ( iquadpt) - 1.0)
      prs ( iquadpt) = AMAX1 ( prs ( iquadpt), prsfloor)

      eng_t ( iquadpt) = prs ( iquadpt)
     1                 / ( gamma_eos_quadpt ( iquadpt) - 1.0)

      u_quadpt ( 5, iquadpt) = eng_t ( iquadpt)
     1                           + tempab + tempac

      END IF
#else
      IF ( prs ( iquadpt) .LE. prsfloor) THEN

      prs ( iquadpt) = prsfloor

      eng_t ( iquadpt) = prs ( iquadpt)
     1                 / ( gamma_eos_quadpt ( iquadpt) - 1.0)

      u_quadpt ( 5, iquadpt) = eng_t ( iquadpt)
     1                           + tempab + tempac

      END IF
#endif

#if ( NFLUID >= 1)
      specie ( 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( 1: NFLUID, iquadpt) * rho ( iquadpt)
#endif

! ----------

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 8.a) Build fluxes ( f_quadpt, g_quadpt, h_quadpt) from the primitive
! variables that have now become available at spatial nodes at t=0.

! Multiply "f_quadpt" by "dtmove/dx" and so on.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

! -----

      tempaa = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      tempab = ( vlx ( iquadpt) * bfx ( iquadpt)
     1         + vly ( iquadpt) * bfy ( iquadpt)
     1         + vlz ( iquadpt) * bfz ( iquadpt)) * ob4pi

! -----

      f_quadpt ( 1, iquadpt) = u_quadpt ( 2, iquadpt)

      f_quadpt ( 2, iquadpt) = u_quadpt ( 2, iquadpt)
     1         * vlx ( iquadpt)
     1         + prs ( iquadpt)
     1         + tempaa - bfx ( iquadpt)**2 * ob4pi

      f_quadpt ( 3, iquadpt) = u_quadpt ( 2, iquadpt)
     1         * vly ( iquadpt)
     1         - bfx ( iquadpt) * bfy ( iquadpt) * ob4pi

      f_quadpt ( 4, iquadpt) = u_quadpt ( 2, iquadpt)
     1         * vlz ( iquadpt)
     1         - bfx ( iquadpt) * bfz ( iquadpt) * ob4pi

      f_quadpt ( 5, iquadpt) = ( u_quadpt ( 5, iquadpt)
     1         + prs ( iquadpt) + tempaa) * vlx ( iquadpt)
     1         - bfx ( iquadpt) * tempab

      f_quadpt ( 6, iquadpt) = 0.0

      f_quadpt ( 7, iquadpt) =
     1           vlx ( iquadpt) * bfy ( iquadpt)
     1         - vly ( iquadpt) * bfx ( iquadpt)

      f_quadpt ( 8, iquadpt) =
     1           vlx ( iquadpt) * bfz ( iquadpt)
     1         - vlz ( iquadpt) * bfx ( iquadpt)

#if ( NFLUID >= 1)
      f_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlx ( iquadpt)
#endif

! -----

      g_quadpt ( 1, iquadpt) = u_quadpt ( 3, iquadpt)

      g_quadpt ( 2, iquadpt) = u_quadpt ( 3, iquadpt)
     1         * vlx ( iquadpt)
     1         - bfy ( iquadpt) * bfx ( iquadpt) * ob4pi

      g_quadpt ( 3, iquadpt) = u_quadpt ( 3, iquadpt)
     1         * vly ( iquadpt)
     1         + prs ( iquadpt)
     1         + tempaa - bfy ( iquadpt)**2 * ob4pi

      g_quadpt ( 4, iquadpt) = u_quadpt ( 3, iquadpt)
     1         * vlz ( iquadpt)
     1         - bfy ( iquadpt) * bfz ( iquadpt) * ob4pi

      g_quadpt ( 5, iquadpt) = ( u_quadpt ( 5, iquadpt)
     1         + prs ( iquadpt) + tempaa) * vly ( iquadpt)
     1                       - bfy ( iquadpt) * tempab

      g_quadpt ( 6, iquadpt) =
     1           vly ( iquadpt) * bfx ( iquadpt)
     1         - vlx ( iquadpt) * bfy ( iquadpt)

      g_quadpt ( 7, iquadpt) = 0.0

      g_quadpt ( 8, iquadpt) =
     1           vly ( iquadpt) * bfz ( iquadpt)
     1         - vlz ( iquadpt) * bfy ( iquadpt)

#if ( NFLUID >= 1)
      g_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vly ( iquadpt)
#endif

! -----

      h_quadpt ( 1, iquadpt) = u_quadpt ( 4, iquadpt)

      h_quadpt ( 2, iquadpt) = u_quadpt ( 4, iquadpt)
     1         * vlx ( iquadpt)
     1         - bfz ( iquadpt) * bfx ( iquadpt) * ob4pi

      h_quadpt ( 3, iquadpt) = u_quadpt ( 4, iquadpt)
     1         * vly ( iquadpt)
     1         - bfz ( iquadpt) * bfy ( iquadpt) * ob4pi

      h_quadpt ( 4, iquadpt) = u_quadpt ( 4, iquadpt)
     1         * vlz ( iquadpt)
     1         + prs ( iquadpt)
     1         + tempaa - bfz ( iquadpt)**2 * ob4pi

      h_quadpt ( 5, iquadpt) = ( u_quadpt ( 5, iquadpt)
     1         + prs ( iquadpt) + tempaa) * vlz ( iquadpt)
     1                       - bfz ( iquadpt) * tempab

      h_quadpt ( 6, iquadpt) =
     1           vlz ( iquadpt) * bfx ( iquadpt)
     1         - vlx ( iquadpt) * bfz ( iquadpt)

      h_quadpt ( 7, iquadpt) =
     1           vlz ( iquadpt) * bfy ( iquadpt)
     1         - vly ( iquadpt) * bfz ( iquadpt)

      h_quadpt ( 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      h_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt (  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlz ( iquadpt)
#endif

! -----

#if ( ISOTHERMAL_EQNS == 1)
      f_quadpt ( 5, iquadpt) = 0.0
      g_quadpt ( 5, iquadpt) = 0.0
      h_quadpt ( 5, iquadpt) = 0.0
#endif

! -----

      END DO
! Finished "DO iquadpt ="

! ------------------------------

! To make cylindrical and spherical geometries accessible to non-rel.
! and relativistic flow, we put this part after either flux construction.

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      x1 = nodal_points ( 1, iquadpt)
      x2 = nodal_points ( 2, iquadpt)
      x3 = nodal_points ( 3, iquadpt)

      obx1 = 1.0 / x1

      IF ( igeom == 2 .or. igeom == 3) THEN

         g_quadpt ( :, iquadpt) = g_quadpt ( :, iquadpt) * obx1

      END IF

      IF ( igeom == 3) THEN

         sin_x2 = 1.0 / ( x1 * sin ( x2))

         h_quadpt ( :, iquadpt) = h_quadpt ( :, iquadpt) * sin_x2

      END IF

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      tempaa = dtmove / dx

      f_quadpt ( :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS) =
     1   tempaa
     1 * f_quadpt ( :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS)

      tempab = dtmove / dy

      g_quadpt ( :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS) =
     1   tempab
     1 * g_quadpt ( :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS)

      tempac = dtmove / dz

      h_quadpt ( :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS) =
     1   tempac
     1 * h_quadpt ( :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS)

! ----------------------------------------------------------------------

! STEP 8.b) Now project the fluxes "( f_quadpt, g_quadpt, h_quadpt)"
! from the t>0 space-time nodal points back to the modal
! space time variables for those variables, i.e. to "( f, g, h)" .

! Note that some of the "f/g/h_quadpt" values with t>0 will be overwritten.
! But that is ok since these are also constructed anew in each ADER iteration.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal

! ----------

!      f ( :, 05) = f_quadpt ( :, 07) - f ( :, 01)

! ----------

!      g ( :, 05) = g_quadpt ( :, 07) - g ( :, 01)

! ----------

!      h ( :, 05) = h_quadpt ( :, 07) - h ( :, 01)

! ----------

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal

! ----------

      f ( :, 13) = 2.0 * ( f_quadpt ( :, 14)
     1 - f_quadpt ( :, 15) - f ( :, 02) )

!      f ( :, 14) = 2.0 * ( f_quadpt ( :, 16)
!     1 - f_quadpt ( :, 17) - f ( :, 03) )

!      f ( :, 15) = 2.0 * ( f_quadpt ( :, 18)
!     1 - f_quadpt ( :, 19) - f ( :, 04) )

!      f ( :, 11) = 2.0 * ( f_quadpt ( :, 14)
!     1 + f_quadpt ( :, 15) + f_quadpt ( :, 16)
!     1 + f_quadpt ( :, 17) + f_quadpt ( :, 18)
!     1 + f_quadpt ( :, 19) - 6.0 * f ( :, 01) ) / 3.0
!     1 - f_quadpt ( :, 20) + f_quadpt ( :, 01)

!      f ( :, 12) = f_quadpt ( :, 20) - f_quadpt ( :, 01)
!     1 - f ( :, 11)

! ----------

!      g ( :, 13) = 2.0 * ( g_quadpt ( :, 14)
!     1 - g_quadpt ( :, 15) - g ( :, 02) )

      g ( :, 14) = 2.0 * ( g_quadpt ( :, 16)
     1 - g_quadpt ( :, 17) - g ( :, 03) )

!      g ( :, 15) = 2.0 * ( g_quadpt ( :, 18)
!     1 - g_quadpt ( :, 19) - g ( :, 04) )

!      g ( :, 11) = 2.0 * ( g_quadpt ( :, 14)
!     1 + g_quadpt ( :, 15) + g_quadpt ( :, 16)
!     1 + g_quadpt ( :, 17) + g_quadpt ( :, 18)
!     1 + g_quadpt ( :, 19) - 6.0 * g ( :, 01) ) / 3.0
!     1 - g_quadpt ( :, 20) + g_quadpt ( :, 01)

!      g ( :, 12) = g_quadpt ( :, 20) - g_quadpt ( :, 01)
!     1 - g ( :, 11)

! ----------

!      h ( :, 13) = 2.0 * ( h_quadpt ( :, 14)
!     1 - h_quadpt ( :, 15) - h ( :, 02) )

!      h ( :, 14) = 2.0 * ( h_quadpt ( :, 16)
!     1 - h_quadpt ( :, 17) - h ( :, 03) )

      h ( :, 15) = 2.0 * ( h_quadpt ( :, 18)
     1 - h_quadpt ( :, 19) - h ( :, 04) )

!      h ( :, 11) = 2.0 * ( h_quadpt ( :, 14)
!     1 + h_quadpt ( :, 15) + h_quadpt ( :, 16)
!     1 + h_quadpt ( :, 17) + h_quadpt ( :, 18)
!     1 + h_quadpt ( :, 19) - 6.0 * h ( :, 01) ) / 3.0
!     1 - h_quadpt ( :, 20) + h_quadpt ( :, 01)

!      h ( :, 12) = h_quadpt ( :, 20) - h_quadpt ( :, 01)
!     1 - h ( :, 11)

! ----------

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal

! ----------

      tf ( :, 07) = f_quadpt ( :, 23) + f_quadpt ( :, 24)
     1 - f_quadpt ( :, 25) - f_quadpt ( :, 26)
     1 - f_quadpt ( :, 27) - f_quadpt ( :, 28)
     1 + f_quadpt ( :, 29) + f_quadpt ( :, 30)

      tf ( :, 08) = f_quadpt ( :, 32) + f_quadpt ( :, 34)
     1 - f_quadpt ( :, 24) - f_quadpt ( :, 28)
     1 - f_quadpt ( :, 31) - f_quadpt ( :, 33)
     1 + f_quadpt ( :, 23) + f_quadpt ( :, 27)

      tf ( :, 09) = f_quadpt ( :, 25) + f_quadpt ( :, 29)
     1 - f_quadpt ( :, 31) - f_quadpt ( :, 32)
     1 - f_quadpt ( :, 26) - f_quadpt ( :, 30)
     1 + f_quadpt ( :, 33) + f_quadpt ( :, 34)

      f ( :, 16) = tf ( :, 07) + 0.5 * ( tf ( :, 01)
     1 - tf ( :, 02) ) - 2.0 * ( f_quadpt ( :, 02)
     1 - f_quadpt ( :, 05) )

      f ( :, 18) = - 2.0 * tf ( :, 07) + f ( :, 16)

      f ( :, 14) = tf ( :, 08) + 0.5 * ( tf ( :, 03)
     1 - tf ( :, 04) ) - 2.0 * ( f_quadpt ( :, 06)
     1 - f_quadpt ( :, 09) )

      f ( :, 19) = - 2.0 * tf ( :, 08) + f ( :, 14)

      f ( :, 15) = tf ( :, 09) + 0.5 * ( tf ( :, 05)
     1 - tf ( :, 06) ) - 2.0 * ( f_quadpt ( :, 10)
     1 - f_quadpt ( :, 13) )

      f ( :, 17) = - 2.0 * tf ( :, 09) + f ( :, 15)


      tf ( :, 10) = f_quadpt ( :, 23) + f_quadpt ( :, 24)
     1 + f_quadpt ( :, 25) + f_quadpt ( :, 26)
     1 - f_quadpt ( :, 27) - f_quadpt ( :, 28)
     1 - f_quadpt ( :, 29) - f_quadpt ( :, 30)

!      tf ( :, 11) = f_quadpt ( :, 23) + f_quadpt ( :, 27)
!     1 + f_quadpt ( :, 31) + f_quadpt ( :, 33)
!     1 - f_quadpt ( :, 24) - f_quadpt ( :, 28)
!     1 - f_quadpt ( :, 32) - f_quadpt ( :, 34)

!      tf ( :, 12) = f_quadpt ( :, 25) + f_quadpt ( :, 29)
!     1 + f_quadpt ( :, 31) + f_quadpt ( :, 32)
!     1 - f_quadpt ( :, 26) - f_quadpt ( :, 30)
!     1 - f_quadpt ( :, 33) - f_quadpt ( :, 34)

      f ( :, 27) = 1.125 * ( tf ( :, 01) - tf ( :, 02) )
     1 + 4.5 * ( f_quadpt ( :, 36) - f_quadpt ( :, 37) )
     1 - 2.25 * tf ( :, 10)

      f ( :, 24) = 3.0 * ( f_quadpt ( :, 36)
     1 - f_quadpt ( :, 37) ) - 0.75 * tf ( :, 10)
     1 + 0.375 * ( f ( :, 16) + f ( :, 18) )
     1 - f ( :, 27)

!      f ( :, 28) = 1.125 * ( tf ( :, 03) - tf ( :, 04) )
!     1 + 4.5 * ( f_quadpt ( :, 38) - f_quadpt ( :, 39) )
!     1 - 2.25 * tf ( :, 11)

!      f ( :, 25) = 3.0 * ( f_quadpt ( :, 38)
!     1 - f_quadpt ( :, 39) ) - 0.75 * tf ( :, 11)
!     1 + 0.375 * ( f ( :, 14) + f ( :, 19) )
!     1 - f ( :, 28)

!      f ( :, 29) = 1.125 * ( tf ( :, 05) - tf ( :, 06) )
!     1 + 4.5 * ( f_quadpt ( :, 40) - f_quadpt ( :, 41) )
!     1 - 2.25 * tf ( :, 12)

!      f ( :, 26) = 3.0 * ( f_quadpt ( :, 40)
!     1 - f_quadpt ( :, 41) ) - 0.75 * tf ( :, 12)
!     1 + 0.375 * ( f ( :, 15) + f ( :, 17) )
!     1 - f ( :, 29)


      f ( :, 30) = 3.0 * ( f_quadpt ( :, 36)
     1 - 2.0 * f_quadpt ( :, 35) + f_quadpt ( :, 37) )
     1 - 1.5 * f ( :, 05)

!      f ( :, 31) = 3.0 * ( f_quadpt ( :, 38)
!     1 - 2.0 * f_quadpt ( :, 35) + f_quadpt ( :, 39) )
!     1 - 1.5 * f ( :, 06)

!      f ( :, 32) = 3.0 * ( f_quadpt ( :, 40)
!     1 - 2.0 * f_quadpt ( :, 35) + f_quadpt ( :, 41) )
!     1 - 1.5 * f ( :, 07)

      f ( :, 33) = 3.0 * ( f_quadpt ( :, 23)
     1 - f_quadpt ( :, 27) - f_quadpt ( :, 24)
     1 + f_quadpt ( :, 28) ) - 3.0 * f ( :, 08)

!      f ( :, 34) = 3.0 * ( f_quadpt ( :, 31)
!     1 - f_quadpt ( :, 33) - f_quadpt ( :, 32)
!     1 + f_quadpt ( :, 34) ) - 3.0 * f ( :, 09)

      f ( :, 35) = 3.0 * ( f_quadpt ( :, 25)
     1 - f_quadpt ( :, 29) - f_quadpt ( :, 26)
     1 + f_quadpt ( :, 30) ) - 3.0 * f ( :, 10)


!      f ( :, 23) = 4.5 * ( f_quadpt ( :, 42)
!     1 - f_quadpt ( :, 01) ) - 13.5 * ( f_quadpt ( :, 35)
!     1 - f_quadpt ( :, 22) )

!      f ( :, 22) = 2.25 * ( f_quadpt ( :, 42)
!     1 + f_quadpt ( :, 01) - f_quadpt ( :, 35)
!     1 - f_quadpt ( :, 22) ) - 1.5 * f ( :, 23)

!      f ( :, 21) = f_quadpt ( :, 42) - f_quadpt ( :, 01)
!     1 - f ( :, 22) - f ( :, 23)
!     1 + ( f ( :, 30) + f ( :, 31)
!     1 + f ( :, 32) ) / 12.0

! ----------

      tg ( :, 07) = g_quadpt ( :, 23) + g_quadpt ( :, 24)
     1 - g_quadpt ( :, 25) - g_quadpt ( :, 26)
     1 - g_quadpt ( :, 27) - g_quadpt ( :, 28)
     1 + g_quadpt ( :, 29) + g_quadpt ( :, 30)

      tg ( :, 08) = g_quadpt ( :, 32) + g_quadpt ( :, 34)
     1 - g_quadpt ( :, 24) - g_quadpt ( :, 28)
     1 - g_quadpt ( :, 31) - g_quadpt ( :, 33)
     1 + g_quadpt ( :, 23) + g_quadpt ( :, 27)

      tg ( :, 09) = g_quadpt ( :, 25) + g_quadpt ( :, 29)
     1 - g_quadpt ( :, 31) - g_quadpt ( :, 32)
     1 - g_quadpt ( :, 26) - g_quadpt ( :, 30)
     1 + g_quadpt ( :, 33) + g_quadpt ( :, 34)

      g ( :, 16) = tg ( :, 07) + 0.5 * ( tg ( :, 01)
     1 - tg ( :, 02) ) - 2.0 * ( g_quadpt ( :, 02)
     1 - g_quadpt ( :, 05) )

      g ( :, 18) = - 2.0 * tg ( :, 07) + g ( :, 16)

      g ( :, 14) = tg ( :, 08) + 0.5 * ( tg ( :, 03)
     1 - tg ( :, 04) ) - 2.0 * ( g_quadpt ( :, 06)
     1 - g_quadpt ( :, 09) )

      g ( :, 19) = - 2.0 * tg ( :, 08) + g ( :, 14)

      g ( :, 15) = tg ( :, 09) + 0.5 * ( tg ( :, 05)
     1 - tg ( :, 06) ) - 2.0 * ( g_quadpt ( :, 10)
     1 - g_quadpt ( :, 13) )

      g ( :, 17) = - 2.0 * tg ( :, 09) + g ( :, 15)


!      tg ( :, 10) = g_quadpt ( :, 23) + g_quadpt ( :, 24)
!     1 + g_quadpt ( :, 25) + g_quadpt ( :, 26)
!     1 - g_quadpt ( :, 27) - g_quadpt ( :, 28)
!     1 - g_quadpt ( :, 29) - g_quadpt ( :, 30)

      tg ( :, 11) = g_quadpt ( :, 23) + g_quadpt ( :, 27)
     1 + g_quadpt ( :, 31) + g_quadpt ( :, 33)
     1 - g_quadpt ( :, 24) - g_quadpt ( :, 28)
     1 - g_quadpt ( :, 32) - g_quadpt ( :, 34)

!      tg ( :, 12) = g_quadpt ( :, 25) + g_quadpt ( :, 29)
!     1 + g_quadpt ( :, 31) + g_quadpt ( :, 32)
!     1 - g_quadpt ( :, 26) - g_quadpt ( :, 30)
!     1 - g_quadpt ( :, 33) - g_quadpt ( :, 34)

!      g ( :, 27) = 1.125 * ( tg ( :, 01) - tg ( :, 02) )
!     1 + 4.5 * ( g_quadpt ( :, 36) - g_quadpt ( :, 37) )
!     1 - 2.25 * tg ( :, 10)

!      g ( :, 24) = 3.0 * ( g_quadpt ( :, 36)
!     1 - g_quadpt ( :, 37) ) - 0.75 * tg ( :, 10)
!     1 + 0.375 * ( g ( :, 16) + g ( :, 18) )
!     1 - g ( :, 27)

      g ( :, 28) = 1.125 * ( tg ( :, 03) - tg ( :, 04) )
     1 + 4.5 * ( g_quadpt ( :, 38) - g_quadpt ( :, 39) )
     1 - 2.25 * tg ( :, 11)

      g ( :, 25) = 3.0 * ( g_quadpt ( :, 38)
     1 - g_quadpt ( :, 39) ) - 0.75 * tg ( :, 11)
     1 + 0.375 * ( g ( :, 14) + g ( :, 19) )
     1 - g ( :, 28)

!      g ( :, 29) = 1.125 * ( tg ( :, 05) - tg ( :, 06) )
!     1 + 4.5 * ( g_quadpt ( :, 40) - g_quadpt ( :, 41) )
!     1 - 2.25 * tg ( :, 12)

!      g ( :, 26) = 3.0 * ( g_quadpt ( :, 40)
!     1 - g_quadpt ( :, 41) ) - 0.75 * tg ( :, 12)
!     1 + 0.375 * ( g ( :, 15) + g ( :, 17) )
!     1 - g ( :, 29)


!      g ( :, 30) = 3.0 * ( g_quadpt ( :, 36)
!     1 - 2.0 * g_quadpt ( :, 35) + g_quadpt ( :, 37) )
!     1 - 1.5 * g ( :, 05)

      g ( :, 31) = 3.0 * ( g_quadpt ( :, 38)
     1 - 2.0 * g_quadpt ( :, 35) + g_quadpt ( :, 39) )
     1 - 1.5 * g ( :, 06)

!      g ( :, 32) = 3.0 * ( g_quadpt ( :, 40)
!     1 - 2.0 * g_quadpt ( :, 35) + g_quadpt ( :, 41) )
!     1 - 1.5 * g ( :, 07)

      g ( :, 33) = 3.0 * ( g_quadpt ( :, 23)
     1 - g_quadpt ( :, 27) - g_quadpt ( :, 24)
     1 + g_quadpt ( :, 28) ) - 3.0 * g ( :, 08)

      g ( :, 34) = 3.0 * ( g_quadpt ( :, 31)
     1 - g_quadpt ( :, 33) - g_quadpt ( :, 32)
     1 + g_quadpt ( :, 34) ) - 3.0 * g ( :, 09)

!      g ( :, 35) = 3.0 * ( g_quadpt ( :, 25)
!     1 - g_quadpt ( :, 29) - g_quadpt ( :, 26)
!     1 + g_quadpt ( :, 30) ) - 3.0 * g ( :, 10)


!      g ( :, 23) = 4.5 * ( g_quadpt ( :, 42)
!     1 - g_quadpt ( :, 01) ) - 13.5 * ( g_quadpt ( :, 35)
!     1 - g_quadpt ( :, 22) )

!      g ( :, 22) = 2.25 * ( g_quadpt ( :, 42)
!     1 + g_quadpt ( :, 01) - g_quadpt ( :, 35)
!     1 - g_quadpt ( :, 22) ) - 1.5 * g ( :, 23)

!      g ( :, 21) = g_quadpt ( :, 42) - g_quadpt ( :, 01)
!     1 - g ( :, 22) - g ( :, 23)
!     1 + ( g ( :, 30) + g ( :, 31)
!     1 + g ( :, 32) ) / 12.0

! ----------

      th ( :, 07) = h_quadpt ( :, 23) + h_quadpt ( :, 24)
     1 - h_quadpt ( :, 25) - h_quadpt ( :, 26)
     1 - h_quadpt ( :, 27) - h_quadpt ( :, 28)
     1 + h_quadpt ( :, 29) + h_quadpt ( :, 30)

      th ( :, 08) = h_quadpt ( :, 32) + h_quadpt ( :, 34)
     1 - h_quadpt ( :, 24) - h_quadpt ( :, 28)
     1 - h_quadpt ( :, 31) - h_quadpt ( :, 33)
     1 + h_quadpt ( :, 23) + h_quadpt ( :, 27)

      th ( :, 09) = h_quadpt ( :, 25) + h_quadpt ( :, 29)
     1 - h_quadpt ( :, 31) - h_quadpt ( :, 32)
     1 - h_quadpt ( :, 26) - h_quadpt ( :, 30)
     1 + h_quadpt ( :, 33) + h_quadpt ( :, 34)

      h ( :, 16) = th ( :, 07) + 0.5 * ( th ( :, 01)
     1 - th ( :, 02) ) - 2.0 * ( h_quadpt ( :, 02)
     1 - h_quadpt ( :, 05) )

      h ( :, 18) = - 2.0 * th ( :, 07) + h ( :, 16)

      h ( :, 14) = th ( :, 08) + 0.5 * ( th ( :, 03)
     1 - th ( :, 04) ) - 2.0 * ( h_quadpt ( :, 06)
     1 - h_quadpt ( :, 09) )

      h ( :, 19) = - 2.0 * th ( :, 08) + h ( :, 14)

      h ( :, 15) = th ( :, 09) + 0.5 * ( th ( :, 05)
     1 - th ( :, 06) ) - 2.0 * ( h_quadpt ( :, 10)
     1 - h_quadpt ( :, 13) )

      h ( :, 17) = - 2.0 * th ( :, 09) + h ( :, 15)


!      th ( :, 10) = h_quadpt ( :, 23) + h_quadpt ( :, 24)
!     1 + h_quadpt ( :, 25) + h_quadpt ( :, 26)
!     1 - h_quadpt ( :, 27) - h_quadpt ( :, 28)
!     1 - h_quadpt ( :, 29) - h_quadpt ( :, 30)

!      th ( :, 11) = h_quadpt ( :, 23) + h_quadpt ( :, 27)
!     1 + h_quadpt ( :, 31) + h_quadpt ( :, 33)
!     1 - h_quadpt ( :, 24) - h_quadpt ( :, 28)
!     1 - h_quadpt ( :, 32) - h_quadpt ( :, 34)

      th ( :, 12) = h_quadpt ( :, 25) + h_quadpt ( :, 29)
     1 + h_quadpt ( :, 31) + h_quadpt ( :, 32)
     1 - h_quadpt ( :, 26) - h_quadpt ( :, 30)
     1 - h_quadpt ( :, 33) - h_quadpt ( :, 34)

!      h ( :, 27) = 1.125 * ( th ( :, 01) - th ( :, 02) )
!     1 + 4.5 * ( h_quadpt ( :, 36) - h_quadpt ( :, 37) )
!     1 - 2.25 * th ( :, 10)

!      h ( :, 24) = 3.0 * ( h_quadpt ( :, 36)
!     1 - h_quadpt ( :, 37) ) - 0.75 * th ( :, 10)
!     1 + 0.375 * ( h ( :, 16) + h ( :, 18) )
!     1 - h ( :, 27)

!      h ( :, 28) = 1.125 * ( th ( :, 03) - th ( :, 04) )
!     1 + 4.5 * ( h_quadpt ( :, 38) - h_quadpt ( :, 39) )
!     1 - 2.25 * th ( :, 11)

!      h ( :, 25) = 3.0 * ( h_quadpt ( :, 38)
!     1 - h_quadpt ( :, 39) ) - 0.75 * th ( :, 11)
!     1 + 0.375 * ( h ( :, 14) + h ( :, 19) )
!     1 - h ( :, 28)

      h ( :, 29) = 1.125 * ( th ( :, 05) - th ( :, 06) )
     1 + 4.5 * ( h_quadpt ( :, 40) - h_quadpt ( :, 41) )
     1 - 2.25 * th ( :, 12)

      h ( :, 26) = 3.0 * ( h_quadpt ( :, 40)
     1 - h_quadpt ( :, 41) ) - 0.75 * th ( :, 12)
     1 + 0.375 * ( h ( :, 15) + h ( :, 17) )
     1 - h ( :, 29)


!      h ( :, 30) = 3.0 * ( h_quadpt ( :, 36)
!     1 - 2.0 * h_quadpt ( :, 35) + h_quadpt ( :, 37) )
!     1 - 1.5 * h ( :, 05)

!      h ( :, 31) = 3.0 * ( h_quadpt ( :, 38)
!     1 - 2.0 * h_quadpt ( :, 35) + h_quadpt ( :, 39) )
!     1 - 1.5 * h ( :, 06)

      h ( :, 32) = 3.0 * ( h_quadpt ( :, 40)
     1 - 2.0 * h_quadpt ( :, 35) + h_quadpt ( :, 41) )
     1 - 1.5 * h ( :, 07)

!      h ( :, 33) = 3.0 * ( h_quadpt ( :, 23)
!     1 - h_quadpt ( :, 27) - h_quadpt ( :, 24)
!     1 + h_quadpt ( :, 28) ) - 3.0 * h ( :, 08)

      h ( :, 34) = 3.0 * ( h_quadpt ( :, 31)
     1 - h_quadpt ( :, 33) - h_quadpt ( :, 32)
     1 + h_quadpt ( :, 34) ) - 3.0 * h ( :, 09)

      h ( :, 35) = 3.0 * ( h_quadpt ( :, 25)
     1 - h_quadpt ( :, 29) - h_quadpt ( :, 26)
     1 + h_quadpt ( :, 30) ) - 3.0 * h ( :, 10)


!      h ( :, 23) = 4.5 * ( h_quadpt ( :, 42)
!     1 - h_quadpt ( :, 01) ) - 13.5 * ( h_quadpt ( :, 35)
!     1 - h_quadpt ( :, 22) )

!      h ( :, 22) = 2.25 * ( h_quadpt ( :, 42)
!     1 + h_quadpt ( :, 01) - h_quadpt ( :, 35)
!     1 - h_quadpt ( :, 22) ) - 1.5 * h ( :, 23)

!      h ( :, 21) = h_quadpt ( :, 42) - h_quadpt ( :, 01)
!     1 - h ( :, 22) - h ( :, 23)
!     1 + ( h ( :, 30) + h ( :, 31)
!     1 + h ( :, 32) ) / 12.0

! ----------

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
! Default is to build source terms at t>0 . In some ADER-CG updates
! (on Cartesian meshes), we may not require any source terms.
#if ( BUILD_SOURCE_TERMS == 1)
! ----------------------------------------------------------------------

! STEP 9.a) Build source terms "s_quadpt" at the t>0 quadrature points.
! Remember that the primitive variables are still available from the previous
! loop at each of the t=0 quadrature points for constructing the source terms.

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      s_quadpt ( :, iquadpt) = 0.0

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

#if ( 0 == 1)

! For single-fluid adiabatic MHD turbulence, we may want to provide a gradual
! (explicit) heating/cooling term which ensures that pressure always
! tends to a mean pressure "prs_bar" in a time "tau_bar" (comparable
! to a turnover time). Such a source term can be implemented explicitly
! and its implementation is given below.
! Remember to declare "prs_bar, tau_bar".

      prs_bar = 1.0
      tau_bar = 1.0

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

! ----------

      s_quadpt ( 5, iquadpt) = - ( prs ( iquadpt) - prs_bar)
     1     / ( tau_bar * ( gamma_eos_quadpt ( iquadpt) - 1.0) )

! ----------

      END DO
! Finished "DO iquadpt ="

      END IF

#endif

! ------------------------------

! Cylindrical Geometry

      IF ( igeom .EQ. 2) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

! ----------

      x1 = nodal_points ( 1, iquadpt)
      x2 = nodal_points ( 2, iquadpt)
      x3 = nodal_points ( 3, iquadpt)

      obx1 = 1.0 / x1

      tempaa = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      tempab = ( vlx ( iquadpt) * bfx ( iquadpt)
     1         + vly ( iquadpt) * bfy ( iquadpt)
     1         + vlz ( iquadpt) * bfz ( iquadpt)) * ob4pi

      tempac = u_quadpt ( 5, iquadpt) + prs ( iquadpt) + tempaa

! ----------

      s_quadpt ( 1, iquadpt) = - rho ( iquadpt)
     1                           * vlx ( iquadpt) * obx1

      s_quadpt ( 2, iquadpt) =
     1 ( rho ( iquadpt) * vly ( iquadpt)**2
     1 - rho ( iquadpt) * vlx ( iquadpt)**2
     1 + bfx ( iquadpt)**2 * ob4pi
     1 - bfy ( iquadpt)**2 * ob4pi) * obx1

      s_quadpt ( 3, iquadpt) =
     1 - 2.0 * ( rho ( iquadpt) * vlx ( iquadpt)
     1                              * vly ( iquadpt)
     1 - bfx ( iquadpt) * bfy ( iquadpt) * ob4pi) * obx1

      s_quadpt ( 4, iquadpt) =
     1 - ( rho ( iquadpt) * vlx ( iquadpt)
     1                        * vlz ( iquadpt)
     1 - bfx ( iquadpt) * bfz ( iquadpt) * ob4pi) * obx1

      s_quadpt ( 5, iquadpt) = - ( tempac * vlx ( iquadpt)
     1 - bfx ( iquadpt) * tempab) * obx1

      s_quadpt ( 6, iquadpt) = 0.0

      s_quadpt ( 7, iquadpt) = 0.0

      s_quadpt ( 8, iquadpt) = - f_quadpt ( 8, iquadpt) * obx1

! ----------

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Spherical Geometry

      IF ( igeom .EQ. 3) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

! ----------

      x1 = nodal_points ( 1, iquadpt)
      x2 = nodal_points ( 2, iquadpt)
      x3 = nodal_points ( 3, iquadpt)

      obx1 = 1.0 / x1

      cos_x2 = cos ( x2)
      sin_x2 = sin ( x2)
      cot_x2 = cos_x2 / sin_x2

      tempaa = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      tempab = ( vlx ( iquadpt) * bfx ( iquadpt)
     1         + vly ( iquadpt) * bfy ( iquadpt)
     1         + vlz ( iquadpt) * bfz ( iquadpt)) * ob4pi

      tempac = u_quadpt ( 5, iquadpt) + prs ( iquadpt) + tempaa

! ----------

      s_quadpt ( 1, iquadpt) = - rho ( iquadpt)
     1 * ( 2.0 * vlx ( iquadpt)
     1 + cot_x2 * vly ( iquadpt) ) * obx1

      s_quadpt ( 2, iquadpt) = rho ( iquadpt)
     1 * ( - 2.0 * vlx ( iquadpt)**2 + vly ( iquadpt)**2
     1   + vlz ( iquadpt)**2) * obx1
     1 - cot_x2
     1 * ( rho ( iquadpt) * vlx ( iquadpt) * vly ( iquadpt)
     1   - bfx ( iquadpt) * bfy ( iquadpt) * ob4pi) * obx1
     1 - ( - 2.0 * bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1   + bfz ( iquadpt)**2) * ob4pi * obx1

      s_quadpt ( 3, iquadpt) = - rho ( iquadpt)
     1 * ( 3.0 * vlx ( iquadpt) * vly ( iquadpt)
     1 + cot_x2
     1 * ( vly ( iquadpt)**2 - vlz ( iquadpt)**2)) * obx1
     1 + ( 3.0 * bfx ( iquadpt) * bfy ( iquadpt)
     1   + cot_x2
     1   * ( bfy ( iquadpt)**2 - bfz ( iquadpt)**2))
     1 * ob4pi * obx1

      s_quadpt ( 4, iquadpt) =
     1 - rho ( iquadpt) * vlz ( iquadpt)
     1 * ( 3.0 * vlx ( iquadpt) + 2.0 * cot_x2
     1   * vly ( iquadpt) ) * obx1
     1 + bfz ( iquadpt)
     1 * ( 3.0 * bfx ( iquadpt) + 2.0 * cot_x2
     1   * bfy ( iquadpt) ) * ob4pi * obx1

      s_quadpt ( 5, iquadpt) =
     1 - 2.0 * ( tempac * vlx ( iquadpt)
     1         - bfx ( iquadpt) * tempab) * obx1
     1 - cot_x2 * ( tempac * vly ( iquadpt)
     1 - bfy ( iquadpt) * tempab) * obx1

      s_quadpt ( 6, iquadpt) = - cot_x2
     1                           * g_quadpt ( 6, iquadpt) * obx1

      s_quadpt ( 7, iquadpt) = g_quadpt ( 6, iquadpt) * obx1

      s_quadpt ( 8, iquadpt) = - f_quadpt ( 8, iquadpt) * obx1

! ----------

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Multiply "s_quadpt" by "dtmove".

      s_quadpt ( :, NUM_ZONECEN_QUADPTS + 1: NUM_SPACETIME_QUADPTS)=
     1 dtmove
     1*s_quadpt ( :, NUM_ZONECEN_QUADPTS + 1: NUM_SPACETIME_QUADPTS)

! ------------------------------

! STEP 9.b) Now project the sources "s_quadpt" from the t>0 space-time nodal
! points back to the modal space time variables for those variables,
! i.e. to "s".

! Note that some of the "s_quadpt" values with t>0 will be overwritten.
! But that is ok since these are also constructed anew in each ADER iteration.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal

! ----------

      s ( :, 05) = s_quadpt ( :, 07) - s ( :, 01)

! ----------

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal

! ----------

      s ( :, 13) = 2.0 * ( s_quadpt ( :, 14)
     1 - s_quadpt ( :, 15) - s ( :, 02) )

      s ( :, 14) = 2.0 * ( s_quadpt ( :, 16)
     1 - s_quadpt ( :, 17) - s ( :, 03) )

      s ( :, 15) = 2.0 * ( s_quadpt ( :, 18)
     1 - s_quadpt ( :, 19) - s ( :, 04) )

      s ( :, 11) = 2.0 * ( s_quadpt ( :, 14)
     1 + s_quadpt ( :, 15) + s_quadpt ( :, 16)
     1 + s_quadpt ( :, 17) + s_quadpt ( :, 18)
     1 + s_quadpt ( :, 19) - 6.0 * s ( :, 01) ) / 3.0
     1 - s_quadpt ( :, 20) + s_quadpt ( :, 01)

      s ( :, 12) = s_quadpt ( :, 20) - s_quadpt ( :, 01)
     1 - s ( :, 11)

! ----------

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal

! ----------

      ts ( :, 07) = s_quadpt ( :, 23) + s_quadpt ( :, 24)
     1 - s_quadpt ( :, 25) - s_quadpt ( :, 26)
     1 - s_quadpt ( :, 27) - s_quadpt ( :, 28)
     1 + s_quadpt ( :, 29) + s_quadpt ( :, 30)

      ts ( :, 08) = s_quadpt ( :, 32) + s_quadpt ( :, 34)
     1 - s_quadpt ( :, 24) - s_quadpt ( :, 28)
     1 - s_quadpt ( :, 31) - s_quadpt ( :, 33)
     1 + s_quadpt ( :, 23) + s_quadpt ( :, 27)

      ts ( :, 09) = s_quadpt ( :, 25) + s_quadpt ( :, 29)
     1 - s_quadpt ( :, 31) - s_quadpt ( :, 32)
     1 - s_quadpt ( :, 26) - s_quadpt ( :, 30)
     1 + s_quadpt ( :, 33) + s_quadpt ( :, 34)

      s ( :, 16) = ts ( :, 07) + 0.5 * ( ts ( :, 01)
     1 - ts ( :, 02) ) - 2.0 * ( s_quadpt ( :, 02)
     1 - s_quadpt ( :, 05) )

      s ( :, 18) = - 2.0 * ts ( :, 07) + s ( :, 16)

      s ( :, 14) = ts ( :, 08) + 0.5 * ( ts ( :, 03)
     1 - ts ( :, 04) ) - 2.0 * ( s_quadpt ( :, 06)
     1 - s_quadpt ( :, 09) )

      s ( :, 19) = - 2.0 * ts ( :, 08) + s ( :, 14)

      s ( :, 15) = ts ( :, 09) + 0.5 * ( ts ( :, 05)
     1 - ts ( :, 06) ) - 2.0 * ( s_quadpt ( :, 10)
     1 - s_quadpt ( :, 13) )

      s ( :, 17) = - 2.0 * ts ( :, 09) + s ( :, 15)


      ts ( :, 10) = s_quadpt ( :, 23) + s_quadpt ( :, 24)
     1 + s_quadpt ( :, 25) + s_quadpt ( :, 26)
     1 - s_quadpt ( :, 27) - s_quadpt ( :, 28)
     1 - s_quadpt ( :, 29) - s_quadpt ( :, 30)

      ts ( :, 11) = s_quadpt ( :, 23) + s_quadpt ( :, 27)
     1 + s_quadpt ( :, 31) + s_quadpt ( :, 33)
     1 - s_quadpt ( :, 24) - s_quadpt ( :, 28)
     1 - s_quadpt ( :, 32) - s_quadpt ( :, 34)

      ts ( :, 12) = s_quadpt ( :, 25) + s_quadpt ( :, 29)
     1 + s_quadpt ( :, 31) + s_quadpt ( :, 32)
     1 - s_quadpt ( :, 26) - s_quadpt ( :, 30)
     1 - s_quadpt ( :, 33) - s_quadpt ( :, 34)

      s ( :, 27) = 1.125 * ( ts ( :, 01) - ts ( :, 02) )
     1 + 4.5 * ( s_quadpt ( :, 36) - s_quadpt ( :, 37) )
     1 - 2.25 * ts ( :, 10)

      s ( :, 24) = 3.0 * ( s_quadpt ( :, 36)
     1 - s_quadpt ( :, 37) ) - 0.75 * ts ( :, 10)
     1 + 0.375 * ( s ( :, 16) + s ( :, 18) )
     1 - s ( :, 27)

      s ( :, 28) = 1.125 * ( ts ( :, 03) - ts ( :, 04) )
     1 + 4.5 * ( s_quadpt ( :, 38) - s_quadpt ( :, 39) )
     1 - 2.25 * ts ( :, 11)

      s ( :, 25) = 3.0 * ( s_quadpt ( :, 38)
     1 - s_quadpt ( :, 39) ) - 0.75 * ts ( :, 11)
     1 + 0.375 * ( s ( :, 14) + s ( :, 19) )
     1 - s ( :, 28)

      s ( :, 29) = 1.125 * ( ts ( :, 05) - ts ( :, 06) )
     1 + 4.5 * ( s_quadpt ( :, 40) - s_quadpt ( :, 41) )
     1 - 2.25 * ts ( :, 12)

      s ( :, 26) = 3.0 * ( s_quadpt ( :, 40)
     1 - s_quadpt ( :, 41) ) - 0.75 * ts ( :, 12)
     1 + 0.375 * ( s ( :, 15) + s ( :, 17) )
     1 - s ( :, 29)


      s ( :, 30) = 3.0 * ( s_quadpt ( :, 36)
     1 - 2.0 * s_quadpt ( :, 35) + s_quadpt ( :, 37) )
     1 - 1.5 * s ( :, 05)

      s ( :, 31) = 3.0 * ( s_quadpt ( :, 38)
     1 - 2.0 * s_quadpt ( :, 35) + s_quadpt ( :, 39) )
     1 - 1.5 * s ( :, 06)

      s ( :, 32) = 3.0 * ( s_quadpt ( :, 40)
     1 - 2.0 * s_quadpt ( :, 35) + s_quadpt ( :, 41) )
     1 - 1.5 * s ( :, 07)

      s ( :, 33) = 3.0 * ( s_quadpt ( :, 23)
     1 - s_quadpt ( :, 27) - s_quadpt ( :, 24)
     1 + s_quadpt ( :, 28) ) - 3.0 * s ( :, 08)

      s ( :, 34) = 3.0 * ( s_quadpt ( :, 31)
     1 - s_quadpt ( :, 33) - s_quadpt ( :, 32)
     1 + s_quadpt ( :, 34) ) - 3.0 * s ( :, 09)

      s ( :, 35) = 3.0 * ( s_quadpt ( :, 25)
     1 - s_quadpt ( :, 29) - s_quadpt ( :, 26)
     1 + s_quadpt ( :, 30) ) - 3.0 * s ( :, 10)


      s ( :, 23) = 4.5 * ( s_quadpt ( :, 42)
     1 - s_quadpt ( :, 01) ) - 13.5 * ( s_quadpt ( :, 35)
     1 - s_quadpt ( :, 22) )

      s ( :, 22) = 2.25 * ( s_quadpt ( :, 42)
     1 + s_quadpt ( :, 01) - s_quadpt ( :, 35)
     1 - s_quadpt ( :, 22) ) - 1.5 * s ( :, 23)

      s ( :, 21) = s_quadpt ( :, 42) - s_quadpt ( :, 01)
     1 - s ( :, 22) - s ( :, 23)
     1 + ( s ( :, 30) + s ( :, 31)
     1 + s ( :, 32) ) / 12.0

! ----------

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
#endif
! End of "BUILD_SOURCE_TERMS == 1"
! ----------------------------------------------------------------------

! STEP 10) By this point, the modal flux coefficients ( f, g, h) and the modal
! source coefficients "s" that enable us to define the fluxes and source terms
! any where in space and time within the zone have already been loaded.
! Now we can use the analytical inversion of the matrix to obtain
! the next iterate for the conserved variables "u" . The coefficients "u"
! enable us to construct the conserved variables at any space-time point in
! the zone.

! Since "u" is made by the end of this step, we also build and return "s_int",
! the space-time integral of the source term.

! ----------------------------------------------------------------------

#if ( NUMORDER == 2)

! ----------

      u ( :, 5) = - f ( :, 2) - g ( :, 3) - h ( :, 4)
     1 + s ( :, 1) + 2.0 * s ( :, 5) / 3.0

! ----------

#endif
! End of "NUMORDER == 2".

! ----------------------------------------------------------------------

#if ( NUMORDER == 3)

! ----------

      u ( :, 11) = - f ( :, 2) - g ( :, 3) - h ( :, 4)
     1 + s ( :, 1) - 0.3 * s ( :, 12)

      u ( :, 12) = - 0.5 * ( f ( :, 13) + g ( :, 14)
     1 + h ( :, 15) - s ( :, 11) ) + 0.6 * s ( :, 12)

      u ( :, 13) = - 2.0 * f ( :, 5) - g ( :, 8)
     1 - h ( :, 10) + s ( :, 2) + 2.0 * s ( :, 13) / 3.0

      u ( :, 14) = - f ( :, 8) - 2.0 * g ( :, 6)
     1 - h ( :, 9) + s ( :, 3) + 2.0 * s ( :, 14) / 3.0

      u ( :, 15) = - f ( :, 10) - g ( :, 9)
     1 - 2.0 * h ( :, 7) + s ( :, 4)
     1 + 2.0 * s ( :, 15) / 3.0

! ----------

#endif
! End of "NUMORDER == 3".

! ----------------------------------------------------------------------

#if ( NUMORDER == 4)

! ----------

      u ( :, 21) =
     1 - 0.1 * ( f ( :, 11) + g ( :, 12) + h ( :, 13) )
     1 - f ( :, 2) - g ( :, 3) - h ( :, 4)
     1 + s ( :, 1) + 8.0 * s ( :, 23) / 70.0

      u ( :, 22) =
     1 - 0.5 * ( f ( :, 24) + g ( :, 25) + h ( :, 26)
     1 - s ( :, 21) ) - 3.0 * s ( :, 23) / 7.0

      u ( :, 23) =
     1 - ( f ( :, 27) + g ( :, 28) + h ( :, 29)
     1 - s ( :, 22) ) / 3.0 + 4.0 * s ( :, 23) / 7.0

      u ( :, 24) =
     1 - 2.0 * f ( :, 5) - g ( :, 8) - h ( :, 10)
     1 + s ( :, 2) - 0.3 * s ( :, 27)

      u ( :, 25) =
     1 - f ( :, 8) - 2.0 * g ( :, 6) - h ( :, 9)
     1 + s ( :, 3) - 0.3 * s ( :, 28)

      u ( :, 26) =
     1 - f ( :, 10) - g ( :, 9) - 2.0 * h ( :, 7)
     1 + s ( :, 4) - 0.3 * s ( :, 29)

      u ( :, 27) =
     1 - 0.5 * ( g ( :, 33) + h ( :, 35) - s ( :, 24) )
     1 - f ( :, 30) + 0.6 * s ( :, 27)

      u ( :, 28) =
     1 - 0.5 * ( f ( :, 33) + h ( :, 34) - s ( :, 25) )
     1 - g ( :, 31) + 0.6 * s ( :, 28)

      u ( :, 29) =
     1 - 0.5 * ( f ( :, 35) + g ( :, 34) - s ( :, 26) )
     1 - h ( :, 32) + 0.6 * s ( :, 29)

      u ( :, 30) =
     1 - 3.0 * f ( :, 11) - g ( :, 14) - h ( :, 15)
     1 + s ( :, 5) + 2.0 * s ( :, 30) / 3.0

      u ( :, 31) =
     1 - f ( :, 16) - 3.0 * g ( :, 12) - h ( :, 17)
     1 + s ( :, 6) + 2.0 * s ( :, 31) / 3.0

      u ( :, 32) =
     1 - f ( :, 18) - g ( :, 19) - 3.0 * h ( :, 13)
     1 + s ( :, 7) + 2.0 * s ( :, 32) / 3.0

      u ( :, 33) =
     1 - 2.0 * f ( :, 14) - 2.0 * g ( :, 16) - h ( :, 20)
     1 + s ( :, 8) + 2.0 * s ( :, 33) / 3.0

      u ( :, 34) =
     1 - f ( :, 20) - 2.0 * g ( :, 17) - 2.0 * h ( :, 19)
     1 + s ( :, 9) + 2.0 * s ( :, 34) / 3.0

      u ( :, 35) =
     1 - 2.0 * f ( :, 15) - g ( :, 20) - 2.0 * h ( :, 18)
     1 + s ( :, 10) + 2.0 * s ( :, 35) / 3.0

! ----------

#endif
! End of "NUMORDER == 4".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO
! End of "DO iter_ader =" loop.

! End the outer ADER iteration here.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Beginning of building all the nodal conserved and primitive variables at the
! t > 0 nodes. This is needed either if we want to rebuild the modal
! fluxes "f,g,h" or if we want to rebuild the modal source term "s".

#if ( BUILD_SOURCE_TERMS == 1)

! ----------------------------------------------------------------------

! STEP 6_Polish-off) Make the conserved variables at the t>0 space-time nodes.
! We put in the flexibility of having different "gamma_eos_spacetime" at each
! space-time node, but that is not used.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : modal -> nodal

! ----------

      u_quadpt ( :, 07) = u ( :, 01) + u ( :, 05)

! ----------

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : modal -> nodal

! ----------

      u_quadpt ( :, 14) = 0.25 * ( u_quadpt ( :, 02)
     1 + u_quadpt ( :, 03) + u_quadpt ( :, 04)
     1 + u_quadpt ( :, 05) + u ( :, 12) + u ( :, 13) )
     1 + 0.5 * u ( :, 11) - 0.125 * ( u ( :, 06)
     1 + u ( :, 07) )

      u_quadpt ( :, 15) = u_quadpt ( :, 14)
     1 - u ( :, 02) - 0.5 * u ( :, 13)

      u_quadpt ( :, 16) = 0.25 * ( u_quadpt ( :, 02)
     1 + u_quadpt ( :, 06) + u_quadpt ( :, 10)
     1 + u_quadpt ( :, 12) + u ( :, 12) + u ( :, 14) )
     1 + 0.5 * u ( :, 11) - 0.125 * ( u ( :, 05)
     1 + u ( :, 07) )

      u_quadpt ( :, 17) = u_quadpt ( :, 16)
     1 - u ( :, 03) - 0.5 * u ( :, 14)

      u_quadpt ( :, 18) = 0.25 * ( u_quadpt ( :, 04)
     1 + u_quadpt ( :, 08) + u_quadpt ( :, 10)
     1 + u_quadpt ( :, 11) + u ( :, 12) + u ( :, 15) )
     1 + 0.5 * u ( :, 11) - 0.125 * ( u ( :, 05)
     1 + u ( :, 06) )

      u_quadpt ( :, 19) = u_quadpt ( :, 18) - u ( :, 04)
     1 - 0.5 * u ( :, 15)

      u_quadpt ( :, 20) = u_quadpt ( :, 01) + u ( :, 11)
     1 + u ( :, 12)

! ----------

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : modal -> nodal

! ----------

      u_quadpt ( :, 22) = u_quadpt ( :, 01)
     1 + ( 36.0 * u ( :, 21) + 12.0 * u ( :, 22)
     1 + 4.0 * u ( :, 23) - 3.0 * u ( :, 30)
     1 - 3.0 * u ( :, 31) - 3.0 * u ( :, 32) ) / 108.0

      u_quadpt ( :, 23) = u_quadpt ( :, 22)
     1 + u_quadpt ( :, 02) + u_quadpt ( :, 06)
     1 - 2.0 * u_quadpt ( :, 01) + ( 18.0 * u ( :, 08)
     1 + 9.0 * u ( :, 14) + 9.0 * u ( :, 16)
     1 + 12.0 * u ( :, 24) + 12.0 * u ( :, 25)
     1 + 4.0 * u ( :, 27) + 4.0 * u ( :, 28)
     1 + 6.0 * u ( :, 30) + 6.0 * u ( :, 31)
     1 + 6.0 * u ( :, 33) ) / 72.0

      u_quadpt ( :, 24) = u_quadpt ( :, 23)
     1 - u_quadpt ( :, 06) + u_quadpt ( :, 09)
     1 - 0.5 * u ( :, 08) - 0.25 * u ( :, 14)
     1 - ( 6.0 * u ( :, 25) + 2.0 * u ( :, 28)
     1 + 3.0 * u ( :, 33) ) / 18.0

      u_quadpt ( :, 25) = 0.5 * ( u_quadpt ( :, 23)
     1 + u_quadpt ( :, 24) + u_quadpt ( :, 14)
     1 + u_quadpt ( :, 15) ) - u_quadpt ( :, 02)
     1 - 0.5 * u ( :, 06) - 0.25 * u ( :, 16)
     1 - 0.125 * u ( :, 17) + ( 6.0 * u ( :, 26)
     1 + 2.0 * u ( :, 29) - 3.0 * u ( :, 31)
     1 + 3.0 * u ( :, 32) + 3.0 * u ( :, 35) ) / 36.0

      u_quadpt ( :, 26) = u_quadpt ( :, 23)
     1 + u_quadpt ( :, 24) - u_quadpt ( :, 25)
     1 - 0.5 * ( u ( :, 06) - u ( :, 07) )
     1 - 0.25 * ( u ( :, 16) - u ( :, 18) )
     1 - ( u ( :, 31) - u ( :, 32) ) / 6.0

      u_quadpt ( :, 27) = u_quadpt ( :, 23)
     1 + u_quadpt ( :, 05) - u_quadpt ( :, 02)
     1 - 0.5 * u ( :, 08) - 0.25 * u ( :, 16)
     1 - ( 6.0 * u ( :, 24) + 2.0 * u ( :, 27)
     1 + 3.0 * u ( :, 33) ) / 18.0

      u_quadpt ( :, 28) = u_quadpt ( :, 27)
     1 - u_quadpt ( :, 06) + u_quadpt ( :, 09)
     1 + 0.5 * u ( :, 08) - 0.25 * u ( :, 14)
     1 - ( 6.0 * u ( :, 25) + 2.0 * u ( :, 28)
     1 - 3.0 * u ( :, 33) ) / 18.0

      u_quadpt ( :, 29) = 0.5 * ( u_quadpt ( :, 27)
     1 + u_quadpt ( :, 28) + u_quadpt ( :, 18)
     1 + u_quadpt ( :, 19) ) - u_quadpt ( :, 05)
     1 - 0.5 * u ( :, 06) + 0.25 * u ( :, 16)
     1 - 0.125 * u ( :, 17) + ( 6.0 * u ( :, 26)
     1 + 2.0 * u ( :, 29) - 3.0 * u ( :, 31)
     1 + 3.0 * u ( :, 32)  - 3.0 * u ( :, 35) ) / 36.0

      u_quadpt ( :, 30) = u_quadpt ( :, 27)
     1 + u_quadpt ( :, 28) - u_quadpt ( :, 29)
     1 - 0.5 * ( u ( :, 06) - u ( :, 07) )
     1 + 0.25 * ( u ( :, 16) - u ( :, 18) )
     1 - ( u ( :, 31) - u ( :, 32) ) / 6.0

      u_quadpt ( :, 31) = 0.5 * ( u_quadpt ( :, 25)
     1 + u_quadpt ( :, 29) + u_quadpt ( :, 14)
     1 + u_quadpt ( :, 18) ) - u_quadpt ( :, 10)
     1 - 0.5 * u ( :, 05) - 0.25 * u ( :, 15)
     1 - 0.125 * u ( :, 14) + ( 6.0 * u ( :, 25)
     1 + 2.0 * u ( :, 28) - 3.0 * u ( :, 30)
     1 + 3.0 * u ( :, 31) + 3.0 * u ( :, 34) ) / 36.0

      u_quadpt ( :, 32) = u_quadpt ( :, 31)
     1 - u_quadpt ( :, 06) + u_quadpt ( :, 09)
     1 - 0.5 * u ( :, 09) - 0.25 * u ( :, 19)
     1 - ( 6.0 * u ( :, 25) + 2.0 * u ( :, 28)
     1 + 3.0 * u ( :, 34) ) / 18.0

      u_quadpt ( :, 33) = u_quadpt ( :, 23)
     1 + u_quadpt ( :, 27) - u_quadpt ( :, 31)
     1 - 0.5 * ( u ( :, 05) - u ( :, 07) )
     1 - 0.25 * ( u ( :, 14) - u ( :, 19) )
     1 - ( u ( :, 30) - u ( :, 32) ) / 6.0

      u_quadpt ( :, 34) = u_quadpt ( :, 24)
     1 + u_quadpt ( :, 28) - u_quadpt ( :, 32)
     1 - 0.5 * ( u ( :, 05) - u ( :, 07) )
     1 + 0.25 * ( u ( :, 14) - u ( :, 19) )
     1 - ( u ( :, 30) - u ( :, 32) ) / 6.0


      tu ( :, 01) = ( 36.0 * u ( :, 21) + 24.0 * u ( :, 22)
     1 + 16.0 * u ( :, 23) ) / 54.0

      tu ( :, 02) = ( 18.0 * u ( :, 24)
     1 + 12.0 * u ( :, 27) ) / 54.0

      tu ( :, 03) = ( 6.0 * u ( :, 30) - 3.0 * u ( :, 31)
     1 - 3.0 * u ( :, 32) ) / 54.0

      tu ( :, 04) = ( 18.0 * u ( :, 25)
     1 + 12.0 * u ( :, 28) ) / 54.0

      tu ( :, 05) = ( - 3.0 * u ( :, 30) + 6.0 * u ( :, 31)
     1 - 3.0 * u ( :, 32) ) / 54.0

      tu ( :, 06) = ( 18.0 * u ( :, 26)
     1 + 12.0 * u ( :, 29) ) / 54.0

      tu ( :, 07) = ( - 3.0 * u ( :, 30) - 3.0 * u ( :, 31)
     1 + 6.0 * u ( :, 32) ) / 54.0

      u_quadpt ( :, 35) = ( 2.0 * u_quadpt ( :, 22)
     1 - u_quadpt ( :, 01) ) + 2.0 * ( u ( :, 22)
     1 + u ( :, 23) ) / 9.0

      u_quadpt ( :, 36) = u_quadpt ( :, 02) + tu ( :, 01)
     1 + tu ( :, 02) + tu ( :, 03)

      u_quadpt ( :, 37) = u_quadpt ( :, 05) + tu ( :, 01)
     1 - tu ( :, 02) + tu ( :, 03)

      u_quadpt ( :, 38) = u_quadpt ( :, 06) + tu ( :, 01)
     1 + tu ( :, 04) + tu ( :, 05)

      u_quadpt ( :, 39) = u_quadpt ( :, 09) + tu ( :, 01)
     1 - tu ( :, 04) + tu ( :, 05)

      u_quadpt ( :, 40) = u_quadpt ( :, 10) + tu ( :, 01)
     1 + tu ( :, 06) + tu ( :, 07)

      u_quadpt ( :, 41) = u_quadpt ( :, 13) + tu ( :, 01)
     1 - tu ( :, 06) + tu ( :, 07)


      u_quadpt ( :, 42) = u_quadpt ( :, 01) + u ( :, 21)
     1 + u ( :, 22) + u ( :, 23) - ( u ( :, 30)
     1 + u ( :, 31) + u ( :, 32) ) / 12.0

! ----------

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 7_Polish-off) Get the primitive variables from the conserved variables
! at the spatial nodes at t>0. The density and pressure variables are
! given positive floor values to ensure that they can be used for physically
! meaningful source terms.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

! ----------

      u_quadpt ( 1, iquadpt) =
     1                   AMAX1 ( u_quadpt ( 1, iquadpt), rhofloor)

      rho ( iquadpt) = u_quadpt ( 1, iquadpt)

      tempaa = 1.0 / rho ( iquadpt)

      vlx ( iquadpt) = u_quadpt ( 2, iquadpt) * tempaa
      vly ( iquadpt) = u_quadpt ( 3, iquadpt) * tempaa
      vlz ( iquadpt) = u_quadpt ( 4, iquadpt) * tempaa

      bfx ( iquadpt) = u_quadpt ( 6, iquadpt)
      bfy ( iquadpt) = u_quadpt ( 7, iquadpt)
      bfz ( iquadpt) = u_quadpt ( 8, iquadpt)

      tempab = 0.5 * rho ( iquadpt) * ( vlx ( iquadpt)**2
     1         + vly ( iquadpt)**2 + vlz ( iquadpt)**2)

      tempac = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      eng_t ( iquadpt) = u_quadpt ( 5, iquadpt)
     1                     - tempab - tempac

      prs ( iquadpt) = eng_t ( iquadpt)
     1                   * ( gamma_eos_quadpt ( iquadpt) - 1.0)

#if ( ENFORCE_POSITIVITY == 1)
      IF ( prs ( iquadpt) .LE. prsfloor) THEN

      prs ( iquadpt) = entropy ( iquadpt)
     1* rho ( iquadpt)**( gamma_eos_quadpt ( iquadpt) - 1.0)
      prs ( iquadpt) = AMAX1 ( prs ( iquadpt), prsfloor)

      eng_t ( iquadpt) = prs ( iquadpt)
     1                 / ( gamma_eos_quadpt ( iquadpt) - 1.0)

      u_quadpt ( 5, iquadpt) = eng_t ( iquadpt)
     1                           + tempab + tempac

      END IF
#else
      IF ( prs ( iquadpt) .LE. prsfloor) THEN

      prs ( iquadpt) = prsfloor

      eng_t ( iquadpt) = prs ( iquadpt)
     1                 / ( gamma_eos_quadpt ( iquadpt) - 1.0)

      u_quadpt ( 5, iquadpt) = eng_t ( iquadpt)
     1                           + tempab + tempac

      END IF
#endif

#if ( NFLUID >= 1)
      specie ( 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( 1: NFLUID, iquadpt) * rho ( iquadpt)
#endif

! ----------

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

#endif
! End of "BUILD_SOURCE_TERMS == 1"

! End of building all the nodal conserved and primitive variables at the
! t > 0 nodes. This is needed either if we want to rebuild the modal
! fluxes "f,g,h" or if we want to rebuild the modal source term "s".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Default is to build source terms at t>0 . In some ADER-CG updates
! (on Cartesian meshes), we may not require any source terms.

#if ( BUILD_SOURCE_TERMS == 1)
! ----------------------------------------------------------------------

! STEP 9.a_Polish-off) Build source terms "s_quadpt" at the t>0 quadrature
! points.
! Remember that the primitive variables are still available from the previous
! loop at each of the t=0 quadrature points for constructing the source terms.

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      s_quadpt ( :, iquadpt) = 0.0

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

#if ( 0 == 1)

! For single-fluid adiabatic MHD turbulence, we may want to provide a gradual
! (explicit) heating/cooling term which ensures that pressure always
! tends to a mean pressure "prs_bar" in a time "tau_bar" (comparable
! to a turnover time). Such a source term can be implemented explicitly
! and its implementation is given below.
! Remember to declare "prs_bar, tau_bar".

      prs_bar = 1.0
      tau_bar = 1.0

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

! ----------

      s_quadpt ( 5, iquadpt) = - ( prs ( iquadpt) - prs_bar)
     1     / ( tau_bar * ( gamma_eos_quadpt ( iquadpt) - 1.0) )

! ----------

      END DO
! Finished "DO iquadpt ="

      END IF

#endif

! ------------------------------

! Cylindrical Geometry

      IF ( igeom .EQ. 2) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

! ----------

      x1 = nodal_points ( 1, iquadpt)
      x2 = nodal_points ( 2, iquadpt)
      x3 = nodal_points ( 3, iquadpt)

      obx1 = 1.0 / x1

      tempaa = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      tempab = ( vlx ( iquadpt) * bfx ( iquadpt)
     1         + vly ( iquadpt) * bfy ( iquadpt)
     1         + vlz ( iquadpt) * bfz ( iquadpt)) * ob4pi

      tempac = u_quadpt ( 5, iquadpt) + prs ( iquadpt) + tempaa

! ----------

      s_quadpt ( 1, iquadpt) = - rho ( iquadpt)
     1                           * vlx ( iquadpt) * obx1

      s_quadpt ( 2, iquadpt) =
     1 ( rho ( iquadpt) * vly ( iquadpt)**2
     1 - rho ( iquadpt) * vlx ( iquadpt)**2
     1 + bfx ( iquadpt)**2 * ob4pi
     1 - bfy ( iquadpt)**2 * ob4pi) * obx1

      s_quadpt ( 3, iquadpt) =
     1 - 2.0 * ( rho ( iquadpt) * vlx ( iquadpt)
     1                              * vly ( iquadpt)
     1 - bfx ( iquadpt) * bfy ( iquadpt) * ob4pi) * obx1

      s_quadpt ( 4, iquadpt) =
     1 - ( rho ( iquadpt) * vlx ( iquadpt)
     1                        * vlz ( iquadpt)
     1 - bfx ( iquadpt) * bfz ( iquadpt) * ob4pi) * obx1

      s_quadpt ( 5, iquadpt) = - ( tempac * vlx ( iquadpt)
     1 - bfx ( iquadpt) * tempab) * obx1

      s_quadpt ( 6, iquadpt) = 0.0

      s_quadpt ( 7, iquadpt) = 0.0

      s_quadpt ( 8, iquadpt) = - f_quadpt ( 8, iquadpt) * obx1

! ----------

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Spherical Geometry

      IF ( igeom .EQ. 3) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

! ----------

      x1 = nodal_points ( 1, iquadpt)
      x2 = nodal_points ( 2, iquadpt)
      x3 = nodal_points ( 3, iquadpt)

      obx1 = 1.0 / x1

      cos_x2 = cos ( x2)
      sin_x2 = sin ( x2)
      cot_x2 = cos_x2 / sin_x2

      tempaa = ( bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1         + bfz ( iquadpt)**2) * ob8pi

      tempab = ( vlx ( iquadpt) * bfx ( iquadpt)
     1         + vly ( iquadpt) * bfy ( iquadpt)
     1         + vlz ( iquadpt) * bfz ( iquadpt)) * ob4pi

      tempac = u_quadpt ( 5, iquadpt) + prs ( iquadpt) + tempaa

! ----------

      s_quadpt ( 1, iquadpt) = - rho ( iquadpt)
     1 * ( 2.0 * vlx ( iquadpt)
     1 + cot_x2 * vly ( iquadpt) ) * obx1

      s_quadpt ( 2, iquadpt) = rho ( iquadpt)
     1 * ( - 2.0 * vlx ( iquadpt)**2 + vly ( iquadpt)**2
     1   + vlz ( iquadpt)**2) * obx1
     1 - cot_x2
     1 * ( rho ( iquadpt) * vlx ( iquadpt) * vly ( iquadpt)
     1   - bfx ( iquadpt) * bfy ( iquadpt) * ob4pi) * obx1
     1 - ( - 2.0 * bfx ( iquadpt)**2 + bfy ( iquadpt)**2
     1   + bfz ( iquadpt)**2) * ob4pi * obx1

      s_quadpt ( 3, iquadpt) = - rho ( iquadpt)
     1 * ( 3.0 * vlx ( iquadpt) * vly ( iquadpt)
     1 + cot_x2
     1 * ( vly ( iquadpt)**2 - vlz ( iquadpt)**2)) * obx1
     1 + ( 3.0 * bfx ( iquadpt) * bfy ( iquadpt)
     1   + cot_x2
     1   * ( bfy ( iquadpt)**2 - bfz ( iquadpt)**2))
     1 * ob4pi * obx1

      s_quadpt ( 4, iquadpt) =
     1 - rho ( iquadpt) * vlz ( iquadpt)
     1 * ( 3.0 * vlx ( iquadpt) + 2.0 * cot_x2
     1   * vly ( iquadpt) ) * obx1
     1 + bfz ( iquadpt)
     1 * ( 3.0 * bfx ( iquadpt) + 2.0 * cot_x2
     1   * bfy ( iquadpt) ) * ob4pi * obx1

      s_quadpt ( 5, iquadpt) =
     1 - 2.0 * ( tempac * vlx ( iquadpt)
     1         - bfx ( iquadpt) * tempab) * obx1
     1 - cot_x2 * ( tempac * vly ( iquadpt)
     1 - bfy ( iquadpt) * tempab) * obx1

      s_quadpt ( 6, iquadpt) = - cot_x2
     1                           * g_quadpt ( 6, iquadpt) * obx1

      s_quadpt ( 7, iquadpt) = g_quadpt ( 6, iquadpt) * obx1

      s_quadpt ( 8, iquadpt) = - f_quadpt ( 8, iquadpt) * obx1

! ----------

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Multiply "s_quadpt" by "dtmove".

      s_quadpt ( :, NUM_ZONECEN_QUADPTS + 1: NUM_SPACETIME_QUADPTS)=
     1 dtmove
     1*s_quadpt ( :, NUM_ZONECEN_QUADPTS + 1: NUM_SPACETIME_QUADPTS)

! ------------------------------

! STEP 9.b_Polish-off) Now project the sources "s_quadpt" from the t>0
! space-time nodal points back to the modal space time variables for
! those variables, i.e. to "s".

! Note that some of the "s_quadpt" values with t>0 will be overwritten.
! But that is ok since these are also constructed anew in each ADER iteration.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal

! ----------

      s ( :, 05) = s_quadpt ( :, 07) - s ( :, 01)

! ----------

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal

! ----------

      s ( :, 13) = 2.0 * ( s_quadpt ( :, 14)
     1 - s_quadpt ( :, 15) - s ( :, 02) )

      s ( :, 14) = 2.0 * ( s_quadpt ( :, 16)
     1 - s_quadpt ( :, 17) - s ( :, 03) )

      s ( :, 15) = 2.0 * ( s_quadpt ( :, 18)
     1 - s_quadpt ( :, 19) - s ( :, 04) )

      s ( :, 11) = 2.0 * ( s_quadpt ( :, 14)
     1 + s_quadpt ( :, 15) + s_quadpt ( :, 16)
     1 + s_quadpt ( :, 17) + s_quadpt ( :, 18)
     1 + s_quadpt ( :, 19) - 6.0 * s ( :, 01) ) / 3.0
     1 - s_quadpt ( :, 20) + s_quadpt ( :, 01)

      s ( :, 12) = s_quadpt ( :, 20) - s_quadpt ( :, 01)
     1 - s ( :, 11)

! ----------

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal

! ----------

      ts ( :, 07) = s_quadpt ( :, 23) + s_quadpt ( :, 24)
     1 - s_quadpt ( :, 25) - s_quadpt ( :, 26)
     1 - s_quadpt ( :, 27) - s_quadpt ( :, 28)
     1 + s_quadpt ( :, 29) + s_quadpt ( :, 30)

      ts ( :, 08) = s_quadpt ( :, 32) + s_quadpt ( :, 34)
     1 - s_quadpt ( :, 24) - s_quadpt ( :, 28)
     1 - s_quadpt ( :, 31) - s_quadpt ( :, 33)
     1 + s_quadpt ( :, 23) + s_quadpt ( :, 27)

      ts ( :, 09) = s_quadpt ( :, 25) + s_quadpt ( :, 29)
     1 - s_quadpt ( :, 31) - s_quadpt ( :, 32)
     1 - s_quadpt ( :, 26) - s_quadpt ( :, 30)
     1 + s_quadpt ( :, 33) + s_quadpt ( :, 34)

      s ( :, 16) = ts ( :, 07) + 0.5 * ( ts ( :, 01)
     1 - ts ( :, 02) ) - 2.0 * ( s_quadpt ( :, 02)
     1 - s_quadpt ( :, 05) )

      s ( :, 18) = - 2.0 * ts ( :, 07) + s ( :, 16)

      s ( :, 14) = ts ( :, 08) + 0.5 * ( ts ( :, 03)
     1 - ts ( :, 04) ) - 2.0 * ( s_quadpt ( :, 06)
     1 - s_quadpt ( :, 09) )

      s ( :, 19) = - 2.0 * ts ( :, 08) + s ( :, 14)

      s ( :, 15) = ts ( :, 09) + 0.5 * ( ts ( :, 05)
     1 - ts ( :, 06) ) - 2.0 * ( s_quadpt ( :, 10)
     1 - s_quadpt ( :, 13) )

      s ( :, 17) = - 2.0 * ts ( :, 09) + s ( :, 15)


      ts ( :, 10) = s_quadpt ( :, 23) + s_quadpt ( :, 24)
     1 + s_quadpt ( :, 25) + s_quadpt ( :, 26)
     1 - s_quadpt ( :, 27) - s_quadpt ( :, 28)
     1 - s_quadpt ( :, 29) - s_quadpt ( :, 30)

      ts ( :, 11) = s_quadpt ( :, 23) + s_quadpt ( :, 27)
     1 + s_quadpt ( :, 31) + s_quadpt ( :, 33)
     1 - s_quadpt ( :, 24) - s_quadpt ( :, 28)
     1 - s_quadpt ( :, 32) - s_quadpt ( :, 34)

      ts ( :, 12) = s_quadpt ( :, 25) + s_quadpt ( :, 29)
     1 + s_quadpt ( :, 31) + s_quadpt ( :, 32)
     1 - s_quadpt ( :, 26) - s_quadpt ( :, 30)
     1 - s_quadpt ( :, 33) - s_quadpt ( :, 34)

      s ( :, 27) = 1.125 * ( ts ( :, 01) - ts ( :, 02) )
     1 + 4.5 * ( s_quadpt ( :, 36) - s_quadpt ( :, 37) )
     1 - 2.25 * ts ( :, 10)

      s ( :, 24) = 3.0 * ( s_quadpt ( :, 36)
     1 - s_quadpt ( :, 37) ) - 0.75 * ts ( :, 10)
     1 + 0.375 * ( s ( :, 16) + s ( :, 18) )
     1 - s ( :, 27)

      s ( :, 28) = 1.125 * ( ts ( :, 03) - ts ( :, 04) )
     1 + 4.5 * ( s_quadpt ( :, 38) - s_quadpt ( :, 39) )
     1 - 2.25 * ts ( :, 11)

      s ( :, 25) = 3.0 * ( s_quadpt ( :, 38)
     1 - s_quadpt ( :, 39) ) - 0.75 * ts ( :, 11)
     1 + 0.375 * ( s ( :, 14) + s ( :, 19) )
     1 - s ( :, 28)

      s ( :, 29) = 1.125 * ( ts ( :, 05) - ts ( :, 06) )
     1 + 4.5 * ( s_quadpt ( :, 40) - s_quadpt ( :, 41) )
     1 - 2.25 * ts ( :, 12)

      s ( :, 26) = 3.0 * ( s_quadpt ( :, 40)
     1 - s_quadpt ( :, 41) ) - 0.75 * ts ( :, 12)
     1 + 0.375 * ( s ( :, 15) + s ( :, 17) )
     1 - s ( :, 29)


      s ( :, 30) = 3.0 * ( s_quadpt ( :, 36)
     1 - 2.0 * s_quadpt ( :, 35) + s_quadpt ( :, 37) )
     1 - 1.5 * s ( :, 05)

      s ( :, 31) = 3.0 * ( s_quadpt ( :, 38)
     1 - 2.0 * s_quadpt ( :, 35) + s_quadpt ( :, 39) )
     1 - 1.5 * s ( :, 06)

      s ( :, 32) = 3.0 * ( s_quadpt ( :, 40)
     1 - 2.0 * s_quadpt ( :, 35) + s_quadpt ( :, 41) )
     1 - 1.5 * s ( :, 07)

      s ( :, 33) = 3.0 * ( s_quadpt ( :, 23)
     1 - s_quadpt ( :, 27) - s_quadpt ( :, 24)
     1 + s_quadpt ( :, 28) ) - 3.0 * s ( :, 08)

      s ( :, 34) = 3.0 * ( s_quadpt ( :, 31)
     1 - s_quadpt ( :, 33) - s_quadpt ( :, 32)
     1 + s_quadpt ( :, 34) ) - 3.0 * s ( :, 09)

      s ( :, 35) = 3.0 * ( s_quadpt ( :, 25)
     1 - s_quadpt ( :, 29) - s_quadpt ( :, 26)
     1 + s_quadpt ( :, 30) ) - 3.0 * s ( :, 10)


      s ( :, 23) = 4.5 * ( s_quadpt ( :, 42)
     1 - s_quadpt ( :, 01) ) - 13.5 * ( s_quadpt ( :, 35)
     1 - s_quadpt ( :, 22) )

      s ( :, 22) = 2.25 * ( s_quadpt ( :, 42)
     1 + s_quadpt ( :, 01) - s_quadpt ( :, 35)
     1 - s_quadpt ( :, 22) ) - 1.5 * s ( :, 23)

      s ( :, 21) = s_quadpt ( :, 42) - s_quadpt ( :, 01)
     1 - s ( :, 22) - s ( :, 23)
     1 + ( s ( :, 30) + s ( :, 31)
     1 + s ( :, 32) ) / 12.0

! ----------

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
#endif
! End of "BUILD_SOURCE_TERMS == 1"

! All the modes of "s" are made up to the desired order by the end
! of this step.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 11) Returning the space-time integral of source terms that are treated
! time-explicitly is easy and is done here. If the treatment is time-implicit,
! the source term "s ( :, 1)" should be rebuilt here using
! "u ( :, :)" from the most recent iteration.

! IF "BUILD_SOURCE_TERMS == 0" then "s" should be zero, so "s_int" will be 0.

! In principle, the factor of "dtmove", that was previously multiplied into
! the source terms, is divided out of "s_int" here. However,
! because the "dudt" terms in the subroutine that calls this one already
! carry the source terms, we leave the "dtmove" multiplied into the
! "s_int" terms.

#if ( NUMORDER == 2)

      s_int ( :) = s ( :, 1) + 0.5 * s ( :, 5)

#endif

#if ( NUMORDER == 3)

      s_int ( :) = s ( :, 1) + 0.5 * s ( :, 11)
     1 + s ( :, 12) / 3.0

#endif

#if ( NUMORDER == 4)

      s_int ( :) = s ( :, 1) + 0.5 * s ( :, 21)
     1 + s ( :, 22) / 3.0 + 0.25 * s ( :, 23)

#endif

! ----------------------------------------------------------------------

! STEP 12) Here we give ourselves the option of trapping a zone that has
! not converged. This may consist of zeroing out higher order moments, i.e.
! we fall back to first order here if we absolutely must.

! This is also where we can use Dumbser's trick of zeroing out the
! higher moments in zones that produce negative pressures or densities.

      corrupted_zone = 0

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

! ----------

        IF ( ( rho ( iquadpt) .LT. rhofloor) .OR.
     1       ( prs ( iquadpt) .LT. prsfloor) )
     1                                      corrupted_zone = 1

! ----------

      END DO

! ------------------------------

! For ADER_CG ( but not for ADER_DG) "f ( :, 1)" is just the
! zone-centered flux at t=0. For ADER_DG, this will have to change.

! ----------

      IF ( corrupted_zone .EQ. 1) THEN

! -----

      u ( :, 1) = w ( :, 1)

      s ( :, 1) = 0.0

      s_int ( :) = 0.0

! -----

      w ( :, 2: NUM_ZONECEN_BASES) = 0.0

      u ( :, 2: NUM_SPACETIME_BASES) = 0.0

      f ( :, 2: NUM_SPACETIME_BASES) = 0.0

      g ( :, 2: NUM_SPACETIME_BASES) = 0.0

      h ( :, 2: NUM_SPACETIME_BASES) = 0.0

      s ( :, 2: NUM_SPACETIME_BASES) = 0.0

! -----

      END IF

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ADER_CG_MODAL_3D_MHD_PTWISE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADER_DG_MODAL_3D_TIME_EVOL_MHD

     1           ( oned_begin, oned_end, oned_off,
     1             dtmove,
     1             igeom, pi, rhofloor, prsfloor, smallnum,
     1             gamma, gamma_euler, gasconst,
     1             ader_tolerance, max_iter_ader_dg,
     1             tol_rootsolve, tol_superluminal,

     1             a_kerr, viscosity_coef, resistivity_coef,
     1           ambipolar_diffusion_drag_coef, molewt_ion, molewt_neut,
     1             isospeed, barotropic_stiffness,

     1             corrupted_zone,

     1             xc, yc, zc, dx, dy, dz,

     1             w, entropy_1d, gamma_eos_quadpt, molewt_quadpt,
     1             lorenz_quadpt, mu_quadpt,

     1             u, f, g, h, s, s_int)

! ----------------------------------------------------------------------

! Written by Dinshaw Balsara

! Incorporates GMRES-free and quadrature-free ADER DG with stiff source terms
! as contained in the paper by Balsara et al 2010. The present version
! accomodated Galerkin elements that are discontinuous in time and is,
! therefore, suitable for use with stiff source terms. It is useful for
! situations where the hyperbolic system is large and/or the Jacobian of
! the flux matrix is difficult to evaluate. Current version is for
! single-fluid MHD with multiple species. It scales linearly as the number of
! species and is, therefore, very useful for reactive flow.

! One MHD fluid and "N_EULER" number of Euler fluids (which can also be 
! repurposed as rediation fields in various frequency bins) are treated here.

! Takes the conserved variables and their spatial moments "w" as input
! and returns the space-time moments of the conserved variables as "u".
! This is all done on a reference space-time element in 3d that spans
! (x,y,z,t) in the range [-.5,.5]X[-.5,.5]X[-.5,.5]X[0,1].

! The exit condition is kept same for ADER CG and ADER DG so that the two
! subroutines can be used interchangably with the same value of
! "ader_tolerance".

! The ADER CG and ADER DG subroutines should always have the same interfaces.

! "dtmove, dx, dy, dz" are timestep and x,y,z zone sizes. The update has
! to be done in physical variables so various combinations of these
! variables are multiplied into the flux and source terms.

! "gamma, pi, rhofloor, prsfloor, smallnum" are the usual inputs.

! gamma law gas is assumed, but allowance is made for more general eos if
! the spatial moments are provided in "gamma_eos_spacetime" with the time terms zeroed.

! "ader_tolerance, max_iter_ader_dg" are the desired accuracy
! tolerance measured in the reference space-time element; the maximum number
! of permitted iterations and the number of iterations that have actually been
! taken to converge.

! "istiff" Set it to 1 for stiff source terms, 0 otherwise.

! "spacetime_testfn_int" gives the values of the space-time test function
! at the quadrature/nodal points in space and time.

! "coef_spacetime_setup" gives normalized values of the space-time test
! function at the quadrature points. Useful for projecting the nodal variables
! on to a modal basis set.

! "f, g, h, s" are fluxes in x,y,z-directions and the source terms respectively
! and are also returned. This is done in a zone-by-zone fashion.

! "s_int" is the space-time integral of the source term which is also returned.

! ----------------------------------------------------------------------

! It helps to know the space-time modal basis set into which the solution is
! being projected. Here it is for orders 2, 3, and 4.

! P0[x_] := 1 ; P1[x_] := x ; P2[x_] := ( x^2 - 1/12) ;
! P3[x_] := ( x^3 - 3 * x / 20)

! Q0[t_] := 1 ; Q1[t_] := 2 * t - 1 ; Q2[t_] := 6 * t^2 - 6 * t + 1 ;
! Q3[t_] := 20 * t^3 - 30 * t^2 + 12 * t - 1 ; DQ0[t_] := D[Q0[t],t] ;
! DQ1[t_] := D[Q1[t],t] ; DQ2[t_] := D[Q2[t],t] ; DQ3[t_] := D[Q3[t],t]

! ------------------------------

! Second order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P0[x] * P0[y] * P0[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
! [0,0,0,1/2], [1/2,0,0,1/2], [-1/2,0,0,1/2], [0,1/2,0,1/2], [0,-1/2,0,1/2],
! [0,0,1/2,1/2], [0,0,-1/2,1/2],
!   [0,0,0,1]

! ------------------------------

! Third order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P2[x] * P0[y] * P0[z] * Q0[t] +
!                   u6 * P0[x] * P2[y] * P0[z] * Q0[t] +
!                   u7 * P0[x] * P0[y] * P2[z] * Q0[t] +
!                   u8 * P1[x] * P1[y] * P0[z] * Q0[t] +
!                   u9 * P0[x] * P1[y] * P1[z] * Q0[t] +
!                   u10 * P1[x] * P0[y] * P1[z] * Q0[t] +
!                         u11 * P0[x] * P0[y] * P0[z] * Q1[t] +
!                         u12 * P0[x] * P0[y] * P0[z] * Q2[t] +
!                   u13 * P1[x] * P0[y] * P0[z] * Q1[t] +
!                   u14 * P0[x] * P1[y] * P0[z] * Q1[t] +
!                   u15 * P0[x] * P0[y] * P1[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
!   [0,0,0,1/3], [1/2,0,0,1/3], [-1/2,0,0,1/3], [0,1/2,0,1/3], [0,-1/2,0,1/3],
!   [0,0,1/2,1/3], [0,0,-1/2,1/3],
! [1/2,1/2,1/2,1/3],[-1/2,1/2,1/2,1/3], [1/2,-1/2,1/2,1/3], [-1/2,-1/2,1/2,1/3],
! [1/2,1/2,-1/2,1/3], [-1/2,1/2,-1/2,1/3], [1/2,-1/2,-1/2,1/3],
! [-1/2,-1/2,-1/2,1/3],
!   [1/2,0,0,2/3], [-1/2,0,0,2/3], [0,1/2,0,2/3], [0,-1/2,0,2/3],
!   [0,0,1/2,2/3], [0,0,-1/2,2/3],
! [0,0,0,1]

! ------------------------------

! Fourth order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P2[x] * P0[y] * P0[z] * Q0[t] +
!                   u6 * P0[x] * P2[y] * P0[z] * Q0[t] +
!                   u7 * P0[x] * P0[y] * P2[z] * Q0[t] +
!                   u8 * P1[x] * P1[y] * P0[z] * Q0[t] +
!                   u9 * P0[x] * P1[y] * P1[z] * Q0[t] +
!                   u10 * P1[x] * P0[y] * P1[z] * Q0[t] +
!                         u11 * P3[x] * P0[y] * P0[z] * Q0[t] +
!                         u12 * P0[x] * P3[y] * P0[z] * Q0[t] +
!                         u13 * P0[x] * P0[y] * P3[z] * Q0[t] +
!                         u14 * P2[x] * P1[y] * P0[z] * Q0[t] +
!                         u15 * P2[x] * P0[y] * P1[z] * Q0[t] +
!                         u16 * P1[x] * P2[y] * P0[z] * Q0[t] +
!                         u17 * P0[x] * P2[y] * P1[z] * Q0[t] +
!                         u18 * P1[x] * P0[y] * P2[z] * Q0[t] +
!                         u19 * P0[x] * P1[y] * P2[z] * Q0[t] +
!                         u20 * P1[x] * P1[y] * P1[z] * Q0[t] +
!                   u21 * P0[x] * P0[y] * P0[z] * Q1[t] +
!                   u22 * P0[x] * P0[y] * P0[z] * Q2[t] +
!                   u23 * P0[x] * P0[y] * P0[z] * Q3[t] +
!                         u24 * P1[x] * P0[y] * P0[z] * Q1[t] +
!                         u25 * P0[x] * P1[y] * P0[z] * Q1[t] +
!                         u26 * P0[x] * P0[y] * P1[z] * Q1[t] +
!                         u27 * P1[x] * P0[y] * P0[z] * Q2[t] +
!                         u28 * P0[x] * P1[y] * P0[z] * Q2[t] +
!                         u29 * P0[x] * P0[y] * P1[z] * Q2[t] +
!                   u30 * P2[x] * P0[y] * P0[z] * Q1[t] +
!                   u31 * P0[x] * P2[y] * P0[z] * Q1[t] +
!                   u32 * P0[x] * P0[y] * P2[z] * Q1[t] +
!                   u33 * P1[x] * P1[y] * P0[z] * Q1[t] +
!                   u34 * P0[x] * P1[y] * P1[z] * Q1[t] +
!                   u35 * P1[x] * P0[y] * P1[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
! [0,0,0,1/4], [1/2,0,0,1/4], [1/4,0,0,1/4], [-1/4,0,0,1/4], [-1/2,0,0,1/4],
! [0,1/2,0,1/4], [0,1/4,0,1/4], [0,-1/4,0,1/4], [0,-1/2,0,1/4],
! [0,0,1/2,1/4], [0,0,1/4,1/4], [0,0,-1/4,1/4], [0,0,-1/2,1/4],
!   [1/2,1/2,1/2,1/4], [0,1/2,1/2,1/4], [-1/2,1/2,1/2,1/4], [1/2,-1/2,1/2,1/4],
!   [0,-1/2,1/2,1/4], [-1/2,-1/2,1/2,1/4], [1/2,1/2,-1/2,1/4], [0,1/2,-1/2,1/4],
!   [-1/2,1/2,-1/2,1/4], [1/2,-1/2,-1/2,1/4], [0,-1/2,-1/2,1/4],[-1/2,-1/2,-1/2,1/4],
! [1/2,0,1/2,1/4], [-1/2,0,1/2,1/4], [1/2,0,-1/2,1/4], [-1/2,0,-1/2,1/4],
! [1/2,1/2,0,1/4], [-1/2,1/2,0,1/4], [1/2,-1/2,0,1/4], [-1/2,-1/2,0,1/4],
!   [0,0,0,1/2], [1/2,0,0,1/2], [-1/2,0,0,1/2], [0,1/2,0,1/2], [0,0,-1/2,1/2],
! [1/2,1/2,1/2,1/2],[-1/2,1/2,1/2,1/2], [1/2,-1/2,1/2,1/2], [-1/2,-1/2,1/2,1/2],
! [1/2,1/2,-1/2,1/2], [-1/2,1/2,-1/2,1/2], [1/2,-1/2,-1/2,1/2],
! [-1/2,-1/2,-1/2,1/2],
!   [1/2,0,0,3/4], [-1/2,0,0,3/4], [0,1/2,0,3/4], [0,-1/2,0,3/4],
!   [0,0,1/2,3/4], [0,0,-1/2,3/4],
! [0,0,0,1]

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER oned_begin, oned_end, oned_off, max_iter_ader_dg, igeom

      REAL dtmove, gamma, gamma_euler, pi, gasconst, ader_tolerance,
     1     tol_rootsolve, tol_superluminal,
     1     rhofloor, prsfloor, smallnum,
     1     a_kerr, viscosity_coef, resistivity_coef,
     1     ambipolar_diffusion_drag_coef, molewt_ion, molewt_neut,
     1     isospeed, barotropic_stiffness

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1                     corrupted_zone

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1                  xc, yc, zc, dx, dy, dz

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_ZONECEN_BASES) :: w

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, 5) :: entropy_1d

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  1: N_EULER + 1, NUM_SPACETIME_QUADPTS) ::
     1                gamma_eos_quadpt, molewt_quadpt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NUM_SPACETIME_QUADPTS) ::
     1                lorenz_quadpt, mu_quadpt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_SPACETIME_BASES) ::
     1                                              u, f, g, h, s

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, NVAR) :: s_int

! ----------------------------------------------------------------------

! Variables declared past this point do not go through the subroutine
! interface.

      INTEGER ix, iter_ader, iquadpt, ibases, ichk, itest, i_euler,
     1        ivar_lo, ivar

      REAL fourpi, ob4pi, ob8pi, ob3, tempaa, tempab, tempac, tempad,
     1     tempae, tempaf,
     1     x1, x2, x3, obx1, cos_x2, sin_x2, cot_x2

! These are the space-time-averaged values around which we construct the
! Jacobian of the source term during each iteration.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar,
     1     bfxbar, bfybar, bfzbar, gammabar, va_bar, vf_bar

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, 0: NFLUID) ::
     1     speciebar

#if ( N_EULER >= 1)
      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, N_EULER) ::
     1     rhobar_euler, prsbar_euler,
     1     vlxbar_euler, vlybar_euler, vlzbar_euler,
     1     gammabar_euler, vf_bar_euler

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1     0: NFLUID_EULER, N_EULER) :: speciebar_euler
#endif

! These are the flow variables (and their spatial derivatives, if needed)
! at the quadrature points.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NUM_SPACETIME_QUADPTS) ::

     1      rho, prs, eng_t, vlx, vly, vlz, bfx, bfy, bfz,
     1      entropy,

     1      drho_dx, dtemperature_dx, dvlx_dx, dvly_dx, dvlz_dx,
     1      dbfx_dx, dbfy_dx, dbfz_dx, dv2_dx, db2_dx,

     1      drho_dy, dtemperature_dy, dvlx_dy, dvly_dy, dvlz_dy,
     1      dbfx_dy, dbfy_dy, dbfz_dy, dv2_dy, db2_dy,

     1      drho_dz, dtemperature_dz, dvlx_dz, dvly_dz, dvlz_dz,
     1      dbfx_dz, dbfy_dz, dbfz_dz, dv2_dz, db2_dz,

     1      div_vl, temperature

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0: NFLUID, NUM_SPACETIME_QUADPTS) :: specie

#if ( N_EULER >= 1)
      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NUM_SPACETIME_QUADPTS, N_EULER) ::
     1      rho_euler, prs_euler, eng_t_euler,
     1      vlx_euler, vly_euler, vlz_euler

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1             0: NFLUID_EULER, NUM_SPACETIME_QUADPTS, N_EULER) ::
     1      specie_euler
#endif

! These are the conserved variables, fluxes and source terms at the
! quadrature points.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_SPACETIME_QUADPTS) ::

     1      u_quadpt, f_quadpt, g_quadpt, h_quadpt, s_quadpt,
     1      tu, tf, tg, th, ts,

     1      du_quadpt_dx, du_quadpt_dy, du_quadpt_dz,
     1      fni_quadpt, gni_quadpt, hni_quadpt,

     1      u_quadpt_prev

! These arrays help measure how well the ADER-DG iteration has converged.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, NVAR) ::
     1      convergence_array, tolerance_array

! List of nodal points in each zone.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  1:4, NUM_SPACETIME_QUADPTS) :: nodal_points

! rhs for the tridiagonal matrix inversion and the incremental change in
! the modal values at each ADER-DG iteration.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_SPACETIME_BASES) :: rhs_um, d_um

! ----------------------------------------------------------------------

! Variables that help out with the linear algebra stemming from
! the block tridiagonal matrices. Their dimensioning matches that
! of the tridiagonal solver.

      INTEGER j, iblock, nblocks_done,
     1        num_1by1_blocks, num_2by2_blocks,
     1        num_3by3_blocks, num_4by4_blocks,
     1        nriembegin, nriemend, mat_storage, mat_rank,
     1        nmax_blocks, n_blocks, new_call

      PARAMETER ( mat_storage = NVAR, mat_rank = NVAR,
     1            nmax_blocks = NUMORDER)

      INTEGER, DIMENSION ( nmax_blocks) :: iblock_index_begin,
     1      iblock_index_end

      REAL, DIMENSION ( nmax_blocks, nmax_blocks) :: a_mat_coef,
     1      b_mat_coef, c_mat_coef

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  mat_storage, mat_storage, nmax_blocks ) ::
     1      a, b, c, dinv, dinv_c

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  mat_storage, nmax_blocks ) :: r

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  mat_storage, mat_storage) :: E1

      REAL, DIMENSION ( mat_storage, mat_storage) :: Identity_Matrix

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Make some sanity checks here.

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_DG_MODAL_3D_TIME_EVOL_MHD, NFI","ELD .NE. 7",
     1              NFIELD
      STOP
      END IF

      IF ( NFIELD_EULER .NE. 5) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_DG_MODAL_3D_TIME_EVOL_MHD, NFI",
     1  "ELD_EULER .NE. 5", NFIELD
      STOP
      END IF

#if ( NUMORDER == 2)
      IF ( NUM_ZONECEN_BASES .NE. 4) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_BASES .NE. 4", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 5) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_BASES .NE. 5", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 6) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 6", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 7", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

#if ( NUMORDER == 3)
      IF ( NUM_ZONECEN_BASES .NE. 10) THEN
      WRITE ( 6, *)
     1  "Stopping in third order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_BASES .NE. 10", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 15) THEN
      WRITE ( 6, *)
     1  "Stopping in third order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_BASES .NE. 15", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 13) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 15", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 20) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 22", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

#if ( NUMORDER == 4)
      IF ( NUM_ZONECEN_BASES .NE. 20) THEN
      WRITE ( 6, *)
     1  "Stopping in fourth order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_BASES .NE. 20", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 35) THEN
      WRITE ( 6, *)
     1  "Stopping in fourth order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_BASES .NE. 35", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 21) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 33", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 42) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 55", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

      ivar =
     1  NFIELD + 1 + NFLUID + N_EULER * (NFIELD_EULER + NFLUID_EULER)

      IF ( NVAR .NE. ivar) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1  "NVAR .NE. ivar","NVA","R = ", NVAR," ivar = ", ivar
      STOP
      END IF

      IF ( max_iter_ader_dg .LT. NUMORDER ) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1  "max_iter_ader_dg .LT. NUMORDER ; max_iter_ader_dg = ",
     1   max_iter_ader_dg, NUMORDER
      STOP
      END IF

! ------------------------------

! Declare the number of blocks that will be used at each order.
! Use that to build the limits of the indices to be used in the
! block tridiagonal solve as well as the coefficients for the
! block matrices.

#if ( NUMORDER == 2)
      num_1by1_blocks = 3
      num_2by2_blocks = 1
      num_3by3_blocks = 0
      num_4by4_blocks = 0
#endif

#if ( NUMORDER == 3)
      num_1by1_blocks = 6
      num_2by2_blocks = 3
      num_3by3_blocks = 1
      num_4by4_blocks = 0
#endif

#if ( NUMORDER == 4)
      num_1by1_blocks = 10
      num_2by2_blocks = 6
      num_3by3_blocks = 3
      num_4by4_blocks = 1
#endif

! ----------

      a_mat_coef ( :, :) = 0.0
      b_mat_coef ( :, :) = 0.0
      c_mat_coef ( :, :) = 0.0

#if ( NUMORDER >= 2)
! Set starting and ending indices and coefficients for 1x1 and 2x2 blocks.

      iblock_index_begin ( 1) = NUM_SPACETIME_BASES
      iblock_index_end ( 1) = iblock_index_begin ( 1)
     1                      - num_1by1_blocks + 1

      iblock_index_begin ( 2) = iblock_index_end ( 1) - 1
      iblock_index_end ( 2) = iblock_index_begin ( 2)
     1                      - 2 * num_2by2_blocks + 1

      a_mat_coef ( 1, 1) = 1.0

      a_mat_coef ( 2, 1) = 0.5
      a_mat_coef ( 2, 2) = 1.0 / 6.0

      b_mat_coef ( 2, 1) = - 1.0 / 6.0

      c_mat_coef ( 2, 2) = 0.5

#endif

#if ( NUMORDER >= 3)
! Set starting and ending indices and coefficients for 3x3 blocks.

      iblock_index_begin ( 3) = iblock_index_end ( 2) - 1
      iblock_index_end ( 3) = iblock_index_begin ( 3)
     1                      - 3 * num_3by3_blocks + 1

      a_mat_coef ( 3, 1) = 0.5
      a_mat_coef ( 3, 2) = 0.0
      a_mat_coef ( 3, 3) = 0.1

      b_mat_coef ( 3, 1) = - 1.0 / 6.0
      b_mat_coef ( 3, 2) = - 0.1

      c_mat_coef ( 3, 2) = 0.5
      c_mat_coef ( 3, 3) = 1.0 / 6.0

#endif

#if ( NUMORDER >= 4)
! Set starting and ending indices and coefficients for 4x4 blocks.

      iblock_index_begin ( 4) = iblock_index_end ( 3) - 1
      iblock_index_end ( 4) = iblock_index_begin ( 4)
     1                      - 4 * num_4by4_blocks + 1

      a_mat_coef ( 4, 1) = 0.5
      a_mat_coef ( 4, 2) = 0.0
      a_mat_coef ( 4, 3) = 0.0
      a_mat_coef ( 4, 4) = 1.0 / 14.0

      b_mat_coef ( 4, 1) = - 1.0 / 6.0
      b_mat_coef ( 4, 2) = - 0.1
      b_mat_coef ( 4, 3) = - 1.0 / 14.0

      c_mat_coef ( 4, 2) = 0.5
      c_mat_coef ( 4, 3) = 1.0 / 6.0
      c_mat_coef ( 4, 4) = 0.1

#endif

! ----------

! Define "Identity_Matrix" here.

      Identity_Matrix ( :, :) = 0.0
      DO j = 1, mat_storage
      Identity_Matrix ( j, j) = 1.0
      END DO

! ----------------------------------------------------------------------

! Set up some constants that are used all over in this subroutine.

      fourpi = 4.0 * pi
      ob4pi = 1.0 / ( 4.0 * pi)
      ob8pi = 1.0 / ( 8.0 * pi)
      ob3 = 1.0 / 3.0


#if ( 0 == 1)


! Figure out what to do with these once we get to linear algebra part.


! Set up dual timestep to trivial initial value. Values greater than 1.0e32
! make the iteration revert to a pure Newton iteration anyway.
      dual_tau = 1.0e32

! Set up tolerance for linear algebra package.
      tolerance_lin_alg = 1.0e-12

! Set up the choice of exit criterion for linear algebra package.
      itol_lin_alg = 1
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 0) Set up the nodal points that are to be used at each of the orders.

! ------------------------------
#if ( NUMORDER == 2)
! cg, 2nd order nodal_points

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1:4, 01) = (/  0.5 ,  0.0 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 02) = (/ -0.5 ,  0.0 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 03) = (/  0.0 ,  0.5 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 04) = (/  0.0 , -0.5 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 05) = (/  0.0 ,  0.0 ,  0.5 , 0.5  /)
      nodal_points ( ix, 1:4, 06) = (/  0.0 ,  0.0 , -0.5 , 0.5  /)

      nodal_points ( ix, 1:4, 07) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /)

      END DO
! End of "DO ix =" loop.

#endif
! End "NUMORDER == 2".
! ------------------------------
#if ( NUMORDER == 3)
! cg, 3rd order nodal_points

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1:4, 01) = (/  0.0 ,  0.0 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 02) = (/  0.5 ,  0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 03) = (/  0.5 , -0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 04) = (/  0.5 ,  0.0 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 05) = (/  0.5 ,  0.0 , -0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 06) = (/ -0.5 ,  0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 07) = (/ -0.5 , -0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 08) = (/ -0.5 ,  0.0 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 09) = (/ -0.5 ,  0.0 , -0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 10) = (/  0.0 ,  0.5 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 11) = (/  0.0 , -0.5 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 12) = (/  0.0 ,  0.5 , -0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 13) = (/  0.0 , -0.5 , -0.5 , 1.0/3.0 /)

      nodal_points ( ix, 1:4, 14) = (/  0.5 ,  0.0 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 15) = (/ -0.5 ,  0.0 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 16) = (/  0.0 ,  0.5 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 17) = (/  0.0 , -0.5 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 18) = (/  0.0 ,  0.0 ,  0.5 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 19) = (/  0.0 ,  0.0 , -0.5 , 2.0/3.0 /)

      nodal_points ( ix, 1:4, 20) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /)

      END DO
! End of "DO ix =" loop.

#endif
! End "NUMORDER == 3".
! ------------------------------
#if ( NUMORDER == 4)
! cg, 4th order nodal_points

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1:4, 01) = (/  0.0 ,  0.0 ,  0.0 , 0.25 /)
      nodal_points ( ix, 1:4, 02) = (/  0.5 ,  0.0 ,  0.0 , 0.25 /)
      nodal_points ( ix, 1:4, 03) = (/  0.25,  0.0 ,  0.0 , 0.25 /)
      nodal_points ( ix, 1:4, 04) = (/ -0.25,  0.0 ,  0.0 , 0.25 /)
      nodal_points ( ix, 1:4, 05) = (/ -0.5 ,  0.0 ,  0.0 , 0.25 /)
      nodal_points ( ix, 1:4, 06) = (/  0.0 ,  0.5 ,  0.0 , 0.25 /)
      nodal_points ( ix, 1:4, 07) = (/  0.0 ,  0.25,  0.0 , 0.25 /)
      nodal_points ( ix, 1:4, 08) = (/  0.0 , -0.25,  0.0 , 0.25 /)
      nodal_points ( ix, 1:4, 09) = (/  0.0 , -0.5 ,  0.0 , 0.25 /)
      nodal_points ( ix, 1:4, 10) = (/  0.0 ,  0.0 ,  0.5 , 0.25 /)
      nodal_points ( ix, 1:4, 11) = (/  0.0 ,  0.0 ,  0.25, 0.25 /)
      nodal_points ( ix, 1:4, 12) = (/  0.0 ,  0.0 , -0.25, 0.25 /)
      nodal_points ( ix, 1:4, 13) = (/  0.0 ,  0.0 , -0.5 , 0.25 /)
      nodal_points ( ix, 1:4, 14) = (/  0.5 ,  0.5 ,  0.5 , 0.25 /)
      nodal_points ( ix, 1:4, 15) = (/  0.5 , -0.5 ,  0.5 , 0.25 /)
      nodal_points ( ix, 1:4, 16) = (/  0.5 ,  0.5 , -0.5 , 0.25 /)
      nodal_points ( ix, 1:4, 17) = (/  0.5 , -0.5 , -0.5 , 0.25 /)
      nodal_points ( ix, 1:4, 18) = (/ -0.5 ,  0.5 ,  0.5 , 0.25 /)
      nodal_points ( ix, 1:4, 19) = (/ -0.5 , -0.5 ,  0.5 , 0.25 /)
      nodal_points ( ix, 1:4, 20) = (/ -0.5 ,  0.5 , -0.5 , 0.25 /)
      nodal_points ( ix, 1:4, 21) = (/ -0.5 , -0.5 , -0.5 , 0.25 /)

      nodal_points ( ix, 1:4, 22) = (/  0.0 ,  0.0 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 23) = (/  0.5 ,  0.5 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 24) = (/  0.5 , -0.5 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 25) = (/  0.5 ,  0.0 ,  0.5 , 0.5  /)
      nodal_points ( ix, 1:4, 26) = (/  0.5 ,  0.0 , -0.5 , 0.5  /)
      nodal_points ( ix, 1:4, 27) = (/ -0.5 ,  0.5 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 28) = (/ -0.5 , -0.5 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 29) = (/ -0.5 ,  0.0 ,  0.5 , 0.5  /)
      nodal_points ( ix, 1:4, 30) = (/ -0.5 ,  0.0 , -0.5 , 0.5  /)
      nodal_points ( ix, 1:4, 31) = (/  0.0 ,  0.5 ,  0.5 , 0.5  /)
      nodal_points ( ix, 1:4, 32) = (/  0.0 , -0.5 ,  0.5 , 0.5  /)
      nodal_points ( ix, 1:4, 33) = (/  0.0 ,  0.5 , -0.5 , 0.5  /)
      nodal_points ( ix, 1:4, 34) = (/  0.0 , -0.5 , -0.5 , 0.5  /)

      nodal_points ( ix, 1:4, 35) = (/  0.0 ,  0.0 ,  0.0 , 0.75 /)
      nodal_points ( ix, 1:4, 36) = (/  0.5 ,  0.0 ,  0.0 , 0.75 /)
      nodal_points ( ix, 1:4, 37) = (/ -0.5 ,  0.0 ,  0.0 , 0.75 /)
      nodal_points ( ix, 1:4, 38) = (/  0.0 ,  0.5 ,  0.0 , 0.75 /)
      nodal_points ( ix, 1:4, 39) = (/  0.0 , -0.5 ,  0.0 , 0.75 /)
      nodal_points ( ix, 1:4, 40) = (/  0.0 ,  0.0 ,  0.5 , 0.75 /)
      nodal_points ( ix, 1:4, 41) = (/  0.0 ,  0.0 , -0.5 , 0.75 /)

      nodal_points ( ix, 1:4, 42) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /)

      END DO
! End of "DO ix =" loop.

#endif
! End "NUMORDER == 4".
! ------------------------------

#if ( ENFORCE_POSITIVITY == 1)

! While "nodal_points" still represents positions in the reference space-time
! hypercube, we use it to build the "entropy" at all the nodes.

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

      entropy ( ix, iquadpt) = entropy_1d ( ix, 01)
     1    + nodal_points ( ix, 1, iquadpt) * entropy_1d ( ix, 02)
     1    + nodal_points ( ix, 2, iquadpt) * entropy_1d ( ix, 03)
     1    + nodal_points ( ix, 3, iquadpt) * entropy_1d ( ix, 04)
     1    + nodal_points ( ix, 4, iquadpt) * entropy_1d ( ix, 05)

      END DO ! End of "DO ix =" loop.

      END DO ! End of "DO iquadpt =" loop.

#endif
! End of "ENFORCE_POSITIVITY == 1".

! ------------------------------

! Set up the nodal points in a global space-time. These are coordinate
! values.

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1, :) = nodal_points ( ix, 1, :) * dx ( ix)
     1                                                    + xc ( ix)
      nodal_points ( ix, 2, :) = nodal_points ( ix, 2, :) * dy ( ix)
     1                                                    + yc ( ix)
      nodal_points ( ix, 3, :) = nodal_points ( ix, 3, :) * dz ( ix)
     1                                                    + zc ( ix)
      nodal_points ( ix, 4, :) = nodal_points ( ix, 4, :) * dtmove

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 1) Set up arrays that control convergence.

! ------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      rhobar ( ix) = AMAX1 ( w ( ix, 1, 1), rhofloor)

      tempaa = 1.0 / rhobar ( ix)
      vlxbar ( ix) = w ( ix, 2, 1) * tempaa
      vlybar ( ix) = w ( ix, 3, 1) * tempaa
      vlzbar ( ix) = w ( ix, 4, 1) * tempaa

      prsbar ( ix) = w ( ix, 5, 1)

      bfxbar ( ix) = w ( ix, 6, 1)
      bfybar ( ix) = w ( ix, 7, 1)
      bfzbar ( ix) = w ( ix, 8, 1)

#if ( NUMORDER >= 3)
      gammabar ( ix) = gamma_eos_quadpt ( ix, 1, 1)
#else
      gammabar ( ix) = ( gamma_eos_quadpt ( ix, 1, 1)
     1                 + gamma_eos_quadpt ( ix, 1, 2)
     1                 + gamma_eos_quadpt ( ix, 1, 3)
     1                 + gamma_eos_quadpt ( ix, 1, 4)
     1                 + gamma_eos_quadpt ( ix, 1, 5)
     1                 + gamma_eos_quadpt ( ix, 1, 6) ) / 6.0
#endif

      prsbar ( ix) = prsbar ( ix) - 0.5 * rhobar ( ix)
     1 * ( vlxbar ( ix)**2 + vlybar ( ix)**2 + vlzbar ( ix)**2)
     1 - ( bfxbar ( ix)**2 + bfybar ( ix)**2 + bfzbar ( ix)**2) * ob8pi
      prsbar ( ix) = prsbar ( ix) * ( gammabar ( ix) - 1.0)

#if ( ISOTHERMAL_EQNS == 1)
      gammabar ( ix) = 1.1
      prsbar ( ix) = isospeed**2 * rhobar ( ix)
     1       + barotropic_stiffness * rhobar ( ix)**gamma
#endif

      va_bar ( ix) = ( bfxbar ( ix)**2 + bfybar ( ix)**2
     1                            + bfzbar ( ix)**2) * ob4pi * tempaa
      vf_bar ( ix) = gammabar ( ix) * prsbar ( ix) * tempaa
     1             + va_bar ( ix)
      va_bar ( ix) = SQRT ( va_bar ( ix))
      vf_bar ( ix) = SQRT ( vf_bar ( ix))

#if ( NFLUID >= 1)
      speciebar ( ix, 1: NFLUID) = 
     1           w ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, 1) * tempaa
#endif

! ----------

      tolerance_array ( ix, 1) = rhobar ( ix) * ader_tolerance
      tolerance_array ( ix, 2) = rhobar ( ix) * vf_bar ( ix)
     1                                           * ader_tolerance
      tolerance_array ( ix, 3) = tolerance_array ( ix, 2)
      tolerance_array ( ix, 4) = tolerance_array ( ix, 2)
      tolerance_array ( ix, 5) = rhobar ( ix) * vf_bar ( ix)**2
     1                                           * ader_tolerance
      tolerance_array ( ix, 6) = va_bar ( ix)
     1            * SQRT (fourpi * rhobar ( ix)) * ader_tolerance
      tolerance_array ( ix, 7) = tolerance_array ( ix, 6)
      tolerance_array ( ix, 8) = tolerance_array ( ix, 6)
#if ( NFLUID >= 1)
      tolerance_array ( ix, NFIELD + 2: NFIELD + 1 + NFLUID) =
     1                tolerance_array ( ix, 1)
#endif

! ----------

      END DO
! End of "DO ix =" loop.

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      rhobar_euler ( ix, i_euler) = 
     1                  AMAX1 ( w ( ix, ivar + 1, 1), rhofloor)

      tempaa = 1.0 / rhobar_euler ( ix, i_euler)
      vlxbar_euler ( ix, i_euler) = w ( ix, ivar + 2, 1) * tempaa
      vlybar_euler ( ix, i_euler) = w ( ix, ivar + 3, 1) * tempaa
      vlzbar_euler ( ix, i_euler) = w ( ix, ivar + 4, 1) * tempaa

      prsbar_euler ( ix, i_euler) = w ( ix, ivar + 5, 1)

#if ( NUMORDER >= 3)
      gammabar_euler ( ix, i_euler) = 
     1                     gamma_eos_quadpt ( ix, i_euler + 1, 1)
#else
      gammabar_euler ( ix, i_euler) =
     1           ( gamma_eos_quadpt ( ix, i_euler + 1, 1)
     1           + gamma_eos_quadpt ( ix, i_euler + 1, 2)
     1           + gamma_eos_quadpt ( ix, i_euler + 1, 3)
     1           + gamma_eos_quadpt ( ix, i_euler + 1, 4)
     1           + gamma_eos_quadpt ( ix, i_euler + 1, 5)
     1           + gamma_eos_quadpt ( ix, i_euler + 1, 6) ) / 6.0
#endif

      prsbar_euler ( ix, i_euler) = prsbar_euler ( ix, i_euler)
     1 - 0.5 * rhobar_euler ( ix, i_euler)
     1       * ( vlxbar_euler ( ix, i_euler)**2
     1         + vlybar_euler ( ix, i_euler)**2
     1         + vlzbar_euler ( ix, i_euler)**2)
      prsbar_euler ( ix, i_euler) = prsbar_euler ( ix, i_euler)
     1                        * ( gammabar_euler ( ix, i_euler) - 1.0)

#if ( ISOTHERMAL_EQNS == 1)
      gammabar_euler ( ix, i_euler) = 1.1
      prsbar_euler ( ix, i_euler) = isospeed**2
     1     * rhobar_euler ( ix, i_euler)
     1     + barotropic_stiffness * rhobar_euler ( ix, i_euler)**gamma
#endif

      vf_bar_euler ( ix, i_euler) = gammabar_euler ( ix, i_euler)
     1                        * prsbar_euler ( ix, i_euler) * tempaa
      vf_bar_euler ( ix, i_euler) = SQRT ( vf_bar_euler ( ix, i_euler))

#if ( NFLUID_EULER >= 1)
      speciebar_euler ( ix, 1: NFLUID_EULER, i_euler) =
     1        w ( ix, ivar + NFIELD_EULER + 1:
     1                ivar + NFIELD_EULER + NFLUID_EULER, 1)
     1                       * tempaa
#endif

! ----------

      tolerance_array ( ix, ivar + 1) = rhobar_euler ( ix, i_euler)
     1                            * ader_tolerance
      tolerance_array ( ix, ivar + 2) = rhobar_euler ( ix, i_euler)
     1                            * vf_bar_euler ( ix, i_euler)
     1                            * ader_tolerance
      tolerance_array ( ix, ivar + 3) = tolerance_array ( ix, ivar + 2)
      tolerance_array ( ix, ivar + 4) = tolerance_array ( ix, ivar + 2)
      tolerance_array ( ix, ivar + 5) = rhobar_euler ( ix, i_euler)
     1                            * vf_bar_euler ( ix, i_euler)**2
     1                            * ader_tolerance

#if ( NFLUID_EULER >= 1)
      tolerance_array ( ix, ivar + NFIELD_EULER + 1:
     1                      ivar + NFIELD_EULER + NFLUID_EULER) =
     1                              tolerance_array ( ix, ivar + 1)
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ------------------------------

      tolerance_array ( oned_begin - oned_off :
     1                       oned_end + oned_off, :) =
     1    tolerance_array ( oned_begin - oned_off :
     1                       oned_end + oned_off, :)**2

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 2) Having come into the subroutine with variables loaded into "w",
! load the appropriate spatial variation of "w" into "u" while setting
! the temporal variation of "u" to zero initially.

! ------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

      u ( ix, :, 1: NUM_ZONECEN_BASES) = 
     1                            w ( ix, :, 1: NUM_ZONECEN_BASES)
      u ( ix, :, NUM_ZONECEN_BASES + 1: NUM_SPACETIME_BASES) = 0.0

      u_quadpt_prev ( ix, :, :) = 0.0

      f ( ix, :, :) = 0.0
      g ( ix, :, :) = 0.0
      h ( ix, :, :) = 0.0
      s ( ix, :, :) = 0.0

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Start the outer ADER iteration here.

      DO iter_ader = 1, max_iter_ader_dg

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 3) Make the conserved variables at the space-time nodes.

! --------------------
#if ( NUMORDER == 2)
! dg, 2nd order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 01) = u ( ix, :, 01) + 0.5 * u ( ix, :, 02)

      u_quadpt ( ix, :, 02) = u ( ix, :, 01) - 0.5 * u ( ix, :, 02)

      u_quadpt ( ix, :, 03) = u ( ix, :, 01) + 0.5 * u ( ix, :, 03)

      u_quadpt ( ix, :, 04) = u ( ix, :, 01) - 0.5 * u ( ix, :, 03)

      u_quadpt ( ix, :, 05) = u ( ix, :, 01) + 0.5 * u ( ix, :, 04)

      u_quadpt ( ix, :, 06) = u ( ix, :, 01) - 0.5 * u ( ix, :, 04)

! -----

      u_quadpt ( ix, :, 07) = u ( ix, :, 01) + u ( ix, :, 05)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
! --------------------
#if ( NUMORDER == 3)
! dg, 3rd order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 01) = u ( ix, :, 01) + ( - u ( ix, :, 05)
     1 - u ( ix, :, 06) - u ( ix, :, 07) - 4.0 * u ( ix, :, 11)
     1 - 4.0 * u ( ix, :, 12) ) / 12.0

      u_quadpt ( ix, :, 02) = u_quadpt ( ix, :, 01)
     1 + ( 6.0 * u ( ix, :, 02) + 6.0 * u ( ix, :, 03)
     1 + 3.0 * u ( ix, :, 05) + 3.0 * u ( ix, :, 06)
     1 + 3.0 * u ( ix, :, 08) - 2.0 * u ( ix, :, 13)
     1 - 2.0 * u ( ix, :, 14) ) / 12.0

      u_quadpt ( ix, :, 03) = u_quadpt ( ix, :, 02)
     1 + ( - 6.0 * u ( ix, :, 03) - 3.0 * u ( ix, :, 08)
     1 + 2.0 * u ( ix, :, 14) ) / 6.0

      u_quadpt ( ix, :, 04) = 0.5 * ( u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 03) ) + ( 6.0 * u ( ix, :, 04)
     1 - 3.0 * u ( ix, :, 06) + 3.0 * u ( ix, :, 07)
     1 + 3.0 * u ( ix, :, 10) - 2.0 * u ( ix, :, 15) ) / 12.0

      u_quadpt ( ix, :, 05) = u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 03) - u_quadpt ( ix, :, 04)
     1 - 0.5 * ( u ( ix, :, 06) - u ( ix, :, 07) )

      u_quadpt ( ix, :, 06) = 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 03) + 0.5 * ( u ( ix, :, 05)
     1 + u ( ix, :, 06) - u ( ix, :, 08) )

      u_quadpt ( ix, :, 07) = 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 02) + 0.5 * ( u ( ix, :, 05)
     1 + u ( ix, :, 06) + u ( ix, :, 08) )

      u_quadpt ( ix, :, 08) = 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 05) + 0.5 * ( u ( ix, :, 05)
     1 + u ( ix, :, 07) - u ( ix, :, 10) )

      u_quadpt ( ix, :, 09) = 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 04) + 0.5 * ( u ( ix, :, 05)
     1 + u ( ix, :, 07) + u ( ix, :, 10) )

      u_quadpt ( ix, :, 10) = 0.5 * ( u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 04) + u_quadpt ( ix, :, 06)
     1 + u_quadpt ( ix, :, 08) ) - u_quadpt ( ix, :, 01)
     1 - 0.5 * u ( ix, :, 05) + 0.25 * u ( ix, :, 09)

      u_quadpt ( ix, :, 11) = u_quadpt ( ix, :, 04)
     1 + u_quadpt ( ix, :, 08) - u_quadpt ( ix, :, 10)
     1 - 0.5 * ( u ( ix, :, 05) - u ( ix, :, 06) )

      u_quadpt ( ix, :, 12) = 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 11) + 0.5 * ( u ( ix, :, 06)
     1 + u ( ix, :, 07) - u ( ix, :, 09) )

      u_quadpt ( ix, :, 13) = 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 10) + 0.5 * ( u ( ix, :, 06)
     1 + u ( ix, :, 07) + u ( ix, :, 09) )

! ----------

      tu ( ix, :, 01) = u_quadpt ( ix, :, 01)
     1 + 2.0 * u ( ix, :, 11) / 3.0

      u_quadpt ( ix, :, 14) = tu ( ix, :, 01) + ( 6.0 * u ( ix, :, 02)
     1 + 3.0 * u ( ix, :, 05) + 2.0 * u ( ix, :, 13) ) / 12.0

      u_quadpt ( ix, :, 15) = 2.0 * tu ( ix, :, 01)
     1 - u_quadpt ( ix, :, 14) + 0.5 * u ( ix, :, 05)

      u_quadpt ( ix, :, 16) = tu ( ix, :, 01) + ( 6.0 * u ( ix, :, 03)
     1 + 3.0 * u ( ix, :, 06) + 2.0 * u ( ix, :, 14) ) / 12.0

      u_quadpt ( ix, :, 17) = 2.0 * tu ( ix, :, 01)
     1 - u_quadpt ( ix, :, 16) + 0.5 * u ( ix, :, 06)

      u_quadpt ( ix, :, 18) = tu ( ix, :, 01) + ( 6.0 * u ( ix, :, 04)
     1 + 3.0 * u ( ix, :, 07) + 2.0 * u ( ix, :, 15) ) / 12.0

      u_quadpt ( ix, :, 19) = 2.0 * tu ( ix, :, 01)
     1 - u_quadpt ( ix, :, 18) + 0.5 * u ( ix, :, 07)

! ----------

      u_quadpt ( ix, :, 20) = u_quadpt ( ix, :, 01)
     1 + 4.0 * ( u ( ix, :, 11) + u ( ix, :, 12) ) / 3.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
! --------------------
#if ( NUMORDER == 4)
! dg, 4th order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 01) = ( 48.0 * u ( ix, :, 01)
     1 - 4.0 * u ( ix, :, 05) - 4.0 * u ( ix, :, 06)
     1 - 4.0 * u ( ix, :, 07) - 24.0 * u ( ix, :, 21)
     1 - 6.0 * u ( ix, :, 22) + 21.0 * u ( ix, :, 23)
     1 + 2.0 * u ( ix, :, 30) + 2.0 * u ( ix, :, 31)
     1 + 2.0 * u ( ix, :, 32) ) / 48.0

      u_quadpt ( ix, :, 02) = u_quadpt ( ix, :, 01)
     1 + ( 120.0 * u ( ix, :, 02) + 60.0 * u ( ix, :, 05)
     1 + 12.0 * u ( ix, :, 11) - 10.0 * u ( ix, :, 16)
     1 - 10.0 * u ( ix, :, 18) - 60.0 * u ( ix, :, 24)
     1 - 15.0 * u ( ix, :, 27) - 30.0 * u ( ix, :, 30) ) / 240.0

      u_quadpt ( ix, :, 03) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 02) ) + 0.015625 * ( - 4.0 * u ( ix, :, 05)
     1 - 3.0 * u ( ix, :, 11) + 2.0 * u ( ix, :, 30) )

      u_quadpt ( ix, :, 04) = 3.0 * ( u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 03) ) + u_quadpt ( ix, :, 02)
     1 - 0.09375 * u ( ix, :, 11)

      u_quadpt ( ix, :, 05) = 3.0 * ( u_quadpt ( ix, :, 04)
     1 - u_quadpt ( ix, :, 01) ) + u_quadpt ( ix, :, 03)
     1 - 0.09375 * u ( ix, :, 11)

      u_quadpt ( ix, :, 06) = u_quadpt ( ix, :, 01)
     1 + ( 120.0 * u ( ix, :, 03) + 60.0 * u ( ix, :, 06)
     1 + 12.0 * u ( ix, :, 12) - 10.0 * u ( ix, :, 14)
     1 - 10.0 * u ( ix, :, 19) - 60.0 * u ( ix, :, 25)
     1 - 15.0 * u ( ix, :, 28) - 30.0 * u ( ix, :, 31) ) / 240.0

      u_quadpt ( ix, :, 07) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 06) ) + 0.015625 * ( - 4.0 * u ( ix, :, 06)
     1 - 3.0 * u ( ix, :, 12) + 2.0 * u ( ix, :, 31) )

      u_quadpt ( ix, :, 08) = 3.0 * ( u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 07) ) + u_quadpt ( ix, :, 06)
     1 - 0.09375 * u ( ix, :, 12)

      u_quadpt ( ix, :, 09) = 3.0 * ( u_quadpt ( ix, :, 08)
     1 - u_quadpt ( ix, :, 01) ) + u_quadpt ( ix, :, 07)
     1 - 0.09375 * u ( ix, :, 12)

      u_quadpt ( ix, :, 10) = u_quadpt ( ix, :, 01)
     1 + ( 120.0 * u ( ix, :, 04) + 60.0 * u ( ix, :, 07)
     1 + 12.0 * u ( ix, :, 13) - 10.0 * u ( ix, :, 15)
     1 - 10.0 * u ( ix, :, 17) - 60.0 * u ( ix, :, 26)
     1 - 15.0 * u ( ix, :, 29) - 30.0 * u ( ix, :, 32) ) / 240.0

      u_quadpt ( ix, :, 11) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 10) ) + 0.015625 * ( - 4.0 * u ( ix, :, 07)
     1 - 3.0 * u ( ix, :, 13) + 2.0 * u ( ix, :, 32) )

      u_quadpt ( ix, :, 12) = 3.0 * ( u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 11) ) + u_quadpt ( ix, :, 10)
     1 - 0.09375 * u ( ix, :, 13)

      u_quadpt ( ix, :, 13) = 3.0 * ( u_quadpt ( ix, :, 12)
     1 - u_quadpt ( ix, :, 01) ) + u_quadpt ( ix, :, 11)
     1 - 0.09375 * u ( ix, :, 13)

      u_quadpt ( ix, :, 14) = u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 06) + u_quadpt ( ix, :, 10)
     1 - 2.0 * u_quadpt ( ix, :, 01) + 0.125 * ( 2.0 * u ( ix, :, 08)
     1 + 2.0 * u ( ix, :, 09) + 2.0 * u ( ix, :, 10) + u ( ix, :, 14)
     1 + u ( ix, :, 15) + u ( ix, :, 16) + u ( ix, :, 17)
     1 + u ( ix, :, 18) + u ( ix, :, 19) + u ( ix, :, 20)
     1 - u ( ix, :, 33) - u ( ix, :, 34) - u ( ix, :, 35) )

      u_quadpt ( ix, :, 15) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 09) - u_quadpt ( ix, :, 06)
     1 - 0.5 * ( u ( ix, :, 08) + u ( ix, :, 09) )
     1 + 0.25 * ( - u ( ix, :, 14) - u ( ix, :, 19)
     1 - u ( ix, :, 20) + u ( ix, :, 33) + u ( ix, :, 34) )

      u_quadpt ( ix, :, 16) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 13) - u_quadpt ( ix, :, 10)
     1 - 0.5 * ( u ( ix, :, 09) + u ( ix, :, 10))
     1 + 0.25 * ( - u ( ix, :, 15) - u ( ix, :, 17)
     1 - u ( ix, :, 20) + u ( ix, :, 34) + u ( ix, :, 35) )

      u_quadpt ( ix, :, 17) = u_quadpt ( ix, :, 15)
     1 + u_quadpt ( ix, :, 16) - u_quadpt ( ix, :, 14)
     1 + u ( ix, :, 09) + 0.5 * ( u ( ix, :, 20) - u ( ix, :, 34) )

      u_quadpt ( ix, :, 18) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 02)
     1 - 0.5 * ( u ( ix, :, 08) + u ( ix, :, 10) )
     1 + 0.25 * ( - u ( ix, :, 16) - u ( ix, :, 18)
     1 - u ( ix, :, 20) + u ( ix, :, 33) + u ( ix, :, 35) )

      u_quadpt ( ix, :, 19) = u_quadpt ( ix, :, 15)
     1 + u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 02)
     1 - 0.5 * ( u ( ix, :, 10) - u ( ix, :, 08) )
     1 + 0.25 * ( - u ( ix, :, 16) - u ( ix, :, 18)
     1 + u ( ix, :, 20) - u ( ix, :, 33) + u ( ix, :, 35) )

      u_quadpt ( ix, :, 20) = u_quadpt ( ix, :, 16)
     1 + u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 02)
     1 - 0.5 * ( u ( ix, :, 08) - u ( ix, :, 10) )
     1 + 0.25 * ( - u ( ix, :, 16) - u ( ix, :, 18)
     1 + u ( ix, :, 20) + u ( ix, :, 33) - u ( ix, :, 35) )

      u_quadpt ( ix, :, 21) = u_quadpt ( ix, :, 19)
     1 + u_quadpt ( ix, :, 20) - u_quadpt ( ix, :, 18)
     1 + u ( ix, :, 09) - 0.5 * ( u ( ix, :, 20) + u ( ix, :, 34) )

! ----------

      u_quadpt ( ix, :, 22) = u_quadpt ( ix, :, 01)
     1 + ( 24.0 * u ( ix, :, 21) - 18.0 * u ( ix, :, 22)
     1 - 21.0 * u ( ix, :, 23) - 2.0 * u ( ix, :, 30)
     1 - 2.0 * u ( ix, :, 31) - 2.0 * u ( ix, :, 32) ) / 48.0

      u_quadpt ( ix, :, 23) = u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 06) + u_quadpt ( ix, :, 22)
     1 - 2.0 * u_quadpt ( ix, :, 01) + 0.0625 * ( 4.0 * u ( ix, :, 08)
     1 + 2.0 * u ( ix, :, 14) + 2.0 * u ( ix, :, 16)
     1 + 4.0 * u ( ix, :, 24) + 4.0 * u ( ix, :, 25)
     1 - 3.0 * u ( ix, :, 27) - 3.0 * u ( ix, :, 28)
     1 + 2.0 * u ( ix, :, 30) + 2.0 * u ( ix, :, 31) )

      u_quadpt ( ix, :, 24) = u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 09) - u_quadpt ( ix, :, 06)
     1 - 0.5 * ( u ( ix, :, 08) + u ( ix, :, 25) )
     1 - 0.25 * u ( ix, :, 14) + 0.375 * u ( ix, :, 28)

      u_quadpt ( ix, :, 25) = ( u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 24) + u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 15)) / 2.0 - u_quadpt ( ix, :, 02)
     1 - 0.5 * u ( ix, :, 06) - 0.25 * ( u ( ix, :, 16)
     1 - u ( ix, :, 26) ) - 0.125 * ( u ( ix, :, 17)
     1 - u ( ix, :, 31) - u ( ix, :, 32) - u ( ix, :, 35) )
     1 - 0.1875 * u ( ix, :, 29)

      u_quadpt ( ix, :, 26) = u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 24) - u_quadpt ( ix, :, 25)
     1 - 0.5 * ( u ( ix, :, 06) - u ( ix, :, 07) )
     1 - 0.25 * ( u ( ix, :, 16) - u ( ix, :, 18) )

      u_quadpt ( ix, :, 27) = 2.0 * u_quadpt ( ix, :, 22)
     1 - u_quadpt ( ix, :, 24) + 0.5 * ( u ( ix, :, 05)
     1 + u ( ix, :, 06) - u ( ix, :, 08) )

      u_quadpt ( ix, :, 28) = 2.0 * u_quadpt ( ix, :, 22)
     1 - u_quadpt ( ix, :, 23) + 0.5 * ( u ( ix, :, 05)
     1 + u ( ix, :, 06) + u ( ix, :, 08) )

      u_quadpt ( ix, :, 29) = 2.0 * u_quadpt ( ix, :, 22)
     1 - u_quadpt ( ix, :, 26) + 0.5 * ( u ( ix, :, 05)
     1 + u ( ix, :, 07) - u ( ix, :, 10) )

      u_quadpt ( ix, :, 30) = 2.0 * u_quadpt ( ix, :, 22)
     1 - u_quadpt ( ix, :, 25) + 0.5 * ( u ( ix, :, 05)
     1 + u ( ix, :, 07) + u ( ix, :, 10) )

      u_quadpt ( ix, :, 31) = 0.5 * ( u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 25) + u_quadpt ( ix, :, 27)
     1 + u_quadpt ( ix, :, 29) ) - u_quadpt ( ix, :, 22)
     1 - 0.5 * u ( ix, :, 05) + 0.25 * u ( ix, :, 09)
     1 - 0.125 * ( u ( ix, :, 14) + u ( ix, :, 15) - u ( ix, :, 17)
     1 - u ( ix, :, 19) )

      u_quadpt ( ix, :, 32) = u_quadpt ( ix, :, 25)
     1 + u_quadpt ( ix, :, 29) - u_quadpt ( ix, :, 31)
     1 - 0.5 * ( u ( ix, :, 05) - u ( ix, :, 06) )
     1 - 0.25 * ( u ( ix, :, 15) - u ( ix, :, 17) )

      u_quadpt ( ix, :, 33) = 2.0 * u_quadpt ( ix, :, 22)
     1 - u_quadpt ( ix, :, 32) + 0.5 * ( u ( ix, :, 06)
     1 + u ( ix, :, 07) - u ( ix, :, 09) )


      u_quadpt ( ix, :, 34) = 2.0 * u_quadpt ( ix, :, 22)
     1 - u_quadpt ( ix, :, 31) + 0.5 * ( u ( ix, :, 06)
     1 + u ( ix, :, 07) + u ( ix, :, 09) )

! ----------

      u_quadpt ( ix, :, 35) = 2.0 * u_quadpt ( ix, :, 22)
     1 - u_quadpt ( ix, :, 01) + 0.75 * u ( ix, :, 22)

      u_quadpt ( ix, :, 36) = u_quadpt ( ix, :, 35)
     1 + u_quadpt ( ix, :, 02) - u_quadpt ( ix, :, 01)
     1 + 0.5 * u ( ix, :, 24) + 0.25 * u ( ix, :, 30)

      u_quadpt ( ix, :, 37) = 2.0 * u_quadpt ( ix, :, 35)
     1 - u_quadpt ( ix, :, 36) + 0.5 * u ( ix, :, 05)
     1 + 0.25 * u ( ix, :, 30)

      u_quadpt ( ix, :, 38) = u_quadpt ( ix, :, 35)
     1 + u_quadpt ( ix, :, 06) - u_quadpt ( ix, :, 01)
     1 + 0.5 * u ( ix, :, 25) + 0.25 * u ( ix, :, 31)

      u_quadpt ( ix, :, 39) = 2.0 * u_quadpt ( ix, :, 35)
     1 - u_quadpt ( ix, :, 38) + 0.5 * u ( ix, :, 06)
     1 + 0.25 * u ( ix, :, 31)

      u_quadpt ( ix, :, 40) = u_quadpt ( ix, :, 35)
     1 + u_quadpt ( ix, :, 10) - u_quadpt ( ix, :, 01)
     1 + 0.5 * u ( ix, :, 26) + 0.25 * u ( ix, :, 32)

      u_quadpt ( ix, :, 41) = 2.0 * u_quadpt ( ix, :, 35)
     1 - u_quadpt ( ix, :, 40) + 0.5 * u ( ix, :, 07)
     1 + 0.25 * u ( ix, :, 32)

! ----------

      u_quadpt ( ix, :, 42) = 3.0 * ( u_quadpt ( ix, :, 35)
     1 - u_quadpt ( ix, :, 22) ) + u_quadpt ( ix, :, 01)
     1 + 1.875 * u ( ix, :, 23)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! --------------------

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 4) Check for convergence and evaluate whether it is acceptable to
! leave the main ADER iteration loop. Do that by comparing convergence
! to the desired level at the quadrature points across iterations.
! We can, therefore, do this step immediately after a new "u_quadpt" is built.

! The exit condition is kept same for ADER CG and ADER DG so that the two
! subroutines can be used interchangably with the same value of
! "ader_tolerance".


      convergence_array ( oned_begin-oned_off: oned_end+oned_off, :)=0.0

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

        DO ix = oned_begin - oned_off, oned_end + oned_off

      convergence_array ( ix, :) = convergence_array ( ix, :)
     1   + ( u_quadpt ( ix, :, iquadpt)
     1     - u_quadpt_prev ( ix, :, iquadpt))**2

        END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

      tempaa = 1.0 / ( NUM_SPACETIME_QUADPTS )

      convergence_array ( oned_begin-oned_off: oned_end+oned_off, :) =
     1 convergence_array ( oned_begin-oned_off: oned_end+oned_off, :)
     1 * tempaa


! Check to see if we are close to convergence. If so, we can exit.
! NOTE : Eventually "itest" can become a 1d array and the evaluation
! of zones can be switched off if it is unity in that zone.

      IF ( iter_ader .GT. 1) THEN

        itest = 1

        DO ichk = 1, NVAR

          DO ix = oned_begin - oned_off, oned_end + oned_off

          IF ( convergence_array ( ix, ichk) .GT.
     1         tolerance_array ( ix, ichk)) itest = 0

          END DO
! End of "DO ix =" loop.

        END DO

        IF ( itest .EQ. 1) EXIT

      END IF


! Now hand back data to the previous set of conserved variables.

      u_quadpt_prev ( oned_begin-oned_off: oned_end+oned_off, :, :) =
     1  u_quadpt ( oned_begin-oned_off: oned_end+oned_off, :, :)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 5) Build space-time-averaged primitive variables, "rhobar, prsbar" etc.,
! that can be used as "good" variables around which to linearize
! implicit source terms during this round of the iteration.

! These variables can even be used to build a new value of "tolerance_array",
! though that has never seemed necessary.

! ------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      rhobar ( ix) = AMAX1 ( u ( ix, 1, 1), rhofloor)

      tempaa = 1.0 / rhobar ( ix)
      vlxbar ( ix) = u ( ix, 2, 1) * tempaa
      vlybar ( ix) = u ( ix, 3, 1) * tempaa
      vlzbar ( ix) = u ( ix, 4, 1) * tempaa

      prsbar ( ix) = u ( ix, 5, 1)

      bfxbar ( ix) = u ( ix, 6, 1)
      bfybar ( ix) = u ( ix, 7, 1)
      bfzbar ( ix) = u ( ix, 8, 1)

#if ( NUMORDER >= 3)
      gammabar ( ix) = gamma_eos_quadpt ( ix, 1, 1)
#else
      gammabar ( ix) = ( gamma_eos_quadpt ( ix, 1, 1)
     1                 + gamma_eos_quadpt ( ix, 1, 2)
     1                 + gamma_eos_quadpt ( ix, 1, 3)
     1                 + gamma_eos_quadpt ( ix, 1, 4)
     1                 + gamma_eos_quadpt ( ix, 1, 5)
     1                 + gamma_eos_quadpt ( ix, 1, 6) ) / 6.0
#endif

      prsbar ( ix) = prsbar ( ix) - 0.5 * rhobar ( ix)
     1 * ( vlxbar ( ix)**2 + vlybar ( ix)**2 + vlzbar ( ix)**2)
     1 - ( bfxbar ( ix)**2 + bfybar ( ix)**2 + bfzbar ( ix)**2) * ob8pi
      prsbar ( ix) = prsbar ( ix) * ( gammabar ( ix) - 1.0)
#if ( ENFORCE_POSITIVITY == 1)
      IF ( prsbar ( ix) .LE. prsfloor) prsbar ( ix) =
     1     entropy_1d ( ix, 01) * rhobar ( ix)**( gammabar ( ix) - 1.0)
#endif
      prsbar ( ix) = AMAX1 ( prsbar ( ix), prsfloor)

#if ( NFLUID >= 1)
      speciebar ( ix, 1: NFLUID) =
     1           u ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, 1) * tempaa
#endif

! ----------

      END DO
! End of "DO ix =" loop.

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      rhobar_euler ( ix, i_euler) =
     1                  AMAX1 ( u ( ix, ivar + 1, 1), rhofloor)

      tempaa = 1.0 / rhobar_euler ( ix, i_euler)
      vlxbar_euler ( ix, i_euler) = u ( ix, ivar + 2, 1) * tempaa
      vlybar_euler ( ix, i_euler) = u ( ix, ivar + 3, 1) * tempaa
      vlzbar_euler ( ix, i_euler) = u ( ix, ivar + 4, 1) * tempaa

      prsbar_euler ( ix, i_euler) = u ( ix, ivar + 5, 1)

#if ( NUMORDER >= 3)
      gammabar_euler ( ix, i_euler) =
     1                     gamma_eos_quadpt ( ix, i_euler + 1, 1)
#else
      gammabar_euler ( ix, i_euler) =
     1           ( gamma_eos_quadpt ( ix, i_euler + 1, 1)
     1           + gamma_eos_quadpt ( ix, i_euler + 1, 2)
     1           + gamma_eos_quadpt ( ix, i_euler + 1, 3)
     1           + gamma_eos_quadpt ( ix, i_euler + 1, 4)
     1           + gamma_eos_quadpt ( ix, i_euler + 1, 5)
     1           + gamma_eos_quadpt ( ix, i_euler + 1, 6) ) / 6.0
#endif

      prsbar_euler ( ix, i_euler) = prsbar_euler ( ix, i_euler)
     1 - 0.5 * rhobar_euler ( ix, i_euler)
     1       * ( vlxbar_euler ( ix, i_euler)**2
     1         + vlybar_euler ( ix, i_euler)**2
     1         + vlzbar_euler ( ix, i_euler)**2)
      prsbar_euler ( ix, i_euler) = prsbar_euler ( ix, i_euler)
     1                        * ( gammabar_euler ( ix, i_euler) - 1.0)

#if ( ISOTHERMAL_EQNS == 1)
      gammabar_euler ( ix, i_euler) = 1.1
      prsbar_euler ( ix, i_euler) = isospeed**2
     1     * rhobar_euler ( ix, i_euler)
     1     + barotropic_stiffness * rhobar_euler ( ix, i_euler)**gamma
#endif

#if ( NFLUID_EULER >= 1)
      speciebar_euler ( ix, 1: NFLUID_EULER, i_euler) =
     1        u ( ix, ivar + NFIELD_EULER + 1:
     1                ivar + NFIELD_EULER + NFLUID_EULER, 1)
     1                       * tempaa
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 6) Get the primitive variables from the conserved variables at
! the nodes. The density and pressure variables are
! given positive floor values to ensure that they can be used for physically
! meaningful source terms.

! ------------------------------

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, 1, iquadpt) =
     1                   AMAX1 ( u_quadpt ( ix, 1, iquadpt), rhofloor)

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempaa = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempaa
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempaa
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempaa

      bfx ( ix, iquadpt) = u_quadpt ( ix, 6, iquadpt)
      bfy ( ix, iquadpt) = u_quadpt ( ix, 7, iquadpt)
      bfz ( ix, iquadpt) = u_quadpt ( ix, 8, iquadpt)

      tempab = 0.5 * rho ( ix, iquadpt) * ( vlx ( ix, iquadpt)**2
     1         + vly ( ix, iquadpt)**2 + vlz ( ix, iquadpt)**2)

      tempac = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      eng_t ( ix, iquadpt) = u_quadpt ( ix, 5, iquadpt)
     1                     - tempab - tempac

      prs ( ix, iquadpt) = eng_t ( ix, iquadpt)
     1                   * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

#if ( ENFORCE_POSITIVITY == 1)
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = entropy ( ix, iquadpt)
     1* rho ( ix, iquadpt)**( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)
      prs ( ix, iquadpt) = AMAX1 ( prs ( ix, iquadpt), prsfloor)

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac 

      END IF
#else
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = prsfloor

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#endif

#if ( NFLUID >= 1)
      specie ( ix, 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( ix, 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( ix, 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( ix, 1: NFLUID, iquadpt) * rho ( ix, iquadpt)
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, ivar + 1, iquadpt) =
     1         AMAX1 ( u_quadpt ( ix, ivar + 1, iquadpt), rhofloor)

      rho_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 1, iquadpt)

      tempaa = 1.0 / rho_euler ( ix, iquadpt, i_euler)

      vlx_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 2, iquadpt) * tempaa
      vly_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 3, iquadpt) * tempaa
      vlz_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 4, iquadpt) * tempaa

      tempab = 0.5 * rho_euler ( ix, iquadpt, i_euler)
     1           * ( vlx_euler ( ix, iquadpt, i_euler)**2
     1             + vly_euler ( ix, iquadpt, i_euler)**2
     1             + vlz_euler ( ix, iquadpt, i_euler)**2)

      prs_euler ( ix, iquadpt, i_euler) =
     1         ( u_quadpt ( ix, ivar + 5, iquadpt) - tempab)
     1       * ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      prs_euler ( ix, iquadpt, i_euler) =
     1    AMAX1 ( prs_euler ( ix, iquadpt, i_euler), prsfloor)

      eng_t_euler ( ix, iquadpt, i_euler) =
     1         prs_euler ( ix, iquadpt, i_euler)
     1       / ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      IF ( prs_euler ( ix, iquadpt, i_euler) .LT. prsfloor * 1.00000001)
     1   u_quadpt ( ix, ivar + 5, iquadpt) =
     1                   eng_t_euler ( ix, iquadpt, i_euler) + tempab

#if ( NFLUID_EULER >= 1)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) = tempaa
     1       * u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1                    ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) =
     1    AMIN1 ( AMAX1 (
     1    specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler),
     1            0.0), 1.0 )

      u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1       specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler)
     1     * rho_euler ( ix, iquadpt, i_euler)

#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 7.a) Build fluxes ( f_quadpt, g_quadpt, h_quadpt) from the primitive
! variables that have now become available at the space-time nodes.

! Multiply "f_quadpt" by "dtmove/dx" and so on.

! ------------------------------

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! -----

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

! -----

      f_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 2, iquadpt)

      f_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfx ( ix, iquadpt)**2 * ob4pi

      f_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vly ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi

      f_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi

      f_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab

      f_quadpt ( ix, 6, iquadpt) = 0.0

      f_quadpt ( ix, 7, iquadpt) =
     1           vlx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         - vly ( ix, iquadpt) * bfx ( ix, iquadpt)

      f_quadpt ( ix, 8, iquadpt) =
     1           vlx ( ix, iquadpt) * bfz ( ix, iquadpt)
     1         - vlz ( ix, iquadpt) * bfx ( ix, iquadpt)

#if ( NFLUID >= 1)
      f_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlx ( ix, iquadpt)
#endif

! -----

      g_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 3, iquadpt)

      g_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         - bfy ( ix, iquadpt) * bfx ( ix, iquadpt) * ob4pi

      g_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vly ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfy ( ix, iquadpt)**2 * ob4pi

      g_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         - bfy ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi

      g_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vly ( ix, iquadpt)
     1                       - bfy ( ix, iquadpt) * tempab

      g_quadpt ( ix, 6, iquadpt) =
     1           vly ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         - vlx ( ix, iquadpt) * bfy ( ix, iquadpt)

      g_quadpt ( ix, 7, iquadpt) = 0.0

      g_quadpt ( ix, 8, iquadpt) =
     1           vly ( ix, iquadpt) * bfz ( ix, iquadpt)
     1         - vlz ( ix, iquadpt) * bfy ( ix, iquadpt)

#if ( NFLUID >= 1)
      g_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vly ( ix, iquadpt)
#endif

! -----

      h_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 4, iquadpt)

      h_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         - bfz ( ix, iquadpt) * bfx ( ix, iquadpt) * ob4pi

      h_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vly ( ix, iquadpt)
     1         - bfz ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi

      h_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfz ( ix, iquadpt)**2 * ob4pi

      h_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vlz ( ix, iquadpt)
     1                       - bfz ( ix, iquadpt) * tempab

      h_quadpt ( ix, 6, iquadpt) =
     1           vlz ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         - vlx ( ix, iquadpt) * bfz ( ix, iquadpt)

      h_quadpt ( ix, 7, iquadpt) =
     1           vlz ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         - vly ( ix, iquadpt) * bfz ( ix, iquadpt)

      h_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      h_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlz ( ix, iquadpt)
#endif

! -----

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! -----

      f_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)

      f_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vlx_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      f_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vlx_euler ( ix, iquadpt, i_euler)
#endif

! -----

      g_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)

      g_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vly_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      g_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vly_euler ( ix, iquadpt, i_euler)
#endif

! -----

      h_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)

      h_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vlz_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      h_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vlz_euler ( ix, iquadpt, i_euler)
#endif

! -----

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ------------------------------

! To make cylindrical and spherical geometries accessible to non-rel.
! and relativistic flow, we put this part after either flux construction.

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      IF ( igeom == 2 .or. igeom == 3) THEN

         g_quadpt ( ix, :, iquadpt) = g_quadpt ( ix, :, iquadpt) * obx1

      END IF

      IF ( igeom == 3) THEN

         sin_x2 = 1.0 / ( x1 * sin ( x2))

         h_quadpt ( ix, :, iquadpt) = h_quadpt ( ix, :, iquadpt)
     1                              * sin_x2

      END IF

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 7.b) Optional step needed only for viscous, resistive or conductive MHD.

! ASSUME that "du_quadpt_dx, du_quadpt_dy and du_quadpt_dz" have been built at
! all orders (NOT DONE YET).
! We set about building the derivatives of the primitive
! variables and various auxiliary variables at the nodal points. These are
! needed only in the presence of viscous, resistive or conduction terms.

! These derivatives are then used to build non-ideal fluxes at the quadrature
! points given by "fni_quadpt, gni_quadpt and hni_quadpt".

! Lastly, the non-ideal fluxes are added into the ideal fluxes

! ----------------------------------------------------------------------

      IF ( ( ( viscosity_coef .GT. 0.0) .OR.
     1       ( resistivity_coef .GT. 0.0)) .AND. ( NUMORDER >= 3)) THEN

! ------------------------------

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = 1.0 / rho ( ix, iquadpt)

      drho_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 1, iquadpt)
      drho_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 1, iquadpt)
      drho_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 1, iquadpt)


      dvlx_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvlx_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvlx_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dvly_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvly_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvly_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dvlz_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvlz_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvlz_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dbfx_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 6, iquadpt)
      dbfx_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 6, iquadpt)
      dbfx_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 6, iquadpt)

      dbfy_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 7, iquadpt)
      dbfy_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 7, iquadpt)
      dbfy_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 7, iquadpt)

      dbfz_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 8, iquadpt)
      dbfz_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 8, iquadpt)
      dbfz_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 8, iquadpt)


      dv2_dx ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dx ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dx ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dx ( ix, iquadpt) )

      dv2_dy ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dy ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dy ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dy ( ix, iquadpt) )

      dv2_dz ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dz ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dz ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dz ( ix, iquadpt) )

      db2_dx ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dx ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dx ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dx ( ix, iquadpt) )

      db2_dy ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dy ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dy ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dy ( ix, iquadpt) )

      db2_dz ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dz ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dz ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dz ( ix, iquadpt) )


      tempac = vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1       + vlz ( ix, iquadpt)**2

      tempab = molewt_quadpt ( ix, 1, iquadpt)
     1       / ( gasconst * rho ( ix, iquadpt) )

      temperature ( ix, iquadpt) = prs ( ix, iquadpt) * tempab

      tempab = tempab * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      dtemperature_dx ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dx ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dx ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dx ( ix, iquadpt)
     1    - ob8pi * db2_dx ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dx ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

      dtemperature_dy ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dy ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dy ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dy ( ix, iquadpt)
     1    - ob8pi * db2_dy ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dy ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

      dtemperature_dz ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dz ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dz ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dz ( ix, iquadpt)
     1    - ob8pi * db2_dz ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dz ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

! The derivatives are within the reference element's local coordinates.
! To get true physical fluxes, we need to divide by the size of the zone here.

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = 1.0 / dx ( ix)
      tempab = 1.0 / dy ( ix)
      tempac = 1.0 / dz ( ix)

      drho_dx ( ix, iquadpt) = tempaa * drho_dx ( ix, iquadpt)
      dtemperature_dx ( ix, iquadpt) =
     1                 tempaa * dtemperature_dx ( ix, iquadpt)
      dvlx_dx ( ix, iquadpt) = tempaa * dvlx_dx ( ix, iquadpt)
      dvly_dx ( ix, iquadpt) = tempaa * dvly_dx ( ix, iquadpt)
      dvlz_dx ( ix, iquadpt) = tempaa * dvlz_dx ( ix, iquadpt)
      dbfx_dx ( ix, iquadpt) = tempaa * dbfx_dx ( ix, iquadpt)
      dbfy_dx ( ix, iquadpt) = tempaa * dbfy_dx ( ix, iquadpt)
      dbfz_dx ( ix, iquadpt) = tempaa * dbfz_dx ( ix, iquadpt)
      dv2_dx ( ix, iquadpt) = tempaa * dv2_dx ( ix, iquadpt)
      db2_dx ( ix, iquadpt) = tempaa * db2_dx ( ix, iquadpt)

      drho_dy ( ix, iquadpt) = tempab * drho_dy ( ix, iquadpt)
      dtemperature_dy ( ix, iquadpt) =
     1                 tempab * dtemperature_dy ( ix, iquadpt)
      dvlx_dy ( ix, iquadpt) = tempab * dvlx_dy ( ix, iquadpt)
      dvly_dy ( ix, iquadpt) = tempab * dvly_dy ( ix, iquadpt)
      dvlz_dy ( ix, iquadpt) = tempab * dvlz_dy ( ix, iquadpt)
      dbfx_dy ( ix, iquadpt) = tempab * dbfx_dy ( ix, iquadpt)
      dbfy_dy ( ix, iquadpt) = tempab * dbfy_dy ( ix, iquadpt)
      dbfz_dy ( ix, iquadpt) = tempab * dbfz_dy ( ix, iquadpt)
      dv2_dy ( ix, iquadpt) = tempab * dv2_dy ( ix, iquadpt)
      db2_dy ( ix, iquadpt) = tempab * db2_dy ( ix, iquadpt)

      drho_dz ( ix, iquadpt) = tempac * drho_dz ( ix, iquadpt)
      dtemperature_dz ( ix, iquadpt) =
     1                 tempac * dtemperature_dz ( ix, iquadpt)
      dvlx_dz ( ix, iquadpt) = tempac * dvlx_dz ( ix, iquadpt)
      dvly_dz ( ix, iquadpt) = tempac * dvly_dz ( ix, iquadpt)
      dvlz_dz ( ix, iquadpt) = tempac * dvlz_dz ( ix, iquadpt)
      dbfx_dz ( ix, iquadpt) = tempac * dbfx_dz ( ix, iquadpt)
      dbfy_dz ( ix, iquadpt) = tempac * dbfy_dz ( ix, iquadpt)
      dbfz_dz ( ix, iquadpt) = tempac * dbfz_dz ( ix, iquadpt)
      dv2_dz ( ix, iquadpt) = tempac * dv2_dz ( ix, iquadpt)
      db2_dz ( ix, iquadpt) = tempac * db2_dz ( ix, iquadpt)

      div_vl ( ix, iquadpt) = dvlx_dx ( ix, iquadpt)
     1                      + dvly_dy ( ix, iquadpt)
     1                      + dvlz_dz ( ix, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = rho ( ix, iquadpt) * viscosity_coef

! ----------

      fni_quadpt ( ix, 1, iquadpt) = 0.0

      fni_quadpt ( ix, 2, iquadpt) =
     1               2.0 * tempaa * ( dvlx_dx ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      fni_quadpt ( ix, 3, iquadpt) =
     1                     tempaa * ( dvlx_dy ( ix, iquadpt)
     1                              + dvly_dx ( ix, iquadpt) )

      fni_quadpt ( ix, 4, iquadpt) =
     1                     tempaa * ( dvlx_dz ( ix, iquadpt)
     1                              + dvlz_dx ( ix, iquadpt) )

      fni_quadpt ( ix, 5, iquadpt) =
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vlx ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvlx_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvlx_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvlx_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dx ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dx ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfx_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfx_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfx_dz ( ix, iquadpt) )

      fni_quadpt ( ix, 6, iquadpt) = 0.0

      fni_quadpt ( ix, 7, iquadpt) =
     1   resistivity_coef * ( dbfy_dx ( ix, iquadpt)
     1                      - dbfx_dy ( ix, iquadpt) )

      fni_quadpt ( ix, 8, iquadpt) =
     1   resistivity_coef * ( dbfz_dx ( ix, iquadpt)
     1                      - dbfx_dz ( ix, iquadpt) )

#if ( NVAR > 8)
      fni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      gni_quadpt ( ix, 1, iquadpt) = 0.0

      gni_quadpt ( ix, 2, iquadpt) =
     1                     tempaa * ( dvly_dx ( ix, iquadpt)
     1                              + dvlx_dy ( ix, iquadpt) )

      gni_quadpt ( ix, 3, iquadpt) =
     1               2.0 * tempaa * ( dvly_dy ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      gni_quadpt ( ix, 4, iquadpt) =
     1                     tempaa * ( dvly_dz ( ix, iquadpt)
     1                              + dvlz_dy ( ix, iquadpt) )

      gni_quadpt ( ix, 5, iquadpt) = 
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vly ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvly_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvly_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvly_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dy ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dy ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfy_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfy_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfy_dz ( ix, iquadpt) )

      gni_quadpt ( ix, 6, iquadpt) =
     1   resistivity_coef * ( dbfx_dy ( ix, iquadpt)
     1                      - dbfy_dx ( ix, iquadpt) )

      gni_quadpt ( ix, 7, iquadpt) = 0.0

      gni_quadpt ( ix, 8, iquadpt) =
     1   resistivity_coef * ( dbfz_dy ( ix, iquadpt)
     1                      - dbfy_dz ( ix, iquadpt) )

#if ( NVAR > 8)
      gni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      hni_quadpt ( ix, 1, iquadpt) = 0.0

      hni_quadpt ( ix, 2, iquadpt) =
     1                     tempaa * ( dvlz_dx ( ix, iquadpt)
     1                              + dvlx_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 3, iquadpt) =
     1                     tempaa * ( dvlz_dy ( ix, iquadpt)
     1                              + dvly_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 4, iquadpt) =
     1               2.0 * tempaa * ( dvlz_dz ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      hni_quadpt ( ix, 5, iquadpt) =
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vlz ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvlz_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvlz_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvlz_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dz ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dz ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfz_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfz_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfz_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 6, iquadpt) =
     1   resistivity_coef * ( dbfx_dz ( ix, iquadpt)
     1                      - dbfz_dx ( ix, iquadpt) )

      hni_quadpt ( ix, 7, iquadpt) =
     1   resistivity_coef * ( dbfy_dz ( ix, iquadpt)
     1                      - dbfz_dy ( ix, iquadpt) )

      hni_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NVAR > 8)
      hni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      f_quadpt ( ix, :, iquadpt) = f_quadpt ( ix, :, iquadpt)
     1                       + fni_quadpt ( ix, :, iquadpt)

      g_quadpt ( ix, :, iquadpt) = g_quadpt ( ix, :, iquadpt)
     1                       + gni_quadpt ( ix, :, iquadpt)

      h_quadpt ( ix, :, iquadpt) = h_quadpt ( ix, :, iquadpt)
     1                       + hni_quadpt ( ix, :, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      END IF
! End of "( viscosity_coef .GT. 0.0) .OR. ( resistivity_coef .GT. 0.0)".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = dtmove / dx ( ix)
      f_quadpt ( ix, :, 1: NUM_SPACETIME_QUADPTS) = tempaa
     1         * f_quadpt ( ix, :, 1: NUM_SPACETIME_QUADPTS)

      tempab = dtmove / dy ( ix)
      g_quadpt ( ix, :, 1: NUM_SPACETIME_QUADPTS) = tempab
     1         * g_quadpt ( ix, :, 1: NUM_SPACETIME_QUADPTS)

      tempac = dtmove / dz ( ix)
      h_quadpt ( ix, :, 1: NUM_SPACETIME_QUADPTS) = tempac
     1         * h_quadpt ( ix, :, 1: NUM_SPACETIME_QUADPTS)

! ----------

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 7.c) Now project the fluxes "( f_quadpt, g_quadpt, h_quadpt)"
! from the space-time nodal points back to the modal
! space time variables for those variables, i.e. to "( f, g, h)" .

! --------------------
#if ( NUMORDER == 2)
! dg, 2nd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

!      f ( ix, :, 01) = ( f_quadpt ( ix, :, 01) + f_quadpt ( ix, :, 02)
!     1         + f_quadpt ( ix, :, 03) + f_quadpt ( ix, :, 04)
!     1         + f_quadpt ( ix, :, 05) + f_quadpt ( ix, :, 06) ) / 6.0

      f ( ix, :, 02) = f_quadpt ( ix, :, 01) - f_quadpt ( ix, :, 02)

!      f ( ix, :, 03) = f_quadpt ( ix, :, 03) - f_quadpt ( ix, :, 04)

!      f ( ix, :, 04) = f_quadpt ( ix, :, 05) - f_quadpt ( ix, :, 06)


!      f ( ix, :, 05) = f_quadpt ( ix, :, 07) - f ( ix, :, 01)

! ----------

!      g ( ix, :, 01) = ( g_quadpt ( ix, :, 01) + g_quadpt ( ix, :, 02)
!     1         + g_quadpt ( ix, :, 03) + g_quadpt ( ix, :, 04)
!     1         + g_quadpt ( ix, :, 05) + g_quadpt ( ix, :, 06) ) / 6.0

!      g ( ix, :, 02) = g_quadpt ( ix, :, 01) - g_quadpt ( ix, :, 02)

      g ( ix, :, 03) = g_quadpt ( ix, :, 03) - g_quadpt ( ix, :, 04)

!      g ( ix, :, 04) = g_quadpt ( ix, :, 05) - g_quadpt ( ix, :, 06)


!      g ( ix, :, 05) = g_quadpt ( ix, :, 07) - g ( ix, :, 01)

! ----------

!      h ( ix, :, 01) = ( h_quadpt ( ix, :, 01) + h_quadpt ( ix, :, 02)
!     1         + h_quadpt ( ix, :, 03) + h_quadpt ( ix, :, 04)
!     1         + h_quadpt ( ix, :, 05) + h_quadpt ( ix, :, 06) ) / 6.0

!      h ( ix, :, 02) = h_quadpt ( ix, :, 01) - h_quadpt ( ix, :, 02)

!      h ( ix, :, 03) = h_quadpt ( ix, :, 03) - h_quadpt ( ix, :, 04)

      h ( ix, :, 04) = h_quadpt ( ix, :, 05) - h_quadpt ( ix, :, 06)


!      h ( ix, :, 05) = h_quadpt ( ix, :, 07) - h ( ix, :, 01)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
! --------------------
#if ( NUMORDER == 3)
! dg, 3rd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tf ( ix, :, 01) = f_quadpt ( ix, :, 02) + f_quadpt ( ix, :, 03)
     1                + f_quadpt ( ix, :, 04) + f_quadpt ( ix, :, 05)

      tf ( ix, :, 02) = f_quadpt ( ix, :, 06) + f_quadpt ( ix, :, 07)
     1                + f_quadpt ( ix, :, 08) + f_quadpt ( ix, :, 09)

      tf ( ix, :, 03) = f_quadpt ( ix, :, 10) + f_quadpt ( ix, :, 11)
     1                + f_quadpt ( ix, :, 12) + f_quadpt ( ix, :, 13)

      tf ( ix, :, 06) = f_quadpt ( ix, :, 04) + f_quadpt ( ix, :, 05)
     1                + f_quadpt ( ix, :, 08) + f_quadpt ( ix, :, 09)

      tf ( ix, :, 09) = f_quadpt ( ix, :, 02) + f_quadpt ( ix, :, 03)
     1                + f_quadpt ( ix, :, 06) + f_quadpt ( ix, :, 07)


      f ( ix, :, 02) = 0.125 * ( tf ( ix, :, 01) - tf ( ix, :, 02) )
     1 + 0.5 * ( f_quadpt ( ix, :, 14) - f_quadpt ( ix, :, 15) )

      f ( ix, :, 05) = 0.5 * ( - tf ( ix, :, 03) + tf ( ix, :, 06)
     1 + tf ( ix, :, 09) ) - 2.0 * f_quadpt ( ix, :, 01)

      f ( ix, :, 08) = f_quadpt ( ix, :, 02) - f_quadpt ( ix, :, 06)
     1 - f_quadpt ( ix, :, 03) + f_quadpt ( ix, :, 07)

      f ( ix, :, 10) = f_quadpt ( ix, :, 04) - f_quadpt ( ix, :, 05)
     1 - f_quadpt ( ix, :, 08) + f_quadpt ( ix, :, 09)

      f ( ix, :, 13) = 1.5 * ( f_quadpt ( ix, :, 14)
     1                       - f_quadpt ( ix, :, 15) )
     1 - 0.375 * ( tf ( ix, :, 01) - tf ( ix, :, 02) )

! ----------

      tg ( ix, :, 03) = g_quadpt ( ix, :, 10) + g_quadpt ( ix, :, 11)
     1                + g_quadpt ( ix, :, 12) + g_quadpt ( ix, :, 13)

      tg ( ix, :, 04) = g_quadpt ( ix, :, 02) + g_quadpt ( ix, :, 06)
     1                + g_quadpt ( ix, :, 10) + g_quadpt ( ix, :, 12)

      tg ( ix, :, 05) = g_quadpt ( ix, :, 03) + g_quadpt ( ix, :, 07)
     1                + g_quadpt ( ix, :, 11) + g_quadpt ( ix, :, 13)

      tg ( ix, :, 06) = g_quadpt ( ix, :, 04) + g_quadpt ( ix, :, 05)
     1                + g_quadpt ( ix, :, 08) + g_quadpt ( ix, :, 09)

      tg ( ix, :, 09) = g_quadpt ( ix, :, 02) + g_quadpt ( ix, :, 03)
     1                + g_quadpt ( ix, :, 06) + g_quadpt ( ix, :, 07)


      g ( ix, :, 03) = 0.125 * ( tg ( ix, :, 04) - tg ( ix, :, 05) )
     1 + 0.5 * ( g_quadpt ( ix, :, 16) - g_quadpt ( ix, :, 17) )

      g ( ix, :, 06) = 0.5 * ( tg ( ix, :, 03) - tg ( ix, :, 06)
     1 + tg ( ix, :, 09) ) - 2.0 * g_quadpt ( ix, :, 01)

      g ( ix, :, 08) = g_quadpt ( ix, :, 02) - g_quadpt ( ix, :, 06)
     1 - g_quadpt ( ix, :, 03) + g_quadpt ( ix, :, 07)

      g ( ix, :, 09) = g_quadpt ( ix, :, 10) - g_quadpt ( ix, :, 11)
     1 - g_quadpt ( ix, :, 12) + g_quadpt ( ix, :, 13)

      g ( ix, :, 14) = 1.5 * ( g_quadpt ( ix, :, 16)
     1 - g_quadpt ( ix, :, 17) ) - 0.375 * ( tg ( ix, :, 04)
     1 - tg ( ix, :, 05) )


! ----------

      th ( ix, :, 03) = h_quadpt ( ix, :, 10) + h_quadpt ( ix, :, 11)
     1                + h_quadpt ( ix, :, 12) + h_quadpt ( ix, :, 13)

      th ( ix, :, 06) = h_quadpt ( ix, :, 04) + h_quadpt ( ix, :, 05)
     1                + h_quadpt ( ix, :, 08) + h_quadpt ( ix, :, 09)

      th ( ix, :, 07) = h_quadpt ( ix, :, 04) + h_quadpt ( ix, :, 08)
     1                + h_quadpt ( ix, :, 10) + h_quadpt ( ix, :, 11)

      th ( ix, :, 08) = h_quadpt ( ix, :, 05) + h_quadpt ( ix, :, 09)
     1                + h_quadpt ( ix, :, 12) + h_quadpt ( ix, :, 13)

      th ( ix, :, 09) = h_quadpt ( ix, :, 02) + h_quadpt ( ix, :, 03)
     1                + h_quadpt ( ix, :, 06) + h_quadpt ( ix, :, 07)


      h ( ix, :, 04) = 0.125 * ( th ( ix, :, 07) - th ( ix, :, 08) )
     1 + 0.5 * ( h_quadpt ( ix, :, 18) - h_quadpt ( ix, :, 19) )

      h ( ix, :, 07) = 0.5 * ( th ( ix, :, 03) + th ( ix, :, 06)
     1 - th ( ix, :, 09) ) - 2.0 * h_quadpt ( ix, :, 01)

      h ( ix, :, 09) = h_quadpt ( ix, :, 10) - h_quadpt ( ix, :, 11)
     1 - h_quadpt ( ix, :, 12) + h_quadpt ( ix, :, 13)

      h ( ix, :, 10) = h_quadpt ( ix, :, 04) - h_quadpt ( ix, :, 05)
     1 - h_quadpt ( ix, :, 08) + h_quadpt ( ix, :, 09)

      h ( ix, :, 15) = 1.5 * ( h_quadpt ( ix, :, 18)
     1 - h_quadpt ( ix, :, 19) ) - 0.375 * ( th ( ix, :, 07)
     1 - th ( ix, :, 08) )

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
! --------------------
#if ( NUMORDER == 4)
! dg, 4th order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tf ( ix, :, 01) = f_quadpt ( ix, :, 02)
     1 - 2.0 * f_quadpt ( ix, :, 01) + f_quadpt ( ix, :, 05)

      tf ( ix, :, 02) = f_quadpt ( ix, :, 36)
     1 - 2.0 * f_quadpt ( ix, :, 35) + f_quadpt ( ix, :, 37)

      tf ( ix, :, 07) = f_quadpt ( ix, :, 14) + f_quadpt ( ix, :, 15)
     1 + f_quadpt ( ix, :, 16) + f_quadpt ( ix, :, 17)

      tf ( ix, :, 08) = f_quadpt ( ix, :, 18) + f_quadpt ( ix, :, 19)
     1 + f_quadpt ( ix, :, 20) + f_quadpt ( ix, :, 21)

      tf ( ix, :, 09) = f_quadpt ( ix, :, 14) + f_quadpt ( ix, :, 16)
     1 + f_quadpt ( ix, :, 18) + f_quadpt ( ix, :, 20)

      tf ( ix, :, 10) = f_quadpt ( ix, :, 15) + f_quadpt ( ix, :, 17)
     1 + f_quadpt ( ix, :, 19) + f_quadpt ( ix, :, 21)

      tf ( ix, :, 11) = f_quadpt ( ix, :, 14) + f_quadpt ( ix, :, 15)
     1 + f_quadpt ( ix, :, 18) + f_quadpt ( ix, :, 19)

      tf ( ix, :, 12) = f_quadpt ( ix, :, 16) + f_quadpt ( ix, :, 17)
     1 + f_quadpt ( ix, :, 20) + f_quadpt ( ix, :, 21)

      tf ( ix, :, 13) = f_quadpt ( ix, :, 23) + f_quadpt ( ix, :, 24)
     1 - f_quadpt ( ix, :, 25) - f_quadpt ( ix, :, 26)
     1 - f_quadpt ( ix, :, 27) - f_quadpt ( ix, :, 28)
     1 + f_quadpt ( ix, :, 29) + f_quadpt ( ix, :, 30)

      tf ( ix, :, 14) = f_quadpt ( ix, :, 32) + f_quadpt ( ix, :, 34)
     1 - f_quadpt ( ix, :, 24) - f_quadpt ( ix, :, 28)
     1 - f_quadpt ( ix, :, 31) - f_quadpt ( ix, :, 33)
     1 + f_quadpt ( ix, :, 23) + f_quadpt ( ix, :, 27)

      tf ( ix, :, 15) = f_quadpt ( ix, :, 25) + f_quadpt ( ix, :, 29)
     1 - f_quadpt ( ix, :, 31) - f_quadpt ( ix, :, 32)
     1 - f_quadpt ( ix, :, 26) - f_quadpt ( ix, :, 30)
     1 + f_quadpt ( ix, :, 33) + f_quadpt ( ix, :, 34)

      tf ( ix, :, 16) = 0.5 * ( f_quadpt ( ix, :, 23)
     1 + f_quadpt ( ix, :, 24) - f_quadpt ( ix, :, 27)
     1 - f_quadpt ( ix, :, 28) )


      f ( ix, :, 05) = tf ( ix, :, 01) + tf ( ix, :, 02)

      f ( ix, :, 08) = f_quadpt ( ix, :, 23) - f_quadpt ( ix, :, 27)
     1 - f_quadpt ( ix, :, 24) + f_quadpt ( ix, :, 28)

      f ( ix, :, 10) = f_quadpt ( ix, :, 25) - f_quadpt ( ix, :, 26)
     1 - f_quadpt ( ix, :, 29) + f_quadpt ( ix, :, 30)

      f ( ix, :, 11) = 16 * ( f_quadpt ( ix, :, 02)
     1 - 2.0 * f_quadpt ( ix, :, 03) + 2.0 * f_quadpt ( ix, :, 04)
     1 - f_quadpt ( ix, :, 05) ) / 3.0

      f ( ix, :, 16) = tf ( ix, :, 13) + 0.5 * ( tf ( ix, :, 07)
     1 - tf ( ix, :, 08) ) - 2.0 * ( f_quadpt ( ix, :, 02)
     1 - f_quadpt ( ix, :, 05) )

      f ( ix, :, 18) = - 2.0 * tf ( ix, :, 13) + f ( ix, :, 16)

      f ( ix, :, 14) = tf ( ix, :, 14) + 0.5 * ( tf ( ix, :, 09)
     1 - tf ( ix, :, 10) ) - 2.0 * ( f_quadpt ( ix, :, 06)
     1 - f_quadpt ( ix, :, 09) )

      f ( ix, :, 19) = - 2.0 * tf ( ix, :, 14) + f ( ix, :, 14)

      f ( ix, :, 15) = tf ( ix, :, 15) + 0.5 * ( tf ( ix, :, 11)
     1 - tf ( ix, :, 12) ) - 2.0 * ( f_quadpt ( ix, :, 10)
     1 - f_quadpt ( ix, :, 13) )

      f ( ix, :, 17) = - 2.0 * tf ( ix, :, 15) + f ( ix, :, 15)

      f ( ix, :, 20) = f_quadpt ( ix, :, 14) - f_quadpt ( ix, :, 18)
     1 - f_quadpt ( ix, :, 15) + f_quadpt ( ix, :, 19)
     1 - f_quadpt ( ix, :, 16) + f_quadpt ( ix, :, 17)
     1 + f_quadpt ( ix, :, 20) - f_quadpt ( ix, :, 21)

      f ( ix, :, 24) = - f_quadpt ( ix, :, 02) + f_quadpt ( ix, :, 05)
     1 + f_quadpt ( ix, :, 36) - f_quadpt ( ix, :, 37)

      f ( ix, :, 02) = ( 4.0 * f_quadpt ( ix, :, 02)
     1 - 4.0 * f_quadpt ( ix, :, 05) + 2.0 * f ( ix, :, 24)
     1 - tf ( ix, :, 16) ) / 3.0 - 0.1 * f ( ix, :, 11)
     1 + ( 2.0 * f ( ix, :, 16) + f ( ix, :, 18) ) / 12.0

      f ( ix, :, 27) = - 2.0 * tf ( ix, :, 16) + 2.0 * f ( ix, :, 02)
     1 + 0.2 * f ( ix, :, 11) + ( 2.0 * f ( ix, :, 16)
     1 - f ( ix, :, 18) ) / 6.0

      f ( ix, :, 30) = 2.0 * ( tf ( ix, :, 02) - tf ( ix, :, 01) )

      f ( ix, :, 33) = 2.0 * f ( ix, :, 08) - f_quadpt ( ix, :, 14)
     1 + f_quadpt ( ix, :, 15) - f_quadpt ( ix, :, 16)
     1 + f_quadpt ( ix, :, 17) + f_quadpt ( ix, :, 18)
     1 - f_quadpt ( ix, :, 19) + f_quadpt ( ix, :, 20)
     1 - f_quadpt ( ix, :, 21)

      f ( ix, :, 35) = 2.0 * f ( ix, :, 10) - f_quadpt ( ix, :, 14)
     1 - f_quadpt ( ix, :, 15) + f_quadpt ( ix, :, 16)
     1 + f_quadpt ( ix, :, 17) + f_quadpt ( ix, :, 18)
     1 + f_quadpt ( ix, :, 19) - f_quadpt ( ix, :, 20)
     1 - f_quadpt ( ix, :, 21)

! ----------

      tg ( ix, :, 03) = g_quadpt ( ix, :, 06)
     1 - 2.0 * g_quadpt ( ix, :, 01) + g_quadpt ( ix, :, 09)

      tg ( ix, :, 04) = g_quadpt ( ix, :, 38)
     1 - 2.0 * g_quadpt ( ix, :, 35) + g_quadpt ( ix, :, 39)

      tg ( ix, :, 07) = g_quadpt ( ix, :, 14) + g_quadpt ( ix, :, 15)
     1 + g_quadpt ( ix, :, 16) + g_quadpt ( ix, :, 17)

      tg ( ix, :, 08) = g_quadpt ( ix, :, 18) + g_quadpt ( ix, :, 19)
     1 + g_quadpt ( ix, :, 20) + g_quadpt ( ix, :, 21)

      tg ( ix, :, 09) = g_quadpt ( ix, :, 14) + g_quadpt ( ix, :, 16)
     1 + g_quadpt ( ix, :, 18) + g_quadpt ( ix, :, 20)

      tg ( ix, :, 10) = g_quadpt ( ix, :, 15) + g_quadpt ( ix, :, 17)
     1 + g_quadpt ( ix, :, 19) + g_quadpt ( ix, :, 21)

      tg ( ix, :, 11) = g_quadpt ( ix, :, 14) + g_quadpt ( ix, :, 15)
     1 + g_quadpt ( ix, :, 18) + g_quadpt ( ix, :, 19)

      tg ( ix, :, 12) = g_quadpt ( ix, :, 16) + g_quadpt ( ix, :, 17)
     1 + g_quadpt ( ix, :, 20) + g_quadpt ( ix, :, 21)

      tg ( ix, :, 13) = g_quadpt ( ix, :, 23) + g_quadpt ( ix, :, 24)
     1 - g_quadpt ( ix, :, 25) - g_quadpt ( ix, :, 26)
     1 - g_quadpt ( ix, :, 27) - g_quadpt ( ix, :, 28)
     1 + g_quadpt ( ix, :, 29) + g_quadpt ( ix, :, 30)

      tg ( ix, :, 14) = g_quadpt ( ix, :, 32) + g_quadpt ( ix, :, 34)
     1 - g_quadpt ( ix, :, 24) - g_quadpt ( ix, :, 28)
     1 - g_quadpt ( ix, :, 31) - g_quadpt ( ix, :, 33)
     1 + g_quadpt ( ix, :, 23) + g_quadpt ( ix, :, 27)

      tg ( ix, :, 15) = g_quadpt ( ix, :, 25) + g_quadpt ( ix, :, 29)
     1 - g_quadpt ( ix, :, 31) - g_quadpt ( ix, :, 32)
     1 - g_quadpt ( ix, :, 26) - g_quadpt ( ix, :, 30)
     1 + g_quadpt ( ix, :, 33) + g_quadpt ( ix, :, 34)

      tg ( ix, :, 17) = 0.5 * ( g_quadpt ( ix, :, 23)
     1 + g_quadpt ( ix, :, 27) - g_quadpt ( ix, :, 24)
     1 - g_quadpt ( ix, :, 28) )


      g ( ix, :, 06) = tg ( ix, :, 03) + tg ( ix, :, 04)

      g ( ix, :, 08) = g_quadpt ( ix, :, 23) - g_quadpt ( ix, :, 27)
     1 - g_quadpt ( ix, :, 24) + g_quadpt ( ix, :, 28)

      g ( ix, :, 09) = g_quadpt ( ix, :, 31) - g_quadpt ( ix, :, 32)
     1 - g_quadpt ( ix, :, 33) + g_quadpt ( ix, :, 34)

      g ( ix, :, 12) = 16 * ( g_quadpt ( ix, :, 06)
     1 - 2.0 * g_quadpt ( ix, :, 07) + 2.0 * g_quadpt ( ix, :, 08)
     1 - g_quadpt ( ix, :, 09) ) / 3.0

      g ( ix, :, 16) = tg ( ix, :, 13) + 0.5 * ( tg ( ix, :, 07)
     1 - tg ( ix, :, 08) ) - 2.0 * ( g_quadpt ( ix, :, 02)
     1 - g_quadpt ( ix, :, 05) )

      g ( ix, :, 18) = - 2.0 * tg ( ix, :, 13) + g ( ix, :, 16)

      g ( ix, :, 14) = tg ( ix, :, 14) + 0.5 * ( tg ( ix, :, 09)
     1 - tg ( ix, :, 10) ) - 2.0 * ( g_quadpt ( ix, :, 06)
     1 - g_quadpt ( ix, :, 09) )

      g ( ix, :, 19) = - 2.0 * tg ( ix, :, 14) + g ( ix, :, 14)

      g ( ix, :, 15) = tg ( ix, :, 15) + 0.5 * ( tg ( ix, :, 11)
     1 - tg ( ix, :, 12) ) - 2.0 * ( g_quadpt ( ix, :, 10)
     1 - g_quadpt ( ix, :, 13) )

      g ( ix, :, 17) = - 2.0 * tg ( ix, :, 15) + g ( ix, :, 15)

      g ( ix, :, 20) = g_quadpt ( ix, :, 14) - g_quadpt ( ix, :, 18)
     1 - g_quadpt ( ix, :, 15) + g_quadpt ( ix, :, 19)
     1 - g_quadpt ( ix, :, 16) + g_quadpt ( ix, :, 17)
     1 + g_quadpt ( ix, :, 20) - g_quadpt ( ix, :, 21)

      g ( ix, :, 25) = - g_quadpt ( ix, :, 06) + g_quadpt ( ix, :, 09)
     1 + g_quadpt ( ix, :, 38) - g_quadpt ( ix, :, 39)

      g ( ix, :, 03) = ( 4.0 * g_quadpt ( ix, :, 06)
     1 - 4.0 * g_quadpt ( ix, :, 09) + 2.0 * g ( ix, :, 25)
     1 - tg ( ix, :, 17) ) / 3.0 - 0.1 * g ( ix, :, 12)
     1 + ( 2.0 * g ( ix, :, 14) + g ( ix, :, 19) ) / 12.0

      g ( ix, :, 28) = - 2.0 * tg ( ix, :, 17) + 2.0 * g ( ix, :, 03)
     1 + 0.2 * g ( ix, :, 12) + ( 2.0 * g ( ix, :, 14)
     1 - g ( ix, :, 19) ) / 6.0

      g ( ix, :, 31) = 2.0 * ( tg ( ix, :, 04) - tg ( ix, :, 03) )

      g ( ix, :, 33) = 2.0 * g ( ix, :, 08) - g_quadpt ( ix, :, 14)
     1 + g_quadpt ( ix, :, 15) - g_quadpt ( ix, :, 16)
     1 + g_quadpt ( ix, :, 17) + g_quadpt ( ix, :, 18)
     1 - g_quadpt ( ix, :, 19) + g_quadpt ( ix, :, 20)
     1 - g_quadpt ( ix, :, 21)

      g ( ix, :, 34) = 2.0 * g ( ix, :, 09) - g_quadpt ( ix, :, 14)
     1 + g_quadpt ( ix, :, 15) + g_quadpt ( ix, :, 16)
     1 - g_quadpt ( ix, :, 17) - g_quadpt ( ix, :, 18)
     1 + g_quadpt ( ix, :, 19) + g_quadpt ( ix, :, 20)
     1 - g_quadpt ( ix, :, 21)

! ----------

      th ( ix, :, 05) = h_quadpt ( ix, :, 10)
     1 - 2.0 * h_quadpt ( ix, :, 01) + h_quadpt ( ix, :, 13)

      th ( ix, :, 06) = h_quadpt ( ix, :, 40)
     1 - 2.0 * h_quadpt ( ix, :, 35) + h_quadpt ( ix, :, 41)

      th ( ix, :, 07) = h_quadpt ( ix, :, 14) + h_quadpt ( ix, :, 15)
     1 + h_quadpt ( ix, :, 16) + h_quadpt ( ix, :, 17)

      th ( ix, :, 08) = h_quadpt ( ix, :, 18) + h_quadpt ( ix, :, 19)
     1 + h_quadpt ( ix, :, 20) + h_quadpt ( ix, :, 21)

      th ( ix, :, 09) = h_quadpt ( ix, :, 14) + h_quadpt ( ix, :, 16)
     1 + h_quadpt ( ix, :, 18) + h_quadpt ( ix, :, 20)

      th ( ix, :, 10) = h_quadpt ( ix, :, 15) + h_quadpt ( ix, :, 17)
     1 + h_quadpt ( ix, :, 19) + h_quadpt ( ix, :, 21)

      th ( ix, :, 11) = h_quadpt ( ix, :, 14) + h_quadpt ( ix, :, 15)
     1 + h_quadpt ( ix, :, 18) + h_quadpt ( ix, :, 19)

      th ( ix, :, 12) = h_quadpt ( ix, :, 16) + h_quadpt ( ix, :, 17)
     1 + h_quadpt ( ix, :, 20) + h_quadpt ( ix, :, 21)

      th ( ix, :, 13) = h_quadpt ( ix, :, 23) + h_quadpt ( ix, :, 24)
     1 - h_quadpt ( ix, :, 25) - h_quadpt ( ix, :, 26)
     1 - h_quadpt ( ix, :, 27) - h_quadpt ( ix, :, 28)
     1 + h_quadpt ( ix, :, 29) + h_quadpt ( ix, :, 30)

      th ( ix, :, 14) = h_quadpt ( ix, :, 32) + h_quadpt ( ix, :, 34)
     1 - h_quadpt ( ix, :, 24) - h_quadpt ( ix, :, 28)
     1 - h_quadpt ( ix, :, 31) - h_quadpt ( ix, :, 33)
     1 + h_quadpt ( ix, :, 23) + h_quadpt ( ix, :, 27)

      th ( ix, :, 15) = h_quadpt ( ix, :, 25) + h_quadpt ( ix, :, 29)
     1 - h_quadpt ( ix, :, 31) - h_quadpt ( ix, :, 32)
     1 - h_quadpt ( ix, :, 26) - h_quadpt ( ix, :, 30)
     1 + h_quadpt ( ix, :, 33) + h_quadpt ( ix, :, 34)

      th ( ix, :, 18) = 0.5 * ( h_quadpt ( ix, :, 31)
     1 + h_quadpt ( ix, :, 32) - h_quadpt ( ix, :, 33)
     1 - h_quadpt ( ix, :, 34) )


      h ( ix, :, 07) = th ( ix, :, 05) + th ( ix, :, 06)

      h ( ix, :, 09) = h_quadpt ( ix, :, 31) - h_quadpt ( ix, :, 32)
     1 - h_quadpt ( ix, :, 33) + h_quadpt ( ix, :, 34)

      h ( ix, :, 10) = h_quadpt ( ix, :, 25) - h_quadpt ( ix, :, 26)
     1 - h_quadpt ( ix, :, 29) + h_quadpt ( ix, :, 30)

      h ( ix, :, 13) = 16 * ( h_quadpt ( ix, :, 10)
     1 - 2.0 * h_quadpt ( ix, :, 11) + 2.0 * h_quadpt ( ix, :, 12)
     1 - h_quadpt ( ix, :, 13) ) / 3.0

      h ( ix, :, 16) = th ( ix, :, 13) + 0.5 * ( th ( ix, :, 07)
     1 - th ( ix, :, 08) ) - 2.0 * ( h_quadpt ( ix, :, 02)
     1 - h_quadpt ( ix, :, 05) )

      h ( ix, :, 18) = - 2.0 * th ( ix, :, 13) + h ( ix, :, 16)

      h ( ix, :, 14) = th ( ix, :, 14) + 0.5 * ( th ( ix, :, 09)
     1 - th ( ix, :, 10) ) - 2.0 * ( h_quadpt ( ix, :, 06)
     1 - h_quadpt ( ix, :, 09) )

      h ( ix, :, 19) = - 2.0 * th ( ix, :, 14) + h ( ix, :, 14)

      h ( ix, :, 15) = th ( ix, :, 15) + 0.5 * ( th ( ix, :, 11)
     1 - th ( ix, :, 12) ) - 2.0 * ( h_quadpt ( ix, :, 10)
     1 - h_quadpt ( ix, :, 13) )

      h ( ix, :, 17) = - 2.0 * th ( ix, :, 15) + h ( ix, :, 15)

      h ( ix, :, 20) = h_quadpt ( ix, :, 14) - h_quadpt ( ix, :, 18)
     1 - h_quadpt ( ix, :, 15) + h_quadpt ( ix, :, 19)
     1 - h_quadpt ( ix, :, 16) + h_quadpt ( ix, :, 17)
     1 + h_quadpt ( ix, :, 20) - h_quadpt ( ix, :, 21)

      h ( ix, :, 26) = - h_quadpt ( ix, :, 10) + h_quadpt ( ix, :, 13)
     1 + h_quadpt ( ix, :, 40) - h_quadpt ( ix, :, 41)

      h ( ix, :, 04) = ( 4.0 * h_quadpt ( ix, :, 10)
     1 - 4.0 * h_quadpt ( ix, :, 13) + 2.0 * h ( ix, :, 26)
     1 - th ( ix, :, 18) ) / 3.0 - 0.1 * h ( ix, :, 13)
     1 + ( h ( ix, :, 15) + 2.0 * h ( ix, :, 17)) / 12.0

      h ( ix, :, 29) = - 2.0 * th ( ix, :, 18) + 2.0 * h ( ix, :, 04)
     1 + 0.2 * h ( ix, :, 13) - ( h ( ix, :, 15)
     1 - 2.0 * h ( ix, :, 17)) / 6.0

      h ( ix, :, 32) = 2.0 * ( th ( ix, :, 06) - th ( ix, :, 05) )

      h ( ix, :, 34) = 2.0 * h ( ix, :, 09) - h_quadpt ( ix, :, 14)
     1 + h_quadpt ( ix, :, 15) + h_quadpt ( ix, :, 16)
     1 - h_quadpt ( ix, :, 17) - h_quadpt ( ix, :, 18)
     1 + h_quadpt ( ix, :, 19) + h_quadpt ( ix, :, 20)
     1 - h_quadpt ( ix, :, 21)

      h ( ix, :, 35) = 2.0 * h ( ix, :, 10) - h_quadpt ( ix, :, 14)
     1 - h_quadpt ( ix, :, 15) + h_quadpt ( ix, :, 16)
     1 + h_quadpt ( ix, :, 17) + h_quadpt ( ix, :, 18)
     1 + h_quadpt ( ix, :, 19) - h_quadpt ( ix, :, 20)
     1 - h_quadpt ( ix, :, 21)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! --------------------

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 8.a) Build source terms "s_quadpt" at the quadrature points.
! Remember that the primitive variables are still available from the previous
! loop at each of the quadrature points for constructing the source terms.

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      s_quadpt ( oned_begin - oned_off: oned_end + oned_off,
     1                                           :, iquadpt) = 0.0

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

#if ( 0 == 1)

! For single-fluid adiabatic MHD turbulence, we may want to provide a gradual
! (explicit) heating/cooling term which ensures that pressure always
! tends to a mean pressure "prs_bar" in a time "tau_bar" (comparable
! to a turnover time). Such a source term can be implemented explicitly
! and its implementation is given below.
! Remember to declare "prs_bar, tau_bar".

      prs_bar = 1.0
      tau_bar = 1.0

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 5, iquadpt) = - ( prs ( ix, iquadpt) - prs_bar)
     1     / ( tau_bar * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0) )

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

#endif

! ------------------------------

! Build two-fluid source terms. In this part, the source terms are
! hard-coded for one family of ionized fluid and one family of neutrals.

#if ( N_EULER == 1)

      IF ( ( igeom .EQ. 1) .AND. ( N_EULER .EQ. 1) .AND.
     1     ( ambipolar_diffusion_drag_coef .GT. 0.0)) THEN

      ivar = NFIELD + 1 + NFLUID

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 1, iquadpt) = 0.0

      tempaa = ambipolar_diffusion_drag_coef
     1       * rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

      s_quadpt ( ix, 2, iquadpt) = tempaa
     1           * ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt))

      s_quadpt ( ix, 3, iquadpt) = tempaa
     1           * ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt))

      s_quadpt ( ix, 4, iquadpt) = tempaa
     1           * ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt))

      s_quadpt ( ix, 5, iquadpt) = ( ambipolar_diffusion_drag_coef
     1   / ( molewt_ion + molewt_neut) )

     1 * ( rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

     1   * ( ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt) )
     1     * ( molewt_neut * vlx_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlx ( ix, iquadpt) )
     1     + ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt) )
     1     * ( molewt_neut * vly_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vly ( ix, iquadpt) )
     1     + ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt) )
     1     * ( molewt_neut * vlz_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlz ( ix, iquadpt) )
     1     )

     1   + 3.0 * ( molewt_neut * prs_euler ( ix, iquadpt, 1)
     1           * rho ( ix, iquadpt)
     1           - molewt_ion * prs ( ix, iquadpt)
     1           * rho_euler ( ix, iquadpt, 1) )
     1   )

      s_quadpt ( ix, 6, iquadpt) = 0.0
      s_quadpt ( ix, 7, iquadpt) = 0.0
      s_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      s_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) = 0.0
#endif


      s_quadpt ( ix, ivar + 1, iquadpt) = 0.0
      s_quadpt ( ix, ivar + 2, iquadpt) = - s_quadpt ( ix, 2, iquadpt)
      s_quadpt ( ix, ivar + 3, iquadpt) = - s_quadpt ( ix, 3, iquadpt)
      s_quadpt ( ix, ivar + 4, iquadpt) = - s_quadpt ( ix, 4, iquadpt)
      s_quadpt ( ix, ivar + 5, iquadpt) = - s_quadpt ( ix, 5, iquadpt)

#if ( NFLUID_EULER >= 1)
      s_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

      END IF
! End of "ambipolar_diffusion_drag_coef .GT. 0.0" part

#endif
! End of "N_EULER == 1"

! ------------------------------

! Cylindrical Geometry

      IF ( igeom .EQ. 2) THEN

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1                           * vlx ( ix, iquadpt) * obx1

      s_quadpt ( ix, 2, iquadpt) =
     1 ( rho ( ix, iquadpt) * vly ( ix, iquadpt)**2
     1 - rho ( ix, iquadpt) * vlx ( ix, iquadpt)**2
     1 + bfx ( ix, iquadpt)**2 * ob4pi
     1 - bfy ( ix, iquadpt)**2 * ob4pi) * obx1

      s_quadpt ( ix, 3, iquadpt) =
     1 - 2.0 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                              * vly ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                        * vlz ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 5, iquadpt) = - ( tempac * vlx ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = 0.0

      s_quadpt ( ix, 7, iquadpt) = 0.0

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Spherical Geometry

      IF ( igeom .EQ. 3) THEN

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      cos_x2 = cos ( x2)
      sin_x2 = sin ( x2)
      cot_x2 = cos_x2 / sin_x2

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 2.0 * vlx ( ix, iquadpt)
     1 + cot_x2 * vly ( ix, iquadpt) ) * obx1

      s_quadpt ( ix, 2, iquadpt) = rho ( ix, iquadpt)
     1 * ( - 2.0 * vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1   + vlz ( ix, iquadpt)**2) * obx1
     1 - cot_x2
     1 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1   - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1
     1 - ( - 2.0 * bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1   + bfz ( ix, iquadpt)**2) * ob4pi * obx1

      s_quadpt ( ix, 3, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1 + cot_x2
     1 * ( vly ( ix, iquadpt)**2 - vlz ( ix, iquadpt)**2)) * obx1
     1 + ( 3.0 * bfx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1   + cot_x2
     1   * ( bfy ( ix, iquadpt)**2 - bfz ( ix, iquadpt)**2))
     1 * ob4pi * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - rho ( ix, iquadpt) * vlz ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) + 2.0 * cot_x2
     1   * vly ( ix, iquadpt) ) * obx1
     1 + bfz ( ix, iquadpt)
     1 * ( 3.0 * bfx ( ix, iquadpt) + 2.0 * cot_x2
     1   * bfy ( ix, iquadpt) ) * ob4pi * obx1

      s_quadpt ( ix, 5, iquadpt) =
     1 - 2.0 * ( tempac * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab) * obx1
     1 - cot_x2 * ( tempac * vly ( ix, iquadpt)
     1 - bfy ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = - cot_x2
     1                           * g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 7, iquadpt) = g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF
! ------------------------------

! Multiply "s_quadpt" by "dtmove".

      DO ix = oned_begin - oned_off, oned_end + oned_off
      s_quadpt ( ix, :, 1: NUM_SPACETIME_QUADPTS) = dtmove
     1         * s_quadpt ( ix, :, 1: NUM_SPACETIME_QUADPTS)
      END DO

! ------------------------------

! STEP 8.b) Now project the sources "s_quadpt" from the space-time nodal
! points back to the modal space time variables for those variables,
! i.e. to "s".

! --------------------
#if ( NUMORDER == 2)
! dg, 2nd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s ( ix, :, 01) = ( s_quadpt ( ix, :, 01) + s_quadpt ( ix, :, 02)
     1 + s_quadpt ( ix, :, 03) + s_quadpt ( ix, :, 04)
     1 + s_quadpt ( ix, :, 05) + s_quadpt ( ix, :, 06) ) / 6.0

      s ( ix, :, 02) = s_quadpt ( ix, :, 01) - s_quadpt ( ix, :, 02)

      s ( ix, :, 03) = s_quadpt ( ix, :, 03) - s_quadpt ( ix, :, 04)

      s ( ix, :, 04) = s_quadpt ( ix, :, 05) - s_quadpt ( ix, :, 06)


      s ( ix, :, 05) = s_quadpt ( ix, :, 07) - s ( ix, :, 01)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
! --------------------
#if ( NUMORDER == 3)
! dg, 3rd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ts ( ix, :, 01) = s_quadpt ( ix, :, 02) + s_quadpt ( ix, :, 03)
     1 + s_quadpt ( ix, :, 04) + s_quadpt ( ix, :, 05)

      ts ( ix, :, 02) = s_quadpt ( ix, :, 06) + s_quadpt ( ix, :, 07)
     1 + s_quadpt ( ix, :, 08) + s_quadpt ( ix, :, 09)

      ts ( ix, :, 03) = s_quadpt ( ix, :, 10) + s_quadpt ( ix, :, 11)
     1 + s_quadpt ( ix, :, 12) + s_quadpt ( ix, :, 13)

      ts ( ix, :, 04) = s_quadpt ( ix, :, 02) + s_quadpt ( ix, :, 06)
     1 + s_quadpt ( ix, :, 10) + s_quadpt ( ix, :, 12)

      ts ( ix, :, 05) = s_quadpt ( ix, :, 03) + s_quadpt ( ix, :, 07)
     1 + s_quadpt ( ix, :, 11) + s_quadpt ( ix, :, 13)

      ts ( ix, :, 06) = s_quadpt ( ix, :, 04) + s_quadpt ( ix, :, 05)
     1 + s_quadpt ( ix, :, 08) + s_quadpt ( ix, :, 09)

      ts ( ix, :, 07) = s_quadpt ( ix, :, 04) + s_quadpt ( ix, :, 08)
     1 + s_quadpt ( ix, :, 10) + s_quadpt ( ix, :, 11)

      ts ( ix, :, 08) = s_quadpt ( ix, :, 05) + s_quadpt ( ix, :, 09)
     1 + s_quadpt ( ix, :, 12) + s_quadpt ( ix, :, 13)

      ts ( ix, :, 09) = s_quadpt ( ix, :, 02) + s_quadpt ( ix, :, 03)
     1 + s_quadpt ( ix, :, 06) + s_quadpt ( ix, :, 07)

      ts ( ix, :, 10) = ( s_quadpt ( ix, :, 14) + s_quadpt ( ix, :, 15)
     1 + s_quadpt ( ix, :, 16) + s_quadpt ( ix, :, 17)
     1 + s_quadpt ( ix, :, 18) + s_quadpt ( ix, :, 19) ) / 6.0


      s ( ix, :, 02) = 0.125 * ( ts ( ix, :, 01) - ts ( ix, :, 02) )
     1 + 0.5 * ( s_quadpt ( ix, :, 14) - s_quadpt ( ix, :, 15) )

      s ( ix, :, 03) = 0.125 * ( ts ( ix, :, 04) - ts ( ix, :, 05) )
     1 + 0.5 * ( s_quadpt ( ix, :, 16) - s_quadpt ( ix, :, 17) )

      s ( ix, :, 04) = 0.125 * ( ts ( ix, :, 07) - ts ( ix, :, 08) )
     1 + 0.5 * ( s_quadpt ( ix, :, 18) - s_quadpt ( ix, :, 19) )

      s ( ix, :, 05) = 0.5 * ( - ts ( ix, :, 03) + ts ( ix, :, 06)
     1 + ts ( ix, :, 09) ) - 2.0 * s_quadpt ( ix, :, 01)

      s ( ix, :, 06) = 0.5 * ( ts ( ix, :, 03) - ts ( ix, :, 06)
     1 + ts ( ix, :, 09) ) - 2.0 * s_quadpt ( ix, :, 01)

      s ( ix, :, 07) = 0.5 * ( ts ( ix, :, 03) + ts ( ix, :, 06)
     1 - ts ( ix, :, 09) ) - 2.0 * s_quadpt ( ix, :, 01)

      s ( ix, :, 08) = s_quadpt ( ix, :, 02) - s_quadpt ( ix, :, 06)
     1 - s_quadpt ( ix, :, 03) + s_quadpt ( ix, :, 07)

      s ( ix, :, 09) = s_quadpt ( ix, :, 10) - s_quadpt ( ix, :, 11)
     1 - s_quadpt ( ix, :, 12) + s_quadpt ( ix, :, 13)

      s ( ix, :, 10) = s_quadpt ( ix, :, 04) - s_quadpt ( ix, :, 05)
     1 - s_quadpt ( ix, :, 08) + s_quadpt ( ix, :, 09)

      s ( ix, :, 13) = 1.5 * ( s_quadpt ( ix, :, 14)
     1 - s_quadpt ( ix, :, 15) ) - 0.375 * ( ts ( ix, :, 01)
     1 - ts ( ix, :, 02) )

      s ( ix, :, 14) = 1.5 * ( s_quadpt ( ix, :, 16)
     1 - s_quadpt ( ix, :, 17) ) - 0.375 * ( ts ( ix, :, 04)
     1 - ts ( ix, :, 05) )

      s ( ix, :, 15) = 1.5 * ( s_quadpt ( ix, :, 18)
     1 - s_quadpt ( ix, :, 19) ) - 0.375 * ( ts ( ix, :, 07)
     1 - ts ( ix, :, 08) )

      ts ( ix, :, 11) = 0.125 * ( s ( ix, :, 05) + s ( ix, :, 06)
     1 + s ( ix, :, 07) )

      s ( ix, :, 11) = 1.5 * ( ts ( ix, :, 10)
     1 - s_quadpt ( ix, :, 01) ) - ts ( ix, :, 11)

      s ( ix, :, 12) = 0.75 * ( s_quadpt ( ix, :, 20)
     1 - 2.0 * ts ( ix, :, 10) + s_quadpt ( ix, :, 01) )
     1 + ts ( ix, :, 11)

      s ( ix, :, 01) = ts ( ix, :, 10) - ( s ( ix, :, 11)
     1 - s ( ix, :, 12) ) / 3.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
! --------------------
#if ( NUMORDER == 4)
! dg, 4th order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ts ( ix, :, 01) = s_quadpt ( ix, :, 02)
     1 - 2.0 * s_quadpt ( ix, :, 01) + s_quadpt ( ix, :, 05)

      ts ( ix, :, 02) = s_quadpt ( ix, :, 36)
     1 - 2.0 * s_quadpt ( ix, :, 35) + s_quadpt ( ix, :, 37)

      ts ( ix, :, 03) = s_quadpt ( ix, :, 06)
     1 - 2.0 * s_quadpt ( ix, :, 01) + s_quadpt ( ix, :, 09)

      ts ( ix, :, 04) = s_quadpt ( ix, :, 38)
     1 - 2.0 * s_quadpt ( ix, :, 35) + s_quadpt ( ix, :, 39)

      ts ( ix, :, 05) = s_quadpt ( ix, :, 10)
     1 - 2.0 * s_quadpt ( ix, :, 01) + s_quadpt ( ix, :, 13)

      ts ( ix, :, 06) = s_quadpt ( ix, :, 40)
     1 - 2.0 * s_quadpt ( ix, :, 35) + s_quadpt ( ix, :, 41)


      ts ( ix, :, 07) = s_quadpt ( ix, :, 14) + s_quadpt ( ix, :, 15)
     1 + s_quadpt ( ix, :, 16) + s_quadpt ( ix, :, 17)

      ts ( ix, :, 08) = s_quadpt ( ix, :, 18) + s_quadpt ( ix, :, 19)
     1 + s_quadpt ( ix, :, 20) + s_quadpt ( ix, :, 21)

      ts ( ix, :, 09) = s_quadpt ( ix, :, 14) + s_quadpt ( ix, :, 16)
     1 + s_quadpt ( ix, :, 18) + s_quadpt ( ix, :, 20)

      ts ( ix, :, 10) = s_quadpt ( ix, :, 15) + s_quadpt ( ix, :, 17)
     1 + s_quadpt ( ix, :, 19) + s_quadpt ( ix, :, 21)

      ts ( ix, :, 11) = s_quadpt ( ix, :, 14) + s_quadpt ( ix, :, 15)
     1 + s_quadpt ( ix, :, 18) + s_quadpt ( ix, :, 19)

      ts ( ix, :, 12) = s_quadpt ( ix, :, 16) + s_quadpt ( ix, :, 17)
     1 + s_quadpt ( ix, :, 20) + s_quadpt ( ix, :, 21)

      ts ( ix, :, 13) = s_quadpt ( ix, :, 23) + s_quadpt ( ix, :, 24)
     1 - s_quadpt ( ix, :, 25) - s_quadpt ( ix, :, 26)
     1 - s_quadpt ( ix, :, 27) - s_quadpt ( ix, :, 28)
     1 + s_quadpt ( ix, :, 29) + s_quadpt ( ix, :, 30)

      ts ( ix, :, 14) = s_quadpt ( ix, :, 32) + s_quadpt ( ix, :, 34)
     1 - s_quadpt ( ix, :, 24) - s_quadpt ( ix, :, 28)
     1 - s_quadpt ( ix, :, 31) - s_quadpt ( ix, :, 33)
     1 + s_quadpt ( ix, :, 23) + s_quadpt ( ix, :, 27)

      ts ( ix, :, 15) = s_quadpt ( ix, :, 25) + s_quadpt ( ix, :, 29)
     1 - s_quadpt ( ix, :, 31) - s_quadpt ( ix, :, 32)
     1 - s_quadpt ( ix, :, 26) - s_quadpt ( ix, :, 30)
     1 + s_quadpt ( ix, :, 33) + s_quadpt ( ix, :, 34)


      ts ( ix, :, 16) = 0.5 * ( s_quadpt ( ix, :, 23)
     1 + s_quadpt ( ix, :, 24) - s_quadpt ( ix, :, 27)
     1 - s_quadpt ( ix, :, 28) )

      ts ( ix, :, 17) = 0.5 * ( s_quadpt ( ix, :, 23)
     1 + s_quadpt ( ix, :, 27) - s_quadpt ( ix, :, 24)
     1 - s_quadpt ( ix, :, 28) ) 

      ts ( ix, :, 18) = 0.5 * ( s_quadpt ( ix, :, 31)
     1 + s_quadpt ( ix, :, 32) - s_quadpt ( ix, :, 33)
     1 - s_quadpt ( ix, :, 34) )


      s ( ix, :, 05) = ts ( ix, :, 01) + ts ( ix, :, 02)

      s ( ix, :, 06) = ts ( ix, :, 03) + ts ( ix, :, 04)

      s ( ix, :, 07) = ts ( ix, :, 05) + ts ( ix, :, 06)

      s ( ix, :, 08) = s_quadpt ( ix, :, 23) - s_quadpt ( ix, :, 27)
     1 - s_quadpt ( ix, :, 24) + s_quadpt ( ix, :, 28)

      s ( ix, :, 09) = s_quadpt ( ix, :, 31) - s_quadpt ( ix, :, 32)
     1 - s_quadpt ( ix, :, 33) + s_quadpt ( ix, :, 34)

      s ( ix, :, 10) = s_quadpt ( ix, :, 25) - s_quadpt ( ix, :, 26)
     1 - s_quadpt ( ix, :, 29) + s_quadpt ( ix, :, 30)


      s ( ix, :, 11) = 16 * ( s_quadpt ( ix, :, 02)
     1 - 2.0 * s_quadpt ( ix, :, 03) + 2.0 * s_quadpt ( ix, :, 04)
     1 - s_quadpt ( ix, :, 05) ) / 3.0

      s ( ix, :, 12) = 16 * ( s_quadpt ( ix, :, 06)
     1 - 2.0 * s_quadpt ( ix, :, 07) + 2.0 * s_quadpt ( ix, :, 08)
     1 - s_quadpt ( ix, :, 09) ) / 3.0

      s ( ix, :, 13) = 16 * ( s_quadpt ( ix, :, 10)
     1 - 2.0 * s_quadpt ( ix, :, 11) + 2.0 * s_quadpt ( ix, :, 12)
     1 - s_quadpt ( ix, :, 13) ) / 3.0

      s ( ix, :, 16) = ts ( ix, :, 13) + 0.5 * ( ts ( ix, :, 07)
     1 - ts ( ix, :, 08) ) - 2.0 * ( s_quadpt ( ix, :, 02)
     1 - s_quadpt ( ix, :, 05) )

      s ( ix, :, 18) = - 2.0 * ts ( ix, :, 13) + s ( ix, :, 16)

      s ( ix, :, 14) = ts ( ix, :, 14)
     1 + 0.5 * ( ts ( ix, :, 09) - ts ( ix, :, 10) )
     1 - 2.0 * ( s_quadpt ( ix, :, 06) - s_quadpt ( ix, :, 09) )

      s ( ix, :, 19) = - 2.0 * ts ( ix, :, 14) + s ( ix, :, 14)

      s ( ix, :, 15) = ts ( ix, :, 15) + 0.5 * ( ts ( ix, :, 11)
     1 - ts ( ix, :, 12) ) - 2.0 * ( s_quadpt ( ix, :, 10)
     1 - s_quadpt ( ix, :, 13) )

      s ( ix, :, 17) = - 2.0 * ts ( ix, :, 15) + s ( ix, :, 15)

      s ( ix, :, 20) = s_quadpt ( ix, :, 14) - s_quadpt ( ix, :, 18)
     1 - s_quadpt ( ix, :, 15) + s_quadpt ( ix, :, 19)
     1 - s_quadpt ( ix, :, 16) + s_quadpt ( ix, :, 17)
     1 + s_quadpt ( ix, :, 20) - s_quadpt ( ix, :, 21)


      s ( ix, :, 30) = 2.0 * ( ts ( ix, :, 02) - ts ( ix, :, 01) )

      s ( ix, :, 31) = 2.0 * ( ts ( ix, :, 04) - ts ( ix, :, 03) )

      s ( ix, :, 32) = 2.0 * ( ts ( ix, :, 06) - ts ( ix, :, 05) )

      s ( ix, :, 33) = 2.0 * s ( ix, :, 08) - s_quadpt ( ix, :, 14)
     1 + s_quadpt ( ix, :, 15) - s_quadpt ( ix, :, 16)
     1 + s_quadpt ( ix, :, 17) + s_quadpt ( ix, :, 18)
     1 - s_quadpt ( ix, :, 19) + s_quadpt ( ix, :, 20)
     1 - s_quadpt ( ix, :, 21)

      s ( ix, :, 34) = 2.0 * s ( ix, :, 09) - s_quadpt ( ix, :, 14)
     1 + s_quadpt ( ix, :, 15) + s_quadpt ( ix, :, 16)
     1 - s_quadpt ( ix, :, 17) - s_quadpt ( ix, :, 18)
     1 + s_quadpt ( ix, :, 19) + s_quadpt ( ix, :, 20)
     1 - s_quadpt ( ix, :, 21)

      s ( ix, :, 35) = 2.0 * s ( ix, :, 10) - s_quadpt ( ix, :, 14)
     1 - s_quadpt ( ix, :, 15) + s_quadpt ( ix, :, 16)
     1 + s_quadpt ( ix, :, 17) + s_quadpt ( ix, :, 18)
     1 + s_quadpt ( ix, :, 19) - s_quadpt ( ix, :, 20)
     1 - s_quadpt ( ix, :, 21)


      s ( ix, :, 22) = 4.0 * ( s_quadpt ( ix, :, 35)
     1 - 2.0 * s_quadpt ( ix, :, 22) + s_quadpt ( ix, :, 01) ) / 3.0

      s ( ix, :, 23) = 8.0 * ( s_quadpt ( ix, :, 42)
     1 - s_quadpt ( ix, :, 01) - 3.0 * s_quadpt ( ix, :, 35)
     1 + 3.0 * s_quadpt ( ix, :, 22) ) / 15.0

      s ( ix, :, 21) = s_quadpt ( ix, :, 35) - s_quadpt ( ix, :, 01)
     1 + 0.875 * s ( ix, :, 23) + ( s ( ix, :, 30) + s ( ix, :, 31)
     1 + s ( ix, :, 32) ) / 12.0

      s ( ix, :, 01) = ( s_quadpt ( ix, :, 01) + s_quadpt ( ix, :, 22)
     1 + s_quadpt ( ix, :, 35) + 0.25 * ( s ( ix, :, 05)
     1 + s ( ix, :, 06) + s ( ix, :, 07) )
     1 + 0.75 * s ( ix, :, 22) ) / 3.0


      s ( ix, :, 24) = - s_quadpt ( ix, :, 02) + s_quadpt ( ix, :, 05)
     1 + s_quadpt ( ix, :, 36) - s_quadpt ( ix, :, 37)

      s ( ix, :, 02) = ( 4.0 * s_quadpt ( ix, :, 02)
     1 - 4.0 * s_quadpt ( ix, :, 05) + 2.0 * s ( ix, :, 24)
     1 - ts ( ix, :, 16) ) / 3.0 - 0.1 * s ( ix, :, 11)
     1 + ( 2.0 * s ( ix, :, 16) + s ( ix, :, 18) ) / 12.0

      s ( ix, :, 27) = - 2.0 * ts ( ix, :, 16) + 2.0 * s ( ix, :, 02)
     1 + 0.2 * s ( ix, :, 11) + ( 2.0 * s ( ix, :, 16)
     1 - s ( ix, :, 18) ) / 6.0


      s ( ix, :, 25) = - s_quadpt ( ix, :, 06) + s_quadpt ( ix, :, 09)
     1 + s_quadpt ( ix, :, 38) - s_quadpt ( ix, :, 39)

      s ( ix, :, 03) = ( 4.0 * s_quadpt ( ix, :, 06)
     1 - 4.0 * s_quadpt ( ix, :, 09) + 2.0 * s ( ix, :, 25)
     1 - ts ( ix, :, 17) ) / 3.0 - 0.1 * s ( ix, :, 12)
     1 + ( 2.0 * s ( ix, :, 14) + s ( ix, :, 19) ) / 12.0

      s ( ix, :, 28) = - 2.0 * ts ( ix, :, 17) + 2.0 * s ( ix, :, 03)
     1 + 0.2 * s ( ix, :, 12) + ( 2.0 * s ( ix, :, 14)
     1 - s ( ix, :, 19) ) / 6.0


      s ( ix, :, 26) = - s_quadpt ( ix, :, 10) + s_quadpt ( ix, :, 13)
     1 + s_quadpt ( ix, :, 40) - s_quadpt ( ix, :, 41)

      s ( ix, :, 04) = ( 4.0 * s_quadpt ( ix, :, 10)
     1 - 4.0 * s_quadpt ( ix, :, 13) + 2.0 * s ( ix, :, 26)
     1 - ts ( ix, :, 18) ) / 3.0 - 0.1 * s ( ix, :, 13)
     1 + ( s ( ix, :, 15) + 2.0 * s ( ix, :, 17)) / 12.0

      s ( ix, :, 29) = - 2.0 * ts ( ix, :, 18) + 2.0 * s ( ix, :, 04)
     1 + 0.2 * s ( ix, :, 13) - ( s ( ix, :, 15)
     1 - 2.0 * s ( ix, :, 17)) / 6.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 9) By this point, the flux coefficients ( f, g, h) and the
! un-iterated source coefficients "s" that enable us to define the
! fluxes and source terms any where in space and time within the
! zone have already been loaded.

! Now we can use the analytical inversion of the matrix to obtain the
! rhs "rhs_um" of the system that has to be iterated to convergence
! when stiff source terms are present. The Jacobian matrix of the source
! terms will be made in STEP 10 as will the block tridiagonal matrices.
! Inversion of the resulting block tri-di matrix in STEP 11 will give
! us the fluctuations that have to be added to the next iterate
! in STEP 12 in order to obtain the updated conserved variables "u" .

! The coefficients "u" will then enable us to construct the conserved
! variables at any space-time point in the zone.

! Note that the matrix inversion in STEP 11 is best accomplished by
! permuting the current order of the modal bases. We therefore make
! that permutation in this step.

! Since "u" is made by the end of the next three steps, we also
! build and return "s_int", the space-time integral of the source term.

! ----------------------------------------------------------------------

#if ( NUMORDER == 2)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

! rhs_um1 <-- d_um1
      rhs_um ( ix, :, 1) = - u ( ix, :, 1)
     1          + ( - 3.0 * f ( ix, :, 2) - 3.0 * g ( ix, :, 3)
     1            - 3.0 * h ( ix, :, 4)
     1            + 3.0 * s ( ix, :, 1) - s ( ix, :, 5)
     1            + 6.0 * w ( ix, :, 1) ) * 0.1666666666666666

! rhs_um3 <-- d_um2
      rhs_um ( ix, :, 3) = - u ( ix, :, 2)
     1     + s ( ix, :, 2) + w ( ix, :, 2)

! rhs_um4 <-- d_um3
      rhs_um ( ix, :, 4) = - u ( ix, :, 3)
     1     + s ( ix, :, 3) + w ( ix, :, 3)

! rhs_um5 <-- d_um4
      rhs_um ( ix, :, 5) = - u ( ix, :, 4)
     1     + s ( ix, :, 4) + w ( ix, :, 4)

! rhs_um2 <-- d_um5
      rhs_um ( ix, :, 2) = - u ( ix, :, 5)
     1          + ( - 3.0 * f ( ix, :, 2) - 3.0 * g ( ix, :, 3)
     1            - 3.0 * h ( ix, :, 4)
     1            + 3.0 * s ( ix, :, 1) + s ( ix, :, 5) )
     1          * 0.1666666666666666

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2".

! ----------------------------------------------------------------------

#if ( NUMORDER == 3)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

! rhs_um1 <-- d_um1
      rhs_um ( ix, :, 1) = - u ( ix, :, 1)
     1           + ( f ( ix, :, 13) - 3.0 * f ( ix, :, 2)
     1             + g ( ix, :, 14) - 3.0 * g ( ix, :, 3)
     1             + h ( ix, :, 15) - 3.0 * h ( ix, :, 4)
     1             + 3.0 * s ( ix, :, 1) - s ( ix, :, 11)
     1             + 6.0 * w ( ix, :, 1) ) * 0.1666666666666666

! rhs_um4 <-- d_um2
      rhs_um ( ix, :, 4) = - u ( ix, :, 2)
     1           + ( - 6.0 * f ( ix, :, 5) - 3.0 * g ( ix, :, 8)
     1             - 3.0 * h ( ix, :, 10)
     1             - s ( ix, :, 13) + 3.0 * s ( ix, :, 2)
     1             + 6.0 * w ( ix, :, 2) ) * 0.1666666666666666

! rhs_um6 <-- d_um3
      rhs_um ( ix, :, 6) = - u ( ix, :, 3)
     1           + ( - 3.0 * f ( ix, :, 8) - 6.0 * g ( ix, :, 6)
     1             - 3.0 * h ( ix, :, 9)
     1             - s ( ix, :, 14) + 3.0 * s ( ix, :, 3)
     1             + 6.0 * w ( ix, :, 3) ) * 0.1666666666666666

! rhs_um8 <-- d_um4
      rhs_um ( ix, :, 8) = - u ( ix, :, 4)
     1           + ( - 3.0 * f ( ix, :, 10) - 3.0 * g ( ix, :, 9)
     1             - 6.0 * h ( ix, :, 7)
     1             - s ( ix, :, 15) + 3.0 * s ( ix, :, 4)
     1             + 6.0 * w ( ix, :, 4) ) * 0.1666666666666666

! rhs_um10 <-- d_um5
      rhs_um ( ix, :, 10) = - u ( ix, :, 5)
     1      + s ( ix, :, 5) + w ( ix, :, 5)

! rhs_um11 <-- d_um6
      rhs_um ( ix, :, 11) = - u ( ix, :, 6)
     1      + s ( ix, :, 6) + w ( ix, :, 6)

! rhs_um12 <-- d_um7
      rhs_um ( ix, :, 12) = - u ( ix, :, 7)
     1      + s ( ix, :, 7) + w ( ix, :, 7)

! rhs_um13 <-- d_um8
      rhs_um ( ix, :, 13) = - u ( ix, :, 8)
     1      + s ( ix, :, 8) + w ( ix, :, 8)

! rhs_um14 <-- d_um9
      rhs_um ( ix, :, 14) = - u ( ix, :, 9)
     1      + s ( ix, :, 9) + w ( ix, :, 9)

! rhs_um15 <-- d_um10
      rhs_um ( ix, :, 15) = - u ( ix, :, 10)
     1     + s ( ix, :, 10) + w ( ix, :, 10)

! rhs_um2 <-- d_um11
      rhs_um ( ix, :, 2) = - u ( ix, :, 11)
     1           + ( - 5.0 * f ( ix, :, 2) - 5.0 * g ( ix, :, 3)
     1             - 5.0 * h ( ix, :, 4) + 5.0 * s ( ix, :, 1)
     1             - s ( ix, :, 12) ) * 0.1

! rhs_um3 <-- d_um12
      rhs_um ( ix, :, 3) = - u ( ix, :, 12)
     1           + ( - 5.0 * f ( ix, :, 13) - 5.0 * g ( ix, :, 14)
     1             - 5.0 * h ( ix, :, 15) + 5.0 * s ( ix, :, 11)
     1             + 3.0 * s ( ix, :, 12) )
     1           * 0.3333333333333333e-1

! rhs_um5 <-- d_um13
      rhs_um ( ix, :, 5) = - u ( ix, :, 13)
     1           + ( - 6.0 * f ( ix, :, 5) - 3.0 * g ( ix, :, 8)
     1             - 3.0 * h ( ix, :, 10)
     1             + s ( ix, :, 13) + 3.0 * s ( ix, :, 2) )
     1           * 0.1666666666666666

! rhs_um7 <-- d_um14
      rhs_um ( ix, :, 7) = - u ( ix, :, 14)
     1           + ( - 3.0 * f ( ix, :, 8) - 6.0 * g ( ix, :, 6)
     1             - 3.0 * h ( ix, :, 9)
     1             + s ( ix, :, 14) + 3.0 * s ( ix, :, 3) )
     1           * 0.1666666666666666

! rhs_um9 <-- d_um15
      rhs_um ( ix, :, 9) = - u ( ix, :, 15)
     1           + ( - 3.0 * f ( ix, :, 10) - 3.0 * g ( ix, :, 9)
     1             - 6.0 * h ( ix, :, 7)
     1             + s ( ix, :, 15) + 3.0 * s ( ix, :, 4) )
     1           * 0.1666666666666666

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3".

! ----------------------------------------------------------------------

#if ( NUMORDER == 4)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

! rhs_um1 <-- d_um1
      rhs_um ( ix, :, 1) = - u ( ix, :, 1)
     1          + ( - 3.0 * f ( ix, :, 11) - 30.0 * f ( ix, :, 2)
     1            + 10.0 * f ( ix, :, 24) - 3.0 * g ( ix, :, 12)
     1            + 10.0 * g ( ix, :, 25) - 30.0 * g ( ix, :, 3)
     1            - 3.0 * h ( ix, :, 13) + 10.0 * h ( ix, :, 26)
     1            - 30.0 * h ( ix, :, 4) + 30.0 * s ( ix, :, 1)
     1            - 10.0 * s ( ix, :, 21) + 60.0 * w ( ix, :, 1) )
     1          * 0.1666666666666666e-1

! rhs_um5 <-- d_um2
      rhs_um ( ix, :, 5) = - u ( ix, :, 2)
     1          + ( 2.0 * f ( ix, :, 30) - 6.0 * f ( ix, :, 5)
     1            + g ( ix, :, 33) - 3.0 * g ( ix, :, 8)
     1            - 3.0 * h ( ix, :, 10) + h ( ix, :, 35)
     1            + 3.0 * s ( ix, :, 2) - s ( ix, :, 24)
     1            + 6.0 * w ( ix, :, 2) ) * 0.1666666666666666

! rhs_um8 <-- d_um3
      rhs_um ( ix, :, 8) = - u ( ix, :, 3)
     1          + ( f ( ix, :, 33) - 3.0 * f ( ix, :, 8)
     1            + 2.0 * g ( ix, :, 31) - 6.0 * g ( ix, :, 6)
     1            + h ( ix, :, 34) - 3.0 * h ( ix, :, 9)
     1            - s ( ix, :, 25) + 3.0 * s ( ix, :, 3)
     1            + 6.0 * w ( ix, :, 3) ) * 0.1666666666666666

! rhs_um11 <-- d_um4
      rhs_um ( ix, :, 11) = - u ( ix, :, 4)
     1          + ( - 3.0 * f ( ix, :, 10) + f ( ix, :, 35)
     1            + g ( ix, :, 34) - 3.0 * g ( ix, :, 9)
     1            + 2.0 * h ( ix, :, 32) - 6.0 * h ( ix, :, 7)
     1            - s ( ix, :, 26) + 3.0 * s ( ix, :, 4)
     1            + 6.0 * w ( ix, :, 4) ) * 0.1666666666666666

! rhs_um14 <-- d_um5
      rhs_um ( ix, :, 14) = - u ( ix, :, 5)
     1          + ( - 9.0 * f ( ix, :, 11) - 3.0 * g ( ix, :, 14)
     1            - 3.0 * h ( ix, :, 15)
     1            - s ( ix, :, 30) + 3.0 * s ( ix, :, 5)
     1            + 6.0 * w ( ix, :, 5) ) * 0.1666666666666666

! rhs_um16 <-- d_um6
      rhs_um ( ix, :, 16) = - u ( ix, :, 6)
     1          + ( - 3.0 * f ( ix, :, 16) - 9.0 * g ( ix, :, 12)
     1            - 3.0 * h ( ix, :, 17)
     1            - s ( ix, :, 31) + 3.0 * s ( ix, :, 6)
     1            + 6.0 * w ( ix, :, 6) ) * 0.1666666666666666

! rhs_um18 <-- d_um7
      rhs_um ( ix, :, 18) = - u ( ix, :, 7)
     1          + ( - 3.0 * f ( ix, :, 18) - 3.0 * g ( ix, :, 19)
     1            - 9.0 * h ( ix, :, 13)
     1            - s ( ix, :, 32) + 3.0 * s ( ix, :, 7)
     1            + 6.0 * w ( ix, :, 7) ) * 0.1666666666666666

! rhs_um20 <-- d_um8
      rhs_um ( ix, :, 20) = - u ( ix, :, 8)
     1          + ( - 6.0 * f ( ix, :, 14) - 6.0 * g ( ix, :, 16)
     1            - 3.0 * h ( ix, :, 20)
     1            - s ( ix, :, 33) + 3.0 * s ( ix, :, 8)
     1            + 6.0 * w ( ix, :, 8) ) * 0.1666666666666666

! rhs_um22 <-- d_um9
      rhs_um ( ix, :, 22) = - u ( ix, :, 9)
     1          + ( - 3.0 * f ( ix, :, 20) - 6.0 * g ( ix, :, 17)
     1            - 6.0 * h ( ix, :, 19)
     1            - s ( ix, :, 34) + 3.0 * s ( ix, :, 9)
     1            + 6.0 * w ( ix, :, 9) ) * 0.1666666666666666

! rhs_um24 <-- d_um10
      rhs_um ( ix, :, 24) = - u ( ix, :, 10)
     1           + ( - 6.0 * f ( ix, :, 15) - 3.0 * g ( ix, :, 20)
     1             - 6.0 * h ( ix, :, 18)
     1             + 3.0 * s ( ix, :, 10) - s ( ix, :, 35)
     1             + 6.0 * w ( ix, :, 10) ) * 0.1666666666666666

! rhs_um32 <-- d_um11
      rhs_um ( ix, :, 32) = - u ( ix, :, 11)
     1     + s ( ix, :, 11) + w ( ix, :, 11)

! rhs_um33 <-- d_um12
      rhs_um ( ix, :, 33) = - u ( ix, :, 12)
     1      + s ( ix, :, 12) + w ( ix, :, 12)

! rhs_um34 <-- d_um13
      rhs_um ( ix, :, 34) = - u ( ix, :, 13)
     1      + s ( ix, :, 13) + w ( ix, :, 13)

! rhs_um28 <-- d_um14
      rhs_um ( ix, :, 28) = - u ( ix, :, 14)
     1      + s ( ix, :, 14) + w ( ix, :, 14)

! rhs_um30 <-- d_um15
      rhs_um ( ix, :, 30) = - u ( ix, :, 15)
     1      + s ( ix, :, 15) + w ( ix, :, 15)

! rhs_um27 <-- d_um16
      rhs_um ( ix, :, 27) = - u ( ix, :, 16)
     1      + s ( ix, :, 16) + w ( ix, :, 16)

! rhs_um31 <-- d_um17
      rhs_um ( ix, :, 31) = - u ( ix, :, 17)
     1      + s ( ix, :, 17) + w ( ix, :, 17)

! rhs_um26 <-- d_um18
      rhs_um ( ix, :, 26) = - u ( ix, :, 18)
     1      + s ( ix, :, 18) + w ( ix, :, 18)

! rhs_um29 <-- d_um19
      rhs_um ( ix, :, 29) = - u ( ix, :, 19)
     1      + s ( ix, :, 19) + w ( ix, :, 19)

! rhs_um35 <-- d_um20
      rhs_um ( ix, :, 35) = - u ( ix, :, 20)
     1      + s ( ix, :, 20) + w ( ix, :, 20)

! rhs_um2 <-- d_um21
      rhs_um ( ix, :, 2) = - u ( ix, :, 21)
     1           + ( - f ( ix, :, 11) - 10.0 * f ( ix, :, 2)
     1             + 2.0 * f ( ix, :, 27)
     1             - g ( ix, :, 12) + 2.0 * g ( ix, :, 28)
     1             - 10.0 * g ( ix, :, 3)
     1             - h ( ix, :, 13) + 2.0 * h ( ix, :, 29)
     1             - 10.0 * h ( ix, :, 4)
     1             + 10.0 * s ( ix, :, 1)
     1             - 2.0 * s ( ix, :, 22) ) * 0.05

! rhs_um3 <-- d_um22
      rhs_um ( ix, :, 3) = - u ( ix, :, 22)
     1           + ( - 7.0 * f ( ix, :, 24) - 7.0 * g ( ix, :, 25)
     1             - 7.0 * h ( ix, :, 26) + 7.0 * s ( ix, :, 21)
     1             - 3.0 * s ( ix, :, 23) ) / 42.0

! rhs_um4 <-- d_um23
      rhs_um ( ix, :, 4) = - u ( ix, :, 23)
     1           + ( - 7.0 * f ( ix, :, 27) - 7.0 * g ( ix, :, 28)
     1             - 7.0 * h ( ix, :, 29) + 7.0 * s ( ix, :, 22)
     1             + 5.0 * s ( ix, :, 23) ) / 70.0

! rhs_um6 <-- d_um24
      rhs_um ( ix, :, 6) = - u ( ix, :, 24)
     1           + ( - 10.0 * f ( ix, :, 5) - 5.0 * g ( ix, :, 8)
     1             - 5.0 * h ( ix, :, 10) + 5.0 * s ( ix, :, 2)
     1             - s ( ix, :, 27) ) * 0.1

! rhs_um9 <-- d_um25
      rhs_um ( ix, :, 9) = - u ( ix, :, 25)
     1           + ( - 5.0 * f ( ix, :, 8) - 10.0 * g ( ix, :, 6)
     1             - 5.0 * h ( ix, :, 9) - s ( ix, :, 28)
     1             + 5.0 * s ( ix, :, 3) ) * 0.1

! rhs_um12 <-- d_um26
      rhs_um ( ix, :, 12) = - u ( ix, :, 26)
     1           + ( - 5.0 * f ( ix, :, 10) - 5.0 * g ( ix, :, 9)
     1             - 10.0 * h ( ix, :, 7) - s ( ix, :, 29)
     1             + 5.0 * s ( ix, :, 4) ) * 0.1

! rhs_um7 <-- d_um27
      rhs_um ( ix, :, 7) = - u ( ix, :, 27)
     1           + ( - 10.0 * f ( ix, :, 30) - 5.0 * g ( ix, :, 33)
     1             - 5.0 * h ( ix, :, 35) + 5.0 * s ( ix, :, 24)
     1             + 3.0 * s ( ix, :, 27) )
     1           * 0.3333333333333333e-1

! rhs_um10 <-- d_um28
      rhs_um ( ix, :, 10) = - u ( ix, :, 28)
     1           + ( - 5.0 * f ( ix, :, 33) - 10.0 * g ( ix, :, 31)
     1             - 5.0 * h ( ix, :, 34) + 5.0 * s ( ix, :, 25)
     1             + 3.0 * s ( ix, :, 28) )
     1           * 0.3333333333333333e-1

! rhs_um13 <-- d_um29
      rhs_um ( ix, :, 13) = - u ( ix, :, 29)
     1           + ( - 5.0 * f ( ix, :, 35) - 5.0 * g ( ix, :, 34)
     1             - 10.0 * h ( ix, :, 32) + 5.0 * s ( ix, :, 26)
     1             + 3.0 * s ( ix, :, 29) )
     1           * 0.3333333333333333e-1

! rhs_um15 <-- d_um30
      rhs_um ( ix, :, 15) = - u ( ix, :, 30)
     1           + ( - 9.0 * f ( ix, :, 11) - 3.0 * g ( ix, :, 14)
     1             - 3.0 * h ( ix, :, 15)
     1             + s ( ix, :, 30) + 3.0 * s ( ix, :, 5) )
     1           * 0.1666666666666666

! rhs_um17 <-- d_um31
      rhs_um ( ix, :, 17) = - u ( ix, :, 31)
     1           + ( - 3.0 * f ( ix, :, 16) - 9.0 * g ( ix, :, 12)
     1             - 3.0 * h ( ix, :, 17)
     1             + s ( ix, :, 31) + 3.0 * s ( ix, :, 6) )
     1           * 0.1666666666666666

! rhs_um19 <-- d_um32
      rhs_um ( ix, :, 19) = - u ( ix, :, 32)
     1           + ( - 3.0 * f ( ix, :, 18) - 3.0 * g ( ix, :, 19)
     1             - 9.0 * h ( ix, :, 13)
     1             + s ( ix, :, 32) + 3.0 * s ( ix, :, 7) )
     1           * 0.1666666666666666

! rhs_um21 <-- d_um33
      rhs_um ( ix, :, 21) = - u ( ix, :, 33)
     1           + ( - 6.0 * f ( ix, :, 14) - 6.0 * g ( ix, :, 16)
     1             - 3.0 * h ( ix, :, 20)
     1             + s ( ix, :, 33) + 3.0 * s ( ix, :, 8) )
     1           * 0.1666666666666666

! rhs_um23 <-- d_um34
      rhs_um ( ix, :, 23) = - u ( ix, :, 34)
     1           + ( - 3.0 * f ( ix, :, 20) - 6.0 * g ( ix, :, 17)
     1             - 6.0 * h ( ix, :, 19)
     1             + s ( ix, :, 34) + 3.0 * s ( ix, :, 9) )
     1           * 0.1666666666666666

! rhs_um25 <-- d_um35
      rhs_um ( ix, :, 25) = - u ( ix, :, 35)
     1           + ( - 6.0 * f ( ix, :, 15) - 3.0 * g ( ix, :, 20)
     1             - 6.0 * h ( ix, :, 18)
     1             + 3.0 * s ( ix, :, 10) + s ( ix, :, 35) )
     1           * 0.1666666666666666

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 10) In STEP 5 we built "good" variables around which to linearize.
! In this step we use those variables to build the lumped in space and time
! Jacobian matrix of the source terms. I.e. E = -dS/dU.

! The base-level function that we provide here is to build the Jacobian
! from the space-time-averaged primitive variables, "rhobar, prsbar" etc.,
! that can be used as "good" variables around which to linearize
! implicit source terms. I.e. this is the lumped in space and time approach.

! Just like we multiplied the source terms in STEP 8.a by "dtmove", here we
! should multiply "E1" by "dtmove".

#if ( N_EULER == 1)

      IF ( ambipolar_diffusion_drag_coef .GT. 0.0) THEN

      i_euler = 1

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      E1 ( ix, :, :) = 0.0

      E1 ( ix, 2, 2) = ambipolar_diffusion_drag_coef
     1           * rhobar_euler ( ix, i_euler)

      E1 ( ix, 3, 3) = E1 ( ix, 2, 2)
      E1 ( ix, 4, 4) = E1 ( ix, 2, 2)
      E1 ( ix, 10, 2) = - E1 ( ix, 2, 2)
      E1 ( ix, 11, 3) = - E1 ( ix, 2, 2)
      E1 ( ix, 12, 4) = - E1 ( ix, 2, 2)

      E1 ( ix, 2, 10) = - ambipolar_diffusion_drag_coef * rhobar ( ix)
      E1 ( ix, 3, 11) = E1 ( ix, 2, 10)
      E1 ( ix, 4, 12) = E1 ( ix, 2, 10)
      E1 ( ix, 10, 10) = - E1 ( ix, 2, 10)
      E1 ( ix, 11, 11) = - E1 ( ix, 2, 10)
      E1 ( ix, 12, 12) = - E1 ( ix, 2, 10)

      E1 ( ix, 2, 1) = - ambipolar_diffusion_drag_coef
     1      * rhobar_euler ( ix, i_euler) * vlxbar_euler ( ix, i_euler)
      E1 ( ix, 3, 1) = - ambipolar_diffusion_drag_coef 
     1      * rhobar_euler ( ix, i_euler) * vlybar_euler ( ix, i_euler)
      E1 ( ix, 4, 1) = - ambipolar_diffusion_drag_coef 
     1      * rhobar_euler ( ix, i_euler) * vlzbar_euler ( ix, i_euler)
      E1 ( ix, 10, 1) = - E1 ( ix, 2, 1)
      E1 ( ix, 11, 1) = - E1 ( ix, 3, 1)
      E1 ( ix, 12, 1) = - E1 ( ix, 4, 1)

      E1 ( ix, 2, 9) = ambipolar_diffusion_drag_coef
     1                                * rhobar ( ix) * vlxbar ( ix)
      E1 ( ix, 3, 9) = ambipolar_diffusion_drag_coef
     1                                * rhobar ( ix) * vlybar ( ix)
      E1 ( ix, 4, 9) = ambipolar_diffusion_drag_coef
     1                                * rhobar ( ix) * vlzbar ( ix)
      E1 ( ix, 10, 9) = - E1 ( ix, 2, 9)
      E1 ( ix, 11, 9) = - E1 ( ix, 3, 9)
      E1 ( ix, 12, 9) = - E1 ( ix, 4, 9)

      tempaa = ambipolar_diffusion_drag_coef
     1              / ( molewt_ion + molewt_neut)
      tempab = molewt_ion * rhobar_euler ( ix, i_euler)
     1    * ( 5.0 - 3.0 * gammabar ( ix) )
      tempac = ( molewt_neut - molewt_ion)
      tempad = molewt_neut * rhobar ( ix)
     1    * ( 5.0 - 3.0 * gammabar_euler ( ix, i_euler) )
      tempae = vlxbar_euler ( ix, i_euler)**2
     1       + vlybar_euler ( ix, i_euler)**2
     1       + vlzbar_euler ( ix, i_euler)**2
      tempaf = vlxbar ( ix)**2 + vlybar ( ix)**2 + vlzbar ( ix)**2


      E1 ( ix, 5, 1) = - tempaa
     1 * (   molewt_neut * ( 3.0 * prsbar_euler ( ix, i_euler)
     1                     + rhobar_euler ( ix, i_euler) * tempae)
     1     + 0.5 * tempab * tempaf )

      E1 ( ix, 5, 2) = tempaa
     1 * ( tempac * rhobar_euler ( ix, i_euler)
     1            * vlxbar_euler ( ix, i_euler)
     1   + tempab * vlxbar ( ix) )

      E1 ( ix, 5, 3) = tempaa
     1 * ( tempac * rhobar_euler ( ix, i_euler)
     1            * vlybar_euler ( ix, i_euler)
     1   + tempab * vlybar ( ix) )

      E1 ( ix, 5, 4) = tempaa
     1 * ( tempac * rhobar_euler ( ix, i_euler)
     1            * vlzbar_euler ( ix, i_euler)
     1   + tempab * vlzbar ( ix) )

      E1 ( ix, 5, 5) = 3.0 * tempaa * molewt_ion
     1       * rhobar_euler ( ix, i_euler) * ( gammabar ( ix) - 1.0)

      E1 ( ix, 5, 6) = - E1 ( ix, 5, 5) * bfxbar ( ix) * ob4pi
      E1 ( ix, 5, 7) = - E1 ( ix, 5, 5) * bfybar ( ix) * ob4pi
      E1 ( ix, 5, 8) = - E1 ( ix, 5, 5) * bfzbar ( ix) * ob4pi

      E1 ( ix, 5, 9) = tempaa
     1 * ( molewt_ion * ( 3.0 * prsbar ( ix) + rhobar ( ix) * tempaf)
     1   + 0.5 * tempad * tempae)

      E1 ( ix, 5, 10) = tempaa * ( tempac * rhobar ( ix) * vlxbar ( ix)
     1                  - tempad * vlxbar_euler ( ix, i_euler) )
      E1 ( ix, 5, 11) = tempaa * ( tempac * rhobar ( ix) * vlybar ( ix)
     1                  - tempad * vlybar_euler ( ix, i_euler) )
      E1 ( ix, 5, 12) = tempaa * ( tempac * rhobar ( ix) * vlzbar ( ix)
     1                  - tempad * vlzbar_euler ( ix, i_euler) )

      E1 ( ix, 5, 13) = - 3.0 * tempaa * molewt_neut
     1        * rhobar ( ix) * ( gammabar_euler ( ix, i_euler) - 1.0)

      E1 ( ix, 13, :) = - E1 ( ix, 5, :)

! ----------

      E1 ( ix, :, :) = dtmove * E1 ( ix, :, :)

! ----------

      END DO
! End of "DO ix =" loop.

      END IF
! End of "IF ( ambipolar_diffusion_drag_coef .GT. 0.0)"

#endif
! End of "N_EULER == 1"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 11) This is where we do the matrix assembly and inversion.

! Assemble the smaller matrix sub-systems for the desired
! order and invert the matrix sub-systems in the most economical
! fashion to get "d_um", the incremental fluctuation that has to
! be added to the modes.

! Note. if we are not solving a stiff system but just want to
! exercise the ADER-DG for whatever other reason, then it is good
! to uncomment the next line.

#if ( SKIP_LINEAR_ALGEBRA == 1)

! When the matrix is the identity matrix, the RHS "rhs_um" is indeed
! the solution "d_um".

      DO ix = oned_begin - oned_off, oned_end + oned_off

      d_um ( ix, :, :) = rhs_um ( ix, :, :)

      END DO ! End of "DO ix =" loop.
      GO TO 1000
#endif

! ----------------------------------------------------------------------

! Now zip thru the next "n_blocks = 1, NUMORDER". In doing so, we
! evaluate "num_1by1_blocks" 1x1 blocks, "num_2by2_blocks" 2x2 blocks,
! "num_3by3_blocks" 3x3 blocks and "num_4by4_blocks" 4x4 blocks.
! But remember that they are one, two, three and four times as large
! as the 1x1 blocks.

! Note though that we work backward from the last to first block for
! any given value of "n_blocks". I.e. the "iblock" loop induction variable
! is decreasing.

      nblocks_done = 0

      DO n_blocks = 1, NUMORDER

      DO iblock = iblock_index_begin ( n_blocks),
     1            iblock_index_end ( n_blocks), - n_blocks

! ------------------------------

      IF ( iblock == iblock_index_begin ( n_blocks) ) THEN

        new_call = 1

! --------------------
! Assemble the block tridiagonal matrices here because it is a fresh
! initialization. Assembling diagonal block matrices here; i.e. the
! "a" matrices.

        DO j = 1, n_blocks

        DO ix = oned_begin - oned_off, oned_end + oned_off

          a ( ix, 1: mat_rank, 1: mat_rank, j) =
     1       Identity_Matrix ( 1: mat_rank, 1: mat_rank)
     1     + a_mat_coef ( n_blocks, j)
     1     * E1 ( ix, 1: mat_rank, 1: mat_rank)

        END DO ! End of "DO ix =" loop.

        END DO ! End of "DO j =" loop.

! --------------------
! Assembling super-diagonal block matrices here; i.e. the "b" matrices.

        DO j = 1, n_blocks - 1

        DO ix = oned_begin - oned_off, oned_end + oned_off

          b ( ix, 1: mat_rank, 1: mat_rank, j) =
     1       b_mat_coef ( n_blocks, j)
     1     * E1 ( ix, 1: mat_rank, 1: mat_rank)

        END DO ! End of "DO ix =" loop.

        END DO ! End of "DO j =" loop.

! --------------------
! Assembling sub-diagonal block matrices here; i.e. the "c" matrices.

        DO j = 2, n_blocks

        DO ix = oned_begin - oned_off, oned_end + oned_off

          c ( ix, 1: mat_rank, 1: mat_rank, j) =
     1       c_mat_coef ( n_blocks, j)
     1     * E1 ( ix, 1: mat_rank, 1: mat_rank)

        END DO ! End of "DO ix =" loop.

        END DO ! End of "DO j =" loop.

! ------------------------------

      ELSE

! ------------------------------

        new_call = 0

! ------------------------------

      END IF

! ------------------------------

! Assemble the right hand sides "r"; this is done for each iteration
! of "iblock".

        DO j = 1, n_blocks

        DO ix = oned_begin - oned_off, oned_end + oned_off

          r ( ix, 1: mat_rank, j) =
     1    rhs_um ( ix, 1: mat_rank, iblock - n_blocks + j)

        END DO ! End of "DO ix =" loop.

        END DO ! End of "DO j =" loop.

! ------------------------------

! And call the "VTHOMAS_BLOCK" subroutine.

        nriembegin = oned_begin - oned_off
        nriemend = oned_end + oned_off

        CALL VTHOMAS_BLOCK
     1       ( nriembegin, nriemend, mat_storage, mat_rank,
     1         nmax_blocks, n_blocks, new_call,
     1         a, b, c, dinv, dinv_c, r)

! ------------------------------

! Hand back the right hand sides "r" to "d_um", i.e. these are the
! incremental changes in the modes.

        DO j = 1, n_blocks

        DO ix = oned_begin - oned_off, oned_end + oned_off

          d_um ( ix, 1: mat_rank, iblock - n_blocks + j) =
     1       r ( ix, 1: mat_rank, j)

        END DO ! End of "DO ix =" loop.

        END DO ! End of "DO j =" loop.

! ------------------------------

      nblocks_done = nblocks_done + n_blocks

! ------------------------------

      END DO ! End of "DO iblock =" for 1x1 through 4x4 blocks.

      END DO ! End of "DO n_blocks =" loop,i.e.loop over "n_block" size blocks.

! ------------------------------

      IF ( nblocks_done .NE. NUM_SPACETIME_BASES) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_DG_MODAL_3D_TIME_EVOL_MHD, ",
     1  "nblocks_done .NE. NUM_","SPACETIME_BASES : nblocks_done = ",
     1   nblocks_done, NUM_SPACETIME_BASES
      STOP
      END IF

! ----------------------------------------------------------------------
! Get here if we want to avoid matrix inversion in the previous step.

1000  CONTINUE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 12) By this point, the matrices for the stiff-source
! term iteration have been iterated. Thus we have "d_um" the
! incremental fluctuation that has to be added to the modes.
! We add that part in here, noting that the permutation that
! was performed in STEP 10 has to be reversed here.

! ----------------------------------------------------------------------

#if ( NUMORDER == 2)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u ( ix, :, 1) = u ( ix, :, 1) + d_um ( ix, :, 1)

      u ( ix, :, 2) = u ( ix, :, 2) + d_um ( ix, :, 3)

      u ( ix, :, 3) = u ( ix, :, 3) + d_um ( ix, :, 4)

      u ( ix, :, 4) = u ( ix, :, 4) + d_um ( ix, :, 5)

      u ( ix, :, 5) = u ( ix, :, 5) + d_um ( ix, :, 2)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2".

! ----------------------------------------------------------------------

#if ( NUMORDER == 3)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u ( ix, :, 1) = u ( ix, :, 1) + d_um ( ix, :, 1)

      u ( ix, :, 2) = u ( ix, :, 2) + d_um ( ix, :, 4)

      u ( ix, :, 3) = u ( ix, :, 3) + d_um ( ix, :, 6)

      u ( ix, :, 4) = u ( ix, :, 4) + d_um ( ix, :, 8)

      u ( ix, :, 5) = u ( ix, :, 5) + d_um ( ix, :, 10)

      u ( ix, :, 6) = u ( ix, :, 6) + d_um ( ix, :, 11)

      u ( ix, :, 7) = u ( ix, :, 7) + d_um ( ix, :, 12)

      u ( ix, :, 8) = u ( ix, :, 8) + d_um ( ix, :, 13)

      u ( ix, :, 9) = u ( ix, :, 9) + d_um ( ix, :, 14)

      u ( ix, :, 10) = u ( ix, :, 10) + d_um ( ix, :, 15)

      u ( ix, :, 11) = u ( ix, :, 11) + d_um ( ix, :, 2)

      u ( ix, :, 12) = u ( ix, :, 12) + d_um ( ix, :, 3)

      u ( ix, :, 13) = u ( ix, :, 13) + d_um ( ix, :, 5)

      u ( ix, :, 14) = u ( ix, :, 14) + d_um ( ix, :, 7)

      u ( ix, :, 15) = u ( ix, :, 15) + d_um ( ix, :, 9)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3".

! ----------------------------------------------------------------------

#if ( NUMORDER == 4)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u ( ix, :, 1) = u ( ix, :, 1) + d_um ( ix, :, 1)

      u ( ix, :, 2) = u ( ix, :, 2) + d_um ( ix, :, 5)

      u ( ix, :, 3) = u ( ix, :, 3) + d_um ( ix, :, 8)

      u ( ix, :, 4) = u ( ix, :, 4) + d_um ( ix, :, 11)

      u ( ix, :, 5) = u ( ix, :, 5) + d_um ( ix, :, 14)

      u ( ix, :, 6) = u ( ix, :, 6) + d_um ( ix, :, 16)

      u ( ix, :, 7) = u ( ix, :, 7) + d_um ( ix, :, 18)

      u ( ix, :, 8) = u ( ix, :, 8) + d_um ( ix, :, 20)

      u ( ix, :, 9) = u ( ix, :, 9) + d_um ( ix, :, 22)

      u ( ix, :, 10) = u ( ix, :, 10) + d_um ( ix, :, 24)

      u ( ix, :, 11) = u ( ix, :, 11) + d_um ( ix, :, 32)

      u ( ix, :, 12) = u ( ix, :, 12) + d_um ( ix, :, 33)

      u ( ix, :, 13) = u ( ix, :, 13) + d_um ( ix, :, 34)

      u ( ix, :, 14) = u ( ix, :, 14) + d_um ( ix, :, 28)

      u ( ix, :, 15) = u ( ix, :, 15) + d_um ( ix, :, 30)

      u ( ix, :, 16) = u ( ix, :, 16) + d_um ( ix, :, 27)

      u ( ix, :, 17) = u ( ix, :, 17) + d_um ( ix, :, 31)

      u ( ix, :, 18) = u ( ix, :, 18) + d_um ( ix, :, 26)

      u ( ix, :, 19) = u ( ix, :, 19) + d_um ( ix, :, 29)

      u ( ix, :, 20) = u ( ix, :, 20) + d_um ( ix, :, 35)

      u ( ix, :, 21) = u ( ix, :, 21) + d_um ( ix, :, 2)

      u ( ix, :, 22) = u ( ix, :, 22) + d_um ( ix, :, 3)

      u ( ix, :, 23) = u ( ix, :, 23) + d_um ( ix, :, 4)

      u ( ix, :, 24) = u ( ix, :, 24) + d_um ( ix, :, 6)

      u ( ix, :, 25) = u ( ix, :, 25) + d_um ( ix, :, 9)

      u ( ix, :, 26) = u ( ix, :, 26) + d_um ( ix, :, 12)

      u ( ix, :, 27) = u ( ix, :, 27) + d_um ( ix, :, 7)

      u ( ix, :, 28) = u ( ix, :, 28) + d_um ( ix, :, 10)

      u ( ix, :, 29) = u ( ix, :, 29) + d_um ( ix, :, 13)

      u ( ix, :, 30) = u ( ix, :, 30) + d_um ( ix, :, 15)

      u ( ix, :, 31) = u ( ix, :, 31) + d_um ( ix, :, 17)

      u ( ix, :, 32) = u ( ix, :, 32) + d_um ( ix, :, 19)

      u ( ix, :, 33) = u ( ix, :, 33) + d_um ( ix, :, 21)

      u ( ix, :, 34) = u ( ix, :, 34) + d_um ( ix, :, 23)

      u ( ix, :, 35) = u ( ix, :, 35) + d_um ( ix, :, 25)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4".

! ----------------------------------------------------------------------

! STEP 13) Returning the space-time integral of source terms that are treated
! time-explicitly is easy and is done here. If the treatment is time-implicit,
! the source term "s ( :, 1)" should be rebuilt here using "u ( :, :)" from
! the most recent iteration.

! If the implicit iteration has converged, then this step hands converged
! source terms that is lagged by one iteration sub-cycle. In most instances
! this might be acceptable, especially if the solution is well-converged.

! In principle, the factor of "dtmove", that was previously multiplied into
! the source terms, is divided out of "s_int" here. However,
! because the "dudt" terms in the subroutine that calls this one already
! carry the source terms, we leave the "dtmove" multiplied into the
! "s_int" terms.

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_int ( ix, :) = s ( ix, :, 1)

! ----------

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO
! End of "DO iter_ader =" loop.

! End the outer ADER iteration here.

! ----------------------------------------------------------------------

! STEP 14) Before exit, we undo all the multiplications by "dtmove / dx,y,z"
! to the fluxes.

      IF ( dtmove .GT. 0.0) THEN

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = dx ( ix) / dtmove

      f ( ix, :, 1: NUM_SPACETIME_BASES) = 
     1               f ( ix, :, 1: NUM_SPACETIME_BASES) * tempaa

      tempaa = dy ( ix) / dtmove

      g ( ix, :, 1: NUM_SPACETIME_BASES) = 
     1               g ( ix, :, 1: NUM_SPACETIME_BASES) * tempaa

      tempaa = dz ( ix) / dtmove

      h ( ix, :, 1: NUM_SPACETIME_BASES) = 
     1               h ( ix, :, 1: NUM_SPACETIME_BASES) * tempaa

! ----------

      END DO
! End of "DO ix =" loop.

      END IF

! ----------------------------------------------------------------------

! STEP 15) Here we give ourselves the option of trapping a zone that has
! not converged. This may consist of zeroing out higher order moments, i.e.
! we fall back to first order here if we absolutely must.

! This is also where we can use Dumbser's trick of zeroing out the
! higher moments in zones that produce negative pressures or densities.

      corrupted_zone ( oned_begin - oned_off : oned_end + oned_off) = 0

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

        IF ( ( rho ( ix, iquadpt) .LT. rhofloor) .OR.
     1       ( prs ( ix, iquadpt) .LT. prsfloor) )
     1                                      corrupted_zone ( ix) = 1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO

! ----------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

        IF ( ( rho_euler ( ix, iquadpt, i_euler) .LT. rhofloor) .OR.
     1       ( prs_euler ( ix, iquadpt, i_euler) .LT. prsfloor) )
     1                                     corrupted_zone ( ix) = 1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ----------

! For ADER_DG "f ( :, 1)" is just replaced by "f_bar", the zone-centered
! flux at t=0. This is done for all zones that are corrupted.

! WARNING : It might also help to reset "s" and "s_int" in the same
! way that "f,g,h" are reset. But that is not done yet.

      DO iquadpt = 2, NUM_SPACETIME_BASES

! ----------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( corrupted_zone ( ix) .EQ. 1) THEN

! -----

      s_int ( ix, :) = 0.0

      w ( ix, :, iquadpt) = 0.0

      u ( ix, :, iquadpt) = 0.0

      f ( ix, :, iquadpt) = 0.0

      g ( ix, :, iquadpt) = 0.0

      h ( ix, :, iquadpt) = 0.0

      s ( ix, :, iquadpt) = 0.0

! -----

      END IF

! ----------

      END DO
! End of "DO ix =" loop.

! ----------

      END DO
! End "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 16) Finally, we transcribe the DG modes to work with the CG flux
! subroutine.  This is required because the time-dependence in the DG modes
! is different than in the CG modes.

      DO ix = oned_begin - oned_off, oned_end + oned_off
      
#if ( NUMORDER == 2 )

      u ( ix, :, 1) = u ( ix, :, 1) - u ( ix, :, 5)

      u ( ix, :, 5) = 2.0 * u ( ix, :, 5)

#elif ( NUMORDER == 3 )

      u ( ix, :, 01) = u ( ix, :, 01) - u ( ix, :, 11) + u ( ix, :, 12)

      u ( ix, :, 02) = u ( ix, :, 02) - u ( ix, :, 13)

      u ( ix, :, 03) = u ( ix, :, 03) - u ( ix, :, 14)

      u ( ix, :, 04) = u ( ix, :, 04) - u ( ix, :, 15)

      u ( ix, :, 11) = 2.0 * u ( ix, :, 11) - 6.0 * u ( ix, :, 12)

      u ( ix, :, 12) = 6.0 * u ( ix, :, 12)

      u ( ix, :, 13) = 2.0 * u ( ix, :, 13)

      u ( ix, :, 14) = 2.0 * u ( ix, :, 14)

      u ( ix, :, 15) = 2.0 * u ( ix, :, 15)

#elif ( NUMORDER == 4 )

      u ( ix, :, 01) = u ( ix, :, 01) - u ( ix, :, 21) + u ( ix, :, 22)
     1                                - u ( ix, :, 23)

      u ( ix, :, 02) = u ( ix, :, 02) - u ( ix, :, 24) + u ( ix, :, 27)

      u ( ix, :, 03) = u ( ix, :, 03) - u ( ix, :, 25) + u ( ix, :, 28)

      u ( ix, :, 04) = u ( ix, :, 04) - u ( ix, :, 26) + u ( ix, :, 29)

      u ( ix, :, 05) = u ( ix, :, 05) - u ( ix, :, 30)

      u ( ix, :, 06) = u ( ix, :, 06) - u ( ix, :, 31)

      u ( ix, :, 07) = u ( ix, :, 07) - u ( ix, :, 32)

      u ( ix, :, 08) = u ( ix, :, 08) - u ( ix, :, 33)

      u ( ix, :, 09) = u ( ix, :, 09) - u ( ix, :, 34)

      u ( ix, :, 10) = u ( ix, :, 10) - u ( ix, :, 35)


      u ( ix, :, 21) = 2.0 * u ( ix, :, 21) - 6.0 * u ( ix, :, 22)
     1              + 12.0 * u ( ix, :, 23)

      u ( ix, :, 22) = 6.0 * u ( ix, :, 22) - 30.0 * u ( ix, :, 23)

      u ( ix, :, 23) = 20.0 * u ( ix, :, 23)

      u ( ix, :, 24) = 2.0 * u ( ix, :, 24) - 6.0 * u ( ix, :, 27)

      u ( ix, :, 25) = 2.0 * u ( ix, :, 25) - 6.0 * u ( ix, :, 28)

      u ( ix, :, 26) = 2.0 * u ( ix, :, 26) - 6.0 * u ( ix, :, 29)

      u ( ix, :, 27) = 6.0 * u ( ix, :, 27)

      u ( ix, :, 28) = 6.0 * u ( ix, :, 28)

      u ( ix, :, 29) = 6.0 * u ( ix, :, 29)

      u ( ix, :, 30) = 2.0 * u ( ix, :, 30)

      u ( ix, :, 31) = 2.0 * u ( ix, :, 31)

      u ( ix, :, 32) = 2.0 * u ( ix, :, 32)

      u ( ix, :, 33) = 2.0 * u ( ix, :, 33)

      u ( ix, :, 34) = 2.0 * u ( ix, :, 34)

      u ( ix, :, 35) = 2.0 * u ( ix, :, 35)

#endif

! ----------

      END DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ADER_DG_MODAL_3D_TIME_EVOL_MHD

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADER_CG_MODAL_3D_TIME_EVOL_MHD

     1           ( oned_begin, oned_end, oned_off,
     1             dtmove,
     1             igeom, pi, rhofloor, prsfloor, smallnum,
     1             gamma, gamma_euler, gasconst,
     1             ader_tolerance, max_iter_ader_cg,
     1             tol_rootsolve, tol_superluminal,

     1             a_kerr, viscosity_coef, resistivity_coef,
     1           ambipolar_diffusion_drag_coef, molewt_ion, molewt_neut,
     1             isospeed, barotropic_stiffness,

     1             corrupted_zone,

     1             xc, yc, zc, dx, dy, dz,

     1             w, entropy_1d, gamma_eos_quadpt, molewt_quadpt,
     1             lorenz_quadpt, mu_quadpt,

     1             u, f, g, h, s, s_int)

! ----------------------------------------------------------------------

! Written by Dinshaw Balsara

! Incorporates ADER CG with source terms as contained
! in the paper by Balsara et al 2009, JCP. The present version
! accomodated Galerkin elements that are continuous in time and is,
! therefore, unsuitable for use with stiff source terms.
! It is substantially faster than ADER DG. It is useful for
! situations where the hyperbolic system is large and/or the Jacobian of
! the flux matrix is difficult to evaluate. Current version is for
! single-fluid MHD with multiple species of 2-fluid MHD. It scales linearly
! as the number of species and is, therefore, very useful for reactive flow.

! One MHD fluid and "N_EULER" number of Euler fluids (which may also be
! repurposed as radiation fields in various frequency bins) are treated here.

! Takes the conserved variables and their spatial moments "w" as input
! and returns the space-time moments of the conserved variables as "u".
! This is all done on a reference space-time element in 3d that spans
! (x,y,z,t) in the range [-.5,.5]X[-.5,.5]X[-.5,.5]X[0,1].

! "dtmove, dx, dy, dz" are timestep and x,y,z zone sizes. The update has
! to be done in physical variables so various combinations of these
! variables are multiplied into the flux and source terms.

! "gamma, pi, rhofloor, prsfloor, smallnum" are the usual inputs.

! gamma law gas is assumed, but allowance is made for more general eos if
! the spatial moments are provided in "gamma_eos_spacetime" with the time terms zeroed.

! "max_iter_ader_cg" the maximum number of permitted iterations.

! "f, g, h, s" are fluxes in x,y,z-directions and the source term respectively.
! They are not stored here. 

! "s_int" is the space-time integral of the source term which is also returned.

! ----------------------------------------------------------------------

! It helps to know the space-time modal basis set into which the solution is
! being projected. Here it is for orders 2, 3, and 4.

! P0[x_] := 1 ; P1[x_] := x ; P2[x_] := ( x^2 - 1/12) ;
! P3[x_] := ( x^3 - 3 * x / 20)

! Q0[t_] := 1 ; Q1[t_] := 2 * t - 1 ; Q2[t_] := 6 * t^2 - 6 * t + 1 ;
! Q3[t_] := 20 * t^3 - 30 * t^2 + 12 * t - 1 ; DQ0[t_] := D[Q0[t],t] ;
! DQ1[t_] := D[Q1[t],t] ; DQ2[t_] := D[Q2[t],t] ; DQ3[t_] := D[Q3[t],t]

! ------------------------------

! Second order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P0[x] * P0[y] * P0[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
! [1/2,0,0,0], [-1/2,0,0,0], [0,1/2,0,0], [0,-1/2,0,0],
! [0,0,1/2,0], [0,0,-1/2,0],
!   [0,0,0,1]

! ------------------------------

! Third order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P2[x] * P0[y] * P0[z] * Q0[t] +
!                   u6 * P0[x] * P2[y] * P0[z] * Q0[t] +
!                   u7 * P0[x] * P0[y] * P2[z] * Q0[t] +
!                   u8 * P1[x] * P1[y] * P0[z] * Q0[t] +
!                   u9 * P0[x] * P1[y] * P1[z] * Q0[t] +
!                   u10 * P1[x] * P0[y] * P1[z] * Q0[t] +
!                         u11 * P0[x] * P0[y] * P0[z] * Q1[t] +
!                         u12 * P0[x] * P0[y] * P0[z] * Q2[t] +
!                   u13 * P1[x] * P0[y] * P0[z] * Q1[t] +
!                   u14 * P0[x] * P1[y] * P0[z] * Q1[t] +
!                   u15 * P0[x] * P0[y] * P1[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
!   [0,0,0,0], [1/2,0,0,0], [-1/2,0,0,0], [0,1/2,0,0], [0,-1/2,0,0],
!   [0,0,1/2,0], [0,0,-1/2,0],
! [1/2,1/2,1/2,0], [-1/2,1/2,1/2,0], [1/2,-1/2,1/2,0], [-1/2,-1/2,1/2,0],
! [1/2,1/2,-1/2,0], [-1/2,1/2,-1/2,0], [1/2,-1/2,-1/2,0], [-1/2,-1/2,-1/2,0],
!   [1/2,0,0,1/2], [-1/2,0,0,1/2], [0,1/2,0,1/2], [0,-1/2,0,1/2],
!   [0,0,1/2,1/2], [0,0,-1/2,1/2],
! [0,0,0,1]

! ------------------------------

! Fourth order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P2[x] * P0[y] * P0[z] * Q0[t] +
!                   u6 * P0[x] * P2[y] * P0[z] * Q0[t] +
!                   u7 * P0[x] * P0[y] * P2[z] * Q0[t] +
!                   u8 * P1[x] * P1[y] * P0[z] * Q0[t] +
!                   u9 * P0[x] * P1[y] * P1[z] * Q0[t] +
!                   u10 * P1[x] * P0[y] * P1[z] * Q0[t] +
!                         u11 * P3[x] * P0[y] * P0[z] * Q0[t] +
!                         u12 * P0[x] * P3[y] * P0[z] * Q0[t] +
!                         u13 * P0[x] * P0[y] * P3[z] * Q0[t] +
!                         u14 * P2[x] * P1[y] * P0[z] * Q0[t] +
!                         u15 * P2[x] * P0[y] * P1[z] * Q0[t] +
!                         u16 * P1[x] * P2[y] * P0[z] * Q0[t] +
!                         u17 * P0[x] * P2[y] * P1[z] * Q0[t] +
!                         u18 * P1[x] * P0[y] * P2[z] * Q0[t] +
!                         u19 * P0[x] * P1[y] * P2[z] * Q0[t] +
!                         u20 * P1[x] * P1[y] * P1[z] * Q0[t] +
!                   u21 * P0[x] * P0[y] * P0[z] * Q1[t] +
!                   u22 * P0[x] * P0[y] * P0[z] * Q2[t] +
!                   u23 * P0[x] * P0[y] * P0[z] * Q3[t] +
!                         u24 * P1[x] * P0[y] * P0[z] * Q1[t] +
!                         u25 * P0[x] * P1[y] * P0[z] * Q1[t] +
!                         u26 * P0[x] * P0[y] * P1[z] * Q1[t] +
!                         u27 * P1[x] * P0[y] * P0[z] * Q2[t] +
!                         u28 * P0[x] * P1[y] * P0[z] * Q2[t] +
!                         u29 * P0[x] * P0[y] * P1[z] * Q2[t] +
!                   u30 * P2[x] * P0[y] * P0[z] * Q1[t] +
!                   u31 * P0[x] * P2[y] * P0[z] * Q1[t] +
!                   u32 * P0[x] * P0[y] * P2[z] * Q1[t] +
!                   u33 * P1[x] * P1[y] * P0[z] * Q1[t] +
!                   u34 * P0[x] * P1[y] * P1[z] * Q1[t] +
!                   u35 * P1[x] * P0[y] * P1[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
! [0,0,0,0], [1/2,0,0,0], [1/4,0,0,0], [-1/4,0,0,0], [-1/2,0,0,0],
! [0,1/2,0,0], [0,1/4,0,0], [0,-1/4,0,0], [0,-1/2,0,0],
! [0,0,1/2,0], [0,0,1/4,0], [0,0,-1/4,0], [0,0,-1/2,0],
!   [1/2,1/2,1/2,0], [0,1/2,1/2,0], [-1/2,1/2,1/2,0], [1/2,-1/2,1/2,0],
!   [0,-1/2,1/2,0], [-1/2,-1/2,1/2,0], [1/2,1/2,-1/2,0], [0,1/2,-1/2,0],
!   [-1/2,1/2,-1/2,0], [1/2,-1/2,-1/2,0], [0,-1/2,-1/2,0], [-1/2,-1/2,-1/2,0],
! [1/2,0,1/2,0], [-1/2,0,1/2,0], [1/2,0,-1/2,0], [-1/2,0,-1/2,0],
! [1/2,1/2,0,0], [-1/2,1/2,0,0], [1/2,-1/2,0,0], [-1/2,-1/2,0,0],
!   [0,0,0,1/3], [1/2,0,0,1/3], [-1/2,0,0,1/3], [0,1/2,0,1/3], [0,-1/2,0,1/3],
!   [0,0,1/2,1/3], [0,0,-1/2,1/3],
! [1/2,1/2,1/2,1/3], [-1/2,1/2,1/2,1/3], [1/2,-1/2,1/2,1/3], [-1/2,-1/2,1/2,1/3],
! [1/2,1/2,-1/2,1/3], [-1/2,1/2,-1/2,1/3], [1/2,-1/2,-1/2,1/3], [-1/2,-1/2,-1/2,1/3],
!   [1/2,0,0,2/3], [-1/2,0,0,2/3], [0,1/2,0,2/3], [0,-1/2,0,2/3],
!   [0,0,1/2,2/3], [0,0,-1/2,2/3],
! [0,0,0,1]

! ----------------------------------------------------------------------

      IMPLICIT NONE

! Variables declared here go through the subroutine interface.

      INTEGER oned_begin, oned_end, oned_off, max_iter_ader_cg, igeom

      REAL dtmove, gamma, gamma_euler, pi, gasconst, ader_tolerance, 
     1     tol_rootsolve, tol_superluminal,
     1     rhofloor, prsfloor, smallnum,
     1     a_kerr, viscosity_coef, resistivity_coef,
     1     ambipolar_diffusion_drag_coef, molewt_ion, molewt_neut,
     1     isospeed, barotropic_stiffness

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1                     corrupted_zone

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1                  xc, yc, zc, dx, dy, dz

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_ZONECEN_BASES) :: w

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, 5) :: entropy_1d

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  1: N_EULER + 1, NUM_SPACETIME_QUADPTS) ::
     1                gamma_eos_quadpt, molewt_quadpt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NUM_SPACETIME_QUADPTS) ::
     1                lorenz_quadpt, mu_quadpt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_SPACETIME_BASES) ::
     1                                              u, f, g, h, s

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, NVAR) :: s_int

! ----------------------------------------------------------------------

! Variables declared past this point don't go through the subroutine interface.

      INTEGER ix, iter_ader, iquadpt, ibases, ichk, itest, i_euler,
     1        ivar_lo, ivar

      REAL fourpi, ob4pi, ob8pi, ob3, tempaa, tempab, tempac,
     1     x1, x2, x3, obx1, cos_x2, sin_x2, cot_x2

! These are the flow variables (and their spatial derivatives, if needed)
! at the quadrature points.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NUM_SPACETIME_QUADPTS) ::

     1      rho, prs, eng_t, vlx, vly, vlz, bfx, bfy, bfz,
     1      entropy,

     1      drho_dx, dtemperature_dx, dvlx_dx, dvly_dx, dvlz_dx,
     1      dbfx_dx, dbfy_dx, dbfz_dx, dv2_dx, db2_dx,

     1      drho_dy, dtemperature_dy, dvlx_dy, dvly_dy, dvlz_dy,
     1      dbfx_dy, dbfy_dy, dbfz_dy, dv2_dy, db2_dy,

     1      drho_dz, dtemperature_dz, dvlx_dz, dvly_dz, dvlz_dz,
     1      dbfx_dz, dbfy_dz, dbfz_dz, dv2_dz, db2_dz,

     1      div_vl, temperature

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0: NFLUID, NUM_SPACETIME_QUADPTS) :: specie

#if ( N_EULER >= 1)
      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NUM_SPACETIME_QUADPTS, N_EULER) ::
     1      rho_euler, prs_euler, eng_t_euler,
     1      vlx_euler, vly_euler, vlz_euler

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1             0: NFLUID_EULER, NUM_SPACETIME_QUADPTS, N_EULER) ::
     1      specie_euler
#endif

! These are the conserved variables, fluxes and source terms at the
! quadrature points.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_SPACETIME_QUADPTS) ::

     1      u_quadpt, f_quadpt, g_quadpt, h_quadpt, s_quadpt,
     1      tu, tf, tg, th, ts,

     1      du_quadpt_dx, du_quadpt_dy, du_quadpt_dz,
     1      fni_quadpt, gni_quadpt, hni_quadpt

! List of nodal points in each zone.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  1:4, NUM_SPACETIME_QUADPTS) :: nodal_points

! ----------------------------------------------------------------------

! Make some sanity checks here.

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_MODAL_3D_TIME_EVOL_MHD, NFI","ELD .NE. 7",
     1              NFIELD
      STOP
      END IF

      IF ( NFIELD_EULER .NE. 5) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_MODAL_3D_TIME_EVOL_MHD, NFI",
     1  "ELD_EULER .NE. 5", NFIELD
      STOP
      END IF

#if ( NUMORDER == 2)
      IF ( NUM_ZONECEN_BASES .NE. 4) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_BASES .NE. 4", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 5) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_BASES .NE. 5", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 6) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 6", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 7", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

#if ( NUMORDER == 3)
      IF ( NUM_ZONECEN_BASES .NE. 10) THEN
      WRITE ( 6, *)
     1  "Stopping in third order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_BASES .NE. 10", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 15) THEN
      WRITE ( 6, *)
     1  "Stopping in third order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_BASES .NE. 15", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 13) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 15", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 20) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 22", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

#if ( NUMORDER == 4)
      IF ( NUM_ZONECEN_BASES .NE. 20) THEN
      WRITE ( 6, *)
     1  "Stopping in fourth order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_BASES .NE. 20", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 35) THEN
      WRITE ( 6, *)
     1  "Stopping in fourth order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_BASES .NE. 35", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 21) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 33", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 42) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 55", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

      ivar =
     1  NFIELD + 1 + NFLUID + N_EULER * (NFIELD_EULER + NFLUID_EULER)

      IF ( NVAR .NE. ivar) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1  "NVAR .NE. ivar","NVA","R = ", NVAR," ivar = ", ivar
      STOP
      END IF

      IF ( max_iter_ader_cg .LT. NUMORDER - 1 ) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_MODAL_3D_TIME_EVOL_MHD, ",
     1  "max_iter_ader_cg .LT. NUMORDER ; max_iter_ader_cg = ",
     1   max_iter_ader_cg, NUMORDER - 1
      STOP
      END IF

! ----------------------------------------------------------------------

! Set up some constants that are used all over in this subroutine.

      fourpi = 4.0 * pi
      ob4pi = 1.0 / ( 4.0 * pi)
      ob8pi = 1.0 / ( 8.0 * pi)
      ob3 = 1.0 / 3.0

! ----------------------------------------------------------------------

! STEP 0) Set up the nodal points that are to be used at each of the orders.

! ------------------------------
#if ( NUMORDER == 2)
! cg, 2nd order nodal_points

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1:4, 01) = (/  0.5 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 02) = (/ -0.5 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 03) = (/  0.0 ,  0.5 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 04) = (/  0.0 , -0.5 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 05) = (/  0.0 ,  0.0 ,  0.5 , 0.0  /) 
      nodal_points ( ix, 1:4, 06) = (/  0.0 ,  0.0 , -0.5 , 0.0  /) 

      nodal_points ( ix, 1:4, 07) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /) 

      END DO
! End of "DO ix =" loop.

#endif
! End "NUMORDER == 2".
! ------------------------------
#if ( NUMORDER == 3)
! cg, 3rd order nodal_points

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1:4, 01) = (/  0.0 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 02) = (/  0.5 ,  0.5 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 03) = (/  0.5 , -0.5 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 04) = (/  0.5 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 05) = (/  0.5 ,  0.0 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 06) = (/ -0.5 ,  0.5 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 07) = (/ -0.5 , -0.5 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 08) = (/ -0.5 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 09) = (/ -0.5 ,  0.0 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 10) = (/  0.0 ,  0.5 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 11) = (/  0.0 , -0.5 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 12) = (/  0.0 ,  0.5 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 13) = (/  0.0 , -0.5 , -0.5 , 0.0  /)

      nodal_points ( ix, 1:4, 14) = (/  0.5 ,  0.0 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 15) = (/ -0.5 ,  0.0 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 16) = (/  0.0 ,  0.5 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 17) = (/  0.0 , -0.5 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 18) = (/  0.0 ,  0.0 ,  0.5 , 0.5  /)
      nodal_points ( ix, 1:4, 19) = (/  0.0 ,  0.0 , -0.5 , 0.5  /)

      nodal_points ( ix, 1:4, 20) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /)

      END DO
! End of "DO ix =" loop.

#endif
! End "NUMORDER == 3".
! ------------------------------
#if ( NUMORDER == 4)
! cg, 4th order nodal_points

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1:4, 01) = (/  0.0 ,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 02) = (/  0.5 ,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 03) = (/  0.25,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 04) = (/ -0.25,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 05) = (/ -0.5 ,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 06) = (/  0.0 ,  0.5 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 07) = (/  0.0 ,  0.25,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 08) = (/  0.0 , -0.25,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 09) = (/  0.0 , -0.5 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 10) = (/  0.0 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 11) = (/  0.0 ,  0.0 ,  0.25, 0.0  /)
      nodal_points ( ix, 1:4, 12) = (/  0.0 ,  0.0 , -0.25, 0.0  /)
      nodal_points ( ix, 1:4, 13) = (/  0.0 ,  0.0 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 14) = (/  0.5 ,  0.5 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 15) = (/  0.5 , -0.5 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 16) = (/  0.5 ,  0.5 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 17) = (/  0.5 , -0.5 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 18) = (/ -0.5 ,  0.5 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 19) = (/ -0.5 , -0.5 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 20) = (/ -0.5 ,  0.5 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 21) = (/ -0.5 , -0.5 , -0.5 , 0.0  /)

      nodal_points ( ix, 1:4, 22) = (/  0.0 ,  0.0 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 23) = (/  0.5 ,  0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 24) = (/  0.5 , -0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 25) = (/  0.5 ,  0.0 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 26) = (/  0.5 ,  0.0 , -0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 27) = (/ -0.5 ,  0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 28) = (/ -0.5 , -0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 29) = (/ -0.5 ,  0.0 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 30) = (/ -0.5 ,  0.0 , -0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 31) = (/  0.0 ,  0.5 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 32) = (/  0.0 , -0.5 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 33) = (/  0.0 ,  0.5 , -0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 34) = (/  0.0 , -0.5 , -0.5 , 1.0/3.0 /)

      nodal_points ( ix, 1:4, 35) = (/  0.0 ,  0.0 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 36) = (/  0.5 ,  0.0 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 37) = (/ -0.5 ,  0.0 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 38) = (/  0.0 ,  0.5 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 39) = (/  0.0 , -0.5 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 40) = (/  0.0 ,  0.0 ,  0.5 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 41) = (/  0.0 ,  0.0 , -0.5 , 2.0/3.0 /)

      nodal_points ( ix, 1:4, 42) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /)

      END DO
! End of "DO ix =" loop.

#endif
! End "NUMORDER == 4".
! ------------------------------

#if ( ENFORCE_POSITIVITY == 1)

! While "nodal_points" still represents positions in the reference space-time
! hypercube, we use it to build the "entropy" at all the nodes.

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

      entropy ( ix, iquadpt) = entropy_1d ( ix, 01)
     1    + nodal_points ( ix, 1, iquadpt) * entropy_1d ( ix, 02)
     1    + nodal_points ( ix, 2, iquadpt) * entropy_1d ( ix, 03)
     1    + nodal_points ( ix, 3, iquadpt) * entropy_1d ( ix, 04)
     1    + nodal_points ( ix, 4, iquadpt) * entropy_1d ( ix, 05)

      END DO ! End of "DO ix =" loop.

      END DO ! End of "DO iquadpt =" loop.

#endif
! End of "ENFORCE_POSITIVITY == 1".

! ------------------------------

! Set up the nodal points in a global space-time. These are coordinate
! values.

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1, :) = nodal_points ( ix, 1, :) * dx ( ix)
     1                                                    + xc ( ix)
      nodal_points ( ix, 2, :) = nodal_points ( ix, 2, :) * dy ( ix)
     1                                                    + yc ( ix)
      nodal_points ( ix, 3, :) = nodal_points ( ix, 3, :) * dz ( ix)
     1                                                    + zc ( ix)
      nodal_points ( ix, 4, :) = nodal_points ( ix, 4, :) * dtmove 

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------

! STEP 1) Having come into the subroutine with variables loaded into "w",
! load the appropriate spatial variation of "w" into "u" while setting
! the temporal variation of "u" as well as "f,g,h & s" to zero initially.

      DO ix = oned_begin - oned_off, oned_end + oned_off

      u ( ix, :, 1: NUM_ZONECEN_BASES)=w ( ix, :, 1: NUM_ZONECEN_BASES)
      u ( ix, :, NUM_ZONECEN_BASES + 1: NUM_SPACETIME_BASES) = 0.0

      f ( ix, :, :) = 0.0
      g ( ix, :, :) = 0.0
      h ( ix, :, :) = 0.0
      s ( ix, :, :) = 0.0

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------

! STEP 2) Use the spatial modes to make
! "u_quadpt ( ix, :, 1: NUM_ZONECEN_QUADPTS)"
! from "u ( ix, :, 1: NUM_ZONECEN_BASES)" at time t = 0.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

      u_quadpt ( ix, :, 01) = u ( ix, :, 01) + 0.5 * u ( ix, :, 02)

      u_quadpt ( ix, :, 02) = u ( ix, :, 01) - 0.5 * u ( ix, :, 02)

      u_quadpt ( ix, :, 03) = u ( ix, :, 01) + 0.5 * u ( ix, :, 03)

      u_quadpt ( ix, :, 04) = u ( ix, :, 01) - 0.5 * u ( ix, :, 03)

      u_quadpt ( ix, :, 05) = u ( ix, :, 01) + 0.5 * u ( ix, :, 04)

      u_quadpt ( ix, :, 06) = u ( ix, :, 01) - 0.5 * u ( ix, :, 04)

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

      u_quadpt ( ix, :, 01) = u ( ix, :, 01) - ( u ( ix, :, 05)
     1 + u ( ix, :, 06) + u ( ix, :, 07) ) / 12.0

      u_quadpt ( ix, :, 02) = u_quadpt ( ix, :, 01)
     1 + 0.5 * ( u ( ix, :, 02) + u ( ix, :, 03) )
     1 + 0.25 * ( u ( ix, :, 05) + u ( ix, :, 06) + u ( ix, :, 08) )

      u_quadpt ( ix, :, 03) = u_quadpt ( ix, :, 02) - u ( ix, :, 03)
     1 - 0.5 * u ( ix, :, 08)

      u_quadpt ( ix, :, 04) = u_quadpt ( ix, :, 01)
     1 + 0.5 * ( u ( ix, :, 02) + u ( ix, :, 04) )
     1 + 0.25 * ( u ( ix, :, 05) + u ( ix, :, 07) + u ( ix, :, 10) )

      u_quadpt ( ix, :, 05) = u_quadpt ( ix, :, 04) - u ( ix, :, 04)
     1 - 0.5 * u ( ix, :, 10)

      u_quadpt ( ix, :, 06) = u_quadpt ( ix, :, 03) - u ( ix, :, 02)
     1 + u ( ix, :, 03)

      u_quadpt ( ix, :, 07) = u_quadpt ( ix, :, 02) - u ( ix, :, 02)
     1 - u ( ix, :, 03)

      u_quadpt ( ix, :, 08) = u_quadpt ( ix, :, 05) - u ( ix, :, 02)
     1 + u ( ix, :, 04)

      u_quadpt ( ix, :, 09) = u_quadpt ( ix, :, 04) - u ( ix, :, 02)
     1 - u ( ix, :, 04)

      u_quadpt ( ix, :, 10) = u_quadpt ( ix, :, 01)
     1 + 0.5 * ( u ( ix, :, 03) + u ( ix, :, 04) )
     1 + 0.25 * ( u ( ix, :, 06) + u ( ix, :, 07) + u ( ix, :, 09) )

      u_quadpt ( ix, :, 11) = u_quadpt ( ix, :, 10) - u ( ix, :, 03)
     1 - 0.5 * u ( ix, :, 09)

      u_quadpt ( ix, :, 12) = u_quadpt ( ix, :, 11) + u ( ix, :, 03)
     1 - u ( ix, :, 04)

      u_quadpt ( ix, :, 13) = u_quadpt ( ix, :, 10) - u ( ix, :, 03)
     1 - u ( ix, :, 04)

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

      u_quadpt ( ix, :, 01) = u ( ix, :, 01) - ( u ( ix, :, 05)
     1 + u ( ix, :, 06) + u ( ix, :, 07) ) / 12.0

      u_quadpt ( ix, :, 02) = u_quadpt ( ix, :, 01)
     1 + ( 60.0 * u ( ix, :, 02) + 30.0 * u ( ix, :, 05)
     1 + 6.0 * u ( ix, :, 11) - 5.0 * u ( ix, :, 16)
     1 - 5.0 * u ( ix, :, 18) ) / 120.0

      u_quadpt ( ix, :, 03) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 02) ) - 0.015625 * ( 4 * u ( ix, :, 05)
     1 + 3.0 * u ( ix, :, 11) )

      u_quadpt ( ix, :, 04) = 3.0 * u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 02) - 3.0 * u_quadpt ( ix, :, 03)
     1 - 0.09375 * u ( ix, :, 11)

      u_quadpt ( ix, :, 05) = 6.0 * u_quadpt ( ix, :, 01)
     1 + 3.0 * u_quadpt ( ix, :, 02) - 8.0 * u_quadpt ( ix, :, 03)
     1 - 0.375 * u ( ix, :, 11)

      u_quadpt ( ix, :, 06) = u_quadpt ( ix, :, 01)
     1 + ( 60.0 * u ( ix, :, 03) + 30.0 * u ( ix, :, 06)
     1 + 6.0 * u ( ix, :, 12) - 5.0 * u ( ix, :, 14)
     1 - 5.0 * u ( ix, :, 19) ) / 120.0

      u_quadpt ( ix, :, 07) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 06) ) - 0.015625 * ( 4.0 * u ( ix, :, 06)
     1 + 3.0 * u ( ix, :, 12) )

      u_quadpt ( ix, :, 08) = 3.0 * u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 06) - 3.0 * u_quadpt ( ix, :, 07)
     1 - 0.09375 * u ( ix, :, 12)

      u_quadpt ( ix, :, 09) = 6.0 * u_quadpt ( ix, :, 01)
     1 + 3.0 * u_quadpt ( ix, :, 06) - 8.0 * u_quadpt ( ix, :, 07)
     1 - 0.375 * u ( ix, :, 12)

      u_quadpt ( ix, :, 10) = u_quadpt ( ix, :, 01)
     1 + ( 60.0 * u ( ix, :, 04) + 30.0 * u ( ix, :, 07)
     1 + 6.0 * u ( ix, :, 13) - 5.0 * u ( ix, :, 15)
     1 - 5.0 * u ( ix, :, 17) ) / 120.0

      u_quadpt ( ix, :, 11) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 10) ) - 0.015625 * ( 4.0 * u ( ix, :, 07)
     1 + 3.0 * u ( ix, :, 13) )

      u_quadpt ( ix, :, 12) = 3.0 * u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 10) - 3.0 * u_quadpt ( ix, :, 11)
     1 - 0.09375 * u ( ix, :, 13)

      u_quadpt ( ix, :, 13) = 6.0 * u_quadpt ( ix, :, 01)
     1 + 3.0 * u_quadpt ( ix, :, 10) - 8.0 * u_quadpt ( ix, :, 11)
     1 - 0.375 * u ( ix, :, 13)

      u_quadpt ( ix, :, 14) = u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 07) + u_quadpt ( ix, :, 10)
     1 - 2.0 * u_quadpt ( ix, :, 01) + 0.25 * ( u ( ix, :, 03)
     1 + u ( ix, :, 08) + u ( ix, :, 09) + u ( ix, :, 10) )
     1 + 0.125 * ( u ( ix, :, 15) + u ( ix, :, 16) + u ( ix, :, 17)
     1 + u ( ix, :, 18) + u ( ix, :, 20) ) + 0.1875 * u ( ix, :, 06)
     1 + 5.0 * ( u ( ix, :, 14) + u ( ix, :, 19) ) / 48.0
     1 + 0.071875 * u ( ix, :, 12)

      u_quadpt ( ix, :, 15) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 09) - u_quadpt ( ix, :, 06)
     1 - 0.5 * ( u ( ix, :, 08) + u ( ix, :, 09) )
     1 - 0.25 * ( u ( ix, :, 14) + u ( ix, :, 19) + u ( ix, :, 20) )

      u_quadpt ( ix, :, 16) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 13) - u_quadpt ( ix, :, 10)
     1 - 0.5 * ( u ( ix, :, 09) + u ( ix, :, 10) )
     1 - 0.25 * ( u ( ix, :, 15) + u ( ix, :, 17) + u ( ix, :, 20) )

      u_quadpt ( ix, :, 17) = u_quadpt ( ix, :, 15)
     1 + u_quadpt ( ix, :, 16) - u_quadpt ( ix, :, 14) + u ( ix, :, 09)
     1 + 0.5 * u ( ix, :, 20)

      u_quadpt ( ix, :, 18) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 02)
     1 - 0.5 * ( u ( ix, :, 08) + u ( ix, :, 10) )
     1 - 0.25 * ( u ( ix, :, 16) + u ( ix, :, 18) + u ( ix, :, 20) )

      u_quadpt ( ix, :, 19) = u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 09) - u_quadpt ( ix, :, 06)
     1 + 0.5 * ( u ( ix, :, 08) - u ( ix, :, 09) )
     1 - 0.25 * ( u ( ix, :, 14) + u ( ix, :, 19) - u ( ix, :, 20) )

      u_quadpt ( ix, :, 20) = u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 13) - u_quadpt ( ix, :, 10)
     1 - 0.5 * ( u ( ix, :, 09) - u ( ix, :, 10) )
     1 - 0.25 * ( u ( ix, :, 15) + u ( ix, :, 17) - u ( ix, :, 20) )

      u_quadpt ( ix, :, 21) = u_quadpt ( ix, :, 19)
     1 + u_quadpt ( ix, :, 20) - u_quadpt ( ix, :, 18) + u ( ix, :, 09)
     1 - 0.5 * u ( ix, :, 20)

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------

! STEP 3) Get the primitive variables from the conserved variables at
! the spatial nodes at t=0. The density and pressure variables are
! given positive floor values to ensure that they can be used for physically
! meaningful source terms.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, 1, iquadpt) =
     1                   AMAX1 ( u_quadpt ( ix, 1, iquadpt), rhofloor)

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempaa = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempaa
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempaa
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempaa

      bfx ( ix, iquadpt) = u_quadpt ( ix, 6, iquadpt)
      bfy ( ix, iquadpt) = u_quadpt ( ix, 7, iquadpt)
      bfz ( ix, iquadpt) = u_quadpt ( ix, 8, iquadpt)

      tempab = 0.5 * rho ( ix, iquadpt) * ( vlx ( ix, iquadpt)**2
     1         + vly ( ix, iquadpt)**2 + vlz ( ix, iquadpt)**2)

      tempac = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      eng_t ( ix, iquadpt) = u_quadpt ( ix, 5, iquadpt)
     1                     - tempab - tempac

      prs ( ix, iquadpt) = eng_t ( ix, iquadpt)
     1                   * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

#if ( ENFORCE_POSITIVITY == 1)
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = entropy ( ix, iquadpt)
     1* rho ( ix, iquadpt)**( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)
      prs ( ix, iquadpt) = AMAX1 ( prs ( ix, iquadpt), prsfloor)

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#else
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = prsfloor

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#endif

#if ( NFLUID >= 1)
      specie ( ix, 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( ix, 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( ix, 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( ix, 1: NFLUID, iquadpt) * rho ( ix, iquadpt)
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, ivar + 1, iquadpt) =
     1         AMAX1 ( u_quadpt ( ix, ivar + 1, iquadpt), rhofloor)

      rho_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 1, iquadpt)

      tempaa = 1.0 / rho_euler ( ix, iquadpt, i_euler)

      vlx_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 2, iquadpt) * tempaa
      vly_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 3, iquadpt) * tempaa
      vlz_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 4, iquadpt) * tempaa

      tempab = 0.5 * rho_euler ( ix, iquadpt, i_euler)
     1           * ( vlx_euler ( ix, iquadpt, i_euler)**2
     1             + vly_euler ( ix, iquadpt, i_euler)**2
     1             + vlz_euler ( ix, iquadpt, i_euler)**2)

      prs_euler ( ix, iquadpt, i_euler) =
     1         ( u_quadpt ( ix, ivar + 5, iquadpt) - tempab)
     1       * ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      prs_euler ( ix, iquadpt, i_euler) =
     1    AMAX1 ( prs_euler ( ix, iquadpt, i_euler), prsfloor)

      eng_t_euler ( ix, iquadpt, i_euler) =
     1         prs_euler ( ix, iquadpt, i_euler)
     1       / ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      IF ( prs_euler ( ix, iquadpt, i_euler) .LT. prsfloor * 1.00000001)
     1   u_quadpt ( ix, ivar + 5, iquadpt) =
     1                   eng_t_euler ( ix, iquadpt, i_euler) + tempab

#if ( NFLUID_EULER >= 1)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) = tempaa
     1       * u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1                    ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) =
     1    AMIN1 ( AMAX1 (
     1    specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler),
     1            0.0), 1.0 )

      u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1       specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler)
     1     * rho_euler ( ix, iquadpt, i_euler)

#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 4.a) Build fluxes ( f_quadpt, g_quadpt, h_quadpt) from the primitive
! variables that have now become available at spatial nodes at t=0.

! Multiply "f_quadpt" by "dtmove/dx" and so on.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! -----

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

! -----

      f_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 2, iquadpt)

      f_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfx ( ix, iquadpt)**2 * ob4pi

      f_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vly ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi

      f_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi

      f_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab

      f_quadpt ( ix, 6, iquadpt) = 0.0

      f_quadpt ( ix, 7, iquadpt) =
     1           vlx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         - vly ( ix, iquadpt) * bfx ( ix, iquadpt)

      f_quadpt ( ix, 8, iquadpt) =
     1           vlx ( ix, iquadpt) * bfz ( ix, iquadpt)
     1         - vlz ( ix, iquadpt) * bfx ( ix, iquadpt)

#if ( NFLUID >= 1)
      f_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlx ( ix, iquadpt)
#endif

! -----

      g_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 3, iquadpt)

      g_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         - bfy ( ix, iquadpt) * bfx ( ix, iquadpt) * ob4pi

      g_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vly ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfy ( ix, iquadpt)**2 * ob4pi

      g_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         - bfy ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi

      g_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vly ( ix, iquadpt)
     1                       - bfy ( ix, iquadpt) * tempab

      g_quadpt ( ix, 6, iquadpt) =
     1           vly ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         - vlx ( ix, iquadpt) * bfy ( ix, iquadpt)

      g_quadpt ( ix, 7, iquadpt) = 0.0

      g_quadpt ( ix, 8, iquadpt) =
     1           vly ( ix, iquadpt) * bfz ( ix, iquadpt)
     1         - vlz ( ix, iquadpt) * bfy ( ix, iquadpt)

#if ( NFLUID >= 1)
      g_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vly ( ix, iquadpt)
#endif

! -----

      h_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 4, iquadpt)

      h_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         - bfz ( ix, iquadpt) * bfx ( ix, iquadpt) * ob4pi

      h_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vly ( ix, iquadpt)
     1         - bfz ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi

      h_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfz ( ix, iquadpt)**2 * ob4pi

      h_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vlz ( ix, iquadpt)
     1                       - bfz ( ix, iquadpt) * tempab

      h_quadpt ( ix, 6, iquadpt) =
     1           vlz ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         - vlx ( ix, iquadpt) * bfz ( ix, iquadpt)

      h_quadpt ( ix, 7, iquadpt) =
     1           vlz ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         - vly ( ix, iquadpt) * bfz ( ix, iquadpt)

      h_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      h_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlz ( ix, iquadpt)
#endif

! -----

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! -----

      f_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)

      f_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vlx_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      f_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vlx_euler ( ix, iquadpt, i_euler)
#endif

! -----

      g_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)

      g_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vly_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      g_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vly_euler ( ix, iquadpt, i_euler)
#endif

! -----

      h_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)

      h_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vlz_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      h_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vlz_euler ( ix, iquadpt, i_euler)
#endif

! -----

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ------------------------------

! To make cylindrical and spherical geometries accessible to non-rel.
! and relativistic flow, we put this part after either flux construction.

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      IF ( igeom == 2 .or. igeom == 3) THEN

         g_quadpt ( ix, :, iquadpt) = g_quadpt ( ix, :, iquadpt) * obx1

      END IF

      IF ( igeom == 3) THEN

         sin_x2 = 1.0 / ( x1 * sin ( x2))

         h_quadpt ( ix, :, iquadpt) = h_quadpt ( ix, :, iquadpt)
     1                              * sin_x2

      END IF

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 4.b) Optional step needed only for viscous, resistive or conductive MHD.

! ASSUME that "du_quadpt_dx, du_quadpt_dy and du_quadpt_dz" have been built at
! all orders (NOT DONE YET).
! We set about building the derivatives of the primitive
! variables and various auxiliary variables at the nodal points. These are
! needed only in the presence of viscous, resistive or conduction terms.

! These derivatives are then used to build non-ideal fluxes at the quadrature
! points given by "fni_quadpt, gni_quadpt and hni_quadpt".

! Lastly, the non-ideal fluxes are added into the ideal fluxes

! ----------------------------------------------------------------------

      IF ( ( ( viscosity_coef .GT. 0.0) .OR.
     1       ( resistivity_coef .GT. 0.0)) .AND. ( NUMORDER >= 3)) THEN

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = 1.0 / rho ( ix, iquadpt)

      drho_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 1, iquadpt)
      drho_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 1, iquadpt)
      drho_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 1, iquadpt)


      dvlx_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvlx_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvlx_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dvly_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvly_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvly_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dvlz_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvlz_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvlz_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dbfx_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 6, iquadpt)
      dbfx_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 6, iquadpt)
      dbfx_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 6, iquadpt)

      dbfy_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 7, iquadpt)
      dbfy_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 7, iquadpt)
      dbfy_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 7, iquadpt)

      dbfz_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 8, iquadpt)
      dbfz_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 8, iquadpt)
      dbfz_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 8, iquadpt)


      dv2_dx ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dx ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dx ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dx ( ix, iquadpt) )

      dv2_dy ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dy ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dy ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dy ( ix, iquadpt) )

      dv2_dz ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dz ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dz ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dz ( ix, iquadpt) )

      db2_dx ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dx ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dx ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dx ( ix, iquadpt) )

      db2_dy ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dy ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dy ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dy ( ix, iquadpt) )

      db2_dz ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dz ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dz ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dz ( ix, iquadpt) )


      tempac = vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1       + vlz ( ix, iquadpt)**2

      tempab = molewt_quadpt ( ix, 1, iquadpt)
     1       / ( gasconst * rho ( ix, iquadpt) )

      temperature ( ix, iquadpt) = prs ( ix, iquadpt) * tempab

      tempab = tempab * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      dtemperature_dx ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dx ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dx ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dx ( ix, iquadpt)
     1    - ob8pi * db2_dx ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dx ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

      dtemperature_dy ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dy ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dy ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dy ( ix, iquadpt)
     1    - ob8pi * db2_dy ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dy ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

      dtemperature_dz ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dz ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dz ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dz ( ix, iquadpt)
     1    - ob8pi * db2_dz ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dz ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

! The derivatives are within the reference element's local coordinates.
! To get true physical fluxes, we need to divide by the size of the zone here.

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = 1.0 / dx ( ix)
      tempab = 1.0 / dy ( ix)
      tempac = 1.0 / dz ( ix)

      drho_dx ( ix, iquadpt) = tempaa * drho_dx ( ix, iquadpt)
      dtemperature_dx ( ix, iquadpt) =
     1                 tempaa * dtemperature_dx ( ix, iquadpt)
      dvlx_dx ( ix, iquadpt) = tempaa * dvlx_dx ( ix, iquadpt)
      dvly_dx ( ix, iquadpt) = tempaa * dvly_dx ( ix, iquadpt)
      dvlz_dx ( ix, iquadpt) = tempaa * dvlz_dx ( ix, iquadpt)
      dbfx_dx ( ix, iquadpt) = tempaa * dbfx_dx ( ix, iquadpt)
      dbfy_dx ( ix, iquadpt) = tempaa * dbfy_dx ( ix, iquadpt)
      dbfz_dx ( ix, iquadpt) = tempaa * dbfz_dx ( ix, iquadpt)
      dv2_dx ( ix, iquadpt) = tempaa * dv2_dx ( ix, iquadpt)
      db2_dx ( ix, iquadpt) = tempaa * db2_dx ( ix, iquadpt)

      drho_dy ( ix, iquadpt) = tempab * drho_dy ( ix, iquadpt)
      dtemperature_dy ( ix, iquadpt) =
     1                 tempab * dtemperature_dy ( ix, iquadpt)
      dvlx_dy ( ix, iquadpt) = tempab * dvlx_dy ( ix, iquadpt)
      dvly_dy ( ix, iquadpt) = tempab * dvly_dy ( ix, iquadpt)
      dvlz_dy ( ix, iquadpt) = tempab * dvlz_dy ( ix, iquadpt)
      dbfx_dy ( ix, iquadpt) = tempab * dbfx_dy ( ix, iquadpt)
      dbfy_dy ( ix, iquadpt) = tempab * dbfy_dy ( ix, iquadpt)
      dbfz_dy ( ix, iquadpt) = tempab * dbfz_dy ( ix, iquadpt)
      dv2_dy ( ix, iquadpt) = tempab * dv2_dy ( ix, iquadpt)
      db2_dy ( ix, iquadpt) = tempab * db2_dy ( ix, iquadpt)

      drho_dz ( ix, iquadpt) = tempac * drho_dz ( ix, iquadpt)
      dtemperature_dz ( ix, iquadpt) =
     1                 tempac * dtemperature_dz ( ix, iquadpt)
      dvlx_dz ( ix, iquadpt) = tempac * dvlx_dz ( ix, iquadpt)
      dvly_dz ( ix, iquadpt) = tempac * dvly_dz ( ix, iquadpt)
      dvlz_dz ( ix, iquadpt) = tempac * dvlz_dz ( ix, iquadpt)
      dbfx_dz ( ix, iquadpt) = tempac * dbfx_dz ( ix, iquadpt)
      dbfy_dz ( ix, iquadpt) = tempac * dbfy_dz ( ix, iquadpt)
      dbfz_dz ( ix, iquadpt) = tempac * dbfz_dz ( ix, iquadpt)
      dv2_dz ( ix, iquadpt) = tempac * dv2_dz ( ix, iquadpt)
      db2_dz ( ix, iquadpt) = tempac * db2_dz ( ix, iquadpt)

      div_vl ( ix, iquadpt) = dvlx_dx ( ix, iquadpt)
     1                      + dvly_dy ( ix, iquadpt)
     1                      + dvlz_dz ( ix, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = rho ( ix, iquadpt) * viscosity_coef

! ----------

      fni_quadpt ( ix, 1, iquadpt) = 0.0

      fni_quadpt ( ix, 2, iquadpt) =
     1               2.0 * tempaa * ( dvlx_dx ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      fni_quadpt ( ix, 3, iquadpt) =
     1                     tempaa * ( dvlx_dy ( ix, iquadpt)
     1                              + dvly_dx ( ix, iquadpt) )

      fni_quadpt ( ix, 4, iquadpt) =
     1                     tempaa * ( dvlx_dz ( ix, iquadpt)
     1                              + dvlz_dx ( ix, iquadpt) )

      fni_quadpt ( ix, 5, iquadpt) =
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vlx ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvlx_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvlx_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvlx_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dx ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dx ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfx_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfx_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfx_dz ( ix, iquadpt) )

      fni_quadpt ( ix, 6, iquadpt) = 0.0

      fni_quadpt ( ix, 7, iquadpt) =
     1   resistivity_coef * ( dbfy_dx ( ix, iquadpt)
     1                      - dbfx_dy ( ix, iquadpt) )

      fni_quadpt ( ix, 8, iquadpt) =
     1   resistivity_coef * ( dbfz_dx ( ix, iquadpt)
     1                      - dbfx_dz ( ix, iquadpt) )

#if ( NVAR > 8)
      fni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      gni_quadpt ( ix, 1, iquadpt) = 0.0

      gni_quadpt ( ix, 2, iquadpt) =
     1                     tempaa * ( dvly_dx ( ix, iquadpt)
     1                              + dvlx_dy ( ix, iquadpt) )

      gni_quadpt ( ix, 3, iquadpt) =
     1               2.0 * tempaa * ( dvly_dy ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      gni_quadpt ( ix, 4, iquadpt) =
     1                     tempaa * ( dvly_dz ( ix, iquadpt)
     1                              + dvlz_dy ( ix, iquadpt) )

      gni_quadpt ( ix, 5, iquadpt) = 
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vly ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvly_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvly_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvly_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dy ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dy ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfy_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfy_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfy_dz ( ix, iquadpt) )

      gni_quadpt ( ix, 6, iquadpt) =
     1   resistivity_coef * ( dbfx_dy ( ix, iquadpt)
     1                      - dbfy_dx ( ix, iquadpt) )

      gni_quadpt ( ix, 7, iquadpt) = 0.0

      gni_quadpt ( ix, 8, iquadpt) =
     1   resistivity_coef * ( dbfz_dy ( ix, iquadpt)
     1                      - dbfy_dz ( ix, iquadpt) )

#if ( NVAR > 8)
      gni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      hni_quadpt ( ix, 1, iquadpt) = 0.0

      hni_quadpt ( ix, 2, iquadpt) =
     1                     tempaa * ( dvlz_dx ( ix, iquadpt)
     1                              + dvlx_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 3, iquadpt) =
     1                     tempaa * ( dvlz_dy ( ix, iquadpt)
     1                              + dvly_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 4, iquadpt) =
     1               2.0 * tempaa * ( dvlz_dz ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      hni_quadpt ( ix, 5, iquadpt) =
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vlz ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvlz_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvlz_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvlz_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dz ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dz ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfz_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfz_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfz_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 6, iquadpt) =
     1   resistivity_coef * ( dbfx_dz ( ix, iquadpt)
     1                      - dbfz_dx ( ix, iquadpt) )

      hni_quadpt ( ix, 7, iquadpt) =
     1   resistivity_coef * ( dbfy_dz ( ix, iquadpt)
     1                      - dbfz_dy ( ix, iquadpt) )

      hni_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NVAR > 8)
      hni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      f_quadpt ( ix, :, iquadpt) = f_quadpt ( ix, :, iquadpt)
     1                       + fni_quadpt ( ix, :, iquadpt)

      g_quadpt ( ix, :, iquadpt) = g_quadpt ( ix, :, iquadpt)
     1                       + gni_quadpt ( ix, :, iquadpt)

      h_quadpt ( ix, :, iquadpt) = h_quadpt ( ix, :, iquadpt)
     1                       + hni_quadpt ( ix, :, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      END IF
! End of "( viscosity_coef .GT. 0.0) .OR. ( resistivity_coef .GT. 0.0)".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = dtmove / dx ( ix)
      f_quadpt ( ix, :, 1: NUM_ZONECEN_QUADPTS) = tempaa
     1         * f_quadpt ( ix, :, 1: NUM_ZONECEN_QUADPTS)

      tempab = dtmove / dy ( ix)
      g_quadpt ( ix, :, 1: NUM_ZONECEN_QUADPTS) = tempab
     1         * g_quadpt ( ix, :, 1: NUM_ZONECEN_QUADPTS)

      tempac = dtmove / dz ( ix)
      h_quadpt ( ix, :, 1: NUM_ZONECEN_QUADPTS) = tempac
     1         * h_quadpt ( ix, :, 1: NUM_ZONECEN_QUADPTS)

! ----------

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 4.c) Now project the fluxes "( f_quadpt, g_quadpt, h_quadpt)"
! from the spatial nodes at t=0 back to the t=0 part of the modal
! space time variables for those variables, i.e. to "( f, g, h)" .

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

!      f ( ix, :, 01) = ( f_quadpt ( ix, :, 01) + f_quadpt ( ix, :, 02)
!     1 + f_quadpt ( ix, :, 03) + f_quadpt ( ix, :, 04)
!     1 + f_quadpt ( ix, :, 05) + f_quadpt ( ix, :, 06)) / 6.0

      f ( ix, :, 02) = f_quadpt ( ix, :, 01) - f_quadpt ( ix, :, 02)

!      f ( ix, :, 03) = f_quadpt ( ix, :, 03) - f_quadpt ( ix, :, 04)

!      f ( ix, :, 04) = f_quadpt ( ix, :, 05) - f_quadpt ( ix, :, 06)

! ----------

!      g ( ix, :, 01) = ( g_quadpt ( ix, :, 01) + g_quadpt ( ix, :, 02)
!     1 + g_quadpt ( ix, :, 03) + g_quadpt ( ix, :, 04)
!     1 + g_quadpt ( ix, :, 05) + g_quadpt ( ix, :, 06)) / 6.0

!      g ( ix, :, 02) = g_quadpt ( ix, :, 01) - g_quadpt ( ix, :, 02)

      g ( ix, :, 03) = g_quadpt ( ix, :, 03) - g_quadpt ( ix, :, 04)

!      g ( ix, :, 04) = g_quadpt ( ix, :, 05) - g_quadpt ( ix, :, 06)

! ----------

!      h ( ix, :, 01) = ( h_quadpt ( ix, :, 01) + h_quadpt ( ix, :, 02)
!     1 + h_quadpt ( ix, :, 03) + h_quadpt ( ix, :, 04)
!     1 + h_quadpt ( ix, :, 05) + h_quadpt ( ix, :, 06)) / 6.0

!      h ( ix, :, 02) = h_quadpt ( ix, :, 01) - h_quadpt ( ix, :, 02)

!      h ( ix, :, 03) = h_quadpt ( ix, :, 03) - h_quadpt ( ix, :, 04)

      h ( ix, :, 04) = h_quadpt ( ix, :, 05) - h_quadpt ( ix, :, 06)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tf ( ix, :, 01) = f_quadpt ( ix, :, 02) + f_quadpt ( ix, :, 03)
     1           + f_quadpt ( ix, :, 04) + f_quadpt ( ix, :, 05)

      tf ( ix, :, 02) = f_quadpt ( ix, :, 06) + f_quadpt ( ix, :, 07)
     1           + f_quadpt ( ix, :, 08) + f_quadpt ( ix, :, 09)

      tf ( ix, :, 03) = f_quadpt ( ix, :, 10) + f_quadpt ( ix, :, 11)
     1           + f_quadpt ( ix, :, 12) + f_quadpt ( ix, :, 13)

!      tf ( ix, :, 04) = f_quadpt ( ix, :, 02) + f_quadpt ( ix, :, 06)
!     1           + f_quadpt ( ix, :, 10) + f_quadpt ( ix, :, 12)

!      tf ( ix, :, 05) = f_quadpt ( ix, :, 03) + f_quadpt ( ix, :, 07)
!     1           + f_quadpt ( ix, :, 11) + f_quadpt ( ix, :, 13)

!      tf ( ix, :, 06) = f_quadpt ( ix, :, 04) + f_quadpt ( ix, :, 05)
!     1           + f_quadpt ( ix, :, 08) + f_quadpt ( ix, :, 09)

!      tf ( ix, :, 07) = f_quadpt ( ix, :, 04) + f_quadpt ( ix, :, 08)
!     1           + f_quadpt ( ix, :, 10) + f_quadpt ( ix, :, 11)

!      tf ( ix, :, 08) = f_quadpt ( ix, :, 05) + f_quadpt ( ix, :, 09)
!     1           + f_quadpt ( ix, :, 12) + f_quadpt ( ix, :, 13)

!      tf ( ix, :, 09) = f_quadpt ( ix, :, 02) + f_quadpt ( ix, :, 03)
!     1           + f_quadpt ( ix, :, 06) + f_quadpt ( ix, :, 07)

      f ( ix, :, 02) = 0.25 * ( tf ( ix, :, 01) - tf ( ix, :, 02) )

!      f ( ix, :, 03) = 0.25 * ( tf ( ix, :, 04) - tf ( ix, :, 05) )

!      f ( ix, :, 04) = 0.25 * ( tf ( ix, :, 07) - tf ( ix, :, 08) )

      f ( ix, :, 05) = 0.5 * ( tf ( ix, :, 01) + tf ( ix, :, 02)
     1 - tf ( ix, :, 03) ) - 2.0 * f_quadpt ( ix, :, 01)

!      f ( ix, :, 06) = 0.5 * ( tf ( ix, :, 04) + tf ( ix, :, 05)
!     1 - tf ( ix, :, 06) ) - 2.0 * f_quadpt ( ix, :, 01)

!      f ( ix, :, 07) = 0.5 * ( tf ( ix, :, 07) + tf ( ix, :, 08)
!     1 - tf ( ix, :, 09) ) - 2.0 * f_quadpt ( ix, :, 01)

      f ( ix, :, 08) = f_quadpt ( ix, :, 02) - f_quadpt ( ix, :, 03)
     1 - f_quadpt ( ix, :, 06) + f_quadpt ( ix, :, 07)

!      f ( ix, :, 09) = f_quadpt ( ix, :, 10) - f_quadpt ( ix, :, 11)
!     1 - f_quadpt ( ix, :, 12) + f_quadpt ( ix, :, 13)

      f ( ix, :, 10) = f_quadpt ( ix, :, 04) - f_quadpt ( ix, :, 05)
     1 - f_quadpt ( ix, :, 08) + f_quadpt ( ix, :, 09)

!      f ( ix, :, 01) = f_quadpt ( ix, :, 01) + ( f ( ix, :, 05)
!     1 + f ( ix, :, 06) + f ( ix, :, 07) ) / 12.0

! ----------

!      tg ( ix, :, 01) = g_quadpt ( ix, :, 02) + g_quadpt ( ix, :, 03)
!     1           + g_quadpt ( ix, :, 04) + g_quadpt ( ix, :, 05)

!      tg ( ix, :, 02) = g_quadpt ( ix, :, 06) + g_quadpt ( ix, :, 07)
!     1           + g_quadpt ( ix, :, 08) + g_quadpt ( ix, :, 09)

!      tg ( ix, :, 03) = g_quadpt ( ix, :, 10) + g_quadpt ( ix, :, 11)
!     1           + g_quadpt ( ix, :, 12) + g_quadpt ( ix, :, 13)

      tg ( ix, :, 04) = g_quadpt ( ix, :, 02) + g_quadpt ( ix, :, 06)
     1           + g_quadpt ( ix, :, 10) + g_quadpt ( ix, :, 12)

      tg ( ix, :, 05) = g_quadpt ( ix, :, 03) + g_quadpt ( ix, :, 07)
     1           + g_quadpt ( ix, :, 11) + g_quadpt ( ix, :, 13)

      tg ( ix, :, 06) = g_quadpt ( ix, :, 04) + g_quadpt ( ix, :, 05)
     1           + g_quadpt ( ix, :, 08) + g_quadpt ( ix, :, 09)

!      tg ( ix, :, 07) = g_quadpt ( ix, :, 04) + g_quadpt ( ix, :, 08)
!     1           + g_quadpt ( ix, :, 10) + g_quadpt ( ix, :, 11)

!      tg ( ix, :, 08) = g_quadpt ( ix, :, 05) + g_quadpt ( ix, :, 09)
!     1           + g_quadpt ( ix, :, 12) + g_quadpt ( ix, :, 13)

!      tg ( ix, :, 09) = g_quadpt ( ix, :, 02) + g_quadpt ( ix, :, 03)
!     1           + g_quadpt ( ix, :, 06) + g_quadpt ( ix, :, 07)

!      g ( ix, :, 02) = 0.25 * ( tg ( ix, :, 01) - tg ( ix, :, 02) )

      g ( ix, :, 03) = 0.25 * ( tg ( ix, :, 04) - tg ( ix, :, 05) )

!      g ( ix, :, 04) = 0.25 * ( tg ( ix, :, 07) - tg ( ix, :, 08) )

!      g ( ix, :, 05) = 0.5 * ( tg ( ix, :, 01) + tg ( ix, :, 02)
!     1 - tg ( ix, :, 03) ) - 2.0 * g_quadpt ( ix, :, 01)

      g ( ix, :, 06) = 0.5 * ( tg ( ix, :, 04) + tg ( ix, :, 05)
     1 - tg ( ix, :, 06) ) - 2.0 * g_quadpt ( ix, :, 01)

!      g ( ix, :, 07) = 0.5 * ( tg ( ix, :, 07) + tg ( ix, :, 08)
!     1 - tg ( ix, :, 09) ) - 2.0 * g_quadpt ( ix, :, 01)

      g ( ix, :, 08) = g_quadpt ( ix, :, 02) - g_quadpt ( ix, :, 03)
     1 - g_quadpt ( ix, :, 06) + g_quadpt ( ix, :, 07)

      g ( ix, :, 09) = g_quadpt ( ix, :, 10) - g_quadpt ( ix, :, 11)
     1 - g_quadpt ( ix, :, 12) + g_quadpt ( ix, :, 13)

!      g ( ix, :, 10) = g_quadpt ( ix, :, 04) - g_quadpt ( ix, :, 05)
!     1 - g_quadpt ( ix, :, 08) + g_quadpt ( ix, :, 09)

!      g ( ix, :, 01) = g_quadpt ( ix, :, 01) + ( g ( ix, :, 05)
!     1 + g ( ix, :, 06) + g ( ix, :, 07) ) / 12.0

! ----------

!      th ( ix, :, 01) = h_quadpt ( ix, :, 02) + h_quadpt ( ix, :, 03)
!     1           + h_quadpt ( ix, :, 04) + h_quadpt ( ix, :, 05)

!      th ( ix, :, 02) = h_quadpt ( ix, :, 06) + h_quadpt ( ix, :, 07)
!     1           + h_quadpt ( ix, :, 08) + h_quadpt ( ix, :, 09)

!      th ( ix, :, 03) = h_quadpt ( ix, :, 10) + h_quadpt ( ix, :, 11)
!     1           + h_quadpt ( ix, :, 12) + h_quadpt ( ix, :, 13)

!      th ( ix, :, 04) = h_quadpt ( ix, :, 02) + h_quadpt ( ix, :, 06)
!     1           + h_quadpt ( ix, :, 10) + h_quadpt ( ix, :, 12)

!      th ( ix, :, 05) = h_quadpt ( ix, :, 03) + h_quadpt ( ix, :, 07)
!     1           + h_quadpt ( ix, :, 11) + h_quadpt ( ix, :, 13)

!      th ( ix, :, 06) = h_quadpt ( ix, :, 04) + h_quadpt ( ix, :, 05)
!     1           + h_quadpt ( ix, :, 08) + h_quadpt ( ix, :, 09)

      th ( ix, :, 07) = h_quadpt ( ix, :, 04) + h_quadpt ( ix, :, 08)
     1           + h_quadpt ( ix, :, 10) + h_quadpt ( ix, :, 11)

      th ( ix, :, 08) = h_quadpt ( ix, :, 05) + h_quadpt ( ix, :, 09)
     1           + h_quadpt ( ix, :, 12) + h_quadpt ( ix, :, 13)

      th ( ix, :, 09) = h_quadpt ( ix, :, 02) + h_quadpt ( ix, :, 03)
     1           + h_quadpt ( ix, :, 06) + h_quadpt ( ix, :, 07)

!      h ( ix, :, 02) = 0.25 * ( th ( ix, :, 01) - th ( ix, :, 02) )

!      h ( ix, :, 03) = 0.25 * ( th ( ix, :, 04) - th ( ix, :, 05) )

      h ( ix, :, 04) = 0.25 * ( th ( ix, :, 07) - th ( ix, :, 08) )

!      h ( ix, :, 05) = 0.5 * ( th ( ix, :, 01) + th ( ix, :, 02)
!     1 - th ( ix, :, 03) ) - 2.0 * h_quadpt ( ix, :, 01)

!      h ( ix, :, 06) = 0.5 * ( th ( ix, :, 04) + th ( ix, :, 05)
!     1 - th ( ix, :, 06) ) - 2.0 * h_quadpt ( ix, :, 01)

      h ( ix, :, 07) = 0.5 * ( th ( ix, :, 07) + th ( ix, :, 08)
     1 - th ( ix, :, 09) ) - 2.0 * h_quadpt ( ix, :, 01)

!      h ( ix, :, 08) = h_quadpt ( ix, :, 02) - h_quadpt ( ix, :, 03)
!     1 - h_quadpt ( ix, :, 06) + h_quadpt ( ix, :, 07)

      h ( ix, :, 09) = h_quadpt ( ix, :, 10) - h_quadpt ( ix, :, 11)
     1 - h_quadpt ( ix, :, 12) + h_quadpt ( ix, :, 13)

      h ( ix, :, 10) = h_quadpt ( ix, :, 04) - h_quadpt ( ix, :, 05)
     1 - h_quadpt ( ix, :, 08) + h_quadpt ( ix, :, 09)

!      h ( ix, :, 01) = h_quadpt ( ix, :, 01) + ( h ( ix, :, 05)
!     1 + h ( ix, :, 06) + h ( ix, :, 07) ) / 12.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      f ( ix, :, 05) = 2.0 * ( f_quadpt ( ix, :, 02)
     1 - 2.0 * f_quadpt ( ix, :, 01) + f_quadpt ( ix, :, 05) )

!      f ( ix, :, 06) = 2.0 * ( f_quadpt ( ix, :, 06)
!     1 - 2.0 * f_quadpt ( ix, :, 01) + f_quadpt ( ix, :, 09) )

!      f ( ix, :, 07) = 2.0 * ( f_quadpt ( ix, :, 10)
!     1 - 2.0 * f_quadpt ( ix, :, 01) + f_quadpt ( ix, :, 13) )

      f ( ix, :, 08) = 0.5 * ( f_quadpt ( ix, :, 14)
     1 - f_quadpt ( ix, :, 15) + f_quadpt ( ix, :, 16)
     1 - f_quadpt ( ix, :, 17) - f_quadpt ( ix, :, 18)
     1 + f_quadpt ( ix, :, 19) - f_quadpt ( ix, :, 20)
     1 + f_quadpt ( ix, :, 21) )

!      f ( ix, :, 09) = 0.5 * ( f_quadpt ( ix, :, 14)
!     1 - f_quadpt ( ix, :, 16) + f_quadpt ( ix, :, 18)
!     1 - f_quadpt ( ix, :, 20) - f_quadpt ( ix, :, 15)
!     1 + f_quadpt ( ix, :, 17) - f_quadpt ( ix, :, 19)
!     1 + f_quadpt ( ix, :, 21) )

      f ( ix, :, 10) = 0.5 * ( f_quadpt ( ix, :, 14)
     1 - f_quadpt ( ix, :, 18) + f_quadpt ( ix, :, 15)
     1 - f_quadpt ( ix, :, 19) - f_quadpt ( ix, :, 16)
     1 + f_quadpt ( ix, :, 20) - f_quadpt ( ix, :, 17)
     1 + f_quadpt ( ix, :, 21) )

      f ( ix, :, 11) = 16 * ( f_quadpt ( ix, :, 02)
     1 - 2.0 * f_quadpt ( ix, :, 03) + 2.0 * f_quadpt ( ix, :, 04)
     1 - f_quadpt ( ix, :, 05) ) / 3.0

!      f ( ix, :, 12) = 16.0 * ( f_quadpt ( ix, :, 06)
!     1 - 2.0 * f_quadpt ( ix, :, 07) + 2.0 * f_quadpt ( ix, :, 08)
!     1 - f_quadpt ( ix, :, 09) ) / 3.0

!      f ( ix, :, 13) = 16.0 * ( f_quadpt ( ix, :, 10)
!     1 - 2.0 * f_quadpt ( ix, :, 11) + 2.0 * f_quadpt ( ix, :, 12)
!     1 - f_quadpt ( ix, :, 13) ) / 3.0

      tf ( ix, :, 01) = f_quadpt ( ix, :, 14) + f_quadpt ( ix, :, 15)
     1 + f_quadpt ( ix, :, 16) + f_quadpt ( ix, :, 17)

      tf ( ix, :, 02) = f_quadpt ( ix, :, 18) + f_quadpt ( ix, :, 19)
     1 + f_quadpt ( ix, :, 20) + f_quadpt ( ix, :, 21)

      tf ( ix, :, 03) = f_quadpt ( ix, :, 14) + f_quadpt ( ix, :, 16)
     1 + f_quadpt ( ix, :, 18) + f_quadpt ( ix, :, 20)

      tf ( ix, :, 04) = f_quadpt ( ix, :, 15) + f_quadpt ( ix, :, 17)
     1 + f_quadpt ( ix, :, 19) + f_quadpt ( ix, :, 21)

      tf ( ix, :, 05) = f_quadpt ( ix, :, 14) + f_quadpt ( ix, :, 15)
     1 + f_quadpt ( ix, :, 18) + f_quadpt ( ix, :, 19)

      tf ( ix, :, 06) = f_quadpt ( ix, :, 16) + f_quadpt ( ix, :, 17)
     1 + f_quadpt ( ix, :, 20) + f_quadpt ( ix, :, 21)

      f ( ix, :, 02) = 2.0 * ( f_quadpt ( ix, :, 02)
     1 - f_quadpt ( ix, :, 05) ) / 3.0 - 0.1 * f ( ix, :, 11)
     1 + 0.25 * ( tf ( ix, :, 01) - tf ( ix, :, 02) ) / 3.0

!      f ( ix, :, 03) = 2.0 * ( f_quadpt ( ix, :, 06)
!     1 - f_quadpt ( ix, :, 09) ) / 3.0 - 0.1 * f ( ix, :, 12)
!     1 + 0.25 * ( tf ( ix, :, 03) - tf ( ix, :, 04) ) / 3.0

!      f ( ix, :, 04) = 2.0 * ( f_quadpt ( ix, :, 10)
!     1 - f_quadpt ( ix, :, 13) ) / 3.0 - 0.1 * f ( ix, :, 13)
!     1 + 0.25 * ( tf ( ix, :, 05) - tf ( ix, :, 06) ) / 3.0

      f ( ix, :, 20) = f_quadpt ( ix, :, 14) - f_quadpt ( ix, :, 18)
     1 - f_quadpt ( ix, :, 16) + f_quadpt ( ix, :, 20)
     1 - f_quadpt ( ix, :, 15) + f_quadpt ( ix, :, 19)
     1 + f_quadpt ( ix, :, 17) - f_quadpt ( ix, :, 21)

!      f ( ix, :, 01) = 0.125 * ( f_quadpt ( ix, :, 14)
!     1 + f_quadpt ( ix, :, 15) + f_quadpt ( ix, :, 16)
!     1 + f_quadpt ( ix, :, 17) + f_quadpt ( ix, :, 18)
!     1 + f_quadpt ( ix, :, 19) + f_quadpt ( ix, :, 20)
!     1 + f_quadpt ( ix, :, 21) ) - 0.5 * ( f ( ix, :, 05)
!     1 + f ( ix, :, 06) + f ( ix, :, 07) ) / 3.0

! ----------

!      g ( ix, :, 05) = 2.0 * ( g_quadpt ( ix, :, 02)
!     1 - 2.0 * g_quadpt ( ix, :, 01) + g_quadpt ( ix, :, 05) )

      g ( ix, :, 06) = 2.0 * ( g_quadpt ( ix, :, 06)
     1 - 2.0 * g_quadpt ( ix, :, 01) + g_quadpt ( ix, :, 09) )

!      g ( ix, :, 07) = 2.0 * ( g_quadpt ( ix, :, 10)
!     1 - 2.0 * g_quadpt ( ix, :, 01) + g_quadpt ( ix, :, 13) )

      g ( ix, :, 08) = 0.5 * ( g_quadpt ( ix, :, 14)
     1 - g_quadpt ( ix, :, 15) + g_quadpt ( ix, :, 16)
     1 - g_quadpt ( ix, :, 17) - g_quadpt ( ix, :, 18)
     1 + g_quadpt ( ix, :, 19) - g_quadpt ( ix, :, 20)
     1 + g_quadpt ( ix, :, 21) )

      g ( ix, :, 09) = 0.5 * ( g_quadpt ( ix, :, 14)
     1 - g_quadpt ( ix, :, 16) + g_quadpt ( ix, :, 18)
     1 - g_quadpt ( ix, :, 20) - g_quadpt ( ix, :, 15)
     1 + g_quadpt ( ix, :, 17) - g_quadpt ( ix, :, 19)
     1 + g_quadpt ( ix, :, 21) )

!      g ( ix, :, 10) = 0.5 * ( g_quadpt ( ix, :, 14)
!     1 - g_quadpt ( ix, :, 18) + g_quadpt ( ix, :, 15)
!     1 - g_quadpt ( ix, :, 19) - g_quadpt ( ix, :, 16)
!     1 + g_quadpt ( ix, :, 20) - g_quadpt ( ix, :, 17)
!     1 + g_quadpt ( ix, :, 21) )

!      g ( ix, :, 11) = 16 * ( g_quadpt ( ix, :, 02)
!     1 - 2.0 * g_quadpt ( ix, :, 03) + 2.0 * g_quadpt ( ix, :, 04)
!     1 - g_quadpt ( ix, :, 05) ) / 3.0

      g ( ix, :, 12) = 16.0 * ( g_quadpt ( ix, :, 06)
     1 - 2.0 * g_quadpt ( ix, :, 07) + 2.0 * g_quadpt ( ix, :, 08)
     1 - g_quadpt ( ix, :, 09) ) / 3.0

!      g ( ix, :, 13) = 16.0 * ( g_quadpt ( ix, :, 10)
!     1 - 2.0 * g_quadpt ( ix, :, 11) + 2.0 * g_quadpt ( ix, :, 12)
!     1 - g_quadpt ( ix, :, 13) ) / 3.0

      tg ( ix, :, 01) = g_quadpt ( ix, :, 14) + g_quadpt ( ix, :, 15)
     1 + g_quadpt ( ix, :, 16) + g_quadpt ( ix, :, 17)

      tg ( ix, :, 02) = g_quadpt ( ix, :, 18) + g_quadpt ( ix, :, 19)
     1 + g_quadpt ( ix, :, 20) + g_quadpt ( ix, :, 21)

      tg ( ix, :, 03) = g_quadpt ( ix, :, 14) + g_quadpt ( ix, :, 16)
     1 + g_quadpt ( ix, :, 18) + g_quadpt ( ix, :, 20)

      tg ( ix, :, 04) = g_quadpt ( ix, :, 15) + g_quadpt ( ix, :, 17)
     1 + g_quadpt ( ix, :, 19) + g_quadpt ( ix, :, 21)

      tg ( ix, :, 05) = g_quadpt ( ix, :, 14) + g_quadpt ( ix, :, 15)
     1 + g_quadpt ( ix, :, 18) + g_quadpt ( ix, :, 19)

      tg ( ix, :, 06) = g_quadpt ( ix, :, 16) + g_quadpt ( ix, :, 17)
     1 + g_quadpt ( ix, :, 20) + g_quadpt ( ix, :, 21)

!      g ( ix, :, 02) = 2.0 * ( g_quadpt ( ix, :, 02)
!     1 - g_quadpt ( ix, :, 05) ) / 3.0 - 0.1 * g ( ix, :, 11)
!     1 + 0.25 * ( tg ( ix, :, 01) - tg ( ix, :, 02) ) / 3.0

      g ( ix, :, 03) = 2.0 * ( g_quadpt ( ix, :, 06)
     1 - g_quadpt ( ix, :, 09) ) / 3.0 - 0.1 * g ( ix, :, 12)
     1 + 0.25 * ( tg ( ix, :, 03) - tg ( ix, :, 04) ) / 3.0

!      g ( ix, :, 04) = 2.0 * ( g_quadpt ( ix, :, 10)
!     1 - g_quadpt ( ix, :, 13) ) / 3.0 - 0.1 * g ( ix, :, 13)
!     1 + 0.25 * ( tg ( ix, :, 05) - tg ( ix, :, 06) ) / 3.0

      g ( ix, :, 20) = g_quadpt ( ix, :, 14) - g_quadpt ( ix, :, 18)
     1 - g_quadpt ( ix, :, 16) + g_quadpt ( ix, :, 20)
     1 - g_quadpt ( ix, :, 15) + g_quadpt ( ix, :, 19)
     1 + g_quadpt ( ix, :, 17) - g_quadpt ( ix, :, 21)

!      g ( ix, :, 01) = 0.125 * ( g_quadpt ( ix, :, 14)
!     1 + g_quadpt ( ix, :, 15) + g_quadpt ( ix, :, 16)
!     1 + g_quadpt ( ix, :, 17) + g_quadpt ( ix, :, 18)
!     1 + g_quadpt ( ix, :, 19) + g_quadpt ( ix, :, 20)
!     1 + g_quadpt ( ix, :, 21) ) - 0.5 * ( g ( ix, :, 05)
!     1 + g ( ix, :, 06) + g ( ix, :, 07) ) / 3.0

! ----------

!      h ( ix, :, 05) = 2.0 * ( h_quadpt ( ix, :, 02)
!     1 - 2.0 * h_quadpt ( ix, :, 01) + h_quadpt ( ix, :, 05) )

!      h ( ix, :, 06) = 2.0 * ( h_quadpt ( ix, :, 06)
!     1 - 2.0 * h_quadpt ( ix, :, 01) + h_quadpt ( ix, :, 09) )

      h ( ix, :, 07) = 2.0 * ( h_quadpt ( ix, :, 10)
     1 - 2.0 * h_quadpt ( ix, :, 01) + h_quadpt ( ix, :, 13) )

!      h ( ix, :, 08) = 0.5 * ( h_quadpt ( ix, :, 14)
!     1 - h_quadpt ( ix, :, 15) + h_quadpt ( ix, :, 16)
!     1 - h_quadpt ( ix, :, 17) - h_quadpt ( ix, :, 18)
!     1 + h_quadpt ( ix, :, 19) - h_quadpt ( ix, :, 20)
!     1 + h_quadpt ( ix, :, 21) )

      h ( ix, :, 09) = 0.5 * ( h_quadpt ( ix, :, 14)
     1 - h_quadpt ( ix, :, 16) + h_quadpt ( ix, :, 18)
     1 - h_quadpt ( ix, :, 20) - h_quadpt ( ix, :, 15)
     1 + h_quadpt ( ix, :, 17) - h_quadpt ( ix, :, 19)
     1 + h_quadpt ( ix, :, 21) )

      h ( ix, :, 10) = 0.5 * ( h_quadpt ( ix, :, 14)
     1 - h_quadpt ( ix, :, 18) + h_quadpt ( ix, :, 15)
     1 - h_quadpt ( ix, :, 19) - h_quadpt ( ix, :, 16)
     1 + h_quadpt ( ix, :, 20) - h_quadpt ( ix, :, 17)
     1 + h_quadpt ( ix, :, 21) )

!      h ( ix, :, 11) = 16 * ( h_quadpt ( ix, :, 02)
!     1 - 2.0 * h_quadpt ( ix, :, 03) + 2.0 * h_quadpt ( ix, :, 04)
!     1 - h_quadpt ( ix, :, 05) ) / 3.0

!      h ( ix, :, 12) = 16.0 * ( h_quadpt ( ix, :, 06)
!     1 - 2.0 * h_quadpt ( ix, :, 07) + 2.0 * h_quadpt ( ix, :, 08)
!     1 - h_quadpt ( ix, :, 09) ) / 3.0

      h ( ix, :, 13) = 16.0 * ( h_quadpt ( ix, :, 10)
     1 - 2.0 * h_quadpt ( ix, :, 11) + 2.0 * h_quadpt ( ix, :, 12)
     1 - h_quadpt ( ix, :, 13) ) / 3.0

      th ( ix, :, 01) = h_quadpt ( ix, :, 14) + h_quadpt ( ix, :, 15)
     1 + h_quadpt ( ix, :, 16) + h_quadpt ( ix, :, 17)

      th ( ix, :, 02) = h_quadpt ( ix, :, 18) + h_quadpt ( ix, :, 19)
     1 + h_quadpt ( ix, :, 20) + h_quadpt ( ix, :, 21)

      th ( ix, :, 03) = h_quadpt ( ix, :, 14) + h_quadpt ( ix, :, 16)
     1 + h_quadpt ( ix, :, 18) + h_quadpt ( ix, :, 20)

      th ( ix, :, 04) = h_quadpt ( ix, :, 15) + h_quadpt ( ix, :, 17)
     1 + h_quadpt ( ix, :, 19) + h_quadpt ( ix, :, 21)

      th ( ix, :, 05) = h_quadpt ( ix, :, 14) + h_quadpt ( ix, :, 15)
     1 + h_quadpt ( ix, :, 18) + h_quadpt ( ix, :, 19)

      th ( ix, :, 06) = h_quadpt ( ix, :, 16) + h_quadpt ( ix, :, 17)
     1 + h_quadpt ( ix, :, 20) + h_quadpt ( ix, :, 21)

!      h ( ix, :, 02) = 2.0 * ( h_quadpt ( ix, :, 02)
!     1 - h_quadpt ( ix, :, 05) ) / 3.0 - 0.1 * h ( ix, :, 11)
!     1 + 0.25 * ( th ( ix, :, 01) - th ( ix, :, 02) ) / 3.0

!      h ( ix, :, 03) = 2.0 * ( h_quadpt ( ix, :, 06)
!     1 - h_quadpt ( ix, :, 09) ) / 3.0 - 0.1 * h ( ix, :, 12)
!     1 + 0.25 * ( th ( ix, :, 03) - th ( ix, :, 04) ) / 3.0

      h ( ix, :, 04) = 2.0 * ( h_quadpt ( ix, :, 10)
     1 - h_quadpt ( ix, :, 13) ) / 3.0 - 0.1 * h ( ix, :, 13)
     1 + 0.25 * ( th ( ix, :, 05) - th ( ix, :, 06) ) / 3.0

      h ( ix, :, 20) = h_quadpt ( ix, :, 14) - h_quadpt ( ix, :, 18)
     1 - h_quadpt ( ix, :, 16) + h_quadpt ( ix, :, 20)
     1 - h_quadpt ( ix, :, 15) + h_quadpt ( ix, :, 19)
     1 + h_quadpt ( ix, :, 17) - h_quadpt ( ix, :, 21)

!      h ( ix, :, 01) = 0.125 * ( h_quadpt ( ix, :, 14)
!     1 + h_quadpt ( ix, :, 15) + h_quadpt ( ix, :, 16)
!     1 + h_quadpt ( ix, :, 17) + h_quadpt ( ix, :, 18)
!     1 + h_quadpt ( ix, :, 19) + h_quadpt ( ix, :, 20)
!     1 + h_quadpt ( ix, :, 21) ) - 0.5 * ( h ( ix, :, 05)
!     1 + h ( ix, :, 06) + h ( ix, :, 07) ) / 3.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Default is to build source terms at t=0 . In some ADER-CG updates
! (on Cartesian meshes), we may not require any source terms.
#if ( BUILD_SOURCE_TERMS == 1)
! ----------------------------------------------------------------------

! STEP 5.a) Build source terms "s_quadpt" at the spatial quadrature points at t=0.
! Remember that the primitive variables are still available from the previous
! loop at each of the t=0 quadrature points for constructing the source terms.

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      s_quadpt ( oned_begin - oned_off: oned_end + oned_off,
     1                                           :, iquadpt) = 0.0

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

#if ( 0 == 1)

! For single-fluid adiabatic MHD turbulence, we may want to provide a gradual
! (explicit) heating/cooling term which ensures that pressure always
! tends to a mean pressure "prs_bar" in a time "tau_bar" (comparable
! to a turnover time). Such a source term can be implemented explicitly
! and its implementation is given below.
! Remember to declare "prs_bar, tau_bar".

      prs_bar = 1.0
      tau_bar = 1.0

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 5, iquadpt) = - ( prs ( ix, iquadpt) - prs_bar)
     1     / ( tau_bar * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0) )

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

#endif

! ------------------------------

! Build two-fluid source terms. In this part, the source terms are
! hard-coded for one family of ionized fluid and one family of neutrals.

#if ( N_EULER == 1) 

      IF ( ( igeom .EQ. 1) .AND. ( N_EULER .EQ. 1) .AND.
     1     ( ambipolar_diffusion_drag_coef .GT. 0.0)) THEN

      ivar = NFIELD + 1 + NFLUID

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 1, iquadpt) = 0.0

      tempaa = ambipolar_diffusion_drag_coef
     1       * rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

      s_quadpt ( ix, 2, iquadpt) = tempaa
     1           * ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt))

      s_quadpt ( ix, 3, iquadpt) = tempaa
     1           * ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt))

      s_quadpt ( ix, 4, iquadpt) = tempaa
     1           * ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt))

      s_quadpt ( ix, 5, iquadpt) = ( ambipolar_diffusion_drag_coef
     1   / ( molewt_ion + molewt_neut) )

     1 * ( rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

     1   * ( ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt) )
     1     * ( molewt_neut * vlx_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlx ( ix, iquadpt) )
     1     + ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt) )
     1     * ( molewt_neut * vly_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vly ( ix, iquadpt) )
     1     + ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt) )
     1     * ( molewt_neut * vlz_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlz ( ix, iquadpt) )
     1     )

     1   + 3.0 * ( molewt_neut * prs_euler ( ix, iquadpt, 1)
     1           * rho ( ix, iquadpt)
     1           - molewt_ion * prs ( ix, iquadpt)
     1           * rho_euler ( ix, iquadpt, 1) )
     1   )

      s_quadpt ( ix, 6, iquadpt) = 0.0
      s_quadpt ( ix, 7, iquadpt) = 0.0
      s_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      s_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) = 0.0
#endif


      s_quadpt ( ix, ivar + 1, iquadpt) = 0.0
      s_quadpt ( ix, ivar + 2, iquadpt) = - s_quadpt ( ix, 2, iquadpt)
      s_quadpt ( ix, ivar + 3, iquadpt) = - s_quadpt ( ix, 3, iquadpt)
      s_quadpt ( ix, ivar + 4, iquadpt) = - s_quadpt ( ix, 4, iquadpt)
      s_quadpt ( ix, ivar + 5, iquadpt) = - s_quadpt ( ix, 5, iquadpt)

#if ( NFLUID_EULER >= 1)
      s_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

      END IF
! End of "ambipolar_diffusion_drag_coef .GT. 0.0" part

#endif
! End of "N_EULER == 1"

! ------------------------------

! Cylindrical Geometry

      IF ( igeom .EQ. 2) THEN

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1                           * vlx ( ix, iquadpt) * obx1

      s_quadpt ( ix, 2, iquadpt) =
     1 ( rho ( ix, iquadpt) * vly ( ix, iquadpt)**2
     1 - rho ( ix, iquadpt) * vlx ( ix, iquadpt)**2
     1 + bfx ( ix, iquadpt)**2 * ob4pi
     1 - bfy ( ix, iquadpt)**2 * ob4pi) * obx1

      s_quadpt ( ix, 3, iquadpt) =
     1 - 2.0 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt) 
     1                              * vly ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - ( rho ( ix, iquadpt) * vlx ( ix, iquadpt) 
     1                        * vlz ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 5, iquadpt) = - ( tempac * vlx ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = 0.0

      s_quadpt ( ix, 7, iquadpt) = 0.0

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Spherical Geometry

      IF ( igeom .EQ. 3) THEN

      DO iquadpt = 1, NUM_ZONECEN_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      cos_x2 = cos ( x2)
      sin_x2 = sin ( x2)
      cot_x2 = cos_x2 / sin_x2

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 2.0 * vlx ( ix, iquadpt)
     1 + cot_x2 * vly ( ix, iquadpt) ) * obx1

      s_quadpt ( ix, 2, iquadpt) = rho ( ix, iquadpt)
     1 * ( - 2.0 * vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1   + vlz ( ix, iquadpt)**2) * obx1
     1 - cot_x2
     1 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1   - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1
     1 - ( - 2.0 * bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1   + bfz ( ix, iquadpt)**2) * ob4pi * obx1

      s_quadpt ( ix, 3, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1 + cot_x2
     1 * ( vly ( ix, iquadpt)**2 - vlz ( ix, iquadpt)**2)) * obx1
     1 + ( 3.0 * bfx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1   + cot_x2
     1   * ( bfy ( ix, iquadpt)**2 - bfz ( ix, iquadpt)**2))
     1 * ob4pi * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - rho ( ix, iquadpt) * vlz ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) + 2.0 * cot_x2
     1   * vly ( ix, iquadpt) ) * obx1
     1 + bfz ( ix, iquadpt)
     1 * ( 3.0 * bfx ( ix, iquadpt) + 2.0 * cot_x2
     1   * bfy ( ix, iquadpt) ) * ob4pi * obx1

      s_quadpt ( ix, 5, iquadpt) =
     1 - 2.0 * ( tempac * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab) * obx1
     1 - cot_x2 * ( tempac * vly ( ix, iquadpt)
     1 - bfy ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = - cot_x2
     1                           * g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 7, iquadpt) = g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Multiply "s_quadpt" by "dtmove".

      DO ix = oned_begin - oned_off, oned_end + oned_off
      s_quadpt ( ix, :, 1: NUM_ZONECEN_QUADPTS) = dtmove
     1         * s_quadpt ( ix, :, 1: NUM_ZONECEN_QUADPTS)
      END DO

! ----------------------------------------------------------------------

! STEP 5.b) Now project the sources "s_quadpt" from the spatial nodes at t=0
! back to the t=0 part of the modal space time variables for those variables,
! i.e. to "s".

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s ( ix, :, 01) = ( s_quadpt ( ix, :, 01) + s_quadpt ( ix, :, 02)
     1 + s_quadpt ( ix, :, 03) + s_quadpt ( ix, :, 04)
     1 + s_quadpt ( ix, :, 05) + s_quadpt ( ix, :, 06)) / 6.0

      s ( ix, :, 02) = s_quadpt ( ix, :, 01) - s_quadpt ( ix, :, 02)

      s ( ix, :, 03) = s_quadpt ( ix, :, 03) - s_quadpt ( ix, :, 04)

      s ( ix, :, 04) = s_quadpt ( ix, :, 05) - s_quadpt ( ix, :, 06)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ts ( ix, :, 01) = s_quadpt ( ix, :, 02) + s_quadpt ( ix, :, 03)
     1           + s_quadpt ( ix, :, 04) + s_quadpt ( ix, :, 05)

      ts ( ix, :, 02) = s_quadpt ( ix, :, 06) + s_quadpt ( ix, :, 07)
     1           + s_quadpt ( ix, :, 08) + s_quadpt ( ix, :, 09)

      ts ( ix, :, 03) = s_quadpt ( ix, :, 10) + s_quadpt ( ix, :, 11)
     1           + s_quadpt ( ix, :, 12) + s_quadpt ( ix, :, 13)

      ts ( ix, :, 04) = s_quadpt ( ix, :, 02) + s_quadpt ( ix, :, 06)
     1           + s_quadpt ( ix, :, 10) + s_quadpt ( ix, :, 12)

      ts ( ix, :, 05) = s_quadpt ( ix, :, 03) + s_quadpt ( ix, :, 07)
     1           + s_quadpt ( ix, :, 11) + s_quadpt ( ix, :, 13)

      ts ( ix, :, 06) = s_quadpt ( ix, :, 04) + s_quadpt ( ix, :, 05)
     1           + s_quadpt ( ix, :, 08) + s_quadpt ( ix, :, 09)

      ts ( ix, :, 07) = s_quadpt ( ix, :, 04) + s_quadpt ( ix, :, 08)
     1           + s_quadpt ( ix, :, 10) + s_quadpt ( ix, :, 11)

      ts ( ix, :, 08) = s_quadpt ( ix, :, 05) + s_quadpt ( ix, :, 09)
     1           + s_quadpt ( ix, :, 12) + s_quadpt ( ix, :, 13)

      ts ( ix, :, 09) = s_quadpt ( ix, :, 02) + s_quadpt ( ix, :, 03)
     1           + s_quadpt ( ix, :, 06) + s_quadpt ( ix, :, 07)

      s ( ix, :, 02) = 0.25 * ( ts ( ix, :, 01) - ts ( ix, :, 02) )

      s ( ix, :, 03) = 0.25 * ( ts ( ix, :, 04) - ts ( ix, :, 05) )

      s ( ix, :, 04) = 0.25 * ( ts ( ix, :, 07) - ts ( ix, :, 08) )

      s ( ix, :, 05) = 0.5 * ( ts ( ix, :, 01) + ts ( ix, :, 02)
     1 - ts ( ix, :, 03) ) - 2.0 * s_quadpt ( ix, :, 01)

      s ( ix, :, 06) = 0.5 * ( ts ( ix, :, 04) + ts ( ix, :, 05)
     1 - ts ( ix, :, 06) ) - 2.0 * s_quadpt ( ix, :, 01)

      s ( ix, :, 07) = 0.5 * ( ts ( ix, :, 07) + ts ( ix, :, 08)
     1 - ts ( ix, :, 09) ) - 2.0 * s_quadpt ( ix, :, 01)

      s ( ix, :, 08) = s_quadpt ( ix, :, 02) - s_quadpt ( ix, :, 03)
     1 - s_quadpt ( ix, :, 06) + s_quadpt ( ix, :, 07)

      s ( ix, :, 09) = s_quadpt ( ix, :, 10) - s_quadpt ( ix, :, 11)
     1 - s_quadpt ( ix, :, 12) + s_quadpt ( ix, :, 13)

      s ( ix, :, 10) = s_quadpt ( ix, :, 04) - s_quadpt ( ix, :, 05)
     1 - s_quadpt ( ix, :, 08) + s_quadpt ( ix, :, 09)

      s ( ix, :, 01) = s_quadpt ( ix, :, 01) + ( s ( ix, :, 05)
     1 + s ( ix, :, 06) + s ( ix, :, 07) ) / 12.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s ( ix, :, 05) = 2.0 * ( s_quadpt ( ix, :, 02)
     1 - 2.0 * s_quadpt ( ix, :, 01) + s_quadpt ( ix, :, 05) )

      s ( ix, :, 06) = 2.0 * ( s_quadpt ( ix, :, 06)
     1 - 2.0 * s_quadpt ( ix, :, 01) + s_quadpt ( ix, :, 09) )

      s ( ix, :, 07) = 2.0 * ( s_quadpt ( ix, :, 10)
     1 - 2.0 * s_quadpt ( ix, :, 01) + s_quadpt ( ix, :, 13) )

      s ( ix, :, 08) = 0.5 * ( s_quadpt ( ix, :, 14)
     1 - s_quadpt ( ix, :, 15) + s_quadpt ( ix, :, 16)
     1 - s_quadpt ( ix, :, 17) - s_quadpt ( ix, :, 18)
     1 + s_quadpt ( ix, :, 19) - s_quadpt ( ix, :, 20)
     1 + s_quadpt ( ix, :, 21) )

      s ( ix, :, 09) = 0.5 * ( s_quadpt ( ix, :, 14)
     1 - s_quadpt ( ix, :, 16) + s_quadpt ( ix, :, 18)
     1 - s_quadpt ( ix, :, 20) - s_quadpt ( ix, :, 15)
     1 + s_quadpt ( ix, :, 17) - s_quadpt ( ix, :, 19)
     1 + s_quadpt ( ix, :, 21) )

      s ( ix, :, 10) = 0.5 * ( s_quadpt ( ix, :, 14)
     1 - s_quadpt ( ix, :, 18) + s_quadpt ( ix, :, 15)
     1 - s_quadpt ( ix, :, 19) - s_quadpt ( ix, :, 16)
     1 + s_quadpt ( ix, :, 20) - s_quadpt ( ix, :, 17)
     1 + s_quadpt ( ix, :, 21) )

      s ( ix, :, 11) = 16 * ( s_quadpt ( ix, :, 02)
     1 - 2.0 * s_quadpt ( ix, :, 03) + 2.0 * s_quadpt ( ix, :, 04)
     1 - s_quadpt ( ix, :, 05) ) / 3.0

      s ( ix, :, 12) = 16.0 * ( s_quadpt ( ix, :, 06)
     1 - 2.0 * s_quadpt ( ix, :, 07) + 2.0 * s_quadpt ( ix, :, 08)
     1 - s_quadpt ( ix, :, 09) ) / 3.0

      s ( ix, :, 13) = 16.0 * ( s_quadpt ( ix, :, 10)
     1 - 2.0 * s_quadpt ( ix, :, 11) + 2.0 * s_quadpt ( ix, :, 12)
     1 - s_quadpt ( ix, :, 13) ) / 3.0

      ts ( ix, :, 01) = s_quadpt ( ix, :, 14) + s_quadpt ( ix, :, 15)
     1 + s_quadpt ( ix, :, 16) + s_quadpt ( ix, :, 17)

      ts ( ix, :, 02) = s_quadpt ( ix, :, 18) + s_quadpt ( ix, :, 19)
     1 + s_quadpt ( ix, :, 20) + s_quadpt ( ix, :, 21)

      ts ( ix, :, 03) = s_quadpt ( ix, :, 14) + s_quadpt ( ix, :, 16)
     1 + s_quadpt ( ix, :, 18) + s_quadpt ( ix, :, 20)

      ts ( ix, :, 04) = s_quadpt ( ix, :, 15) + s_quadpt ( ix, :, 17)
     1 + s_quadpt ( ix, :, 19) + s_quadpt ( ix, :, 21)

      ts ( ix, :, 05) = s_quadpt ( ix, :, 14) + s_quadpt ( ix, :, 15)
     1 + s_quadpt ( ix, :, 18) + s_quadpt ( ix, :, 19)

      ts ( ix, :, 06) = s_quadpt ( ix, :, 16) + s_quadpt ( ix, :, 17)
     1 + s_quadpt ( ix, :, 20) + s_quadpt ( ix, :, 21)

      s ( ix, :, 02) = 2.0 * ( s_quadpt ( ix, :, 02)
     1 - s_quadpt ( ix, :, 05) ) / 3.0 - 0.1 * s ( ix, :, 11)
     1 + 0.25 * ( ts ( ix, :, 01) - ts ( ix, :, 02) ) / 3.0

      s ( ix, :, 03) = 2.0 * ( s_quadpt ( ix, :, 06)
     1 - s_quadpt ( ix, :, 09) ) / 3.0 - 0.1 * s ( ix, :, 12)
     1 + 0.25 * ( ts ( ix, :, 03) - ts ( ix, :, 04) ) / 3.0

      s ( ix, :, 04) = 2.0 * ( s_quadpt ( ix, :, 10)
     1 - s_quadpt ( ix, :, 13) ) / 3.0 - 0.1 * s ( ix, :, 13)
     1 + 0.25 * ( ts ( ix, :, 05) - ts ( ix, :, 06) ) / 3.0

      s ( ix, :, 20) = s_quadpt ( ix, :, 14) - s_quadpt ( ix, :, 18)
     1 - s_quadpt ( ix, :, 16) + s_quadpt ( ix, :, 20)
     1 - s_quadpt ( ix, :, 15) + s_quadpt ( ix, :, 19)
     1 + s_quadpt ( ix, :, 17) - s_quadpt ( ix, :, 21)

      s ( ix, :, 01) = 0.125 * ( s_quadpt ( ix, :, 14)
     1 + s_quadpt ( ix, :, 15) + s_quadpt ( ix, :, 16)
     1 + s_quadpt ( ix, :, 17) + s_quadpt ( ix, :, 18)
     1 + s_quadpt ( ix, :, 19) + s_quadpt ( ix, :, 20)
     1 + s_quadpt ( ix, :, 21) ) - 0.5 * ( s ( ix, :, 05)
     1 + s ( ix, :, 06) + s ( ix, :, 07) ) / 3.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
#endif
! End of "BUILD_SOURCE_TERMS == 1"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Start the outer ADER iteration here.

      DO iter_ader = 1, max_iter_ader_cg

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 6) Make the conserved variables at the t>0 space-time nodes.
! We put in the flexibility of having different "gamma_eos_spacetime" at each
! space-time node, but that is not used.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 07) = u ( ix, :, 01) + u ( ix, :, 05)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 14) = 0.25 * ( u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 03) + u_quadpt ( ix, :, 04)
     1 + u_quadpt ( ix, :, 05) + u ( ix, :, 12) + u ( ix, :, 13) )
     1 + 0.5 * u ( ix, :, 11) - 0.125 * ( u ( ix, :, 06)
     1 + u ( ix, :, 07) )

      u_quadpt ( ix, :, 15) = u_quadpt ( ix, :, 14)
     1 - u ( ix, :, 02) - 0.5 * u ( ix, :, 13)

      u_quadpt ( ix, :, 16) = 0.25 * ( u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 06) + u_quadpt ( ix, :, 10)
     1 + u_quadpt ( ix, :, 12) + u ( ix, :, 12) + u ( ix, :, 14) )
     1 + 0.5 * u ( ix, :, 11) - 0.125 * ( u ( ix, :, 05)
     1 + u ( ix, :, 07) )

      u_quadpt ( ix, :, 17) = u_quadpt ( ix, :, 16)
     1 - u ( ix, :, 03) - 0.5 * u ( ix, :, 14)

      u_quadpt ( ix, :, 18) = 0.25 * ( u_quadpt ( ix, :, 04)
     1 + u_quadpt ( ix, :, 08) + u_quadpt ( ix, :, 10)
     1 + u_quadpt ( ix, :, 11) + u ( ix, :, 12) + u ( ix, :, 15) )
     1 + 0.5 * u ( ix, :, 11) - 0.125 * ( u ( ix, :, 05)
     1 + u ( ix, :, 06) )

      u_quadpt ( ix, :, 19) = u_quadpt ( ix, :, 18) - u ( ix, :, 04)
     1 - 0.5 * u ( ix, :, 15)

      u_quadpt ( ix, :, 20) = u_quadpt ( ix, :, 01) + u ( ix, :, 11)
     1 + u ( ix, :, 12)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 22) = u_quadpt ( ix, :, 01)
     1 + ( 36.0 * u ( ix, :, 21) + 12.0 * u ( ix, :, 22)
     1 + 4.0 * u ( ix, :, 23) - 3.0 * u ( ix, :, 30)
     1 - 3.0 * u ( ix, :, 31) - 3.0 * u ( ix, :, 32) ) / 108.0

      u_quadpt ( ix, :, 23) = u_quadpt ( ix, :, 22)
     1 + u_quadpt ( ix, :, 02) + u_quadpt ( ix, :, 06)
     1 - 2.0 * u_quadpt ( ix, :, 01) + ( 18.0 * u ( ix, :, 08)
     1 + 9.0 * u ( ix, :, 14) + 9.0 * u ( ix, :, 16)
     1 + 12.0 * u ( ix, :, 24) + 12.0 * u ( ix, :, 25)
     1 + 4.0 * u ( ix, :, 27) + 4.0 * u ( ix, :, 28)
     1 + 6.0 * u ( ix, :, 30) + 6.0 * u ( ix, :, 31)
     1 + 6.0 * u ( ix, :, 33) ) / 72.0

      u_quadpt ( ix, :, 24) = u_quadpt ( ix, :, 23)
     1 - u_quadpt ( ix, :, 06) + u_quadpt ( ix, :, 09)
     1 - 0.5 * u ( ix, :, 08) - 0.25 * u ( ix, :, 14)
     1 - ( 6.0 * u ( ix, :, 25) + 2.0 * u ( ix, :, 28)
     1 + 3.0 * u ( ix, :, 33) ) / 18.0

      u_quadpt ( ix, :, 25) = 0.5 * ( u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 24) + u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 15) ) - u_quadpt ( ix, :, 02)
     1 - 0.5 * u ( ix, :, 06) - 0.25 * u ( ix, :, 16)
     1 - 0.125 * u ( ix, :, 17) + ( 6.0 * u ( ix, :, 26)
     1 + 2.0 * u ( ix, :, 29) - 3.0 * u ( ix, :, 31)
     1 + 3.0 * u ( ix, :, 32) + 3.0 * u ( ix, :, 35) ) / 36.0

      u_quadpt ( ix, :, 26) = u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 24) - u_quadpt ( ix, :, 25)
     1 - 0.5 * ( u ( ix, :, 06) - u ( ix, :, 07) )
     1 - 0.25 * ( u ( ix, :, 16) - u ( ix, :, 18) )
     1 - ( u ( ix, :, 31) - u ( ix, :, 32) ) / 6.0

      u_quadpt ( ix, :, 27) = u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 02)
     1 - 0.5 * u ( ix, :, 08) - 0.25 * u ( ix, :, 16)
     1 - ( 6.0 * u ( ix, :, 24) + 2.0 * u ( ix, :, 27)
     1 + 3.0 * u ( ix, :, 33) ) / 18.0

      u_quadpt ( ix, :, 28) = u_quadpt ( ix, :, 27)
     1 - u_quadpt ( ix, :, 06) + u_quadpt ( ix, :, 09)
     1 + 0.5 * u ( ix, :, 08) - 0.25 * u ( ix, :, 14)
     1 - ( 6.0 * u ( ix, :, 25) + 2.0 * u ( ix, :, 28)
     1 - 3.0 * u ( ix, :, 33) ) / 18.0

      u_quadpt ( ix, :, 29) = 0.5 * ( u_quadpt ( ix, :, 27)
     1 + u_quadpt ( ix, :, 28) + u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 19) ) - u_quadpt ( ix, :, 05)
     1 - 0.5 * u ( ix, :, 06) + 0.25 * u ( ix, :, 16)
     1 - 0.125 * u ( ix, :, 17) + ( 6.0 * u ( ix, :, 26)
     1 + 2.0 * u ( ix, :, 29) - 3.0 * u ( ix, :, 31)
     1 + 3.0 * u ( ix, :, 32)  - 3.0 * u ( ix, :, 35) ) / 36.0

      u_quadpt ( ix, :, 30) = u_quadpt ( ix, :, 27)
     1 + u_quadpt ( ix, :, 28) - u_quadpt ( ix, :, 29)
     1 - 0.5 * ( u ( ix, :, 06) - u ( ix, :, 07) )
     1 + 0.25 * ( u ( ix, :, 16) - u ( ix, :, 18) )
     1 - ( u ( ix, :, 31) - u ( ix, :, 32) ) / 6.0

      u_quadpt ( ix, :, 31) = 0.5 * ( u_quadpt ( ix, :, 25)
     1 + u_quadpt ( ix, :, 29) + u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 18) ) - u_quadpt ( ix, :, 10)
     1 - 0.5 * u ( ix, :, 05) - 0.25 * u ( ix, :, 15)
     1 - 0.125 * u ( ix, :, 14) + ( 6.0 * u ( ix, :, 25)
     1 + 2.0 * u ( ix, :, 28) - 3.0 * u ( ix, :, 30)
     1 + 3.0 * u ( ix, :, 31) + 3.0 * u ( ix, :, 34) ) / 36.0

      u_quadpt ( ix, :, 32) = u_quadpt ( ix, :, 31)
     1 - u_quadpt ( ix, :, 06) + u_quadpt ( ix, :, 09)
     1 - 0.5 * u ( ix, :, 09) - 0.25 * u ( ix, :, 19)
     1 - ( 6.0 * u ( ix, :, 25) + 2.0 * u ( ix, :, 28)
     1 + 3.0 * u ( ix, :, 34) ) / 18.0

      u_quadpt ( ix, :, 33) = u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 27) - u_quadpt ( ix, :, 31)
     1 - 0.5 * ( u ( ix, :, 05) - u ( ix, :, 07) )
     1 - 0.25 * ( u ( ix, :, 14) - u ( ix, :, 19) )
     1 - ( u ( ix, :, 30) - u ( ix, :, 32) ) / 6.0

      u_quadpt ( ix, :, 34) = u_quadpt ( ix, :, 24)
     1 + u_quadpt ( ix, :, 28) - u_quadpt ( ix, :, 32)
     1 - 0.5 * ( u ( ix, :, 05) - u ( ix, :, 07) )
     1 + 0.25 * ( u ( ix, :, 14) - u ( ix, :, 19) )
     1 - ( u ( ix, :, 30) - u ( ix, :, 32) ) / 6.0


      tu ( ix, :, 01) = ( 36.0 * u ( ix, :, 21) + 24.0 * u ( ix, :, 22)
     1 + 16.0 * u ( ix, :, 23) ) / 54.0

      tu ( ix, :, 02) = ( 18.0 * u ( ix, :, 24)
     1 + 12.0 * u ( ix, :, 27) ) / 54.0

      tu ( ix, :, 03) = ( 6.0 * u ( ix, :, 30) - 3.0 * u ( ix, :, 31)
     1 - 3.0 * u ( ix, :, 32) ) / 54.0

      tu ( ix, :, 04) = ( 18.0 * u ( ix, :, 25)
     1 + 12.0 * u ( ix, :, 28) ) / 54.0

      tu ( ix, :, 05) = ( - 3.0 * u ( ix, :, 30) + 6.0 * u ( ix, :, 31)
     1 - 3.0 * u ( ix, :, 32) ) / 54.0

      tu ( ix, :, 06) = ( 18.0 * u ( ix, :, 26)
     1 + 12.0 * u ( ix, :, 29) ) / 54.0

      tu ( ix, :, 07) = ( - 3.0 * u ( ix, :, 30) - 3.0 * u ( ix, :, 31)
     1 + 6.0 * u ( ix, :, 32) ) / 54.0

      u_quadpt ( ix, :, 35) = ( 2.0 * u_quadpt ( ix, :, 22)
     1 - u_quadpt ( ix, :, 01) ) + 2.0 * ( u ( ix, :, 22)
     1 + u ( ix, :, 23) ) / 9.0

      u_quadpt ( ix, :, 36) = u_quadpt ( ix, :, 02) + tu ( ix, :, 01)
     1 + tu ( ix, :, 02) + tu ( ix, :, 03)

      u_quadpt ( ix, :, 37) = u_quadpt ( ix, :, 05) + tu ( ix, :, 01)
     1 - tu ( ix, :, 02) + tu ( ix, :, 03)

      u_quadpt ( ix, :, 38) = u_quadpt ( ix, :, 06) + tu ( ix, :, 01)
     1 + tu ( ix, :, 04) + tu ( ix, :, 05)

      u_quadpt ( ix, :, 39) = u_quadpt ( ix, :, 09) + tu ( ix, :, 01)
     1 - tu ( ix, :, 04) + tu ( ix, :, 05)

      u_quadpt ( ix, :, 40) = u_quadpt ( ix, :, 10) + tu ( ix, :, 01)
     1 + tu ( ix, :, 06) + tu ( ix, :, 07)

      u_quadpt ( ix, :, 41) = u_quadpt ( ix, :, 13) + tu ( ix, :, 01)
     1 - tu ( ix, :, 06) + tu ( ix, :, 07)


      u_quadpt ( ix, :, 42) = u_quadpt ( ix, :, 01) + u ( ix, :, 21)
     1 + u ( ix, :, 22) + u ( ix, :, 23) - ( u ( ix, :, 30)
     1 + u ( ix, :, 31) + u ( ix, :, 32) ) / 12.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 7) Get the primitive variables from the conserved variables at
! the spatial nodes at t>0. The density and pressure variables are
! given positive floor values to ensure that they can be used for physically
! meaningful source terms.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, 1, iquadpt) =
     1                   AMAX1 ( u_quadpt ( ix, 1, iquadpt), rhofloor)

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempaa = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempaa
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempaa
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempaa

      bfx ( ix, iquadpt) = u_quadpt ( ix, 6, iquadpt)
      bfy ( ix, iquadpt) = u_quadpt ( ix, 7, iquadpt)
      bfz ( ix, iquadpt) = u_quadpt ( ix, 8, iquadpt)

      tempab = 0.5 * rho ( ix, iquadpt) * ( vlx ( ix, iquadpt)**2
     1         + vly ( ix, iquadpt)**2 + vlz ( ix, iquadpt)**2)

      tempac = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      eng_t ( ix, iquadpt) = u_quadpt ( ix, 5, iquadpt)
     1                     - tempab - tempac

      prs ( ix, iquadpt) = eng_t ( ix, iquadpt)
     1                   * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

#if ( ENFORCE_POSITIVITY == 1)
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = entropy ( ix, iquadpt)
     1* rho ( ix, iquadpt)**( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)
      prs ( ix, iquadpt) = AMAX1 ( prs ( ix, iquadpt), prsfloor)

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#else
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = prsfloor

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#endif

#if ( NFLUID >= 1)
      specie ( ix, 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( ix, 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( ix, 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( ix, 1: NFLUID, iquadpt) * rho ( ix, iquadpt)
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, ivar + 1, iquadpt) =
     1         AMAX1 ( u_quadpt ( ix, ivar + 1, iquadpt), rhofloor)

      rho_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 1, iquadpt)

      tempaa = 1.0 / rho_euler ( ix, iquadpt, i_euler)

      vlx_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 2, iquadpt) * tempaa
      vly_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 3, iquadpt) * tempaa
      vlz_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 4, iquadpt) * tempaa

      tempab = 0.5 * rho_euler ( ix, iquadpt, i_euler)
     1           * ( vlx_euler ( ix, iquadpt, i_euler)**2
     1             + vly_euler ( ix, iquadpt, i_euler)**2
     1             + vlz_euler ( ix, iquadpt, i_euler)**2)

      prs_euler ( ix, iquadpt, i_euler) =
     1         ( u_quadpt ( ix, ivar + 5, iquadpt) - tempab)
     1       * ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      prs_euler ( ix, iquadpt, i_euler) =
     1    AMAX1 ( prs_euler ( ix, iquadpt, i_euler), prsfloor)

#if ( ISOTHERMAL_EQNS == 1)
      prs_euler ( ix, iquadpt, i_euler) =
     1     isospeed**2 * rho_euler ( ix, iquadpt, i_euler)
     1   + barotropic_stiffness
     1   * rho_euler ( ix, iquadpt, i_euler)**gamma
#endif

      eng_t_euler ( ix, iquadpt, i_euler) =
     1         prs_euler ( ix, iquadpt, i_euler)
     1       / ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      IF ( prs_euler ( ix, iquadpt, i_euler) .LT. prsfloor * 1.00000001)
     1   u_quadpt ( ix, ivar + 5, iquadpt) =
     1                   eng_t_euler ( ix, iquadpt, i_euler) + tempab

#if ( NFLUID_EULER >= 1)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) = tempaa
     1       * u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1                    ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) =
     1    AMIN1 ( AMAX1 (
     1    specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler),
     1            0.0), 1.0 )

      u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1       specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler)
     1     * rho_euler ( ix, iquadpt, i_euler)

#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 8.a) Build fluxes ( f_quadpt, g_quadpt, h_quadpt) from the primitive
! variables that have now become available at spatial nodes at t=0.

! Multiply "f_quadpt" by "dtmove/dx" and so on.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! -----

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

! -----

      f_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 2, iquadpt)

      f_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfx ( ix, iquadpt)**2 * ob4pi

      f_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vly ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi

      f_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi

      f_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab

      f_quadpt ( ix, 6, iquadpt) = 0.0

      f_quadpt ( ix, 7, iquadpt) =
     1           vlx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         - vly ( ix, iquadpt) * bfx ( ix, iquadpt)

      f_quadpt ( ix, 8, iquadpt) =
     1           vlx ( ix, iquadpt) * bfz ( ix, iquadpt)
     1         - vlz ( ix, iquadpt) * bfx ( ix, iquadpt)

#if ( NFLUID >= 1)
      f_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlx ( ix, iquadpt)
#endif

! -----

      g_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 3, iquadpt)

      g_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         - bfy ( ix, iquadpt) * bfx ( ix, iquadpt) * ob4pi

      g_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vly ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfy ( ix, iquadpt)**2 * ob4pi

      g_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         - bfy ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi

      g_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vly ( ix, iquadpt)
     1                       - bfy ( ix, iquadpt) * tempab

      g_quadpt ( ix, 6, iquadpt) =
     1           vly ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         - vlx ( ix, iquadpt) * bfy ( ix, iquadpt)

      g_quadpt ( ix, 7, iquadpt) = 0.0

      g_quadpt ( ix, 8, iquadpt) =
     1           vly ( ix, iquadpt) * bfz ( ix, iquadpt)
     1         - vlz ( ix, iquadpt) * bfy ( ix, iquadpt)

#if ( NFLUID >= 1)
      g_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vly ( ix, iquadpt)
#endif

! -----

      h_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 4, iquadpt)

      h_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         - bfz ( ix, iquadpt) * bfx ( ix, iquadpt) * ob4pi

      h_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vly ( ix, iquadpt)
     1         - bfz ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi

      h_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfz ( ix, iquadpt)**2 * ob4pi

      h_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vlz ( ix, iquadpt)
     1                       - bfz ( ix, iquadpt) * tempab

      h_quadpt ( ix, 6, iquadpt) =
     1           vlz ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         - vlx ( ix, iquadpt) * bfz ( ix, iquadpt)

      h_quadpt ( ix, 7, iquadpt) =
     1           vlz ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         - vly ( ix, iquadpt) * bfz ( ix, iquadpt)

      h_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      h_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlz ( ix, iquadpt)
#endif

! -----

#if ( ISOTHERMAL_EQNS == 1)
      f_quadpt ( ix, 5, iquadpt) = 0.0
      g_quadpt ( ix, 5, iquadpt) = 0.0
      h_quadpt ( ix, 5, iquadpt) = 0.0
#endif

! -----

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! -----

      f_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)

      f_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vlx_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      f_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vlx_euler ( ix, iquadpt, i_euler)
#endif

! -----

      g_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)

      g_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vly_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      g_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vly_euler ( ix, iquadpt, i_euler)
#endif

! -----

      h_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)

      h_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vlz_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      h_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vlz_euler ( ix, iquadpt, i_euler)
#endif

! -----

#if ( ISOTHERMAL_EQNS == 1)
      f_quadpt ( ix, ivar + 5, iquadpt) = 0.0
      g_quadpt ( ix, ivar + 5, iquadpt) = 0.0
      h_quadpt ( ix, ivar + 5, iquadpt) = 0.0
#endif

! -----

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ------------------------------

! To make cylindrical and spherical geometries accessible to non-rel.
! and relativistic flow, we put this part after either flux construction.

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      IF ( igeom == 2 .or. igeom == 3) THEN

         g_quadpt ( ix, :, iquadpt) = g_quadpt ( ix, :, iquadpt) * obx1

      END IF

      IF ( igeom == 3) THEN

         sin_x2 = 1.0 / ( x1 * sin ( x2))

         h_quadpt ( ix, :, iquadpt) = h_quadpt ( ix, :, iquadpt)
     1                              * sin_x2

      END IF

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 8.b) Optional step needed only for viscous, resistive or conductive MHD.

! ASSUME that "du_quadpt_dx, du_quadpt_dy and du_quadpt_dz" have been built at
! all orders (NOT DONE YET).
! We set about building the derivatives of the primitive
! variables and various auxiliary variables at the nodal points. These are
! needed only in the presence of viscous, resistive or conduction terms.

! These derivatives are then used to build non-ideal fluxes at the quadrature
! points given by "fni_quadpt, gni_quadpt and hni_quadpt".

! Lastly, the non-ideal fluxes are added into the ideal fluxes

! ----------------------------------------------------------------------

      IF ( ( ( viscosity_coef .GT. 0.0) .OR.
     1       ( resistivity_coef .GT. 0.0)) .AND. ( NUMORDER >= 3)) THEN

! ------------------------------

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = 1.0 / rho ( ix, iquadpt)

      drho_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 1, iquadpt)
      drho_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 1, iquadpt)
      drho_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 1, iquadpt)


      dvlx_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvlx_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvlx_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dvly_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvly_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvly_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dvlz_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvlz_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvlz_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dbfx_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 6, iquadpt)
      dbfx_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 6, iquadpt)
      dbfx_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 6, iquadpt)

      dbfy_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 7, iquadpt)
      dbfy_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 7, iquadpt)
      dbfy_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 7, iquadpt)

      dbfz_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 8, iquadpt)
      dbfz_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 8, iquadpt)
      dbfz_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 8, iquadpt)


      dv2_dx ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dx ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dx ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dx ( ix, iquadpt) )

      dv2_dy ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dy ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dy ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dy ( ix, iquadpt) )

      dv2_dz ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dz ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dz ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dz ( ix, iquadpt) )

      db2_dx ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dx ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dx ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dx ( ix, iquadpt) )

      db2_dy ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dy ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dy ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dy ( ix, iquadpt) )

      db2_dz ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dz ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dz ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dz ( ix, iquadpt) )


      tempac = vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1       + vlz ( ix, iquadpt)**2

      tempab = molewt_quadpt ( ix, 1, iquadpt)
     1       / ( gasconst * rho ( ix, iquadpt) )

      temperature ( ix, iquadpt) = prs ( ix, iquadpt) * tempab

      tempab = tempab * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      dtemperature_dx ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dx ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dx ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dx ( ix, iquadpt)
     1    - ob8pi * db2_dx ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dx ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

      dtemperature_dy ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dy ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dy ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dy ( ix, iquadpt)
     1    - ob8pi * db2_dy ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dy ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

      dtemperature_dz ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dz ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dz ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dz ( ix, iquadpt)
     1    - ob8pi * db2_dz ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dz ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

! The derivatives are within the reference element's local coordinates.
! To get true physical fluxes, we need to divide by the size of the zone here.

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = 1.0 / dx ( ix)
      tempab = 1.0 / dy ( ix)
      tempac = 1.0 / dz ( ix)

      drho_dx ( ix, iquadpt) = tempaa * drho_dx ( ix, iquadpt)
      dtemperature_dx ( ix, iquadpt) =
     1                 tempaa * dtemperature_dx ( ix, iquadpt)
      dvlx_dx ( ix, iquadpt) = tempaa * dvlx_dx ( ix, iquadpt)
      dvly_dx ( ix, iquadpt) = tempaa * dvly_dx ( ix, iquadpt)
      dvlz_dx ( ix, iquadpt) = tempaa * dvlz_dx ( ix, iquadpt)
      dbfx_dx ( ix, iquadpt) = tempaa * dbfx_dx ( ix, iquadpt)
      dbfy_dx ( ix, iquadpt) = tempaa * dbfy_dx ( ix, iquadpt)
      dbfz_dx ( ix, iquadpt) = tempaa * dbfz_dx ( ix, iquadpt)
      dv2_dx ( ix, iquadpt) = tempaa * dv2_dx ( ix, iquadpt)
      db2_dx ( ix, iquadpt) = tempaa * db2_dx ( ix, iquadpt)

      drho_dy ( ix, iquadpt) = tempab * drho_dy ( ix, iquadpt)
      dtemperature_dy ( ix, iquadpt) =
     1                 tempab * dtemperature_dy ( ix, iquadpt)
      dvlx_dy ( ix, iquadpt) = tempab * dvlx_dy ( ix, iquadpt)
      dvly_dy ( ix, iquadpt) = tempab * dvly_dy ( ix, iquadpt)
      dvlz_dy ( ix, iquadpt) = tempab * dvlz_dy ( ix, iquadpt)
      dbfx_dy ( ix, iquadpt) = tempab * dbfx_dy ( ix, iquadpt)
      dbfy_dy ( ix, iquadpt) = tempab * dbfy_dy ( ix, iquadpt)
      dbfz_dy ( ix, iquadpt) = tempab * dbfz_dy ( ix, iquadpt)
      dv2_dy ( ix, iquadpt) = tempab * dv2_dy ( ix, iquadpt)
      db2_dy ( ix, iquadpt) = tempab * db2_dy ( ix, iquadpt)

      drho_dz ( ix, iquadpt) = tempac * drho_dz ( ix, iquadpt)
      dtemperature_dz ( ix, iquadpt) =
     1                 tempac * dtemperature_dz ( ix, iquadpt)
      dvlx_dz ( ix, iquadpt) = tempac * dvlx_dz ( ix, iquadpt)
      dvly_dz ( ix, iquadpt) = tempac * dvly_dz ( ix, iquadpt)
      dvlz_dz ( ix, iquadpt) = tempac * dvlz_dz ( ix, iquadpt)
      dbfx_dz ( ix, iquadpt) = tempac * dbfx_dz ( ix, iquadpt)
      dbfy_dz ( ix, iquadpt) = tempac * dbfy_dz ( ix, iquadpt)
      dbfz_dz ( ix, iquadpt) = tempac * dbfz_dz ( ix, iquadpt)
      dv2_dz ( ix, iquadpt) = tempac * dv2_dz ( ix, iquadpt)
      db2_dz ( ix, iquadpt) = tempac * db2_dz ( ix, iquadpt)

      div_vl ( ix, iquadpt) = dvlx_dx ( ix, iquadpt)
     1                      + dvly_dy ( ix, iquadpt)
     1                      + dvlz_dz ( ix, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = rho ( ix, iquadpt) * viscosity_coef

! ----------

      fni_quadpt ( ix, 1, iquadpt) = 0.0

      fni_quadpt ( ix, 2, iquadpt) =
     1               2.0 * tempaa * ( dvlx_dx ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      fni_quadpt ( ix, 3, iquadpt) =
     1                     tempaa * ( dvlx_dy ( ix, iquadpt)
     1                              + dvly_dx ( ix, iquadpt) )

      fni_quadpt ( ix, 4, iquadpt) =
     1                     tempaa * ( dvlx_dz ( ix, iquadpt)
     1                              + dvlz_dx ( ix, iquadpt) )

      fni_quadpt ( ix, 5, iquadpt) =
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vlx ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvlx_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvlx_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvlx_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dx ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dx ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfx_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfx_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfx_dz ( ix, iquadpt) )

      fni_quadpt ( ix, 6, iquadpt) = 0.0

      fni_quadpt ( ix, 7, iquadpt) =
     1   resistivity_coef * ( dbfy_dx ( ix, iquadpt)
     1                      - dbfx_dy ( ix, iquadpt) )

      fni_quadpt ( ix, 8, iquadpt) =
     1   resistivity_coef * ( dbfz_dx ( ix, iquadpt)
     1                      - dbfx_dz ( ix, iquadpt) )

#if ( NVAR > 8)
      fni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      gni_quadpt ( ix, 1, iquadpt) = 0.0

      gni_quadpt ( ix, 2, iquadpt) =
     1                     tempaa * ( dvly_dx ( ix, iquadpt)
     1                              + dvlx_dy ( ix, iquadpt) )

      gni_quadpt ( ix, 3, iquadpt) =
     1               2.0 * tempaa * ( dvly_dy ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      gni_quadpt ( ix, 4, iquadpt) =
     1                     tempaa * ( dvly_dz ( ix, iquadpt)
     1                              + dvlz_dy ( ix, iquadpt) )

      gni_quadpt ( ix, 5, iquadpt) = 
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vly ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvly_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvly_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvly_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dy ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dy ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfy_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfy_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfy_dz ( ix, iquadpt) )

      gni_quadpt ( ix, 6, iquadpt) =
     1   resistivity_coef * ( dbfx_dy ( ix, iquadpt)
     1                      - dbfy_dx ( ix, iquadpt) )

      gni_quadpt ( ix, 7, iquadpt) = 0.0

      gni_quadpt ( ix, 8, iquadpt) =
     1   resistivity_coef * ( dbfz_dy ( ix, iquadpt)
     1                      - dbfy_dz ( ix, iquadpt) )

#if ( NVAR > 8)
      gni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      hni_quadpt ( ix, 1, iquadpt) = 0.0

      hni_quadpt ( ix, 2, iquadpt) =
     1                     tempaa * ( dvlz_dx ( ix, iquadpt)
     1                              + dvlx_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 3, iquadpt) =
     1                     tempaa * ( dvlz_dy ( ix, iquadpt)
     1                              + dvly_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 4, iquadpt) =
     1               2.0 * tempaa * ( dvlz_dz ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      hni_quadpt ( ix, 5, iquadpt) =
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vlz ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvlz_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvlz_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvlz_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dz ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dz ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfz_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfz_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfz_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 6, iquadpt) =
     1   resistivity_coef * ( dbfx_dz ( ix, iquadpt)
     1                      - dbfz_dx ( ix, iquadpt) )

      hni_quadpt ( ix, 7, iquadpt) =
     1   resistivity_coef * ( dbfy_dz ( ix, iquadpt)
     1                      - dbfz_dy ( ix, iquadpt) )

      hni_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NVAR > 8)
      hni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      f_quadpt ( ix, :, iquadpt) = f_quadpt ( ix, :, iquadpt)
     1                       + fni_quadpt ( ix, :, iquadpt)

      g_quadpt ( ix, :, iquadpt) = g_quadpt ( ix, :, iquadpt)
     1                       + gni_quadpt ( ix, :, iquadpt)

      h_quadpt ( ix, :, iquadpt) = h_quadpt ( ix, :, iquadpt)
     1                       + hni_quadpt ( ix, :, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      END IF
! End of "( viscosity_coef .GT. 0.0) .OR. ( resistivity_coef .GT. 0.0)".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = dtmove / dx ( ix)
      f_quadpt ( ix, :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS) =
     1   tempaa
     1 * f_quadpt ( ix, :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS)

      tempab = dtmove / dy ( ix)
      g_quadpt ( ix, :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS) =
     1   tempab
     1 * g_quadpt ( ix, :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS)

      tempac = dtmove / dz ( ix)
      h_quadpt ( ix, :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS) =
     1   tempac
     1 * h_quadpt ( ix, :, NUM_ZONECEN_QUADPTS+1: NUM_SPACETIME_QUADPTS)

! ----------

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------

! STEP 8.c) Now project the fluxes "( f_quadpt, g_quadpt, h_quadpt)"
! from the t>0 space-time nodal points back to the modal
! space time variables for those variables, i.e. to "( f, g, h)" .

! Note that some of the "f/g/h_quadpt" values with t>0 will be overwritten.
! But that is ok since these are also constructed anew in each ADER iteration.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

!      f ( ix, :, 05) = f_quadpt ( ix, :, 07) - f ( ix, :, 01)

! ----------

!      g ( ix, :, 05) = g_quadpt ( ix, :, 07) - g ( ix, :, 01)

! ----------

!      h ( ix, :, 05) = h_quadpt ( ix, :, 07) - h ( ix, :, 01)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      f ( ix, :, 13) = 2.0 * ( f_quadpt ( ix, :, 14)
     1 - f_quadpt ( ix, :, 15) - f ( ix, :, 02) )

!      f ( ix, :, 14) = 2.0 * ( f_quadpt ( ix, :, 16)
!     1 - f_quadpt ( ix, :, 17) - f ( ix, :, 03) )

!      f ( ix, :, 15) = 2.0 * ( f_quadpt ( ix, :, 18)
!     1 - f_quadpt ( ix, :, 19) - f ( ix, :, 04) )

!      f ( ix, :, 11) = 2.0 * ( f_quadpt ( ix, :, 14)
!     1 + f_quadpt ( ix, :, 15) + f_quadpt ( ix, :, 16)
!     1 + f_quadpt ( ix, :, 17) + f_quadpt ( ix, :, 18)
!     1 + f_quadpt ( ix, :, 19) - 6.0 * f ( ix, :, 01) ) / 3.0
!     1 - f_quadpt ( ix, :, 20) + f_quadpt ( ix, :, 01)

!      f ( ix, :, 12) = f_quadpt ( ix, :, 20) - f_quadpt ( ix, :, 01)
!     1 - f ( ix, :, 11)

! ----------

!      g ( ix, :, 13) = 2.0 * ( g_quadpt ( ix, :, 14)
!     1 - g_quadpt ( ix, :, 15) - g ( ix, :, 02) )

      g ( ix, :, 14) = 2.0 * ( g_quadpt ( ix, :, 16)
     1 - g_quadpt ( ix, :, 17) - g ( ix, :, 03) )

!      g ( ix, :, 15) = 2.0 * ( g_quadpt ( ix, :, 18)
!     1 - g_quadpt ( ix, :, 19) - g ( ix, :, 04) )

!      g ( ix, :, 11) = 2.0 * ( g_quadpt ( ix, :, 14)
!     1 + g_quadpt ( ix, :, 15) + g_quadpt ( ix, :, 16)
!     1 + g_quadpt ( ix, :, 17) + g_quadpt ( ix, :, 18)
!     1 + g_quadpt ( ix, :, 19) - 6.0 * g ( ix, :, 01) ) / 3.0
!     1 - g_quadpt ( ix, :, 20) + g_quadpt ( ix, :, 01)

!      g ( ix, :, 12) = g_quadpt ( ix, :, 20) - g_quadpt ( ix, :, 01)
!     1 - g ( ix, :, 11)

! ----------

!      h ( ix, :, 13) = 2.0 * ( h_quadpt ( ix, :, 14)
!     1 - h_quadpt ( ix, :, 15) - h ( ix, :, 02) )

!      h ( ix, :, 14) = 2.0 * ( h_quadpt ( ix, :, 16)
!     1 - h_quadpt ( ix, :, 17) - h ( ix, :, 03) )

      h ( ix, :, 15) = 2.0 * ( h_quadpt ( ix, :, 18)
     1 - h_quadpt ( ix, :, 19) - h ( ix, :, 04) )

!      h ( ix, :, 11) = 2.0 * ( h_quadpt ( ix, :, 14)
!     1 + h_quadpt ( ix, :, 15) + h_quadpt ( ix, :, 16)
!     1 + h_quadpt ( ix, :, 17) + h_quadpt ( ix, :, 18)
!     1 + h_quadpt ( ix, :, 19) - 6.0 * h ( ix, :, 01) ) / 3.0
!     1 - h_quadpt ( ix, :, 20) + h_quadpt ( ix, :, 01)

!      h ( ix, :, 12) = h_quadpt ( ix, :, 20) - h_quadpt ( ix, :, 01)
!     1 - h ( ix, :, 11)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tf ( ix, :, 07) = f_quadpt ( ix, :, 23) + f_quadpt ( ix, :, 24)
     1 - f_quadpt ( ix, :, 25) - f_quadpt ( ix, :, 26)
     1 - f_quadpt ( ix, :, 27) - f_quadpt ( ix, :, 28)
     1 + f_quadpt ( ix, :, 29) + f_quadpt ( ix, :, 30)

      tf ( ix, :, 08) = f_quadpt ( ix, :, 32) + f_quadpt ( ix, :, 34)
     1 - f_quadpt ( ix, :, 24) - f_quadpt ( ix, :, 28)
     1 - f_quadpt ( ix, :, 31) - f_quadpt ( ix, :, 33)
     1 + f_quadpt ( ix, :, 23) + f_quadpt ( ix, :, 27)

      tf ( ix, :, 09) = f_quadpt ( ix, :, 25) + f_quadpt ( ix, :, 29)
     1 - f_quadpt ( ix, :, 31) - f_quadpt ( ix, :, 32)
     1 - f_quadpt ( ix, :, 26) - f_quadpt ( ix, :, 30)
     1 + f_quadpt ( ix, :, 33) + f_quadpt ( ix, :, 34)

      f ( ix, :, 16) = tf ( ix, :, 07) + 0.5 * ( tf ( ix, :, 01)
     1 - tf ( ix, :, 02) ) - 2.0 * ( f_quadpt ( ix, :, 02)
     1 - f_quadpt ( ix, :, 05) )

      f ( ix, :, 18) = - 2.0 * tf ( ix, :, 07) + f ( ix, :, 16)

      f ( ix, :, 14) = tf ( ix, :, 08) + 0.5 * ( tf ( ix, :, 03)
     1 - tf ( ix, :, 04) ) - 2.0 * ( f_quadpt ( ix, :, 06)
     1 - f_quadpt ( ix, :, 09) )

      f ( ix, :, 19) = - 2.0 * tf ( ix, :, 08) + f ( ix, :, 14)

      f ( ix, :, 15) = tf ( ix, :, 09) + 0.5 * ( tf ( ix, :, 05)
     1 - tf ( ix, :, 06) ) - 2.0 * ( f_quadpt ( ix, :, 10)
     1 - f_quadpt ( ix, :, 13) )

      f ( ix, :, 17) = - 2.0 * tf ( ix, :, 09) + f ( ix, :, 15)


      tf ( ix, :, 10) = f_quadpt ( ix, :, 23) + f_quadpt ( ix, :, 24)
     1 + f_quadpt ( ix, :, 25) + f_quadpt ( ix, :, 26)
     1 - f_quadpt ( ix, :, 27) - f_quadpt ( ix, :, 28)
     1 - f_quadpt ( ix, :, 29) - f_quadpt ( ix, :, 30)

!      tf ( ix, :, 11) = f_quadpt ( ix, :, 23) + f_quadpt ( ix, :, 27)
!     1 + f_quadpt ( ix, :, 31) + f_quadpt ( ix, :, 33)
!     1 - f_quadpt ( ix, :, 24) - f_quadpt ( ix, :, 28)
!     1 - f_quadpt ( ix, :, 32) - f_quadpt ( ix, :, 34)

!      tf ( ix, :, 12) = f_quadpt ( ix, :, 25) + f_quadpt ( ix, :, 29)
!     1 + f_quadpt ( ix, :, 31) + f_quadpt ( ix, :, 32)
!     1 - f_quadpt ( ix, :, 26) - f_quadpt ( ix, :, 30)
!     1 - f_quadpt ( ix, :, 33) - f_quadpt ( ix, :, 34)

      f ( ix, :, 27) = 1.125 * ( tf ( ix, :, 01) - tf ( ix, :, 02) )
     1 + 4.5 * ( f_quadpt ( ix, :, 36) - f_quadpt ( ix, :, 37) )
     1 - 2.25 * tf ( ix, :, 10)

      f ( ix, :, 24) = 3.0 * ( f_quadpt ( ix, :, 36)
     1 - f_quadpt ( ix, :, 37) ) - 0.75 * tf ( ix, :, 10)
     1 + 0.375 * ( f ( ix, :, 16) + f ( ix, :, 18) )
     1 - f ( ix, :, 27)

!      f ( ix, :, 28) = 1.125 * ( tf ( ix, :, 03) - tf ( ix, :, 04) )
!     1 + 4.5 * ( f_quadpt ( ix, :, 38) - f_quadpt ( ix, :, 39) )
!     1 - 2.25 * tf ( ix, :, 11)

!      f ( ix, :, 25) = 3.0 * ( f_quadpt ( ix, :, 38)
!     1 - f_quadpt ( ix, :, 39) ) - 0.75 * tf ( ix, :, 11)
!     1 + 0.375 * ( f ( ix, :, 14) + f ( ix, :, 19) )
!     1 - f ( ix, :, 28)

!      f ( ix, :, 29) = 1.125 * ( tf ( ix, :, 05) - tf ( ix, :, 06) )
!     1 + 4.5 * ( f_quadpt ( ix, :, 40) - f_quadpt ( ix, :, 41) )
!     1 - 2.25 * tf ( ix, :, 12)

!      f ( ix, :, 26) = 3.0 * ( f_quadpt ( ix, :, 40)
!     1 - f_quadpt ( ix, :, 41) ) - 0.75 * tf ( ix, :, 12)
!     1 + 0.375 * ( f ( ix, :, 15) + f ( ix, :, 17) )
!     1 - f ( ix, :, 29)


      f ( ix, :, 30) = 3.0 * ( f_quadpt ( ix, :, 36)
     1 - 2.0 * f_quadpt ( ix, :, 35) + f_quadpt ( ix, :, 37) )
     1 - 1.5 * f ( ix, :, 05)

!      f ( ix, :, 31) = 3.0 * ( f_quadpt ( ix, :, 38)
!     1 - 2.0 * f_quadpt ( ix, :, 35) + f_quadpt ( ix, :, 39) )
!     1 - 1.5 * f ( ix, :, 06)

!      f ( ix, :, 32) = 3.0 * ( f_quadpt ( ix, :, 40)
!     1 - 2.0 * f_quadpt ( ix, :, 35) + f_quadpt ( ix, :, 41) )
!     1 - 1.5 * f ( ix, :, 07)

      f ( ix, :, 33) = 3.0 * ( f_quadpt ( ix, :, 23)
     1 - f_quadpt ( ix, :, 27) - f_quadpt ( ix, :, 24)
     1 + f_quadpt ( ix, :, 28) ) - 3.0 * f ( ix, :, 08)

!      f ( ix, :, 34) = 3.0 * ( f_quadpt ( ix, :, 31)
!     1 - f_quadpt ( ix, :, 33) - f_quadpt ( ix, :, 32)
!     1 + f_quadpt ( ix, :, 34) ) - 3.0 * f ( ix, :, 09)

      f ( ix, :, 35) = 3.0 * ( f_quadpt ( ix, :, 25)
     1 - f_quadpt ( ix, :, 29) - f_quadpt ( ix, :, 26)
     1 + f_quadpt ( ix, :, 30) ) - 3.0 * f ( ix, :, 10)


!      f ( ix, :, 23) = 4.5 * ( f_quadpt ( ix, :, 42)
!     1 - f_quadpt ( ix, :, 01) ) - 13.5 * ( f_quadpt ( ix, :, 35)
!     1 - f_quadpt ( ix, :, 22) )

!      f ( ix, :, 22) = 2.25 * ( f_quadpt ( ix, :, 42)
!     1 + f_quadpt ( ix, :, 01) - f_quadpt ( ix, :, 35)
!     1 - f_quadpt ( ix, :, 22) ) - 1.5 * f ( ix, :, 23)

!      f ( ix, :, 21) = f_quadpt ( ix, :, 42) - f_quadpt ( ix, :, 01)
!     1 - f ( ix, :, 22) - f ( ix, :, 23)
!     1 + ( f ( ix, :, 30) + f ( ix, :, 31)
!     1 + f ( ix, :, 32) ) / 12.0

! ----------

      tg ( ix, :, 07) = g_quadpt ( ix, :, 23) + g_quadpt ( ix, :, 24)
     1 - g_quadpt ( ix, :, 25) - g_quadpt ( ix, :, 26)
     1 - g_quadpt ( ix, :, 27) - g_quadpt ( ix, :, 28)
     1 + g_quadpt ( ix, :, 29) + g_quadpt ( ix, :, 30)

      tg ( ix, :, 08) = g_quadpt ( ix, :, 32) + g_quadpt ( ix, :, 34)
     1 - g_quadpt ( ix, :, 24) - g_quadpt ( ix, :, 28)
     1 - g_quadpt ( ix, :, 31) - g_quadpt ( ix, :, 33)
     1 + g_quadpt ( ix, :, 23) + g_quadpt ( ix, :, 27)

      tg ( ix, :, 09) = g_quadpt ( ix, :, 25) + g_quadpt ( ix, :, 29)
     1 - g_quadpt ( ix, :, 31) - g_quadpt ( ix, :, 32)
     1 - g_quadpt ( ix, :, 26) - g_quadpt ( ix, :, 30)
     1 + g_quadpt ( ix, :, 33) + g_quadpt ( ix, :, 34)

      g ( ix, :, 16) = tg ( ix, :, 07) + 0.5 * ( tg ( ix, :, 01)
     1 - tg ( ix, :, 02) ) - 2.0 * ( g_quadpt ( ix, :, 02)
     1 - g_quadpt ( ix, :, 05) )

      g ( ix, :, 18) = - 2.0 * tg ( ix, :, 07) + g ( ix, :, 16)

      g ( ix, :, 14) = tg ( ix, :, 08) + 0.5 * ( tg ( ix, :, 03)
     1 - tg ( ix, :, 04) ) - 2.0 * ( g_quadpt ( ix, :, 06)
     1 - g_quadpt ( ix, :, 09) )

      g ( ix, :, 19) = - 2.0 * tg ( ix, :, 08) + g ( ix, :, 14)

      g ( ix, :, 15) = tg ( ix, :, 09) + 0.5 * ( tg ( ix, :, 05)
     1 - tg ( ix, :, 06) ) - 2.0 * ( g_quadpt ( ix, :, 10)
     1 - g_quadpt ( ix, :, 13) )

      g ( ix, :, 17) = - 2.0 * tg ( ix, :, 09) + g ( ix, :, 15)


!      tg ( ix, :, 10) = g_quadpt ( ix, :, 23) + g_quadpt ( ix, :, 24)
!     1 + g_quadpt ( ix, :, 25) + g_quadpt ( ix, :, 26)
!     1 - g_quadpt ( ix, :, 27) - g_quadpt ( ix, :, 28)
!     1 - g_quadpt ( ix, :, 29) - g_quadpt ( ix, :, 30)

      tg ( ix, :, 11) = g_quadpt ( ix, :, 23) + g_quadpt ( ix, :, 27)
     1 + g_quadpt ( ix, :, 31) + g_quadpt ( ix, :, 33)
     1 - g_quadpt ( ix, :, 24) - g_quadpt ( ix, :, 28)
     1 - g_quadpt ( ix, :, 32) - g_quadpt ( ix, :, 34)

!      tg ( ix, :, 12) = g_quadpt ( ix, :, 25) + g_quadpt ( ix, :, 29)
!     1 + g_quadpt ( ix, :, 31) + g_quadpt ( ix, :, 32)
!     1 - g_quadpt ( ix, :, 26) - g_quadpt ( ix, :, 30)
!     1 - g_quadpt ( ix, :, 33) - g_quadpt ( ix, :, 34)

!      g ( ix, :, 27) = 1.125 * ( tg ( ix, :, 01) - tg ( ix, :, 02) )
!     1 + 4.5 * ( g_quadpt ( ix, :, 36) - g_quadpt ( ix, :, 37) )
!     1 - 2.25 * tg ( ix, :, 10)

!      g ( ix, :, 24) = 3.0 * ( g_quadpt ( ix, :, 36)
!     1 - g_quadpt ( ix, :, 37) ) - 0.75 * tg ( ix, :, 10)
!     1 + 0.375 * ( g ( ix, :, 16) + g ( ix, :, 18) )
!     1 - g ( ix, :, 27)

      g ( ix, :, 28) = 1.125 * ( tg ( ix, :, 03) - tg ( ix, :, 04) )
     1 + 4.5 * ( g_quadpt ( ix, :, 38) - g_quadpt ( ix, :, 39) )
     1 - 2.25 * tg ( ix, :, 11)

      g ( ix, :, 25) = 3.0 * ( g_quadpt ( ix, :, 38)
     1 - g_quadpt ( ix, :, 39) ) - 0.75 * tg ( ix, :, 11)
     1 + 0.375 * ( g ( ix, :, 14) + g ( ix, :, 19) )
     1 - g ( ix, :, 28)

!      g ( ix, :, 29) = 1.125 * ( tg ( ix, :, 05) - tg ( ix, :, 06) )
!     1 + 4.5 * ( g_quadpt ( ix, :, 40) - g_quadpt ( ix, :, 41) )
!     1 - 2.25 * tg ( ix, :, 12)

!      g ( ix, :, 26) = 3.0 * ( g_quadpt ( ix, :, 40)
!     1 - g_quadpt ( ix, :, 41) ) - 0.75 * tg ( ix, :, 12)
!     1 + 0.375 * ( g ( ix, :, 15) + g ( ix, :, 17) )
!     1 - g ( ix, :, 29)


!      g ( ix, :, 30) = 3.0 * ( g_quadpt ( ix, :, 36)
!     1 - 2.0 * g_quadpt ( ix, :, 35) + g_quadpt ( ix, :, 37) )
!     1 - 1.5 * g ( ix, :, 05)

      g ( ix, :, 31) = 3.0 * ( g_quadpt ( ix, :, 38)
     1 - 2.0 * g_quadpt ( ix, :, 35) + g_quadpt ( ix, :, 39) )
     1 - 1.5 * g ( ix, :, 06)

!      g ( ix, :, 32) = 3.0 * ( g_quadpt ( ix, :, 40)
!     1 - 2.0 * g_quadpt ( ix, :, 35) + g_quadpt ( ix, :, 41) )
!     1 - 1.5 * g ( ix, :, 07)

      g ( ix, :, 33) = 3.0 * ( g_quadpt ( ix, :, 23)
     1 - g_quadpt ( ix, :, 27) - g_quadpt ( ix, :, 24)
     1 + g_quadpt ( ix, :, 28) ) - 3.0 * g ( ix, :, 08)

      g ( ix, :, 34) = 3.0 * ( g_quadpt ( ix, :, 31)
     1 - g_quadpt ( ix, :, 33) - g_quadpt ( ix, :, 32)
     1 + g_quadpt ( ix, :, 34) ) - 3.0 * g ( ix, :, 09)

!      g ( ix, :, 35) = 3.0 * ( g_quadpt ( ix, :, 25)
!     1 - g_quadpt ( ix, :, 29) - g_quadpt ( ix, :, 26)
!     1 + g_quadpt ( ix, :, 30) ) - 3.0 * g ( ix, :, 10)


!      g ( ix, :, 23) = 4.5 * ( g_quadpt ( ix, :, 42)
!     1 - g_quadpt ( ix, :, 01) ) - 13.5 * ( g_quadpt ( ix, :, 35)
!     1 - g_quadpt ( ix, :, 22) )

!      g ( ix, :, 22) = 2.25 * ( g_quadpt ( ix, :, 42)
!     1 + g_quadpt ( ix, :, 01) - g_quadpt ( ix, :, 35)
!     1 - g_quadpt ( ix, :, 22) ) - 1.5 * g ( ix, :, 23)

!      g ( ix, :, 21) = g_quadpt ( ix, :, 42) - g_quadpt ( ix, :, 01)
!     1 - g ( ix, :, 22) - g ( ix, :, 23)
!     1 + ( g ( ix, :, 30) + g ( ix, :, 31)
!     1 + g ( ix, :, 32) ) / 12.0

! ----------

      th ( ix, :, 07) = h_quadpt ( ix, :, 23) + h_quadpt ( ix, :, 24)
     1 - h_quadpt ( ix, :, 25) - h_quadpt ( ix, :, 26)
     1 - h_quadpt ( ix, :, 27) - h_quadpt ( ix, :, 28)
     1 + h_quadpt ( ix, :, 29) + h_quadpt ( ix, :, 30)

      th ( ix, :, 08) = h_quadpt ( ix, :, 32) + h_quadpt ( ix, :, 34)
     1 - h_quadpt ( ix, :, 24) - h_quadpt ( ix, :, 28)
     1 - h_quadpt ( ix, :, 31) - h_quadpt ( ix, :, 33)
     1 + h_quadpt ( ix, :, 23) + h_quadpt ( ix, :, 27)

      th ( ix, :, 09) = h_quadpt ( ix, :, 25) + h_quadpt ( ix, :, 29)
     1 - h_quadpt ( ix, :, 31) - h_quadpt ( ix, :, 32)
     1 - h_quadpt ( ix, :, 26) - h_quadpt ( ix, :, 30)
     1 + h_quadpt ( ix, :, 33) + h_quadpt ( ix, :, 34)

      h ( ix, :, 16) = th ( ix, :, 07) + 0.5 * ( th ( ix, :, 01)
     1 - th ( ix, :, 02) ) - 2.0 * ( h_quadpt ( ix, :, 02)
     1 - h_quadpt ( ix, :, 05) )

      h ( ix, :, 18) = - 2.0 * th ( ix, :, 07) + h ( ix, :, 16)

      h ( ix, :, 14) = th ( ix, :, 08) + 0.5 * ( th ( ix, :, 03)
     1 - th ( ix, :, 04) ) - 2.0 * ( h_quadpt ( ix, :, 06)
     1 - h_quadpt ( ix, :, 09) )

      h ( ix, :, 19) = - 2.0 * th ( ix, :, 08) + h ( ix, :, 14)

      h ( ix, :, 15) = th ( ix, :, 09) + 0.5 * ( th ( ix, :, 05)
     1 - th ( ix, :, 06) ) - 2.0 * ( h_quadpt ( ix, :, 10)
     1 - h_quadpt ( ix, :, 13) )

      h ( ix, :, 17) = - 2.0 * th ( ix, :, 09) + h ( ix, :, 15)


!      th ( ix, :, 10) = h_quadpt ( ix, :, 23) + h_quadpt ( ix, :, 24)
!     1 + h_quadpt ( ix, :, 25) + h_quadpt ( ix, :, 26)
!     1 - h_quadpt ( ix, :, 27) - h_quadpt ( ix, :, 28)
!     1 - h_quadpt ( ix, :, 29) - h_quadpt ( ix, :, 30)

!      th ( ix, :, 11) = h_quadpt ( ix, :, 23) + h_quadpt ( ix, :, 27)
!     1 + h_quadpt ( ix, :, 31) + h_quadpt ( ix, :, 33)
!     1 - h_quadpt ( ix, :, 24) - h_quadpt ( ix, :, 28)
!     1 - h_quadpt ( ix, :, 32) - h_quadpt ( ix, :, 34)

      th ( ix, :, 12) = h_quadpt ( ix, :, 25) + h_quadpt ( ix, :, 29)
     1 + h_quadpt ( ix, :, 31) + h_quadpt ( ix, :, 32)
     1 - h_quadpt ( ix, :, 26) - h_quadpt ( ix, :, 30)
     1 - h_quadpt ( ix, :, 33) - h_quadpt ( ix, :, 34)

!      h ( ix, :, 27) = 1.125 * ( th ( ix, :, 01) - th ( ix, :, 02) )
!     1 + 4.5 * ( h_quadpt ( ix, :, 36) - h_quadpt ( ix, :, 37) )
!     1 - 2.25 * th ( ix, :, 10)

!      h ( ix, :, 24) = 3.0 * ( h_quadpt ( ix, :, 36)
!     1 - h_quadpt ( ix, :, 37) ) - 0.75 * th ( ix, :, 10)
!     1 + 0.375 * ( h ( ix, :, 16) + h ( ix, :, 18) )
!     1 - h ( ix, :, 27)

!      h ( ix, :, 28) = 1.125 * ( th ( ix, :, 03) - th ( ix, :, 04) )
!     1 + 4.5 * ( h_quadpt ( ix, :, 38) - h_quadpt ( ix, :, 39) )
!     1 - 2.25 * th ( ix, :, 11)

!      h ( ix, :, 25) = 3.0 * ( h_quadpt ( ix, :, 38)
!     1 - h_quadpt ( ix, :, 39) ) - 0.75 * th ( ix, :, 11)
!     1 + 0.375 * ( h ( ix, :, 14) + h ( ix, :, 19) )
!     1 - h ( ix, :, 28)

      h ( ix, :, 29) = 1.125 * ( th ( ix, :, 05) - th ( ix, :, 06) )
     1 + 4.5 * ( h_quadpt ( ix, :, 40) - h_quadpt ( ix, :, 41) )
     1 - 2.25 * th ( ix, :, 12)

      h ( ix, :, 26) = 3.0 * ( h_quadpt ( ix, :, 40)
     1 - h_quadpt ( ix, :, 41) ) - 0.75 * th ( ix, :, 12)
     1 + 0.375 * ( h ( ix, :, 15) + h ( ix, :, 17) )
     1 - h ( ix, :, 29)


!      h ( ix, :, 30) = 3.0 * ( h_quadpt ( ix, :, 36)
!     1 - 2.0 * h_quadpt ( ix, :, 35) + h_quadpt ( ix, :, 37) )
!     1 - 1.5 * h ( ix, :, 05)

!      h ( ix, :, 31) = 3.0 * ( h_quadpt ( ix, :, 38)
!     1 - 2.0 * h_quadpt ( ix, :, 35) + h_quadpt ( ix, :, 39) )
!     1 - 1.5 * h ( ix, :, 06)

      h ( ix, :, 32) = 3.0 * ( h_quadpt ( ix, :, 40)
     1 - 2.0 * h_quadpt ( ix, :, 35) + h_quadpt ( ix, :, 41) )
     1 - 1.5 * h ( ix, :, 07)

!      h ( ix, :, 33) = 3.0 * ( h_quadpt ( ix, :, 23)
!     1 - h_quadpt ( ix, :, 27) - h_quadpt ( ix, :, 24)
!     1 + h_quadpt ( ix, :, 28) ) - 3.0 * h ( ix, :, 08)

      h ( ix, :, 34) = 3.0 * ( h_quadpt ( ix, :, 31)
     1 - h_quadpt ( ix, :, 33) - h_quadpt ( ix, :, 32)
     1 + h_quadpt ( ix, :, 34) ) - 3.0 * h ( ix, :, 09)

      h ( ix, :, 35) = 3.0 * ( h_quadpt ( ix, :, 25)
     1 - h_quadpt ( ix, :, 29) - h_quadpt ( ix, :, 26)
     1 + h_quadpt ( ix, :, 30) ) - 3.0 * h ( ix, :, 10)


!      h ( ix, :, 23) = 4.5 * ( h_quadpt ( ix, :, 42)
!     1 - h_quadpt ( ix, :, 01) ) - 13.5 * ( h_quadpt ( ix, :, 35)
!     1 - h_quadpt ( ix, :, 22) )

!      h ( ix, :, 22) = 2.25 * ( h_quadpt ( ix, :, 42)
!     1 + h_quadpt ( ix, :, 01) - h_quadpt ( ix, :, 35)
!     1 - h_quadpt ( ix, :, 22) ) - 1.5 * h ( ix, :, 23)

!      h ( ix, :, 21) = h_quadpt ( ix, :, 42) - h_quadpt ( ix, :, 01)
!     1 - h ( ix, :, 22) - h ( ix, :, 23)
!     1 + ( h ( ix, :, 30) + h ( ix, :, 31)
!     1 + h ( ix, :, 32) ) / 12.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
! Default is to build source terms at t>0 . In some ADER-CG updates
! (on Cartesian meshes), we may not require any source terms.
#if ( BUILD_SOURCE_TERMS == 1)
! ----------------------------------------------------------------------

! STEP 9.a) Build source terms "s_quadpt" at the t>0 quadrature points.
! Remember that the primitive variables are still available from the previous
! loop at each of the t=0 quadrature points for constructing the source terms.

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      s_quadpt ( oned_begin - oned_off: oned_end + oned_off,
     1                                           :, iquadpt) = 0.0

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

#if ( 0 == 1)

! For single-fluid adiabatic MHD turbulence, we may want to provide a gradual
! (explicit) heating/cooling term which ensures that pressure always
! tends to a mean pressure "prs_bar" in a time "tau_bar" (comparable
! to a turnover time). Such a source term can be implemented explicitly
! and its implementation is given below.
! Remember to declare "prs_bar, tau_bar".

      prs_bar = 1.0
      tau_bar = 1.0

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 5, iquadpt) = - ( prs ( ix, iquadpt) - prs_bar)
     1     / ( tau_bar * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0) )

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

#endif

! ------------------------------

! Build two-fluid source terms. In this part, the source terms are
! hard-coded for one family of ionized fluid and one family of neutrals.

#if ( N_EULER == 1) 

      IF ( ( igeom .EQ. 1) .AND. ( N_EULER .EQ. 1) .AND.
     1     ( ambipolar_diffusion_drag_coef .GT. 0.0)) THEN

      ivar = NFIELD + 1 + NFLUID

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 1, iquadpt) = 0.0

      tempaa = ambipolar_diffusion_drag_coef
     1       * rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

      s_quadpt ( ix, 2, iquadpt) = tempaa
     1           * ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt))

      s_quadpt ( ix, 3, iquadpt) = tempaa
     1           * ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt))

      s_quadpt ( ix, 4, iquadpt) = tempaa
     1           * ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt))

      s_quadpt ( ix, 5, iquadpt) = ( ambipolar_diffusion_drag_coef
     1   / ( molewt_ion + molewt_neut) )

     1 * ( rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

     1   * ( ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt) )
     1     * ( molewt_neut * vlx_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlx ( ix, iquadpt) )
     1     + ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt) )
     1     * ( molewt_neut * vly_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vly ( ix, iquadpt) )
     1     + ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt) )
     1     * ( molewt_neut * vlz_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlz ( ix, iquadpt) )
     1     )

     1   + 3.0 * ( molewt_neut * prs_euler ( ix, iquadpt, 1)
     1           * rho ( ix, iquadpt)
     1           - molewt_ion * prs ( ix, iquadpt)
     1           * rho_euler ( ix, iquadpt, 1) )
     1   )

      s_quadpt ( ix, 6, iquadpt) = 0.0
      s_quadpt ( ix, 7, iquadpt) = 0.0
      s_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      s_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) = 0.0
#endif


      s_quadpt ( ix, ivar + 1, iquadpt) = 0.0
      s_quadpt ( ix, ivar + 2, iquadpt) = - s_quadpt ( ix, 2, iquadpt)
      s_quadpt ( ix, ivar + 3, iquadpt) = - s_quadpt ( ix, 3, iquadpt)
      s_quadpt ( ix, ivar + 4, iquadpt) = - s_quadpt ( ix, 4, iquadpt)
      s_quadpt ( ix, ivar + 5, iquadpt) = - s_quadpt ( ix, 5, iquadpt)

#if ( NFLUID_EULER >= 1)
      s_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

      END IF
! End of "ambipolar_diffusion_drag_coef .GT. 0.0" part

#endif
! End of "N_EULER == 1"

! ------------------------------

! Cylindrical Geometry

      IF ( igeom .EQ. 2) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1                           * vlx ( ix, iquadpt) * obx1

      s_quadpt ( ix, 2, iquadpt) =
     1 ( rho ( ix, iquadpt) * vly ( ix, iquadpt)**2
     1 - rho ( ix, iquadpt) * vlx ( ix, iquadpt)**2
     1 + bfx ( ix, iquadpt)**2 * ob4pi
     1 - bfy ( ix, iquadpt)**2 * ob4pi) * obx1

      s_quadpt ( ix, 3, iquadpt) =
     1 - 2.0 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                              * vly ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                        * vlz ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 5, iquadpt) = - ( tempac * vlx ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = 0.0

      s_quadpt ( ix, 7, iquadpt) = 0.0

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Spherical Geometry

      IF ( igeom .EQ. 3) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      cos_x2 = cos ( x2)
      sin_x2 = sin ( x2)
      cot_x2 = cos_x2 / sin_x2

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 2.0 * vlx ( ix, iquadpt)
     1 + cot_x2 * vly ( ix, iquadpt) ) * obx1

      s_quadpt ( ix, 2, iquadpt) = rho ( ix, iquadpt)
     1 * ( - 2.0 * vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1   + vlz ( ix, iquadpt)**2) * obx1
     1 - cot_x2
     1 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1   - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1
     1 - ( - 2.0 * bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1   + bfz ( ix, iquadpt)**2) * ob4pi * obx1

      s_quadpt ( ix, 3, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1 + cot_x2
     1 * ( vly ( ix, iquadpt)**2 - vlz ( ix, iquadpt)**2)) * obx1
     1 + ( 3.0 * bfx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1   + cot_x2
     1   * ( bfy ( ix, iquadpt)**2 - bfz ( ix, iquadpt)**2))
     1 * ob4pi * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - rho ( ix, iquadpt) * vlz ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) + 2.0 * cot_x2
     1   * vly ( ix, iquadpt) ) * obx1
     1 + bfz ( ix, iquadpt)
     1 * ( 3.0 * bfx ( ix, iquadpt) + 2.0 * cot_x2
     1   * bfy ( ix, iquadpt) ) * ob4pi * obx1

      s_quadpt ( ix, 5, iquadpt) =
     1 - 2.0 * ( tempac * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab) * obx1
     1 - cot_x2 * ( tempac * vly ( ix, iquadpt)
     1 - bfy ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = - cot_x2
     1                           * g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 7, iquadpt) = g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Multiply "s_quadpt" by "dtmove".

      DO ix = oned_begin - oned_off, oned_end + oned_off
      s_quadpt ( ix, :, NUM_ZONECEN_QUADPTS + 1: NUM_SPACETIME_QUADPTS)=
     1 dtmove
     1*s_quadpt ( ix, :, NUM_ZONECEN_QUADPTS + 1: NUM_SPACETIME_QUADPTS)
      END DO

! ------------------------------

! STEP 9.b) Now project the sources "s_quadpt" from the t>0 space-time nodal
! points back to the modal space time variables for those variables,
! i.e. to "s".

! Note that some of the "s_quadpt" values with t>0 will be overwritten.
! But that is ok since these are also constructed anew in each ADER iteration.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s ( ix, :, 05) = s_quadpt ( ix, :, 07) - s ( ix, :, 01)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s ( ix, :, 13) = 2.0 * ( s_quadpt ( ix, :, 14)
     1 - s_quadpt ( ix, :, 15) - s ( ix, :, 02) )

      s ( ix, :, 14) = 2.0 * ( s_quadpt ( ix, :, 16)
     1 - s_quadpt ( ix, :, 17) - s ( ix, :, 03) )

      s ( ix, :, 15) = 2.0 * ( s_quadpt ( ix, :, 18)
     1 - s_quadpt ( ix, :, 19) - s ( ix, :, 04) )

      s ( ix, :, 11) = 2.0 * ( s_quadpt ( ix, :, 14)
     1 + s_quadpt ( ix, :, 15) + s_quadpt ( ix, :, 16)
     1 + s_quadpt ( ix, :, 17) + s_quadpt ( ix, :, 18)
     1 + s_quadpt ( ix, :, 19) - 6.0 * s ( ix, :, 01) ) / 3.0
     1 - s_quadpt ( ix, :, 20) + s_quadpt ( ix, :, 01)

      s ( ix, :, 12) = s_quadpt ( ix, :, 20) - s_quadpt ( ix, :, 01)
     1 - s ( ix, :, 11)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ts ( ix, :, 07) = s_quadpt ( ix, :, 23) + s_quadpt ( ix, :, 24)
     1 - s_quadpt ( ix, :, 25) - s_quadpt ( ix, :, 26)
     1 - s_quadpt ( ix, :, 27) - s_quadpt ( ix, :, 28)
     1 + s_quadpt ( ix, :, 29) + s_quadpt ( ix, :, 30)

      ts ( ix, :, 08) = s_quadpt ( ix, :, 32) + s_quadpt ( ix, :, 34)
     1 - s_quadpt ( ix, :, 24) - s_quadpt ( ix, :, 28)
     1 - s_quadpt ( ix, :, 31) - s_quadpt ( ix, :, 33)
     1 + s_quadpt ( ix, :, 23) + s_quadpt ( ix, :, 27)

      ts ( ix, :, 09) = s_quadpt ( ix, :, 25) + s_quadpt ( ix, :, 29)
     1 - s_quadpt ( ix, :, 31) - s_quadpt ( ix, :, 32)
     1 - s_quadpt ( ix, :, 26) - s_quadpt ( ix, :, 30)
     1 + s_quadpt ( ix, :, 33) + s_quadpt ( ix, :, 34)

      s ( ix, :, 16) = ts ( ix, :, 07) + 0.5 * ( ts ( ix, :, 01)
     1 - ts ( ix, :, 02) ) - 2.0 * ( s_quadpt ( ix, :, 02)
     1 - s_quadpt ( ix, :, 05) )

      s ( ix, :, 18) = - 2.0 * ts ( ix, :, 07) + s ( ix, :, 16)

      s ( ix, :, 14) = ts ( ix, :, 08) + 0.5 * ( ts ( ix, :, 03)
     1 - ts ( ix, :, 04) ) - 2.0 * ( s_quadpt ( ix, :, 06)
     1 - s_quadpt ( ix, :, 09) )

      s ( ix, :, 19) = - 2.0 * ts ( ix, :, 08) + s ( ix, :, 14)

      s ( ix, :, 15) = ts ( ix, :, 09) + 0.5 * ( ts ( ix, :, 05)
     1 - ts ( ix, :, 06) ) - 2.0 * ( s_quadpt ( ix, :, 10)
     1 - s_quadpt ( ix, :, 13) )

      s ( ix, :, 17) = - 2.0 * ts ( ix, :, 09) + s ( ix, :, 15)


      ts ( ix, :, 10) = s_quadpt ( ix, :, 23) + s_quadpt ( ix, :, 24)
     1 + s_quadpt ( ix, :, 25) + s_quadpt ( ix, :, 26)
     1 - s_quadpt ( ix, :, 27) - s_quadpt ( ix, :, 28)
     1 - s_quadpt ( ix, :, 29) - s_quadpt ( ix, :, 30)

      ts ( ix, :, 11) = s_quadpt ( ix, :, 23) + s_quadpt ( ix, :, 27)
     1 + s_quadpt ( ix, :, 31) + s_quadpt ( ix, :, 33)
     1 - s_quadpt ( ix, :, 24) - s_quadpt ( ix, :, 28)
     1 - s_quadpt ( ix, :, 32) - s_quadpt ( ix, :, 34)

      ts ( ix, :, 12) = s_quadpt ( ix, :, 25) + s_quadpt ( ix, :, 29)
     1 + s_quadpt ( ix, :, 31) + s_quadpt ( ix, :, 32)
     1 - s_quadpt ( ix, :, 26) - s_quadpt ( ix, :, 30)
     1 - s_quadpt ( ix, :, 33) - s_quadpt ( ix, :, 34)

      s ( ix, :, 27) = 1.125 * ( ts ( ix, :, 01) - ts ( ix, :, 02) )
     1 + 4.5 * ( s_quadpt ( ix, :, 36) - s_quadpt ( ix, :, 37) )
     1 - 2.25 * ts ( ix, :, 10)

      s ( ix, :, 24) = 3.0 * ( s_quadpt ( ix, :, 36)
     1 - s_quadpt ( ix, :, 37) ) - 0.75 * ts ( ix, :, 10)
     1 + 0.375 * ( s ( ix, :, 16) + s ( ix, :, 18) )
     1 - s ( ix, :, 27)

      s ( ix, :, 28) = 1.125 * ( ts ( ix, :, 03) - ts ( ix, :, 04) )
     1 + 4.5 * ( s_quadpt ( ix, :, 38) - s_quadpt ( ix, :, 39) )
     1 - 2.25 * ts ( ix, :, 11)

      s ( ix, :, 25) = 3.0 * ( s_quadpt ( ix, :, 38)
     1 - s_quadpt ( ix, :, 39) ) - 0.75 * ts ( ix, :, 11)
     1 + 0.375 * ( s ( ix, :, 14) + s ( ix, :, 19) )
     1 - s ( ix, :, 28)

      s ( ix, :, 29) = 1.125 * ( ts ( ix, :, 05) - ts ( ix, :, 06) )
     1 + 4.5 * ( s_quadpt ( ix, :, 40) - s_quadpt ( ix, :, 41) )
     1 - 2.25 * ts ( ix, :, 12)

      s ( ix, :, 26) = 3.0 * ( s_quadpt ( ix, :, 40)
     1 - s_quadpt ( ix, :, 41) ) - 0.75 * ts ( ix, :, 12)
     1 + 0.375 * ( s ( ix, :, 15) + s ( ix, :, 17) )
     1 - s ( ix, :, 29)


      s ( ix, :, 30) = 3.0 * ( s_quadpt ( ix, :, 36)
     1 - 2.0 * s_quadpt ( ix, :, 35) + s_quadpt ( ix, :, 37) )
     1 - 1.5 * s ( ix, :, 05)

      s ( ix, :, 31) = 3.0 * ( s_quadpt ( ix, :, 38)
     1 - 2.0 * s_quadpt ( ix, :, 35) + s_quadpt ( ix, :, 39) )
     1 - 1.5 * s ( ix, :, 06)

      s ( ix, :, 32) = 3.0 * ( s_quadpt ( ix, :, 40)
     1 - 2.0 * s_quadpt ( ix, :, 35) + s_quadpt ( ix, :, 41) )
     1 - 1.5 * s ( ix, :, 07)

      s ( ix, :, 33) = 3.0 * ( s_quadpt ( ix, :, 23)
     1 - s_quadpt ( ix, :, 27) - s_quadpt ( ix, :, 24)
     1 + s_quadpt ( ix, :, 28) ) - 3.0 * s ( ix, :, 08)

      s ( ix, :, 34) = 3.0 * ( s_quadpt ( ix, :, 31)
     1 - s_quadpt ( ix, :, 33) - s_quadpt ( ix, :, 32)
     1 + s_quadpt ( ix, :, 34) ) - 3.0 * s ( ix, :, 09)

      s ( ix, :, 35) = 3.0 * ( s_quadpt ( ix, :, 25)
     1 - s_quadpt ( ix, :, 29) - s_quadpt ( ix, :, 26)
     1 + s_quadpt ( ix, :, 30) ) - 3.0 * s ( ix, :, 10)


      s ( ix, :, 23) = 4.5 * ( s_quadpt ( ix, :, 42)
     1 - s_quadpt ( ix, :, 01) ) - 13.5 * ( s_quadpt ( ix, :, 35)
     1 - s_quadpt ( ix, :, 22) )

      s ( ix, :, 22) = 2.25 * ( s_quadpt ( ix, :, 42)
     1 + s_quadpt ( ix, :, 01) - s_quadpt ( ix, :, 35)
     1 - s_quadpt ( ix, :, 22) ) - 1.5 * s ( ix, :, 23)

      s ( ix, :, 21) = s_quadpt ( ix, :, 42) - s_quadpt ( ix, :, 01)
     1 - s ( ix, :, 22) - s ( ix, :, 23)
     1 + ( s ( ix, :, 30) + s ( ix, :, 31)
     1 + s ( ix, :, 32) ) / 12.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
#endif
! End of "BUILD_SOURCE_TERMS == 1"
! ----------------------------------------------------------------------

! STEP 10) By this point, the modal flux coefficients ( f, g, h) and the modal
! source coefficients "s" that enable us to define the fluxes and source terms
! any where in space and time within the zone have already been loaded.
! Now we can use the analytical inversion of the matrix to obtain
! the next iterate for the conserved variables "u" . The coefficients "u"
! enable us to construct the conserved variables at any space-time point in
! the zone.

! Since "u" is made by the end of this step, we also build and return "s_int",
! the space-time integral of the source term.

! ----------------------------------------------------------------------

#if ( NUMORDER == 2)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u ( ix, :, 5) = - f ( ix, :, 2) - g ( ix, :, 3) - h ( ix, :, 4)
     1 + s ( ix, :, 1) + 2.0 * s ( ix, :, 5) / 3.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2".

! ----------------------------------------------------------------------

#if ( NUMORDER == 3)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u ( ix, :, 11) = - f ( ix, :, 2) - g ( ix, :, 3) - h ( ix, :, 4)
     1 + s ( ix, :, 1) - 0.3 * s ( ix, :, 12)

      u ( ix, :, 12) = - 0.5 * ( f ( ix, :, 13) + g ( ix, :, 14)
     1 + h ( ix, :, 15) - s ( ix, :, 11) ) + 0.6 * s ( ix, :, 12)

      u ( ix, :, 13) = - 2.0 * f ( ix, :, 5) - g ( ix, :, 8)
     1 - h ( ix, :, 10) + s ( ix, :, 2) + 2.0 * s ( ix, :, 13) / 3.0

      u ( ix, :, 14) = - f ( ix, :, 8) - 2.0 * g ( ix, :, 6)
     1 - h ( ix, :, 9) + s ( ix, :, 3) + 2.0 * s ( ix, :, 14) / 3.0

      u ( ix, :, 15) = - f ( ix, :, 10) - g ( ix, :, 9)
     1 - 2.0 * h ( ix, :, 7) + s ( ix, :, 4)
     1 + 2.0 * s ( ix, :, 15) / 3.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3".

! ----------------------------------------------------------------------

#if ( NUMORDER == 4)

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u ( ix, :, 21) =
     1 - 0.1 * ( f ( ix, :, 11) + g ( ix, :, 12) + h ( ix, :, 13) )
     1 - f ( ix, :, 2) - g ( ix, :, 3) - h ( ix, :, 4)
     1 + s ( ix, :, 1) + 8.0 * s ( ix, :, 23) / 70.0

      u ( ix, :, 22) =
     1 - 0.5 * ( f ( ix, :, 24) + g ( ix, :, 25) + h ( ix, :, 26)
     1 - s ( ix, :, 21) ) - 3.0 * s ( ix, :, 23) / 7.0

      u ( ix, :, 23) =
     1 - ( f ( ix, :, 27) + g ( ix, :, 28) + h ( ix, :, 29)
     1 - s ( ix, :, 22) ) / 3.0 + 4.0 * s ( ix, :, 23) / 7.0

      u ( ix, :, 24) =
     1 - 2.0 * f ( ix, :, 5) - g ( ix, :, 8) - h ( ix, :, 10)
     1 + s ( ix, :, 2) - 0.3 * s ( ix, :, 27)

      u ( ix, :, 25) =
     1 - f ( ix, :, 8) - 2.0 * g ( ix, :, 6) - h ( ix, :, 9)
     1 + s ( ix, :, 3) - 0.3 * s ( ix, :, 28)

      u ( ix, :, 26) =
     1 - f ( ix, :, 10) - g ( ix, :, 9) - 2.0 * h ( ix, :, 7)
     1 + s ( ix, :, 4) - 0.3 * s ( ix, :, 29)

      u ( ix, :, 27) =
     1 - 0.5 * ( g ( ix, :, 33) + h ( ix, :, 35) - s ( ix, :, 24) )
     1 - f ( ix, :, 30) + 0.6 * s ( ix, :, 27)

      u ( ix, :, 28) =
     1 - 0.5 * ( f ( ix, :, 33) + h ( ix, :, 34) - s ( ix, :, 25) )
     1 - g ( ix, :, 31) + 0.6 * s ( ix, :, 28)

      u ( ix, :, 29) =
     1 - 0.5 * ( f ( ix, :, 35) + g ( ix, :, 34) - s ( ix, :, 26) )
     1 - h ( ix, :, 32) + 0.6 * s ( ix, :, 29)

      u ( ix, :, 30) =
     1 - 3.0 * f ( ix, :, 11) - g ( ix, :, 14) - h ( ix, :, 15)
     1 + s ( ix, :, 5) + 2.0 * s ( ix, :, 30) / 3.0

      u ( ix, :, 31) =
     1 - f ( ix, :, 16) - 3.0 * g ( ix, :, 12) - h ( ix, :, 17)
     1 + s ( ix, :, 6) + 2.0 * s ( ix, :, 31) / 3.0

      u ( ix, :, 32) =
     1 - f ( ix, :, 18) - g ( ix, :, 19) - 3.0 * h ( ix, :, 13)
     1 + s ( ix, :, 7) + 2.0 * s ( ix, :, 32) / 3.0

      u ( ix, :, 33) =
     1 - 2.0 * f ( ix, :, 14) - 2.0 * g ( ix, :, 16) - h ( ix, :, 20)
     1 + s ( ix, :, 8) + 2.0 * s ( ix, :, 33) / 3.0

      u ( ix, :, 34) =
     1 - f ( ix, :, 20) - 2.0 * g ( ix, :, 17) - 2.0 * h ( ix, :, 19)
     1 + s ( ix, :, 9) + 2.0 * s ( ix, :, 34) / 3.0

      u ( ix, :, 35) =
     1 - 2.0 * f ( ix, :, 15) - g ( ix, :, 20) - 2.0 * h ( ix, :, 18)
     1 + s ( ix, :, 10) + 2.0 * s ( ix, :, 35) / 3.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO
! End of "DO iter_ader =" loop.

! End the outer ADER iteration here.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Beginning of building all the nodal conserved and primitive variables at the
! t > 0 nodes. This is needed either if we want to rebuild the modal
! fluxes "f,g,h" or if we want to rebuild the modal source term "s".

#if ( BUILD_SOURCE_TERMS == 1)

! ----------------------------------------------------------------------

! STEP 6_Polish-off) Make the conserved variables at the t>0 space-time nodes.
! We put in the flexibility of having different "gamma_eos_spacetime" at each
! space-time node, but that is not used.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 07) = u ( ix, :, 01) + u ( ix, :, 05)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 14) = 0.25 * ( u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 03) + u_quadpt ( ix, :, 04)
     1 + u_quadpt ( ix, :, 05) + u ( ix, :, 12) + u ( ix, :, 13) )
     1 + 0.5 * u ( ix, :, 11) - 0.125 * ( u ( ix, :, 06)
     1 + u ( ix, :, 07) )

      u_quadpt ( ix, :, 15) = u_quadpt ( ix, :, 14)
     1 - u ( ix, :, 02) - 0.5 * u ( ix, :, 13)

      u_quadpt ( ix, :, 16) = 0.25 * ( u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 06) + u_quadpt ( ix, :, 10)
     1 + u_quadpt ( ix, :, 12) + u ( ix, :, 12) + u ( ix, :, 14) )
     1 + 0.5 * u ( ix, :, 11) - 0.125 * ( u ( ix, :, 05)
     1 + u ( ix, :, 07) )

      u_quadpt ( ix, :, 17) = u_quadpt ( ix, :, 16)
     1 - u ( ix, :, 03) - 0.5 * u ( ix, :, 14)

      u_quadpt ( ix, :, 18) = 0.25 * ( u_quadpt ( ix, :, 04)
     1 + u_quadpt ( ix, :, 08) + u_quadpt ( ix, :, 10)
     1 + u_quadpt ( ix, :, 11) + u ( ix, :, 12) + u ( ix, :, 15) )
     1 + 0.5 * u ( ix, :, 11) - 0.125 * ( u ( ix, :, 05)
     1 + u ( ix, :, 06) )

      u_quadpt ( ix, :, 19) = u_quadpt ( ix, :, 18) - u ( ix, :, 04)
     1 - 0.5 * u ( ix, :, 15)

      u_quadpt ( ix, :, 20) = u_quadpt ( ix, :, 01) + u ( ix, :, 11)
     1 + u ( ix, :, 12)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 22) = u_quadpt ( ix, :, 01)
     1 + ( 36.0 * u ( ix, :, 21) + 12.0 * u ( ix, :, 22)
     1 + 4.0 * u ( ix, :, 23) - 3.0 * u ( ix, :, 30)
     1 - 3.0 * u ( ix, :, 31) - 3.0 * u ( ix, :, 32) ) / 108.0

      u_quadpt ( ix, :, 23) = u_quadpt ( ix, :, 22)
     1 + u_quadpt ( ix, :, 02) + u_quadpt ( ix, :, 06)
     1 - 2.0 * u_quadpt ( ix, :, 01) + ( 18.0 * u ( ix, :, 08)
     1 + 9.0 * u ( ix, :, 14) + 9.0 * u ( ix, :, 16)
     1 + 12.0 * u ( ix, :, 24) + 12.0 * u ( ix, :, 25)
     1 + 4.0 * u ( ix, :, 27) + 4.0 * u ( ix, :, 28)
     1 + 6.0 * u ( ix, :, 30) + 6.0 * u ( ix, :, 31)
     1 + 6.0 * u ( ix, :, 33) ) / 72.0

      u_quadpt ( ix, :, 24) = u_quadpt ( ix, :, 23)
     1 - u_quadpt ( ix, :, 06) + u_quadpt ( ix, :, 09)
     1 - 0.5 * u ( ix, :, 08) - 0.25 * u ( ix, :, 14)
     1 - ( 6.0 * u ( ix, :, 25) + 2.0 * u ( ix, :, 28)
     1 + 3.0 * u ( ix, :, 33) ) / 18.0

      u_quadpt ( ix, :, 25) = 0.5 * ( u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 24) + u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 15) ) - u_quadpt ( ix, :, 02)
     1 - 0.5 * u ( ix, :, 06) - 0.25 * u ( ix, :, 16)
     1 - 0.125 * u ( ix, :, 17) + ( 6.0 * u ( ix, :, 26)
     1 + 2.0 * u ( ix, :, 29) - 3.0 * u ( ix, :, 31)
     1 + 3.0 * u ( ix, :, 32) + 3.0 * u ( ix, :, 35) ) / 36.0

      u_quadpt ( ix, :, 26) = u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 24) - u_quadpt ( ix, :, 25)
     1 - 0.5 * ( u ( ix, :, 06) - u ( ix, :, 07) )
     1 - 0.25 * ( u ( ix, :, 16) - u ( ix, :, 18) )
     1 - ( u ( ix, :, 31) - u ( ix, :, 32) ) / 6.0

      u_quadpt ( ix, :, 27) = u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 02)
     1 - 0.5 * u ( ix, :, 08) - 0.25 * u ( ix, :, 16)
     1 - ( 6.0 * u ( ix, :, 24) + 2.0 * u ( ix, :, 27)
     1 + 3.0 * u ( ix, :, 33) ) / 18.0

      u_quadpt ( ix, :, 28) = u_quadpt ( ix, :, 27)
     1 - u_quadpt ( ix, :, 06) + u_quadpt ( ix, :, 09)
     1 + 0.5 * u ( ix, :, 08) - 0.25 * u ( ix, :, 14)
     1 - ( 6.0 * u ( ix, :, 25) + 2.0 * u ( ix, :, 28)
     1 - 3.0 * u ( ix, :, 33) ) / 18.0

      u_quadpt ( ix, :, 29) = 0.5 * ( u_quadpt ( ix, :, 27)
     1 + u_quadpt ( ix, :, 28) + u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 19) ) - u_quadpt ( ix, :, 05)
     1 - 0.5 * u ( ix, :, 06) + 0.25 * u ( ix, :, 16)
     1 - 0.125 * u ( ix, :, 17) + ( 6.0 * u ( ix, :, 26)
     1 + 2.0 * u ( ix, :, 29) - 3.0 * u ( ix, :, 31)
     1 + 3.0 * u ( ix, :, 32)  - 3.0 * u ( ix, :, 35) ) / 36.0

      u_quadpt ( ix, :, 30) = u_quadpt ( ix, :, 27)
     1 + u_quadpt ( ix, :, 28) - u_quadpt ( ix, :, 29)
     1 - 0.5 * ( u ( ix, :, 06) - u ( ix, :, 07) )
     1 + 0.25 * ( u ( ix, :, 16) - u ( ix, :, 18) )
     1 - ( u ( ix, :, 31) - u ( ix, :, 32) ) / 6.0

      u_quadpt ( ix, :, 31) = 0.5 * ( u_quadpt ( ix, :, 25)
     1 + u_quadpt ( ix, :, 29) + u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 18) ) - u_quadpt ( ix, :, 10)
     1 - 0.5 * u ( ix, :, 05) - 0.25 * u ( ix, :, 15)
     1 - 0.125 * u ( ix, :, 14) + ( 6.0 * u ( ix, :, 25)
     1 + 2.0 * u ( ix, :, 28) - 3.0 * u ( ix, :, 30)
     1 + 3.0 * u ( ix, :, 31) + 3.0 * u ( ix, :, 34) ) / 36.0

      u_quadpt ( ix, :, 32) = u_quadpt ( ix, :, 31)
     1 - u_quadpt ( ix, :, 06) + u_quadpt ( ix, :, 09)
     1 - 0.5 * u ( ix, :, 09) - 0.25 * u ( ix, :, 19)
     1 - ( 6.0 * u ( ix, :, 25) + 2.0 * u ( ix, :, 28)
     1 + 3.0 * u ( ix, :, 34) ) / 18.0

      u_quadpt ( ix, :, 33) = u_quadpt ( ix, :, 23)
     1 + u_quadpt ( ix, :, 27) - u_quadpt ( ix, :, 31)
     1 - 0.5 * ( u ( ix, :, 05) - u ( ix, :, 07) )
     1 - 0.25 * ( u ( ix, :, 14) - u ( ix, :, 19) )
     1 - ( u ( ix, :, 30) - u ( ix, :, 32) ) / 6.0

      u_quadpt ( ix, :, 34) = u_quadpt ( ix, :, 24)
     1 + u_quadpt ( ix, :, 28) - u_quadpt ( ix, :, 32)
     1 - 0.5 * ( u ( ix, :, 05) - u ( ix, :, 07) )
     1 + 0.25 * ( u ( ix, :, 14) - u ( ix, :, 19) )
     1 - ( u ( ix, :, 30) - u ( ix, :, 32) ) / 6.0


      tu ( ix, :, 01) = ( 36.0 * u ( ix, :, 21) + 24.0 * u ( ix, :, 22)
     1 + 16.0 * u ( ix, :, 23) ) / 54.0

      tu ( ix, :, 02) = ( 18.0 * u ( ix, :, 24)
     1 + 12.0 * u ( ix, :, 27) ) / 54.0

      tu ( ix, :, 03) = ( 6.0 * u ( ix, :, 30) - 3.0 * u ( ix, :, 31)
     1 - 3.0 * u ( ix, :, 32) ) / 54.0

      tu ( ix, :, 04) = ( 18.0 * u ( ix, :, 25)
     1 + 12.0 * u ( ix, :, 28) ) / 54.0

      tu ( ix, :, 05) = ( - 3.0 * u ( ix, :, 30) + 6.0 * u ( ix, :, 31)
     1 - 3.0 * u ( ix, :, 32) ) / 54.0

      tu ( ix, :, 06) = ( 18.0 * u ( ix, :, 26)
     1 + 12.0 * u ( ix, :, 29) ) / 54.0

      tu ( ix, :, 07) = ( - 3.0 * u ( ix, :, 30) - 3.0 * u ( ix, :, 31)
     1 + 6.0 * u ( ix, :, 32) ) / 54.0

      u_quadpt ( ix, :, 35) = ( 2.0 * u_quadpt ( ix, :, 22)
     1 - u_quadpt ( ix, :, 01) ) + 2.0 * ( u ( ix, :, 22)
     1 + u ( ix, :, 23) ) / 9.0

      u_quadpt ( ix, :, 36) = u_quadpt ( ix, :, 02) + tu ( ix, :, 01)
     1 + tu ( ix, :, 02) + tu ( ix, :, 03)

      u_quadpt ( ix, :, 37) = u_quadpt ( ix, :, 05) + tu ( ix, :, 01)
     1 - tu ( ix, :, 02) + tu ( ix, :, 03)

      u_quadpt ( ix, :, 38) = u_quadpt ( ix, :, 06) + tu ( ix, :, 01)
     1 + tu ( ix, :, 04) + tu ( ix, :, 05)

      u_quadpt ( ix, :, 39) = u_quadpt ( ix, :, 09) + tu ( ix, :, 01)
     1 - tu ( ix, :, 04) + tu ( ix, :, 05)

      u_quadpt ( ix, :, 40) = u_quadpt ( ix, :, 10) + tu ( ix, :, 01)
     1 + tu ( ix, :, 06) + tu ( ix, :, 07)

      u_quadpt ( ix, :, 41) = u_quadpt ( ix, :, 13) + tu ( ix, :, 01)
     1 - tu ( ix, :, 06) + tu ( ix, :, 07)


      u_quadpt ( ix, :, 42) = u_quadpt ( ix, :, 01) + u ( ix, :, 21)
     1 + u ( ix, :, 22) + u ( ix, :, 23) - ( u ( ix, :, 30)
     1 + u ( ix, :, 31) + u ( ix, :, 32) ) / 12.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 7_Polish-off) Get the primitive variables from the conserved variables
! at the spatial nodes at t>0. The density and pressure variables are
! given positive floor values to ensure that they can be used for physically
! meaningful source terms.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, 1, iquadpt) =
     1                   AMAX1 ( u_quadpt ( ix, 1, iquadpt), rhofloor)

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempaa = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempaa
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempaa
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempaa

      bfx ( ix, iquadpt) = u_quadpt ( ix, 6, iquadpt)
      bfy ( ix, iquadpt) = u_quadpt ( ix, 7, iquadpt)
      bfz ( ix, iquadpt) = u_quadpt ( ix, 8, iquadpt)

      tempab = 0.5 * rho ( ix, iquadpt) * ( vlx ( ix, iquadpt)**2
     1         + vly ( ix, iquadpt)**2 + vlz ( ix, iquadpt)**2)

      tempac = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      eng_t ( ix, iquadpt) = u_quadpt ( ix, 5, iquadpt)
     1                     - tempab - tempac

      prs ( ix, iquadpt) = eng_t ( ix, iquadpt)
     1                   * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

#if ( ENFORCE_POSITIVITY == 1)
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = entropy ( ix, iquadpt)
     1* rho ( ix, iquadpt)**( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)
      prs ( ix, iquadpt) = AMAX1 ( prs ( ix, iquadpt), prsfloor)

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#else
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = prsfloor

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#endif

#if ( NFLUID >= 1)
      specie ( ix, 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( ix, 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( ix, 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( ix, 1: NFLUID, iquadpt) * rho ( ix, iquadpt)
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, ivar + 1, iquadpt) =
     1         AMAX1 ( u_quadpt ( ix, ivar + 1, iquadpt), rhofloor)

      rho_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 1, iquadpt)

      tempaa = 1.0 / rho_euler ( ix, iquadpt, i_euler)

      vlx_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 2, iquadpt) * tempaa
      vly_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 3, iquadpt) * tempaa
      vlz_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 4, iquadpt) * tempaa

      tempab = 0.5 * rho_euler ( ix, iquadpt, i_euler)
     1           * ( vlx_euler ( ix, iquadpt, i_euler)**2
     1             + vly_euler ( ix, iquadpt, i_euler)**2
     1             + vlz_euler ( ix, iquadpt, i_euler)**2)

      prs_euler ( ix, iquadpt, i_euler) =
     1         ( u_quadpt ( ix, ivar + 5, iquadpt) - tempab)
     1       * ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      prs_euler ( ix, iquadpt, i_euler) =
     1    AMAX1 ( prs_euler ( ix, iquadpt, i_euler), prsfloor)

#if ( ISOTHERMAL_EQNS == 1)
      prs_euler ( ix, iquadpt, i_euler) =
     1     isospeed**2 * rho_euler ( ix, iquadpt, i_euler)
     1   + barotropic_stiffness
     1   * rho_euler ( ix, iquadpt, i_euler)**gamma
#endif

      eng_t_euler ( ix, iquadpt, i_euler) =
     1         prs_euler ( ix, iquadpt, i_euler)
     1       / ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      IF ( prs_euler ( ix, iquadpt, i_euler) .LT. prsfloor * 1.00000001)
     1   u_quadpt ( ix, ivar + 5, iquadpt) =
     1                   eng_t_euler ( ix, iquadpt, i_euler) + tempab

#if ( NFLUID_EULER >= 1)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) = tempaa
     1       * u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1                    ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) =
     1    AMIN1 ( AMAX1 (
     1    specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler),
     1            0.0), 1.0 )

      u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1       specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler)
     1     * rho_euler ( ix, iquadpt, i_euler)

#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ----------------------------------------------------------------------

#endif
! End of "BUILD_SOURCE_TERMS == 1"

! End of building all the nodal conserved and primitive variables at the
! t > 0 nodes. This is needed either if we want to rebuild the modal
! fluxes "f,g,h" or if we want to rebuild the modal source term "s".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Default is to build source terms at t>0 . In some ADER-CG updates
! (on Cartesian meshes), we may not require any source terms.

#if ( BUILD_SOURCE_TERMS == 1)
! ----------------------------------------------------------------------

! STEP 9.a_Polish-off) Build source terms "s_quadpt" at the t>0 quadrature
! points.
! Remember that the primitive variables are still available from the previous
! loop at each of the t=0 quadrature points for constructing the source terms.

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      s_quadpt ( oned_begin - oned_off: oned_end + oned_off,
     1                                           :, iquadpt) = 0.0

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

#if ( 0 == 1)

! For single-fluid adiabatic MHD turbulence, we may want to provide a gradual
! (explicit) heating/cooling term which ensures that pressure always
! tends to a mean pressure "prs_bar" in a time "tau_bar" (comparable
! to a turnover time). Such a source term can be implemented explicitly
! and its implementation is given below.
! Remember to declare "prs_bar, tau_bar".

      prs_bar = 1.0
      tau_bar = 1.0

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 5, iquadpt) = - ( prs ( ix, iquadpt) - prs_bar)
     1     / ( tau_bar * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0) )

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

#endif

! ------------------------------

! Build two-fluid source terms. In this part, the source terms are
! hard-coded for one family of ionized fluid and one family of neutrals.

#if ( N_EULER == 1) 

      IF ( ( igeom .EQ. 1) .AND. ( N_EULER .EQ. 1) .AND.
     1     ( ambipolar_diffusion_drag_coef .GT. 0.0)) THEN

      ivar = NFIELD + 1 + NFLUID

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 1, iquadpt) = 0.0

      tempaa = ambipolar_diffusion_drag_coef
     1       * rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

      s_quadpt ( ix, 2, iquadpt) = tempaa
     1           * ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt))

      s_quadpt ( ix, 3, iquadpt) = tempaa
     1           * ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt))

      s_quadpt ( ix, 4, iquadpt) = tempaa
     1           * ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt))

      s_quadpt ( ix, 5, iquadpt) = ( ambipolar_diffusion_drag_coef
     1   / ( molewt_ion + molewt_neut) )

     1 * ( rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

     1   * ( ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt) )
     1     * ( molewt_neut * vlx_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlx ( ix, iquadpt) )
     1     + ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt) )
     1     * ( molewt_neut * vly_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vly ( ix, iquadpt) )
     1     + ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt) )
     1     * ( molewt_neut * vlz_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlz ( ix, iquadpt) )
     1     )

     1   + 3.0 * ( molewt_neut * prs_euler ( ix, iquadpt, 1)
     1           * rho ( ix, iquadpt)
     1           - molewt_ion * prs ( ix, iquadpt)
     1           * rho_euler ( ix, iquadpt, 1) )
     1   )

      s_quadpt ( ix, 6, iquadpt) = 0.0
      s_quadpt ( ix, 7, iquadpt) = 0.0
      s_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      s_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) = 0.0
#endif


      s_quadpt ( ix, ivar + 1, iquadpt) = 0.0
      s_quadpt ( ix, ivar + 2, iquadpt) = - s_quadpt ( ix, 2, iquadpt)
      s_quadpt ( ix, ivar + 3, iquadpt) = - s_quadpt ( ix, 3, iquadpt)
      s_quadpt ( ix, ivar + 4, iquadpt) = - s_quadpt ( ix, 4, iquadpt)
      s_quadpt ( ix, ivar + 5, iquadpt) = - s_quadpt ( ix, 5, iquadpt)

#if ( NFLUID_EULER >= 1)
      s_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

      END IF
! End of "ambipolar_diffusion_drag_coef .GT. 0.0" part

#endif
! End of "N_EULER == 1"

! ------------------------------

! Cylindrical Geometry

      IF ( igeom .EQ. 2) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1                           * vlx ( ix, iquadpt) * obx1

      s_quadpt ( ix, 2, iquadpt) =
     1 ( rho ( ix, iquadpt) * vly ( ix, iquadpt)**2
     1 - rho ( ix, iquadpt) * vlx ( ix, iquadpt)**2
     1 + bfx ( ix, iquadpt)**2 * ob4pi
     1 - bfy ( ix, iquadpt)**2 * ob4pi) * obx1

      s_quadpt ( ix, 3, iquadpt) =
     1 - 2.0 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                              * vly ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                        * vlz ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 5, iquadpt) = - ( tempac * vlx ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = 0.0

      s_quadpt ( ix, 7, iquadpt) = 0.0

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Spherical Geometry

      IF ( igeom .EQ. 3) THEN

      DO iquadpt = NUM_ZONECEN_QUADPTS + 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      cos_x2 = cos ( x2)
      sin_x2 = sin ( x2)
      cot_x2 = cos_x2 / sin_x2

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 2.0 * vlx ( ix, iquadpt)
     1 + cot_x2 * vly ( ix, iquadpt) ) * obx1

      s_quadpt ( ix, 2, iquadpt) = rho ( ix, iquadpt)
     1 * ( - 2.0 * vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1   + vlz ( ix, iquadpt)**2) * obx1
     1 - cot_x2
     1 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1   - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1
     1 - ( - 2.0 * bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1   + bfz ( ix, iquadpt)**2) * ob4pi * obx1

      s_quadpt ( ix, 3, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1 + cot_x2
     1 * ( vly ( ix, iquadpt)**2 - vlz ( ix, iquadpt)**2)) * obx1
     1 + ( 3.0 * bfx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1   + cot_x2
     1   * ( bfy ( ix, iquadpt)**2 - bfz ( ix, iquadpt)**2))
     1 * ob4pi * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - rho ( ix, iquadpt) * vlz ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) + 2.0 * cot_x2
     1   * vly ( ix, iquadpt) ) * obx1
     1 + bfz ( ix, iquadpt)
     1 * ( 3.0 * bfx ( ix, iquadpt) + 2.0 * cot_x2
     1   * bfy ( ix, iquadpt) ) * ob4pi * obx1

      s_quadpt ( ix, 5, iquadpt) =
     1 - 2.0 * ( tempac * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab) * obx1
     1 - cot_x2 * ( tempac * vly ( ix, iquadpt)
     1 - bfy ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = - cot_x2
     1                           * g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 7, iquadpt) = g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Multiply "s_quadpt" by "dtmove".

      DO ix = oned_begin - oned_off, oned_end + oned_off
      s_quadpt ( ix, :, NUM_ZONECEN_QUADPTS + 1: NUM_SPACETIME_QUADPTS)=
     1 dtmove
     1*s_quadpt ( ix, :, NUM_ZONECEN_QUADPTS + 1: NUM_SPACETIME_QUADPTS)
      END DO

! ------------------------------

! STEP 9.b_Polish-off) Now project the sources "s_quadpt" from the t>0
! space-time nodal points back to the modal space time variables for
! those variables, i.e. to "s".

! Note that some of the "s_quadpt" values with t>0 will be overwritten.
! But that is ok since these are also constructed anew in each ADER iteration.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s ( ix, :, 05) = s_quadpt ( ix, :, 07) - s ( ix, :, 01)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s ( ix, :, 13) = 2.0 * ( s_quadpt ( ix, :, 14)
     1 - s_quadpt ( ix, :, 15) - s ( ix, :, 02) )

      s ( ix, :, 14) = 2.0 * ( s_quadpt ( ix, :, 16)
     1 - s_quadpt ( ix, :, 17) - s ( ix, :, 03) )

      s ( ix, :, 15) = 2.0 * ( s_quadpt ( ix, :, 18)
     1 - s_quadpt ( ix, :, 19) - s ( ix, :, 04) )

      s ( ix, :, 11) = 2.0 * ( s_quadpt ( ix, :, 14)
     1 + s_quadpt ( ix, :, 15) + s_quadpt ( ix, :, 16)
     1 + s_quadpt ( ix, :, 17) + s_quadpt ( ix, :, 18)
     1 + s_quadpt ( ix, :, 19) - 6.0 * s ( ix, :, 01) ) / 3.0
     1 - s_quadpt ( ix, :, 20) + s_quadpt ( ix, :, 01)

      s ( ix, :, 12) = s_quadpt ( ix, :, 20) - s_quadpt ( ix, :, 01)
     1 - s ( ix, :, 11)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ts ( ix, :, 07) = s_quadpt ( ix, :, 23) + s_quadpt ( ix, :, 24)
     1 - s_quadpt ( ix, :, 25) - s_quadpt ( ix, :, 26)
     1 - s_quadpt ( ix, :, 27) - s_quadpt ( ix, :, 28)
     1 + s_quadpt ( ix, :, 29) + s_quadpt ( ix, :, 30)

      ts ( ix, :, 08) = s_quadpt ( ix, :, 32) + s_quadpt ( ix, :, 34)
     1 - s_quadpt ( ix, :, 24) - s_quadpt ( ix, :, 28)
     1 - s_quadpt ( ix, :, 31) - s_quadpt ( ix, :, 33)
     1 + s_quadpt ( ix, :, 23) + s_quadpt ( ix, :, 27)

      ts ( ix, :, 09) = s_quadpt ( ix, :, 25) + s_quadpt ( ix, :, 29)
     1 - s_quadpt ( ix, :, 31) - s_quadpt ( ix, :, 32)
     1 - s_quadpt ( ix, :, 26) - s_quadpt ( ix, :, 30)
     1 + s_quadpt ( ix, :, 33) + s_quadpt ( ix, :, 34)

      s ( ix, :, 16) = ts ( ix, :, 07) + 0.5 * ( ts ( ix, :, 01)
     1 - ts ( ix, :, 02) ) - 2.0 * ( s_quadpt ( ix, :, 02)
     1 - s_quadpt ( ix, :, 05) )

      s ( ix, :, 18) = - 2.0 * ts ( ix, :, 07) + s ( ix, :, 16)

      s ( ix, :, 14) = ts ( ix, :, 08) + 0.5 * ( ts ( ix, :, 03)
     1 - ts ( ix, :, 04) ) - 2.0 * ( s_quadpt ( ix, :, 06)
     1 - s_quadpt ( ix, :, 09) )

      s ( ix, :, 19) = - 2.0 * ts ( ix, :, 08) + s ( ix, :, 14)

      s ( ix, :, 15) = ts ( ix, :, 09) + 0.5 * ( ts ( ix, :, 05)
     1 - ts ( ix, :, 06) ) - 2.0 * ( s_quadpt ( ix, :, 10)
     1 - s_quadpt ( ix, :, 13) )

      s ( ix, :, 17) = - 2.0 * ts ( ix, :, 09) + s ( ix, :, 15)


      ts ( ix, :, 10) = s_quadpt ( ix, :, 23) + s_quadpt ( ix, :, 24)
     1 + s_quadpt ( ix, :, 25) + s_quadpt ( ix, :, 26)
     1 - s_quadpt ( ix, :, 27) - s_quadpt ( ix, :, 28)
     1 - s_quadpt ( ix, :, 29) - s_quadpt ( ix, :, 30)

      ts ( ix, :, 11) = s_quadpt ( ix, :, 23) + s_quadpt ( ix, :, 27)
     1 + s_quadpt ( ix, :, 31) + s_quadpt ( ix, :, 33)
     1 - s_quadpt ( ix, :, 24) - s_quadpt ( ix, :, 28)
     1 - s_quadpt ( ix, :, 32) - s_quadpt ( ix, :, 34)

      ts ( ix, :, 12) = s_quadpt ( ix, :, 25) + s_quadpt ( ix, :, 29)
     1 + s_quadpt ( ix, :, 31) + s_quadpt ( ix, :, 32)
     1 - s_quadpt ( ix, :, 26) - s_quadpt ( ix, :, 30)
     1 - s_quadpt ( ix, :, 33) - s_quadpt ( ix, :, 34)

      s ( ix, :, 27) = 1.125 * ( ts ( ix, :, 01) - ts ( ix, :, 02) )
     1 + 4.5 * ( s_quadpt ( ix, :, 36) - s_quadpt ( ix, :, 37) )
     1 - 2.25 * ts ( ix, :, 10)

      s ( ix, :, 24) = 3.0 * ( s_quadpt ( ix, :, 36)
     1 - s_quadpt ( ix, :, 37) ) - 0.75 * ts ( ix, :, 10)
     1 + 0.375 * ( s ( ix, :, 16) + s ( ix, :, 18) )
     1 - s ( ix, :, 27)

      s ( ix, :, 28) = 1.125 * ( ts ( ix, :, 03) - ts ( ix, :, 04) )
     1 + 4.5 * ( s_quadpt ( ix, :, 38) - s_quadpt ( ix, :, 39) )
     1 - 2.25 * ts ( ix, :, 11)

      s ( ix, :, 25) = 3.0 * ( s_quadpt ( ix, :, 38)
     1 - s_quadpt ( ix, :, 39) ) - 0.75 * ts ( ix, :, 11)
     1 + 0.375 * ( s ( ix, :, 14) + s ( ix, :, 19) )
     1 - s ( ix, :, 28)

      s ( ix, :, 29) = 1.125 * ( ts ( ix, :, 05) - ts ( ix, :, 06) )
     1 + 4.5 * ( s_quadpt ( ix, :, 40) - s_quadpt ( ix, :, 41) )
     1 - 2.25 * ts ( ix, :, 12)

      s ( ix, :, 26) = 3.0 * ( s_quadpt ( ix, :, 40)
     1 - s_quadpt ( ix, :, 41) ) - 0.75 * ts ( ix, :, 12)
     1 + 0.375 * ( s ( ix, :, 15) + s ( ix, :, 17) )
     1 - s ( ix, :, 29)


      s ( ix, :, 30) = 3.0 * ( s_quadpt ( ix, :, 36)
     1 - 2.0 * s_quadpt ( ix, :, 35) + s_quadpt ( ix, :, 37) )
     1 - 1.5 * s ( ix, :, 05)

      s ( ix, :, 31) = 3.0 * ( s_quadpt ( ix, :, 38)
     1 - 2.0 * s_quadpt ( ix, :, 35) + s_quadpt ( ix, :, 39) )
     1 - 1.5 * s ( ix, :, 06)

      s ( ix, :, 32) = 3.0 * ( s_quadpt ( ix, :, 40)
     1 - 2.0 * s_quadpt ( ix, :, 35) + s_quadpt ( ix, :, 41) )
     1 - 1.5 * s ( ix, :, 07)

      s ( ix, :, 33) = 3.0 * ( s_quadpt ( ix, :, 23)
     1 - s_quadpt ( ix, :, 27) - s_quadpt ( ix, :, 24)
     1 + s_quadpt ( ix, :, 28) ) - 3.0 * s ( ix, :, 08)

      s ( ix, :, 34) = 3.0 * ( s_quadpt ( ix, :, 31)
     1 - s_quadpt ( ix, :, 33) - s_quadpt ( ix, :, 32)
     1 + s_quadpt ( ix, :, 34) ) - 3.0 * s ( ix, :, 09)

      s ( ix, :, 35) = 3.0 * ( s_quadpt ( ix, :, 25)
     1 - s_quadpt ( ix, :, 29) - s_quadpt ( ix, :, 26)
     1 + s_quadpt ( ix, :, 30) ) - 3.0 * s ( ix, :, 10)


      s ( ix, :, 23) = 4.5 * ( s_quadpt ( ix, :, 42)
     1 - s_quadpt ( ix, :, 01) ) - 13.5 * ( s_quadpt ( ix, :, 35)
     1 - s_quadpt ( ix, :, 22) )

      s ( ix, :, 22) = 2.25 * ( s_quadpt ( ix, :, 42)
     1 + s_quadpt ( ix, :, 01) - s_quadpt ( ix, :, 35)
     1 - s_quadpt ( ix, :, 22) ) - 1.5 * s ( ix, :, 23)

      s ( ix, :, 21) = s_quadpt ( ix, :, 42) - s_quadpt ( ix, :, 01)
     1 - s ( ix, :, 22) - s ( ix, :, 23)
     1 + ( s ( ix, :, 30) + s ( ix, :, 31)
     1 + s ( ix, :, 32) ) / 12.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
#endif
! End of "BUILD_SOURCE_TERMS == 1"

! All the modes of "s" are made up to the desired order by the end
! of this step.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 11) Returning the space-time integral of source terms that are treated
! time-explicitly is easy and is done here. If the treatment is time-implicit,
! the source term "s ( ix, :, 1)" should be rebuilt here using
! "u ( ix, :, :)" from the most recent iteration.

! IF "BUILD_SOURCE_TERMS == 0" then "s" should be zero, so "s_int" will be 0.

! In principle, the factor of "dtmove", that was previously multiplied into
! the source terms, is divided out of "s_int" here. However,
! because the "dudt" terms in the subroutine that calls this one already
! carry the source terms, we leave the "dtmove" multiplied into the
! "s_int" terms.

#if ( NUMORDER == 2)
      DO ix = oned_begin - oned_off, oned_end + oned_off
      s_int ( ix, :) = s ( ix, :, 1) + 0.5 * s ( ix, :, 5)
      END DO
#endif

#if ( NUMORDER == 3)
      DO ix = oned_begin - oned_off, oned_end + oned_off
      s_int ( ix, :) = s ( ix, :, 1) + 0.5 * s ( ix, :, 11)
     1 + s ( ix, :, 12) / 3.0
      END DO
#endif

#if ( NUMORDER == 4)
      DO ix = oned_begin - oned_off, oned_end + oned_off
      s_int ( ix, :) = s ( ix, :, 1) + 0.5 * s ( ix, :, 21)
     1 + s ( ix, :, 22) / 3.0 + 0.25 * s ( ix, :, 23)
      END DO
#endif

! ----------------------------------------------------------------------

! STEP 12) Here we give ourselves the option of trapping a zone that has
! not converged. This may consist of zeroing out higher order moments, i.e.
! we fall back to first order here if we absolutely must.

! This is also where we can use Dumbser's trick of zeroing out the
! higher moments in zones that produce negative pressures or densities.

      corrupted_zone ( oned_begin - oned_off: oned_end + oned_off) = 0

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

        IF ( ( rho ( ix, iquadpt) .LT. rhofloor) .OR.
     1       ( prs ( ix, iquadpt) .LT. prsfloor) )
     1                                      corrupted_zone ( ix) = 1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

        IF ( ( rho_euler ( ix, iquadpt, i_euler) .LT. rhofloor) .OR.
     1       ( prs_euler ( ix, iquadpt, i_euler) .LT. prsfloor) )
     1                                      corrupted_zone ( ix) = 1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ------------------------------

! For ADER_CG ( but not for ADER_DG) "f ( ix, :, 1)" is just the
! zone-centered flux at t=0. For ADER_DG, this will have to change.

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( corrupted_zone ( ix) .EQ. 1) THEN

! -----

      u ( ix, :, 1) = w ( ix, :, 1)

      s ( ix, :, 1) = 0.0

      s_int ( ix, :) = 0.0

! -----

      w ( ix, :, 2: NUM_ZONECEN_BASES) = 0.0

      u ( ix, :, 2: NUM_SPACETIME_BASES) = 0.0

      f ( ix, :, 2: NUM_SPACETIME_BASES) = 0.0

      g ( ix, :, 2: NUM_SPACETIME_BASES) = 0.0

      h ( ix, :, 2: NUM_SPACETIME_BASES) = 0.0

      s ( ix, :, 2: NUM_SPACETIME_BASES) = 0.0

! -----

      END IF

! ----------

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ADER_CG_MODAL_3D_TIME_EVOL_MHD

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADER_CG_NODAL_3D_TIME_EVOL_MHD

     1           ( oned_begin, oned_end, oned_off,
     1             dtmove,
     1             igeom, pi, rhofloor, prsfloor, smallnum,
     1             gamma, gamma_euler, gasconst,
     1             ader_tolerance, max_iter_ader_cg,
     1             tol_rootsolve, tol_superluminal,

     1             a_kerr, viscosity_coef, resistivity_coef,
     1           ambipolar_diffusion_drag_coef, molewt_ion, molewt_neut,
     1             isospeed, barotropic_stiffness,

     1             corrupted_zone,

     1             xc, yc, zc, dx, dy, dz,

     1             w, entropy_1d, gamma_eos_quadpt, molewt_quadpt,
     1             lorenz_quadpt, mu_quadpt,

     1             u, f, g, h, s, s_int)

! ----------------------------------------------------------------------

! Written by Dinshaw Balsara

! Incorporates nodal ADER CG with source terms as contained
! in the paper by Balsara et al 2010, JCP. The present version
! accomodated Galerkin elements that are continuous in time and is,
! therefore, unsuitable for use with stiff source terms.
! It is substantially faster than ADER DG. It is useful for
! situations where the hyperbolic system is large and/or the Jacobian of
! the flux matrix is difficult to evaluate. Current version is for
! single-fluid MHD with multiple species or 2-fluid MHD. It scales
! linearly as the number of species and is, therefore, very useful
! for reactive flow.

! One MHD fluid and "N_EULER" number of Euler fluids (which may also be
! repurposed as radiation fields in various frequency bins) are treated here.

! Takes the conserved variables and their spatial moments "w" as input
! and returns the space-time moments of the conserved variables as "u".
! This is all done on a reference space-time element in 3d that spans
! (x,y,z,t) in the range [-.5,.5]X[-.5,.5]X[-.5,.5]X[0,1].

! "dtmove, dx, dy, dz" are timestep and x,y,z zone sizes. The update has
! to be done in physical variables so various combinations of these
! variables are multiplied into the flux and source terms.

! "gamma, pi, rhofloor, prsfloor, smallnum" are the usual inputs.

! gamma law gas is assumed, but allowance is made for more general eos if
! the spatial moments are provided in "gamma_eos_spacetime" with the time terms zeroed.

! "max_iter_ader_cg" the maximum number of permitted iterations.

! "f, g, h, s" are fluxes in x,y,z-directions and the source term respectively.
! They are not built here.

! "s_int" is the space-time integral of the source term which is also returned.

! ----------------------------------------------------------------------

! It helps to know the space-time modal basis set into which the solution is
! being projected. Here it is for orders 2, 3, and 4.

! P0[x_] := 1 ; P1[x_] := x ; P2[x_] := ( x^2 - 1/12) ;
! P3[x_] := ( x^3 - 3 * x / 20)

! Q0[t_] := 1 ; Q1[t_] := 2 * t - 1 ; Q2[t_] := 6 * t^2 - 6 * t + 1 ;
! Q3[t_] := 20 * t^3 - 30 * t^2 + 12 * t - 1 ; DQ0[t_] := D[Q0[t],t] ;
! DQ1[t_] := D[Q1[t],t] ; DQ2[t_] := D[Q2[t],t] ; DQ3[t_] := D[Q3[t],t]

! ------------------------------

! Second order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P0[x] * P0[y] * P0[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
! [0,0,0,0], [1/2,0,0,0], [0,1/2,0,0], [0,0,1/2,0], [0,0,0,1]

! ------------------------------

! Third order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P2[x] * P0[y] * P0[z] * Q0[t] +
!                   u6 * P0[x] * P2[y] * P0[z] * Q0[t] +
!                   u7 * P0[x] * P0[y] * P2[z] * Q0[t] +
!                   u8 * P1[x] * P1[y] * P0[z] * Q0[t] +
!                   u9 * P0[x] * P1[y] * P1[z] * Q0[t] +
!                   u10 * P1[x] * P0[y] * P1[z] * Q0[t] +
!                         u11 * P0[x] * P0[y] * P0[z] * Q1[t] +
!                         u12 * P0[x] * P0[y] * P0[z] * Q2[t] +
!                   u13 * P1[x] * P0[y] * P0[z] * Q1[t] +
!                   u14 * P0[x] * P1[y] * P0[z] * Q1[t] +
!                   u15 * P0[x] * P0[y] * P1[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
!   [0,0,0,0], [1/2,0,0,0], [-1/2,0,0,0], [0,1/2,0,0], [0,-1/2,0,0],
!   [0,0,1/2,0], [0,0,-1/2,0], [0,1/2,1/2,0], [1/2,0,1/2,0], [1/2,1/2,0,0],
!   [0,0,0,1/2], [1/2,0,0,1/2], [0,1/2,0,1/2], [0,0,1/2,1/2], [0,0,0,1]

! ------------------------------

! Fourth order:

! u[x_,y_,z_,t_] := u1 * P0[x] * P0[y] * P0[z] * Q0[t] +
!                         u2 * P1[x] * P0[y] * P0[z] * Q0[t] +
!                         u3 * P0[x] * P1[y] * P0[z] * Q0[t] +
!                         u4 * P0[x] * P0[y] * P1[z] * Q0[t] +
!                   u5 * P2[x] * P0[y] * P0[z] * Q0[t] +
!                   u6 * P0[x] * P2[y] * P0[z] * Q0[t] +
!                   u7 * P0[x] * P0[y] * P2[z] * Q0[t] +
!                   u8 * P1[x] * P1[y] * P0[z] * Q0[t] +
!                   u9 * P0[x] * P1[y] * P1[z] * Q0[t] +
!                   u10 * P1[x] * P0[y] * P1[z] * Q0[t] +
!                         u11 * P3[x] * P0[y] * P0[z] * Q0[t] +
!                         u12 * P0[x] * P3[y] * P0[z] * Q0[t] +
!                         u13 * P0[x] * P0[y] * P3[z] * Q0[t] +
!                         u14 * P2[x] * P1[y] * P0[z] * Q0[t] +
!                         u15 * P2[x] * P0[y] * P1[z] * Q0[t] +
!                         u16 * P1[x] * P2[y] * P0[z] * Q0[t] +
!                         u17 * P0[x] * P2[y] * P1[z] * Q0[t] +
!                         u18 * P1[x] * P0[y] * P2[z] * Q0[t] +
!                         u19 * P0[x] * P1[y] * P2[z] * Q0[t] +
!                         u20 * P1[x] * P1[y] * P1[z] * Q0[t] +
!                   u21 * P0[x] * P0[y] * P0[z] * Q1[t] +
!                   u22 * P0[x] * P0[y] * P0[z] * Q2[t] +
!                   u23 * P0[x] * P0[y] * P0[z] * Q3[t] +
!                         u24 * P1[x] * P0[y] * P0[z] * Q1[t] +
!                         u25 * P0[x] * P1[y] * P0[z] * Q1[t] +
!                         u26 * P0[x] * P0[y] * P1[z] * Q1[t] +
!                         u27 * P1[x] * P0[y] * P0[z] * Q2[t] +
!                         u28 * P0[x] * P1[y] * P0[z] * Q2[t] +
!                         u29 * P0[x] * P0[y] * P1[z] * Q2[t] +
!                   u30 * P2[x] * P0[y] * P0[z] * Q1[t] +
!                   u31 * P0[x] * P2[y] * P0[z] * Q1[t] +
!                   u32 * P0[x] * P0[y] * P2[z] * Q1[t] +
!                   u33 * P1[x] * P1[y] * P0[z] * Q1[t] +
!                   u34 * P0[x] * P1[y] * P1[z] * Q1[t] +
!                   u35 * P1[x] * P0[y] * P1[z] * Q1[t]

! We use an ordered set of space-time nodal points at :
! [0,0,0,0], [1/2,1/2,1/2,0], [1/2,0,0,0], [-1/2,0,0,0], [0,1/2,0,0],
! [0,-1/2,0,0], [0,0,1/2,0], [0,0,-1/2,0], [0,1/2,1/2,0], [-1/2,0,1/2,0],
! [1/2,0,1/2,0], [1/2,-1/2,0,0], [1/2,1/2,0,0], [0,1/2,-1/2,0],
! [-1/4,-1/2,-1/2,0], [1/4,-1/2,-1/2,0], [-1/2,-1/4,-1/2,0],
! [-1/2,1/4,-1/2,0], [-1/2,-1/2,-1/4,0], [-1/2,-1/2,1/4,0],
! [0,0,0,1/3], [1/2,0,0,1/3], [-1/2,0,0,1/3], [0,1/2,0,1/3], [0,-1/2,0,1/3],
! [0,0,1/2,1/3], [0,0,-1/2,1/3], [0,1/2,1/2,1/3], [1/2,0,1/2,1/3],
! [1/2,1/2,0,1/3],
! [0,0,0,2/3], [1/2,0,0,2/3], [0,1/2,0,2/3], [0,0,1/2,2/3], [0,0,0,1]

! ----------------------------------------------------------------------

      IMPLICIT NONE

! Variables declared here go through the subroutine interface.

      INTEGER oned_begin, oned_end, oned_off, max_iter_ader_cg, igeom

      REAL dtmove, gamma, gamma_euler, pi, gasconst, ader_tolerance, 
     1     tol_rootsolve, tol_superluminal,
     1     rhofloor, prsfloor, smallnum,
     1     a_kerr, viscosity_coef, resistivity_coef,
     1     ambipolar_diffusion_drag_coef, molewt_ion, molewt_neut,
     1     isospeed, barotropic_stiffness

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1                     corrupted_zone

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1                  xc, yc, zc, dx, dy, dz

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_ZONECEN_BASES) :: w

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, 5) :: entropy_1d

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  1: N_EULER + 1, NUM_SPACETIME_QUADPTS) ::
     1                gamma_eos_quadpt, molewt_quadpt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NUM_SPACETIME_QUADPTS) ::
     1                lorenz_quadpt, mu_quadpt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_SPACETIME_BASES) ::
     1                                              u, f, g, h, s

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, NVAR) :: s_int

! ----------------------------------------------------------------------

! Variables declared past this point don't go through the subroutine interface.

      INTEGER ix, iter_ader, iquadpt, ibases, ichk, itest, i_euler,
     1        ivar_lo, ivar

      REAL fourpi, ob4pi, ob8pi, ob3, tempaa, tempab, tempac,
     1     x1, x2, x3, obx1, cos_x2, sin_x2, cot_x2

! These are the flow variables (and their spatial derivatives, if needed)
! at the quadrature points.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NUM_SPACETIME_BASES) ::

     1      rho, prs, eng_t, vlx, vly, vlz, bfx, bfy, bfz,
     1      entropy,

     1      drho_dx, dtemperature_dx, dvlx_dx, dvly_dx, dvlz_dx,
     1      dbfx_dx, dbfy_dx, dbfz_dx, dv2_dx, db2_dx,

     1      drho_dy, dtemperature_dy, dvlx_dy, dvly_dy, dvlz_dy,
     1      dbfx_dy, dbfy_dy, dbfz_dy, dv2_dy, db2_dy,

     1      drho_dz, dtemperature_dz, dvlx_dz, dvly_dz, dvlz_dz,
     1      dbfx_dz, dbfy_dz, dbfz_dz, dv2_dz, db2_dz,

     1      div_vl, temperature

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0: NFLUID, NUM_SPACETIME_BASES) :: specie

#if ( N_EULER >= 1)
      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NUM_SPACETIME_BASES, N_EULER) ::
     1      rho_euler, prs_euler, eng_t_euler,
     1      vlx_euler, vly_euler, vlz_euler

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1             0: NFLUID_EULER, NUM_SPACETIME_BASES, N_EULER) ::
     1      specie_euler
#endif

! These are the conserved variables, fluxes and source terms at the
! quadrature points.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_SPACETIME_BASES) ::

     1      u_quadpt, f_quadpt, g_quadpt, h_quadpt, s_quadpt,
     1      tu, tf, tg, th, ts,

     1      du_quadpt_dx, du_quadpt_dy, du_quadpt_dz,
     1      fni_quadpt, gni_quadpt, hni_quadpt

! List of nodal points in each zone.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  1:4, NUM_SPACETIME_BASES) :: nodal_points

! ----------------------------------------------------------------------

! Make some sanity checks here.

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_NODAL_3D_TIME_EVOL_MHD, NFI","ELD .NE. 7",
     1              NFIELD
      STOP
      END IF

      IF ( NFIELD_EULER .NE. 5) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_NODAL_3D_TIME_EVOL_MHD, NFI",
     1  "ELD_EULER .NE. 5", NFIELD
      STOP
      END IF

#if ( NUMORDER == 2)
      IF ( NUM_ZONECEN_BASES .NE. 4) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_BASES .NE. 4", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 5) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_BASES .NE. 5", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 6) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 6", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 7", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

#if ( NUMORDER == 3)
      IF ( NUM_ZONECEN_BASES .NE. 10) THEN
      WRITE ( 6, *)
     1  "Stopping in third order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_BASES .NE. 10", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 15) THEN
      WRITE ( 6, *)
     1  "Stopping in third order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_BASES .NE. 15", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 13) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 15", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 20) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 22", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

#if ( NUMORDER == 4)
      IF ( NUM_ZONECEN_BASES .NE. 20) THEN
      WRITE ( 6, *)
     1  "Stopping in fourth order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_BASES .NE. 20", NUM_ZONECEN_BASES
      STOP
      END IF

      IF ( NUM_SPACETIME_BASES .NE. 35) THEN
      WRITE ( 6, *)
     1  "Stopping in fourth order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_BASES .NE. 35", NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( NUM_ZONECEN_QUADPTS .NE. 21) THEN
      WRITE ( 6, *)
     1   "Stopping in second order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_ZONECEN","_QUADPTS .NE. 33", NUM_ZONECEN_QUADPTS
      STOP
      END IF

      IF ( NUM_SPACETIME_QUADPTS .NE. 42) THEN
      WRITE ( 6, *)
     1  "Stopping in second order ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1     "NUM_SPACETIME","_QUADPTS .NE. 55", NUM_SPACETIME_QUADPTS
      STOP
      END IF
#endif

      ivar =
     1  NFIELD + 1 + NFLUID + N_EULER * (NFIELD_EULER + NFLUID_EULER)

      IF ( NVAR .NE. ivar) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1  "NVAR .NE. ivar","NVA","R = ", NVAR," ivar = ", ivar
      STOP
      END IF

      IF ( max_iter_ader_cg .LT. NUMORDER - 1 ) THEN
      WRITE ( 6, *)
     1  "Stopping in ADER_CG_NODAL_3D_TIME_EVOL_MHD, ",
     1  "max_iter_ader_cg .LT. NUMORDER ; max_iter_ader_cg = ",
     1   max_iter_ader_cg, NUMORDER - 1
      STOP
      END IF

! ----------------------------------------------------------------------

! Set up some constants that are used all over in this subroutine.

      fourpi = 4.0 * pi
      ob4pi = 1.0 / ( 4.0 * pi)
      ob8pi = 1.0 / ( 8.0 * pi)
      ob3 = 1.0 / 3.0

! ----------------------------------------------------------------------

! STEP 0) Set up the nodal points that are to be used at each of the orders.

! ------------------------------
#if ( NUMORDER == 2)
! cg, 2nd order nodal_points

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1:4, 01) = (/  0.0 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 02) = (/  0.5 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 03) = (/  0.0 ,  0.5 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 04) = (/  0.0 ,  0.0 ,  0.5 , 0.0  /) 

      nodal_points ( ix, 1:4, 05) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /) 

      END DO
! End of "DO ix =" loop.

#endif
! End "NUMORDER == 2".
! ------------------------------
#if ( NUMORDER == 3)
! cg, 3rd order nodal_points

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1:4, 01) = (/  0.0 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 02) = (/  0.5 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 03) = (/ -0.5 ,  0.0 ,  0.0 , 0.0  /) 
      nodal_points ( ix, 1:4, 04) = (/  0.0 ,  0.5 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 05) = (/  0.5 , -0.5 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 06) = (/  0.0 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 07) = (/  0.0 ,  0.0 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 08) = (/  0.0 ,  0.5 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 09) = (/  0.5 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 10) = (/  0.5 ,  0.5 ,  0.0 , 0.0  /)

      nodal_points ( ix, 1:4, 11) = (/  0.0 ,  0.0 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 12) = (/  0.5 ,  0.0 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 13) = (/  0.0 ,  0.5 ,  0.0 , 0.5  /)
      nodal_points ( ix, 1:4, 14) = (/  0.0 ,  0.0 ,  0.5 , 0.5  /)

      nodal_points ( ix, 1:4, 15) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /)

      END DO
! End of "DO ix =" loop.

#endif
! End "NUMORDER == 3".
! ------------------------------
#if ( NUMORDER == 4)
! cg, 4th order nodal_points

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1:4, 01) = (/  0.0 ,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 02) = (/  0.5 ,  0.5 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 03) = (/  0.5 ,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 04) = (/ -0.5 ,  0.0 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 05) = (/  0.0 ,  0.5 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 06) = (/  0.0 , -0.5 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 07) = (/  0.0 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 08) = (/  0.0 ,  0.0 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 09) = (/  0.0 ,  0.5 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 10) = (/ -0.5 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 11) = (/  0.5 ,  0.0 ,  0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 12) = (/  0.5 , -0.5 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 13) = (/  0.5 ,  0.5 ,  0.0 , 0.0  /)
      nodal_points ( ix, 1:4, 14) = (/  0.0 ,  0.5 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 15) = (/ -0.25, -0.5 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 16) = (/  0.25, -0.5 , -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 17) = (/ -0.5 , -0.25, -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 18) = (/ -0.5 ,  0.25, -0.5 , 0.0  /)
      nodal_points ( ix, 1:4, 19) = (/ -0.5 , -0.5 , -0.25, 0.0  /)
      nodal_points ( ix, 1:4, 20) = (/ -0.5 , -0.5 ,  0.25, 0.0  /)

      nodal_points ( ix, 1:4, 21) = (/  0.0 ,  0.0 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 22) = (/  0.5 ,  0.0 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 23) = (/ -0.5 ,  0.0 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 24) = (/  0.0 ,  0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 25) = (/  0.0 , -0.5 ,  0.0 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 26) = (/  0.0 ,  0.0 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 27) = (/  0.0 ,  0.0 , -0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 28) = (/  0.0 ,  0.5 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 29) = (/  0.5 ,  0.0 ,  0.5 , 1.0/3.0 /)
      nodal_points ( ix, 1:4, 30) = (/  0.5 ,  0.5 ,  0.0 , 1.0/3.0 /)

      nodal_points ( ix, 1:4, 31) = (/  0.0 ,  0.0 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 32) = (/  0.5 ,  0.0 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 33) = (/  0.0 ,  0.5 ,  0.0 , 2.0/3.0 /)
      nodal_points ( ix, 1:4, 34) = (/  0.0 ,  0.0 ,  0.5 , 2.0/3.0 /)

      nodal_points ( ix, 1:4, 35) = (/  0.0 ,  0.0 ,  0.0 , 1.0  /)

      END DO
! End of "DO ix =" loop.

#endif
! End "NUMORDER == 4".
! ------------------------------

#if ( ENFORCE_POSITIVITY == 1)

! While "nodal_points" still represents positions in the reference space-time
! hypercube, we use it to build the "entropy" at all the nodes.

      DO iquadpt = 1, NUM_SPACETIME_QUADPTS

      DO ix = oned_begin - oned_off, oned_end + oned_off

      entropy ( ix, iquadpt) = entropy_1d ( ix, 01)
     1    + nodal_points ( ix, 1, iquadpt) * entropy_1d ( ix, 02)
     1    + nodal_points ( ix, 2, iquadpt) * entropy_1d ( ix, 03)
     1    + nodal_points ( ix, 3, iquadpt) * entropy_1d ( ix, 04)
     1    + nodal_points ( ix, 4, iquadpt) * entropy_1d ( ix, 05)

      END DO ! End of "DO ix =" loop.

      END DO ! End of "DO iquadpt =" loop.

#endif
! End of "ENFORCE_POSITIVITY == 1".

! ------------------------------

! Set up the nodal points in a global space-time. These are coordinate
! values.

      DO ix = oned_begin - oned_off, oned_end + oned_off

      nodal_points ( ix, 1, :) = nodal_points ( ix, 1, :) * dx ( ix)
     1                                                    + xc ( ix)
      nodal_points ( ix, 2, :) = nodal_points ( ix, 2, :) * dy ( ix)
     1                                                    + yc ( ix)
      nodal_points ( ix, 3, :) = nodal_points ( ix, 3, :) * dz ( ix)
     1                                                    + zc ( ix)
      nodal_points ( ix, 4, :) = nodal_points ( ix, 4, :) * dtmove 

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------

! STEP 1) Having come into the subroutine with variables loaded into "w",
! load the appropriate spatial variation of "w" into "u" while setting
! the temporal variation of "u" as well as "f,g,h & s" to zero initially.
! Also, zero out "s_quadpt" so that it doesn't interfere negatively if
! the source terms are not built later on.

! ------------------------------
#if ( ISOTHERMAL_EQNS == 1)
! Zero out the fifth component, i.e. energy density, for isothermal flow.

      DO ix = oned_begin - oned_off, oned_end + oned_off
      w ( ix, 5, 1: NUM_ZONECEN_BASES) = 0.0
      END DO

#if ( N_EULER >= 1)
      ivar_lo = NFIELD + 1 + NFLUID
      DO i_euler = 1, N_EULER
      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)
      DO ix = oned_begin - oned_off, oned_end + oned_off
      w ( ix, ivar + 5, 1: NUM_ZONECEN_BASES) = 0.0
      END DO
      END DO
#endif
! End of "N_EULER >= 1"

#endif
! End of "ISOTHERMAL_EQNS == 1".
! ------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

      u ( ix, :, 1: NUM_ZONECEN_BASES)=w ( ix, :, 1: NUM_ZONECEN_BASES)
      u ( ix, :, NUM_ZONECEN_BASES + 1: NUM_SPACETIME_BASES) = 0.0

      f ( ix, :, :) = 0.0
      g ( ix, :, :) = 0.0
      h ( ix, :, :) = 0.0
      s ( ix, :, :) = 0.0

      s_quadpt ( ix, :, :) = 0.0

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------

! STEP 2) Use the spatial modes "u ( :, 1: NUM_ZONECEN_BASES)" at time t = 0
! to make the nodes "u_quadpt ( :, 1: NUM_ZONECEN_BASES)" at t = 0.
! It is also a good idea to initialize the t > 0 nodes here so that they
! start with reasonably good values when the iteration begins.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

      u_quadpt ( ix, :, 01) = u ( ix, :, 01)

      u_quadpt ( ix, :, 02) = u ( ix, :, 01) + 0.5 * u ( ix, :, 02)

      u_quadpt ( ix, :, 03) = u ( ix, :, 01) + 0.5 * u ( ix, :, 03)

      u_quadpt ( ix, :, 04) = u ( ix, :, 01) + 0.5 * u ( ix, :, 04)

      u_quadpt ( ix, :, 05) = u_quadpt ( ix, :, 01)

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

      u_quadpt ( ix, :, 01) = u ( ix, :, 01) - ( u ( ix, :, 05)
     1 + u ( ix, :, 06) + u ( ix, :, 07) ) / 12.0

      u_quadpt ( ix, :, 02) = u_quadpt ( ix, :, 01)
     1 + 0.5 * u ( ix, :, 02) + 0.25 * u ( ix, :, 05)

      u_quadpt ( ix, :, 03) = 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 02) + 0.5 * u ( ix, :, 05)

      u_quadpt ( ix, :, 04) = u_quadpt ( ix, :, 01)
     1 + 0.5 * u ( ix, :, 03) + 0.25 * u ( ix, :, 06)

      u_quadpt ( ix, :, 05) = 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 04) + 0.5 * u ( ix, :, 06)

      u_quadpt ( ix, :, 06) = u_quadpt ( ix, :, 01)
     1 + 0.5 * u ( ix, :, 04) + 0.25 * u ( ix, :, 07)

      u_quadpt ( ix, :, 07) = 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 06) + 0.5 * u ( ix, :, 07)

      u_quadpt ( ix, :, 08) = u_quadpt ( ix, :, 04)
     1 + u_quadpt ( ix, :, 06) - u_quadpt ( ix, :, 01)
     1 + 0.25 * u ( ix, :, 09)

      u_quadpt ( ix, :, 09) = u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 06) - u_quadpt ( ix, :, 01)
     1 + 0.25 * u ( ix, :, 10)

      u_quadpt ( ix, :, 10) = u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 04) - u_quadpt ( ix, :, 01)
     1 + 0.25 * u ( ix, :, 08)

      u_quadpt ( ix, :, 11) = u_quadpt ( ix, :, 01)

      u_quadpt ( ix, :, 12) = u_quadpt ( ix, :, 02)

      u_quadpt ( ix, :, 13) = u_quadpt ( ix, :, 04)

      u_quadpt ( ix, :, 14) = u_quadpt ( ix, :, 06)

      u_quadpt ( ix, :, 15) = u_quadpt ( ix, :, 01)

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : modal -> nodal

      DO ix = oned_begin - oned_off, oned_end + oned_off

      tu ( ix, :, 01) = 0.25 * u ( ix, :, 02) - 0.125 * ( u ( ix, :, 08)
     1 + u ( ix, :, 10) ) + ( - 21.0 * u ( ix, :, 11)
     1 + 40.0 * u ( ix, :, 16) + 40.0 * u ( ix, :, 18)
     1 + 60.0 * u ( ix, :, 20) ) / 960.0

      tu ( ix, :, 02) = 0.25 * u ( ix, :, 03) - 0.125 * ( u ( ix, :, 08)
     1 + u ( ix, :, 09) ) + ( - 21.0 * u ( ix, :, 12)
     1 + 40.0 * u ( ix, :, 14) + 40.0 * u ( ix, :, 19)
     1 + 60.0 * u ( ix, :, 20) ) / 960.0

      tu ( ix, :, 03) = 0.25 * u ( ix, :, 04) - 0.125 * ( u ( ix, :, 09)
     1 + u ( ix, :, 10) ) + ( - 21.0 * u ( ix, :, 13)
     1 + 40.0 * u ( ix, :, 15) + 40.0 * u ( ix, :, 17)
     1 + 60.0 * u ( ix, :, 20) ) / 960.0

      u_quadpt ( ix, :, 01) = u ( ix, :, 01) - ( u ( ix, :, 05)
     1 + u ( ix, :, 06) + u ( ix, :, 07) ) / 12.0

      u_quadpt ( ix, :, 03) = u_quadpt ( ix, :, 01)
     1 + 0.5 * u ( ix, :, 02) + 0.25 * u ( ix, :, 05)
     1 + 0.05 * u ( ix, :, 11) - ( u ( ix, :, 16)
     1 + u ( ix, :, 18) ) / 24.0

      u_quadpt ( ix, :, 04) = 2 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 03) + 0.5 * u ( ix, :, 05)

      u_quadpt ( ix, :, 05) = u_quadpt ( ix, :, 01)
     1 + 0.5 * u ( ix, :, 03) + 0.25 * u ( ix, :, 06)
     1 + 0.05 * u ( ix, :, 12) - ( u ( ix, :, 14)
     1 + u ( ix, :, 19) ) / 24.0

      u_quadpt ( ix, :, 06) = 2 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 05) + 0.5 * u ( ix, :, 06)

      u_quadpt ( ix, :, 07) = u_quadpt ( ix, :, 01)
     1 + 0.5 * u ( ix, :, 04) + 0.25 * u ( ix, :, 07)
     1 + 0.05 * u ( ix, :, 13) - ( u ( ix, :, 15)
     1 + u ( ix, :, 17) ) / 24.0

      u_quadpt ( ix, :, 08) = 2 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 07) + 0.5 * u ( ix, :, 07)

      u_quadpt ( ix, :, 09) = u_quadpt ( ix, :, 05)
     1 + u_quadpt ( ix, :, 07) - u_quadpt ( ix, :, 01)
     1 + 0.25 * u ( ix, :, 09) + 0.125 * ( u ( ix, :, 17)
     1 + u ( ix, :, 19) )

      u_quadpt ( ix, :, 11) = u_quadpt ( ix, :, 03)
     1 + u_quadpt ( ix, :, 07) - u_quadpt ( ix, :, 01)
     1 + 0.25 * u ( ix, :, 10) + 0.125 * ( u ( ix, :, 15)
     1 + u ( ix, :, 18) )

      u_quadpt ( ix, :, 13) = u_quadpt ( ix, :, 03)
     1 + u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 01)
     1 + 0.25 * u ( ix, :, 08) + 0.125 * ( u ( ix, :, 14)
     1 + u ( ix, :, 16) )

      u_quadpt ( ix, :, 10) = 2 * u_quadpt ( ix, :, 07)
     1 - u_quadpt ( ix, :, 11) + 0.5 * u ( ix, :, 05)
     1 + 0.25 * u ( ix, :, 15)

      u_quadpt ( ix, :, 12) = 2 * u_quadpt ( ix, :, 03)
     1 - u_quadpt ( ix, :, 13) + 0.5 * u ( ix, :, 06)
     1 + 0.25 * u ( ix, :, 16)

      u_quadpt ( ix, :, 14) = 2 * u_quadpt ( ix, :, 05)
     1 - u_quadpt ( ix, :, 09) + 0.5 * u ( ix, :, 07)
     1 + 0.25 * u ( ix, :, 19)

      u_quadpt ( ix, :, 02) = u_quadpt ( ix, :, 11)
     1 + u_quadpt ( ix, :, 09) - u_quadpt ( ix, :, 07)
     1 + 0.25 * u ( ix, :, 08) + 0.125 * ( u ( ix, :, 14)
     1 + u ( ix, :, 16) + u ( ix, :, 20) )

      u_quadpt ( ix, :, 15) = 2 * u_quadpt ( ix, :, 08)
     1 - u_quadpt ( ix, :, 14) - tu ( ix, :, 01)
     1 + 0.0625 * u ( ix, :, 05) + 0.5 * u ( ix, :, 06)
     1 - 0.03125 * ( u ( ix, :, 14) + u ( ix, :, 15) )
     1 - 0.25 * u ( ix, :, 17)

      u_quadpt ( ix, :, 16) = u_quadpt ( ix, :, 15)
     1 + 2.0 * tu ( ix, :, 01)

      u_quadpt ( ix, :, 17) = 2.0 * u_quadpt ( ix, :, 04)
     1 - u_quadpt ( ix, :, 10) - tu ( ix, :, 02)
     1 + 0.0625 * u ( ix, :, 06) + 0.5 * u ( ix, :, 07)
     1 - 0.03125 * ( u ( ix, :, 16) + u ( ix, :, 17) )
     1 - 0.25 * u ( ix, :, 18)

      u_quadpt ( ix, :, 18) = u_quadpt ( ix, :, 17)
     1 + 2.0 * tu ( ix, :, 02)

      u_quadpt ( ix, :, 19) = 2.0 * u_quadpt ( ix, :, 06)
     1 - u_quadpt ( ix, :, 12) - tu ( ix, :, 03)
     1 + 0.5 * u ( ix, :, 05) + 0.0625 * u ( ix, :, 07)
     1 - 0.25 * u ( ix, :, 14)  - 0.03125 * ( u ( ix, :, 18)
     1 + u ( ix, :, 19) )

      u_quadpt ( ix, :, 20) = u_quadpt ( ix, :, 19)
     1 + 2.0 * tu ( ix, :, 03)

      u_quadpt ( ix, :, 21) = u_quadpt ( ix, :, 01)

      u_quadpt ( ix, :, 22) = u_quadpt ( ix, :, 03)

      u_quadpt ( ix, :, 23) = u_quadpt ( ix, :, 04)

      u_quadpt ( ix, :, 24) = u_quadpt ( ix, :, 05)

      u_quadpt ( ix, :, 25) = u_quadpt ( ix, :, 06)

      u_quadpt ( ix, :, 26) = u_quadpt ( ix, :, 07)

      u_quadpt ( ix, :, 27) = u_quadpt ( ix, :, 08)

      u_quadpt ( ix, :, 28) = u_quadpt ( ix, :, 09)

      u_quadpt ( ix, :, 29) = u_quadpt ( ix, :, 11)

      u_quadpt ( ix, :, 30) = u_quadpt ( ix, :, 13)

      u_quadpt ( ix, :, 31) = u_quadpt ( ix, :, 01)

      u_quadpt ( ix, :, 32) = u_quadpt ( ix, :, 03)

      u_quadpt ( ix, :, 33) = u_quadpt ( ix, :, 05)

      u_quadpt ( ix, :, 34) = u_quadpt ( ix, :, 07)

      u_quadpt ( ix, :, 35) = u_quadpt ( ix, :, 01)

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------

! STEP 3) Get the primitive variables from the conserved variables at
! the spatial nodes at t=0. The density and pressure variables are
! given positive floor values to ensure that they can be used for physically
! meaningful source terms.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, 1, iquadpt) =
     1                   AMAX1 ( u_quadpt ( ix, 1, iquadpt), rhofloor)

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempaa = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempaa
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempaa
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempaa

      bfx ( ix, iquadpt) = u_quadpt ( ix, 6, iquadpt)
      bfy ( ix, iquadpt) = u_quadpt ( ix, 7, iquadpt)
      bfz ( ix, iquadpt) = u_quadpt ( ix, 8, iquadpt)

      tempab = 0.5 * rho ( ix, iquadpt) * ( vlx ( ix, iquadpt)**2
     1         + vly ( ix, iquadpt)**2 + vlz ( ix, iquadpt)**2)

      tempac = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      eng_t ( ix, iquadpt) = u_quadpt ( ix, 5, iquadpt)
     1                     - tempab - tempac

      prs ( ix, iquadpt) = eng_t ( ix, iquadpt)
     1                   * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

#if ( ENFORCE_POSITIVITY == 1)
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = entropy ( ix, iquadpt)
     1* rho ( ix, iquadpt)**( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)
      prs ( ix, iquadpt) = AMAX1 ( prs ( ix, iquadpt), prsfloor)

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#else
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = prsfloor

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#endif

#if ( NFLUID >= 1)
      specie ( ix, 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( ix, 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( ix, 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( ix, 1: NFLUID, iquadpt) * rho ( ix, iquadpt)
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, ivar + 1, iquadpt) =
     1         AMAX1 ( u_quadpt ( ix, ivar + 1, iquadpt), rhofloor)

      rho_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 1, iquadpt)

      tempaa = 1.0 / rho_euler ( ix, iquadpt, i_euler)

      vlx_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 2, iquadpt) * tempaa
      vly_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 3, iquadpt) * tempaa
      vlz_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 4, iquadpt) * tempaa

      tempab = 0.5 * rho_euler ( ix, iquadpt, i_euler)
     1           * ( vlx_euler ( ix, iquadpt, i_euler)**2
     1             + vly_euler ( ix, iquadpt, i_euler)**2
     1             + vlz_euler ( ix, iquadpt, i_euler)**2)

      prs_euler ( ix, iquadpt, i_euler) =
     1         ( u_quadpt ( ix, ivar + 5, iquadpt) - tempab)
     1       * ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      prs_euler ( ix, iquadpt, i_euler) =
     1    AMAX1 ( prs_euler ( ix, iquadpt, i_euler), prsfloor)

#if ( ISOTHERMAL_EQNS == 1)
      prs_euler ( ix, iquadpt, i_euler) =
     1     isospeed**2 * rho_euler ( ix, iquadpt, i_euler)
     1   + barotropic_stiffness
     1   * rho_euler ( ix, iquadpt, i_euler)**gamma
#endif

      eng_t_euler ( ix, iquadpt, i_euler) =
     1         prs_euler ( ix, iquadpt, i_euler)
     1       / ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      IF ( prs_euler ( ix, iquadpt, i_euler) .LT. prsfloor * 1.00000001)
     1   u_quadpt ( ix, ivar + 5, iquadpt) =
     1                   eng_t_euler ( ix, iquadpt, i_euler) + tempab

#if ( NFLUID_EULER >= 1)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) = tempaa
     1       * u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1                    ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) =
     1    AMIN1 ( AMAX1 (
     1    specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler),
     1            0.0), 1.0 )

      u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1       specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler)
     1     * rho_euler ( ix, iquadpt, i_euler)

#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 4.a) Build fluxes ( f_quadpt, g_quadpt, h_quadpt) from the primitive
! variables that have now become available at spatial nodes at t=0.

! Multiply "f_quadpt" by "dtmove/dx" and so on.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! -----

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

! -----

      f_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 2, iquadpt)

      f_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfx ( ix, iquadpt)**2 * ob4pi

      f_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vly ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi

      f_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi

      f_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab

      f_quadpt ( ix, 6, iquadpt) = 0.0

      f_quadpt ( ix, 7, iquadpt) =
     1           vlx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         - vly ( ix, iquadpt) * bfx ( ix, iquadpt)

      f_quadpt ( ix, 8, iquadpt) =
     1           vlx ( ix, iquadpt) * bfz ( ix, iquadpt)
     1         - vlz ( ix, iquadpt) * bfx ( ix, iquadpt)

#if ( NFLUID >= 1)
      f_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlx ( ix, iquadpt)
#endif

! -----

      g_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 3, iquadpt)

      g_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         - bfy ( ix, iquadpt) * bfx ( ix, iquadpt) * ob4pi

      g_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vly ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfy ( ix, iquadpt)**2 * ob4pi

      g_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         - bfy ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi

      g_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vly ( ix, iquadpt)
     1                       - bfy ( ix, iquadpt) * tempab

      g_quadpt ( ix, 6, iquadpt) =
     1           vly ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         - vlx ( ix, iquadpt) * bfy ( ix, iquadpt)

      g_quadpt ( ix, 7, iquadpt) = 0.0

      g_quadpt ( ix, 8, iquadpt) =
     1           vly ( ix, iquadpt) * bfz ( ix, iquadpt)
     1         - vlz ( ix, iquadpt) * bfy ( ix, iquadpt)

#if ( NFLUID >= 1)
      g_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vly ( ix, iquadpt)
#endif

! -----

      h_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 4, iquadpt)

      h_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         - bfz ( ix, iquadpt) * bfx ( ix, iquadpt) * ob4pi

      h_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vly ( ix, iquadpt)
     1         - bfz ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi

      h_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfz ( ix, iquadpt)**2 * ob4pi

      h_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vlz ( ix, iquadpt)
     1                       - bfz ( ix, iquadpt) * tempab

      h_quadpt ( ix, 6, iquadpt) =
     1           vlz ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         - vlx ( ix, iquadpt) * bfz ( ix, iquadpt)

      h_quadpt ( ix, 7, iquadpt) =
     1           vlz ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         - vly ( ix, iquadpt) * bfz ( ix, iquadpt)

      h_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      h_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlz ( ix, iquadpt)
#endif

! -----

#if ( ISOTHERMAL_EQNS == 1)
      f_quadpt ( ix, 5, iquadpt) = 0.0
      g_quadpt ( ix, 5, iquadpt) = 0.0
      h_quadpt ( ix, 5, iquadpt) = 0.0
#endif

! -----

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! -----

      f_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)

      f_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vlx_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      f_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vlx_euler ( ix, iquadpt, i_euler)
#endif

! -----

      g_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)

      g_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vly_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      g_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vly_euler ( ix, iquadpt, i_euler)
#endif

! -----

      h_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)

      h_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vlz_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      h_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vlz_euler ( ix, iquadpt, i_euler)
#endif

! -----

#if ( ISOTHERMAL_EQNS == 1)
      f_quadpt ( ix, ivar + 5, iquadpt) = 0.0
      g_quadpt ( ix, ivar + 5, iquadpt) = 0.0
      h_quadpt ( ix, ivar + 5, iquadpt) = 0.0
#endif

! -----

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ------------------------------

! To make cylindrical and spherical geometries accessible to non-rel.
! and relativistic flow, we put this part after either flux construction.

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      IF ( igeom == 2 .or. igeom == 3) THEN

         g_quadpt ( ix, :, iquadpt) = g_quadpt ( ix, :, iquadpt) * obx1

      END IF

      IF ( igeom == 3) THEN

         sin_x2 = 1.0 / ( x1 * sin ( x2))

         h_quadpt ( ix, :, iquadpt) = h_quadpt ( ix, :, iquadpt)
     1                              * sin_x2

      END IF

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 4.b) Optional step needed only for viscous, resistive or conductive MHD.

! ASSUME that "du_quadpt_dx, du_quadpt_dy and du_quadpt_dz" have been built at
! all orders (NOT DONE YET).
! We set about building the derivatives of the primitive
! variables and various auxiliary variables at the nodal points. These are
! needed only in the presence of viscous, resistive or conduction terms.

! These derivatives are then used to build non-ideal fluxes at the quadrature
! points given by "fni_quadpt, gni_quadpt and hni_quadpt".

! Lastly, the non-ideal fluxes are added into the ideal fluxes

! ----------------------------------------------------------------------

      IF ( ( ( viscosity_coef .GT. 0.0) .OR.
     1       ( resistivity_coef .GT. 0.0)) .AND. ( NUMORDER >= 3)) THEN

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = 1.0 / rho ( ix, iquadpt)

      drho_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 1, iquadpt)
      drho_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 1, iquadpt)
      drho_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 1, iquadpt)


      dvlx_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvlx_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvlx_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dvly_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvly_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvly_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dvlz_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvlz_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvlz_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dbfx_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 6, iquadpt)
      dbfx_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 6, iquadpt)
      dbfx_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 6, iquadpt)

      dbfy_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 7, iquadpt)
      dbfy_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 7, iquadpt)
      dbfy_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 7, iquadpt)

      dbfz_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 8, iquadpt)
      dbfz_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 8, iquadpt)
      dbfz_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 8, iquadpt)


      dv2_dx ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dx ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dx ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dx ( ix, iquadpt) )

      dv2_dy ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dy ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dy ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dy ( ix, iquadpt) )

      dv2_dz ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dz ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dz ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dz ( ix, iquadpt) )

      db2_dx ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dx ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dx ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dx ( ix, iquadpt) )

      db2_dy ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dy ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dy ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dy ( ix, iquadpt) )

      db2_dz ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dz ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dz ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dz ( ix, iquadpt) )


      tempac = vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1       + vlz ( ix, iquadpt)**2

      tempab = molewt_quadpt ( ix, 1, iquadpt)
     1       / ( gasconst * rho ( ix, iquadpt) )

      temperature ( ix, iquadpt) = prs ( ix, iquadpt) * tempab

      tempab = tempab * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      dtemperature_dx ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dx ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dx ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dx ( ix, iquadpt)
     1    - ob8pi * db2_dx ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dx ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

      dtemperature_dy ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dy ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dy ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dy ( ix, iquadpt)
     1    - ob8pi * db2_dy ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dy ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

      dtemperature_dz ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dz ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dz ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dz ( ix, iquadpt)
     1    - ob8pi * db2_dz ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dz ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

! The derivatives are within the reference element's local coordinates.
! To get true physical fluxes, we need to divide by the size of the zone here.

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = 1.0 / dx ( ix)
      tempab = 1.0 / dy ( ix)
      tempac = 1.0 / dz ( ix)

      drho_dx ( ix, iquadpt) = tempaa * drho_dx ( ix, iquadpt)
      dtemperature_dx ( ix, iquadpt) =
     1                 tempaa * dtemperature_dx ( ix, iquadpt)
      dvlx_dx ( ix, iquadpt) = tempaa * dvlx_dx ( ix, iquadpt)
      dvly_dx ( ix, iquadpt) = tempaa * dvly_dx ( ix, iquadpt)
      dvlz_dx ( ix, iquadpt) = tempaa * dvlz_dx ( ix, iquadpt)
      dbfx_dx ( ix, iquadpt) = tempaa * dbfx_dx ( ix, iquadpt)
      dbfy_dx ( ix, iquadpt) = tempaa * dbfy_dx ( ix, iquadpt)
      dbfz_dx ( ix, iquadpt) = tempaa * dbfz_dx ( ix, iquadpt)
      dv2_dx ( ix, iquadpt) = tempaa * dv2_dx ( ix, iquadpt)
      db2_dx ( ix, iquadpt) = tempaa * db2_dx ( ix, iquadpt)

      drho_dy ( ix, iquadpt) = tempab * drho_dy ( ix, iquadpt)
      dtemperature_dy ( ix, iquadpt) =
     1                 tempab * dtemperature_dy ( ix, iquadpt)
      dvlx_dy ( ix, iquadpt) = tempab * dvlx_dy ( ix, iquadpt)
      dvly_dy ( ix, iquadpt) = tempab * dvly_dy ( ix, iquadpt)
      dvlz_dy ( ix, iquadpt) = tempab * dvlz_dy ( ix, iquadpt)
      dbfx_dy ( ix, iquadpt) = tempab * dbfx_dy ( ix, iquadpt)
      dbfy_dy ( ix, iquadpt) = tempab * dbfy_dy ( ix, iquadpt)
      dbfz_dy ( ix, iquadpt) = tempab * dbfz_dy ( ix, iquadpt)
      dv2_dy ( ix, iquadpt) = tempab * dv2_dy ( ix, iquadpt)
      db2_dy ( ix, iquadpt) = tempab * db2_dy ( ix, iquadpt)

      drho_dz ( ix, iquadpt) = tempac * drho_dz ( ix, iquadpt)
      dtemperature_dz ( ix, iquadpt) =
     1                 tempac * dtemperature_dz ( ix, iquadpt)
      dvlx_dz ( ix, iquadpt) = tempac * dvlx_dz ( ix, iquadpt)
      dvly_dz ( ix, iquadpt) = tempac * dvly_dz ( ix, iquadpt)
      dvlz_dz ( ix, iquadpt) = tempac * dvlz_dz ( ix, iquadpt)
      dbfx_dz ( ix, iquadpt) = tempac * dbfx_dz ( ix, iquadpt)
      dbfy_dz ( ix, iquadpt) = tempac * dbfy_dz ( ix, iquadpt)
      dbfz_dz ( ix, iquadpt) = tempac * dbfz_dz ( ix, iquadpt)
      dv2_dz ( ix, iquadpt) = tempac * dv2_dz ( ix, iquadpt)
      db2_dz ( ix, iquadpt) = tempac * db2_dz ( ix, iquadpt)

      div_vl ( ix, iquadpt) = dvlx_dx ( ix, iquadpt)
     1                      + dvly_dy ( ix, iquadpt)
     1                      + dvlz_dz ( ix, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = rho ( ix, iquadpt) * viscosity_coef

! ----------

      fni_quadpt ( ix, 1, iquadpt) = 0.0

      fni_quadpt ( ix, 2, iquadpt) =
     1               2.0 * tempaa * ( dvlx_dx ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      fni_quadpt ( ix, 3, iquadpt) =
     1                     tempaa * ( dvlx_dy ( ix, iquadpt)
     1                              + dvly_dx ( ix, iquadpt) )

      fni_quadpt ( ix, 4, iquadpt) =
     1                     tempaa * ( dvlx_dz ( ix, iquadpt)
     1                              + dvlz_dx ( ix, iquadpt) )

      fni_quadpt ( ix, 5, iquadpt) =
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vlx ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvlx_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvlx_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvlx_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dx ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dx ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfx_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfx_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfx_dz ( ix, iquadpt) )

      fni_quadpt ( ix, 6, iquadpt) = 0.0

      fni_quadpt ( ix, 7, iquadpt) =
     1   resistivity_coef * ( dbfy_dx ( ix, iquadpt)
     1                      - dbfx_dy ( ix, iquadpt) )

      fni_quadpt ( ix, 8, iquadpt) =
     1   resistivity_coef * ( dbfz_dx ( ix, iquadpt)
     1                      - dbfx_dz ( ix, iquadpt) )

#if ( NVAR > 8)
      fni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      gni_quadpt ( ix, 1, iquadpt) = 0.0

      gni_quadpt ( ix, 2, iquadpt) =
     1                     tempaa * ( dvly_dx ( ix, iquadpt)
     1                              + dvlx_dy ( ix, iquadpt) )

      gni_quadpt ( ix, 3, iquadpt) =
     1               2.0 * tempaa * ( dvly_dy ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      gni_quadpt ( ix, 4, iquadpt) =
     1                     tempaa * ( dvly_dz ( ix, iquadpt)
     1                              + dvlz_dy ( ix, iquadpt) )

      gni_quadpt ( ix, 5, iquadpt) = 
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vly ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvly_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvly_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvly_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dy ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dy ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfy_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfy_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfy_dz ( ix, iquadpt) )

      gni_quadpt ( ix, 6, iquadpt) =
     1   resistivity_coef * ( dbfx_dy ( ix, iquadpt)
     1                      - dbfy_dx ( ix, iquadpt) )

      gni_quadpt ( ix, 7, iquadpt) = 0.0

      gni_quadpt ( ix, 8, iquadpt) =
     1   resistivity_coef * ( dbfz_dy ( ix, iquadpt)
     1                      - dbfy_dz ( ix, iquadpt) )

#if ( NVAR > 8)
      gni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      hni_quadpt ( ix, 1, iquadpt) = 0.0

      hni_quadpt ( ix, 2, iquadpt) =
     1                     tempaa * ( dvlz_dx ( ix, iquadpt)
     1                              + dvlx_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 3, iquadpt) =
     1                     tempaa * ( dvlz_dy ( ix, iquadpt)
     1                              + dvly_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 4, iquadpt) =
     1               2.0 * tempaa * ( dvlz_dz ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      hni_quadpt ( ix, 5, iquadpt) =
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vlz ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvlz_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvlz_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvlz_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dz ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dz ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfz_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfz_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfz_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 6, iquadpt) =
     1   resistivity_coef * ( dbfx_dz ( ix, iquadpt)
     1                      - dbfz_dx ( ix, iquadpt) )

      hni_quadpt ( ix, 7, iquadpt) =
     1   resistivity_coef * ( dbfy_dz ( ix, iquadpt)
     1                      - dbfz_dy ( ix, iquadpt) )

      hni_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NVAR > 8)
      hni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      f_quadpt ( ix, :, iquadpt) = f_quadpt ( ix, :, iquadpt)
     1                       + fni_quadpt ( ix, :, iquadpt)

      g_quadpt ( ix, :, iquadpt) = g_quadpt ( ix, :, iquadpt)
     1                       + gni_quadpt ( ix, :, iquadpt)

      h_quadpt ( ix, :, iquadpt) = h_quadpt ( ix, :, iquadpt)
     1                       + hni_quadpt ( ix, :, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      END IF
! End of "( viscosity_coef .GT. 0.0) .OR. ( resistivity_coef .GT. 0.0)".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = dtmove / dx ( ix)
      f_quadpt ( ix, :, 1: NUM_ZONECEN_BASES) = tempaa
     1         * f_quadpt ( ix, :, 1: NUM_ZONECEN_BASES)

      tempab = dtmove / dy ( ix)
      g_quadpt ( ix, :, 1: NUM_ZONECEN_BASES) = tempab
     1         * g_quadpt ( ix, :, 1: NUM_ZONECEN_BASES)

      tempac = dtmove / dz ( ix)
      h_quadpt ( ix, :, 1: NUM_ZONECEN_BASES) = tempac
     1         * h_quadpt ( ix, :, 1: NUM_ZONECEN_BASES)

! ----------

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------
! Default is to build source terms at t=0 . In some ADER-CG updates
! (on Cartesian meshes), we may not require any source terms.
#if ( BUILD_SOURCE_TERMS == 1)
! ----------------------------------------------------------------------

! STEP 5.a) Build source terms "s_quadpt" at the spatial quadrature points at t=0.
! Remember that the primitive variables are still available from the previous
! loop at each of the t=0 quadrature points for constructing the source terms.

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = 1, NUM_ZONECEN_BASES

      s_quadpt ( oned_begin - oned_off: oned_end + oned_off,
     1                                           :, iquadpt) = 0.0

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

#if ( 0 == 1)

! For single-fluid adiabatic MHD turbulence, we may want to provide a gradual
! (explicit) heating/cooling term which ensures that pressure always
! tends to a mean pressure "prs_bar" in a time "tau_bar" (comparable
! to a turnover time). Such a source term can be implemented explicitly
! and its implementation is given below.
! Remember to declare "prs_bar, tau_bar".

      prs_bar = 1.0
      tau_bar = 1.0

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 5, iquadpt) = - ( prs ( ix, iquadpt) - prs_bar)
     1     / ( tau_bar * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0) )

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

#endif

! ------------------------------

! Build two-fluid source terms. In this part, the source terms are
! hard-coded for one family of ionized fluid and one family of neutrals.

#if ( N_EULER == 1) 

      IF ( ( igeom .EQ. 1) .AND. ( N_EULER .EQ. 1) .AND.
     1     ( ambipolar_diffusion_drag_coef .GT. 0.0)) THEN

      ivar = NFIELD + 1 + NFLUID

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 1, iquadpt) = 0.0

      tempaa = ambipolar_diffusion_drag_coef
     1       * rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

      s_quadpt ( ix, 2, iquadpt) = tempaa
     1           * ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt))

      s_quadpt ( ix, 3, iquadpt) = tempaa
     1           * ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt))

      s_quadpt ( ix, 4, iquadpt) = tempaa
     1           * ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt))

      s_quadpt ( ix, 5, iquadpt) = ( ambipolar_diffusion_drag_coef
     1   / ( molewt_ion + molewt_neut) )

     1 * ( rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

     1   * ( ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt) )
     1     * ( molewt_neut * vlx_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlx ( ix, iquadpt) )
     1     + ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt) )
     1     * ( molewt_neut * vly_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vly ( ix, iquadpt) )
     1     + ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt) )
     1     * ( molewt_neut * vlz_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlz ( ix, iquadpt) )
     1     )

     1   + 3.0 * ( molewt_neut * prs_euler ( ix, iquadpt, 1)
     1           * rho ( ix, iquadpt)
     1           - molewt_ion * prs ( ix, iquadpt)
     1           * rho_euler ( ix, iquadpt, 1) )
     1   )

      s_quadpt ( ix, 6, iquadpt) = 0.0
      s_quadpt ( ix, 7, iquadpt) = 0.0
      s_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      s_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) = 0.0
#endif


      s_quadpt ( ix, ivar + 1, iquadpt) = 0.0
      s_quadpt ( ix, ivar + 2, iquadpt) = - s_quadpt ( ix, 2, iquadpt)
      s_quadpt ( ix, ivar + 3, iquadpt) = - s_quadpt ( ix, 3, iquadpt)
      s_quadpt ( ix, ivar + 4, iquadpt) = - s_quadpt ( ix, 4, iquadpt)
      s_quadpt ( ix, ivar + 5, iquadpt) = - s_quadpt ( ix, 5, iquadpt)

#if ( NFLUID_EULER >= 1)
      s_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

      END IF
! End of "ambipolar_diffusion_drag_coef .GT. 0.0" part

#endif
! End of "N_EULER == 1"

! ------------------------------

! Cylindrical Geometry

      IF ( igeom .EQ. 2) THEN

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1                           * vlx ( ix, iquadpt) * obx1

      s_quadpt ( ix, 2, iquadpt) =
     1 ( rho ( ix, iquadpt) * vly ( ix, iquadpt)**2
     1 - rho ( ix, iquadpt) * vlx ( ix, iquadpt)**2
     1 + bfx ( ix, iquadpt)**2 * ob4pi
     1 - bfy ( ix, iquadpt)**2 * ob4pi) * obx1

      s_quadpt ( ix, 3, iquadpt) =
     1 - 2.0 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                              * vly ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                        * vlz ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 5, iquadpt) = - ( tempac * vlx ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = 0.0

      s_quadpt ( ix, 7, iquadpt) = 0.0

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Spherical Geometry

      IF ( igeom .EQ. 3) THEN

      DO iquadpt = 1, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      cos_x2 = cos ( x2)
      sin_x2 = sin ( x2)
      cot_x2 = cos_x2 / sin_x2

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 2.0 * vlx ( ix, iquadpt)
     1 + cot_x2 * vly ( ix, iquadpt) ) * obx1

      s_quadpt ( ix, 2, iquadpt) = rho ( ix, iquadpt)
     1 * ( - 2.0 * vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1   + vlz ( ix, iquadpt)**2) * obx1
     1 - cot_x2
     1 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1   - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1
     1 - ( - 2.0 * bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1   + bfz ( ix, iquadpt)**2) * ob4pi * obx1

      s_quadpt ( ix, 3, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1 + cot_x2
     1 * ( vly ( ix, iquadpt)**2 - vlz ( ix, iquadpt)**2)) * obx1
     1 + ( 3.0 * bfx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1   + cot_x2
     1   * ( bfy ( ix, iquadpt)**2 - bfz ( ix, iquadpt)**2))
     1 * ob4pi * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - rho ( ix, iquadpt) * vlz ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) + 2.0 * cot_x2
     1   * vly ( ix, iquadpt) ) * obx1
     1 + bfz ( ix, iquadpt)
     1 * ( 3.0 * bfx ( ix, iquadpt) + 2.0 * cot_x2
     1   * bfy ( ix, iquadpt) ) * ob4pi * obx1

      s_quadpt ( ix, 5, iquadpt) =
     1 - 2.0 * ( tempac * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab) * obx1
     1 - cot_x2 * ( tempac * vly ( ix, iquadpt)
     1 - bfy ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = - cot_x2
     1                           * g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 7, iquadpt) = g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Multiply "s_quadpt" by "dtmove".

      DO ix = oned_begin - oned_off, oned_end + oned_off
      s_quadpt ( ix, :, 1: NUM_ZONECEN_BASES) = dtmove
     1         * s_quadpt ( ix, :, 1: NUM_ZONECEN_BASES)
      END DO

! ----------------------------------------------------------------------
#endif
! End of "BUILD_SOURCE_TERMS == 1"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 6) Make the contributions from the t=0 level to the evolutionary
! equations and lump them all together here in "tu". This ensures that the
! calculations in this step are not repeated for each of the ADER iterations.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : Build "tu" at t=0

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tu ( ix, :, 01) = u_quadpt ( ix, :, 01)
     1 - 2.0 * ( f_quadpt ( ix, :, 02) - f_quadpt ( ix, :, 01) )
     1 - 2.0 * ( g_quadpt ( ix, :, 03) - g_quadpt ( ix, :, 01) )
     1 - 2.0 * ( h_quadpt ( ix, :, 04) - h_quadpt ( ix, :, 01) )
     1 + s_quadpt ( ix, :, 01) / 3.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : Build "tu" at t=0

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tu ( ix, :, 11) = u_quadpt ( ix, :, 01)
     1 + 0.5 * ( f_quadpt ( ix, :, 03) - f_quadpt ( ix, :, 01) )
     1 + 0.5 * ( g_quadpt ( ix, :, 05) - g_quadpt ( ix, :, 01) )
     1 + 0.5 * ( h_quadpt ( ix, :, 07) - h_quadpt ( ix, :, 01) )
     1 + ( - 9.0 * s_quadpt ( ix, :, 01) + 2.0 * s_quadpt ( ix, :, 02)
     1   + 2.0 * s_quadpt ( ix, :, 03) + 2.0 * s_quadpt ( ix, :, 04)
     1   + 2.0 * s_quadpt ( ix, :, 05) + 2.0 * s_quadpt ( ix, :, 06)
     1   + 2.0 * s_quadpt ( ix, :, 07) ) / 24.0

      tu ( ix, :, 12) = u_quadpt ( ix, :, 02)
     1 + ( 36.0 * f_quadpt ( ix, :, 01) - 24.0 * f_quadpt ( ix, :, 02)
     1 - 12.0 * f_quadpt ( ix, :, 03) - 36.0 * g_quadpt ( ix, :, 01)
     1 + 24.0 * g_quadpt ( ix, :, 02) + 24.0 * g_quadpt ( ix, :, 04)
     1 + 12.0 * g_quadpt ( ix, :, 05) - 24.0 * g_quadpt ( ix, :, 10)
     1 - 36.0 * h_quadpt ( ix, :, 01) + 24.0 * h_quadpt ( ix, :, 02)
     1 + 24 * h_quadpt ( ix, :, 06) + 12.0 * h_quadpt ( ix, :, 07)
     1 - 24.0 * h_quadpt ( ix, :, 09) + 7.0 * s_quadpt ( ix, :, 01)
     1 - 8.0 * s_quadpt ( ix, :, 02) - 4.0 * s_quadpt ( ix, :, 03)
     1 + 2.0 * s_quadpt ( ix, :, 04) + 2.0 * s_quadpt ( ix, :, 05)
     1 + 2.0 * s_quadpt ( ix, :, 06) + 2.0 * s_quadpt ( ix, :, 07)
     1   ) / 24.0

      tu ( ix, :, 13) = u_quadpt ( ix, :, 04)
     1 + ( -36.0 * f_quadpt ( ix, :, 01) + 24.0 * f_quadpt ( ix, :, 02)
     1 + 12.0 * f_quadpt ( ix, :, 03) + 24.0 * f_quadpt ( ix, :, 04)
     1 - 24.0 * f_quadpt ( ix, :, 10) + 36.0 * g_quadpt ( ix, :, 01)
     1 - 24.0 * g_quadpt ( ix, :, 04) - 12.0 * g_quadpt ( ix, :, 05)
     1 - 36.0 * h_quadpt ( ix, :, 01) + 24.0 * h_quadpt ( ix, :, 04)
     1 + 24.0 * h_quadpt ( ix, :, 06) + 12.0 * h_quadpt ( ix, :, 07)
     1 - 24.0 * h_quadpt ( ix, :, 08) + 7.0 * s_quadpt ( ix, :, 01)
     1 + 2.0 * s_quadpt ( ix, :, 02) + 2.0 * s_quadpt ( ix, :, 03)
     1 - 8.0 * s_quadpt ( ix, :, 04) - 4.0 * s_quadpt ( ix, :, 05)
     1 + 2.0 * s_quadpt ( ix, :, 06) + 2.0 * s_quadpt ( ix, :, 07)
     1   ) / 24.0

      tu ( ix, :, 14) = u_quadpt ( ix, :, 06)
     1 + ( -36.0 * f_quadpt ( ix, :, 01) + 24.0 * f_quadpt ( ix, :, 02)
     1 + 12.0 * f_quadpt ( ix, :, 03) + 24.0 * f_quadpt ( ix, :, 06)
     1 - 24.0 * f_quadpt ( ix, :, 09) - 36.0 * g_quadpt ( ix, :, 01)
     1 + 24.0 * g_quadpt ( ix, :, 04) + 12.0 * g_quadpt ( ix, :, 05)
     1 + 24.0 * g_quadpt ( ix, :, 06) - 24.0 * g_quadpt ( ix, :, 08)
     1 + 36.0 * h_quadpt ( ix, :, 01) - 24.0 * h_quadpt ( ix, :, 06)
     1 - 12.0 * h_quadpt ( ix, :, 07) + 7.0 * s_quadpt ( ix, :, 01)
     1 + 2.0 * s_quadpt ( ix, :, 02) + 2.0 * s_quadpt ( ix, :, 03)
     1 + 2.0 * s_quadpt ( ix, :, 04) + 2.0 * s_quadpt ( ix, :, 05)
     1 - 8.0 * s_quadpt ( ix, :, 06) - 4.0 * s_quadpt ( ix, :, 07)
     1   ) / 24.0

      tu ( ix, :, 15) = u_quadpt ( ix, :, 01)
     1 + ( -60.0 * f_quadpt ( ix, :, 01) + 30.0 * f_quadpt ( ix, :, 02)
     1 + 30.0 * f_quadpt ( ix, :, 03) - 60.0 * g_quadpt ( ix, :, 01)
     1 + 30.0 * g_quadpt ( ix, :, 04) + 30.0 * g_quadpt ( ix, :, 05)
     1 - 60.0 * h_quadpt ( ix, :, 01) + 30.0 * h_quadpt ( ix, :, 06)
     1 + 30.0 * h_quadpt ( ix, :, 07) - 27.0 * s_quadpt ( ix, :, 01)
     1 + 5.0 * s_quadpt ( ix, :, 02) + 5.0 * s_quadpt ( ix, :, 03)
     1 + 5.0 * s_quadpt ( ix, :, 04) + 5.0 * s_quadpt ( ix, :, 05)
     1 + 5.0 * s_quadpt ( ix, :, 06) + 5.0 * s_quadpt ( ix, :, 07)
     1   ) / 30.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : Build "tu" at t=0

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tu ( ix, :, 21) = u_quadpt ( ix, :, 01)

     1 + ( - 650370 * f_quadpt ( ix, :, 01)
     1 + 191520 * f_quadpt ( ix, :, 02)
     1 + 305340 * f_quadpt ( ix, :, 03)
     1 + 295750 * f_quadpt ( ix, :, 04)
     1 + 238560 * f_quadpt ( ix, :, 05)
     1 + 28000 * f_quadpt ( ix, :, 06)
     1 + 389760 * f_quadpt ( ix, :, 07)
     1 + 159040 * f_quadpt ( ix, :, 08)
     1 - 203840 * f_quadpt ( ix, :, 09)
     1 - 80640 * f_quadpt ( ix, :, 10)
     1 - 284480 * f_quadpt ( ix, :, 11)
     1 + 70560 * f_quadpt ( ix, :, 12)
     1 - 213920 * f_quadpt ( ix, :, 13)
     1 + 10080 * f_quadpt ( ix, :, 14)
     1 + 202720 * f_quadpt ( ix, :, 15)
     1 - 180320 * f_quadpt ( ix, :, 16)
     1 - 89600 * f_quadpt ( ix, :, 17)
     1 - 89600 * f_quadpt ( ix, :, 18)
     1 - 49280 * f_quadpt ( ix, :, 19)
     1 - 49280 * f_quadpt ( ix, :, 20)

     1 - 650370 * g_quadpt ( ix, :, 01)
     1 + 191520 * g_quadpt ( ix, :, 02)
     1 + 389760 * g_quadpt ( ix, :, 03)
     1 + 159040 * g_quadpt ( ix, :, 04)
     1 + 305340 * g_quadpt ( ix, :, 05)
     1 + 295750 * g_quadpt ( ix, :, 06)
     1 + 238560 * g_quadpt ( ix, :, 07)
     1 + 28000 * g_quadpt ( ix, :, 08)
     1 - 213920 * g_quadpt ( ix, :, 09)
     1 + 10080 * g_quadpt ( ix, :, 10)
     1 - 203840 * g_quadpt ( ix, :, 11)
     1 - 80640 * g_quadpt ( ix, :, 12)
     1 - 284480 * g_quadpt ( ix, :, 13)
     1 + 70560 * g_quadpt ( ix, :, 14)
     1 - 49280 * g_quadpt ( ix, :, 15)
     1 - 49280 * g_quadpt ( ix, :, 16)
     1 + 202720 * g_quadpt ( ix, :, 17)
     1 - 180320 * g_quadpt ( ix, :, 18)
     1 - 89600 * g_quadpt ( ix, :, 19)
     1 - 89600 * g_quadpt ( ix, :, 20)

     1 - 650370 * h_quadpt ( ix, :, 01)
     1 + 191520 * h_quadpt ( ix, :, 02)
     1 + 238560 * h_quadpt ( ix, :, 03)
     1 + 28000 * h_quadpt ( ix, :, 04)
     1 + 389760 * h_quadpt ( ix, :, 05)
     1 + 159040 * h_quadpt ( ix, :, 06)
     1 + 305340 * h_quadpt ( ix, :, 07)
     1 + 295750 * h_quadpt ( ix, :, 08)
     1 - 284480 * h_quadpt ( ix, :, 09)
     1 + 70560 * h_quadpt ( ix, :, 10)
     1 - 213920 * h_quadpt ( ix, :, 11)
     1 + 10080 * h_quadpt ( ix, :, 12)
     1 - 203840 * h_quadpt ( ix, :, 13)
     1 - 80640 * h_quadpt ( ix, :, 14)
     1 - 89600 * h_quadpt ( ix, :, 15)
     1 - 89600 * h_quadpt ( ix, :, 16)
     1 - 49280 * h_quadpt ( ix, :, 17)
     1 - 49280 * h_quadpt ( ix, :, 18)
     1 + 202720 * h_quadpt ( ix, :, 19)
     1 - 180320 * h_quadpt ( ix, :, 20)

     1 - 89091 * s_quadpt ( ix, :, 01)
     1 + 17955 * s_quadpt ( ix, :, 03)
     1 + 17955 * s_quadpt ( ix, :, 04)
     1 + 17955 * s_quadpt ( ix, :, 05)
     1 + 17955 * s_quadpt ( ix, :, 06)
     1 + 17955 * s_quadpt ( ix, :, 07)
     1 + 17955 * s_quadpt ( ix, :, 08)  ) / 215460.0

      tu ( ix, :, 22) = u_quadpt ( ix, :, 03) 

     1 + ( 4225410 * f_quadpt ( ix, :, 01) 
     1 - 1149120 * f_quadpt ( ix, :, 02)
     1 - 2262960 * f_quadpt ( ix, :, 03)
     1 - 1666770 * f_quadpt ( ix, :, 04) 
     1 - 1431360 * f_quadpt ( ix, :, 05) 
     1 - 168000 * f_quadpt ( ix, :, 06) 
     1 - 2338560 * f_quadpt ( ix, :, 07) 
     1 - 954240 * f_quadpt ( ix, :, 08) 
     1 + 1223040 * f_quadpt ( ix, :, 09)
     1 + 483840 * f_quadpt ( ix, :, 10) 
     1 + 1706880 * f_quadpt ( ix, :, 11) 
     1 - 423360 * f_quadpt ( ix, :, 12)
     1 + 1283520 * f_quadpt ( ix, :, 13) 
     1 - 60480 * f_quadpt ( ix, :, 14) 
     1 - 1216320 * f_quadpt ( ix, :, 15)
     1 + 1081920 * f_quadpt ( ix, :, 16)
     1 + 537600 * f_quadpt ( ix, :, 17)
     1 + 537600 * f_quadpt ( ix, :, 18)
     1 + 295680 * f_quadpt ( ix, :, 19)
     1 + 295680 * f_quadpt ( ix, :, 20)

     1 - 1735650 * g_quadpt ( ix, :, 01)
     1 + 574560 * g_quadpt ( ix, :, 02)
     1 + 953820 * g_quadpt ( ix, :, 03)
     1 + 477120 * g_quadpt ( ix, :, 04)
     1 + 916020 * g_quadpt ( ix, :, 05)
     1 + 671790 * g_quadpt ( ix, :, 06)
     1 + 715680 * g_quadpt ( ix, :, 07)
     1 + 84000 * g_quadpt ( ix, :, 08)
     1 - 641760 * g_quadpt ( ix, :, 09)
     1 + 30240 * g_quadpt ( ix, :, 10)
     1 - 611520 * g_quadpt ( ix, :, 11)
     1 - 26460 * g_quadpt ( ix, :, 12)
     1 - 853440 * g_quadpt ( ix, :, 13)
     1 + 211680 * g_quadpt ( ix, :, 14)
     1 - 147840 * g_quadpt ( ix, :, 15)
     1 - 147840 * g_quadpt ( ix, :, 16)
     1 + 608160 * g_quadpt ( ix, :, 17)
     1 - 540960 * g_quadpt ( ix, :, 18)
     1 - 268800 * g_quadpt ( ix, :, 19)
     1 - 268800 * g_quadpt ( ix, :, 20)

     1 - 2645370 * h_quadpt ( ix, :, 01)
     1 + 574560 * h_quadpt ( ix, :, 02)
     1 + 981540 * h_quadpt ( ix, :, 03)
     1 + 510720 * h_quadpt ( ix, :, 04)
     1 + 1197000 * h_quadpt ( ix, :, 05)
     1 + 510720 * h_quadpt ( ix, :, 06)
     1 + 1340640 * h_quadpt ( ix, :, 07)
     1 + 1049370 * h_quadpt ( ix, :, 08)
     1 - 853860 * h_quadpt ( ix, :, 09)
     1 - 853860 * h_quadpt ( ix, :, 11)
     1 - 638400 * h_quadpt ( ix, :, 13)
     1 - 215460 * h_quadpt ( ix, :, 14)
     1 - 255360 * h_quadpt ( ix, :, 15)
     1 - 255360 * h_quadpt ( ix, :, 16)
     1 - 255360 * h_quadpt ( ix, :, 17)
     1 - 255360 * h_quadpt ( ix, :, 18) 
     1 + 606480 * h_quadpt ( ix, :, 19)
     1 - 542640 * h_quadpt ( ix, :, 20)

     1 + 474069 * s_quadpt ( ix, :, 01)
     1 - 114912 * s_quadpt ( ix, :, 02)
     1 - 411810 * s_quadpt ( ix, :, 03)
     1 - 192514 * s_quadpt ( ix, :, 04)
     1 - 120561 * s_quadpt ( ix, :, 05)
     1 + 6755 * s_quadpt ( ix, :, 06)
     1 - 145131 * s_quadpt ( ix, :, 07)
     1 - 14539 * s_quadpt ( ix, :, 08)
     1 + 122234 * s_quadpt ( ix, :, 09)
     1 + 13104 * s_quadpt ( ix, :, 10)
     1 + 171248 * s_quadpt ( ix, :, 11)
     1 - 11466 * s_quadpt ( ix, :, 12)
     1 + 159782 * s_quadpt ( ix, :, 13)
     1 - 1638 * s_quadpt ( ix, :, 14)
     1 - 119392 * s_quadpt ( ix, :, 15)
     1 + 110432 * s_quadpt ( ix, :, 16)
     1 + 35840 * s_quadpt ( ix, :, 17)
     1 + 35840 * s_quadpt ( ix, :, 18)
     1 + 29288 * s_quadpt ( ix, :, 19)
     1 + 29288 * s_quadpt ( ix, :, 20)  ) / 646380.0

      tu ( ix, :, 23) = u_quadpt ( ix, :, 04)

     1 + ( 3363570 * f_quadpt ( ix, :, 01)
     1 - 1149120 * f_quadpt ( ix, :, 02)
     1 - 1832040 * f_quadpt ( ix, :, 03)
     1 - 1235850 * f_quadpt ( ix, :, 04)
     1 - 1431360 * f_quadpt ( ix, :, 05)
     1 - 168000 * f_quadpt ( ix, :, 06)
     1 - 2338560 * f_quadpt ( ix, :, 07)
     1 - 954240 * f_quadpt ( ix, :, 08)
     1 + 1223040 * f_quadpt ( ix, :, 09)
     1 + 483840 * f_quadpt ( ix, :, 10)
     1 + 1706880 * f_quadpt ( ix, :, 11)
     1 - 423360 * f_quadpt ( ix, :, 12)
     1 + 1283520 * f_quadpt ( ix, :, 13)
     1 - 60480 * f_quadpt ( ix, :, 14)
     1 - 1216320 * f_quadpt ( ix, :, 15)
     1 + 1081920 * f_quadpt ( ix, :, 16)
     1 + 537600 * f_quadpt ( ix, :, 17)
     1 + 537600 * f_quadpt ( ix, :, 18)
     1 + 295680 * f_quadpt ( ix, :, 19)
     1 + 295680 * f_quadpt ( ix, :, 20)

     1 - 1208970 * g_quadpt ( ix, :, 01)
     1 + 574560 * g_quadpt ( ix, :, 02)
     1 + 966420 * g_quadpt ( ix, :, 03)
     1 + 152880 * g_quadpt ( ix, :, 04)
     1 + 815220 * g_quadpt ( ix, :, 05)
     1 + 249270 * g_quadpt ( ix, :, 06)
     1 + 660240 * g_quadpt ( ix, :, 07)
     1 + 16800 * g_quadpt ( ix, :, 08)
     1 - 588000 * g_quadpt ( ix, :, 09)
     1 - 22680 * g_quadpt ( ix, :, 10)
     1 - 610680 * g_quadpt ( ix, :, 11)
     1 - 34020 * g_quadpt ( ix, :, 12)
     1 - 860160 * g_quadpt ( ix, :, 13)
     1 + 272160 * g_quadpt ( ix, :, 14)
     1 - 144480 * g_quadpt ( ix, :, 15)
     1 - 144480 * g_quadpt ( ix, :, 16)
     1 + 581280 * g_quadpt ( ix, :, 17)
     1 - 567840 * g_quadpt ( ix, :, 18)
     1 - 53760 * g_quadpt ( ix, :, 19)
     1 - 53760 * g_quadpt ( ix, :, 20)

     1 - 2118690 * h_quadpt ( ix, :, 01)
     1 + 574560 * h_quadpt ( ix, :, 02)
     1 + 880740 * h_quadpt ( ix, :, 03)
     1 + 88200 * h_quadpt ( ix, :, 04)
     1 + 1141560 * h_quadpt ( ix, :, 05)
     1 + 443520 * h_quadpt ( ix, :, 06)
     1 + 1353240 * h_quadpt ( ix, :, 07)
     1 + 725130 * h_quadpt ( ix, :, 08)
     1 - 853020 * h_quadpt ( ix, :, 09)
     1 - 7560 * h_quadpt ( ix, :, 10)
     1 - 860580 * h_quadpt ( ix, :, 11)
     1 + 60480 * h_quadpt ( ix, :, 12)
     1 - 584640 * h_quadpt ( ix, :, 13)
     1 - 268380 * h_quadpt ( ix, :, 14)
     1 - 282240 * h_quadpt ( ix, :, 15)
     1 - 282240 * h_quadpt ( ix, :, 16)
     1 - 40320 * h_quadpt ( ix, :, 17)
     1 - 40320 * h_quadpt ( ix, :, 18)
     1 + 609840 * h_quadpt ( ix, :, 19)
     1 - 539280 * h_quadpt ( ix, :, 20)

     1 - 577695 * s_quadpt ( ix, :, 01)
     1 + 114912 * s_quadpt ( ix, :, 02)
     1 + 304080 * s_quadpt ( ix, :, 03)
     1 + 84784 * s_quadpt ( ix, :, 04)
     1 + 228291 * s_quadpt ( ix, :, 05)
     1 + 100975 * s_quadpt ( ix, :, 06)
     1 + 252861 * s_quadpt ( ix, :, 07)
     1 + 122269 * s_quadpt ( ix, :, 08)
     1 - 122234 * s_quadpt ( ix, :, 09)
     1 - 13104 * s_quadpt ( ix, :, 10)
     1 - 171248 * s_quadpt ( ix, :, 11)
     1 + 11466 * s_quadpt ( ix, :, 12)
     1 - 159782 * s_quadpt ( ix, :, 13)
     1 + 1638 * s_quadpt ( ix, :, 14)
     1 + 119392 * s_quadpt ( ix, :, 15)
     1 - 110432 * s_quadpt ( ix, :, 16)
     1 - 35840 * s_quadpt ( ix, :, 17)
     1 - 35840 * s_quadpt ( ix, :, 18)
     1 - 29288 * s_quadpt ( ix, :, 19)
     1 - 29288 * s_quadpt ( ix, :, 20)  ) / 646380.0

      tu ( ix, :, 24) = u_quadpt ( ix, :, 05)

     1 + ( - 2645370 * f_quadpt ( ix, :, 01)
     1 + 574560 * f_quadpt ( ix, :, 02)
     1 + 1340640 * f_quadpt ( ix, :, 03)
     1 + 1049370 * f_quadpt ( ix, :, 04)
     1 + 981540 * f_quadpt ( ix, :, 05)
     1 + 510720 * f_quadpt ( ix, :, 06)
     1 + 1197000 * f_quadpt ( ix, :, 07)
     1 + 510720 * f_quadpt ( ix, :, 08)
     1 - 638400 * f_quadpt ( ix, :, 09)
     1 - 215460 * f_quadpt ( ix, :, 10)
     1 - 853860 * f_quadpt ( ix, :, 11)
     1 - 853860 * f_quadpt ( ix, :, 13)
     1 + 606480 * f_quadpt ( ix, :, 15)
     1 - 542640 * f_quadpt ( ix, :, 16)
     1 - 255360 * f_quadpt ( ix, :, 17)
     1 - 255360 * f_quadpt ( ix, :, 18)
     1 - 255360 * f_quadpt ( ix, :, 19)
     1 - 255360 * f_quadpt ( ix, :, 20)

     1 + 4225410 * g_quadpt ( ix, :, 01)
     1 - 1149120 * g_quadpt ( ix, :, 02)
     1 - 2338560 * g_quadpt ( ix, :, 03)
     1 - 954240 * g_quadpt ( ix, :, 04)
     1 - 2262960 * g_quadpt ( ix, :, 05)
     1 - 1666770 * g_quadpt ( ix, :, 06)
     1 - 1431360 * g_quadpt ( ix, :, 07)
     1 - 168000 * g_quadpt ( ix, :, 08)
     1 + 1283520 * g_quadpt ( ix, :, 09)
     1 - 60480 * g_quadpt ( ix, :, 10)
     1 + 1223040 * g_quadpt ( ix, :, 11)
     1 + 483840 * g_quadpt ( ix, :, 12)
     1 + 1706880 * g_quadpt ( ix, :, 13)
     1 - 423360 * g_quadpt ( ix, :, 14)
     1 + 295680 * g_quadpt ( ix, :, 15)
     1 + 295680 * g_quadpt ( ix, :, 16)
     1 - 1216320 * g_quadpt ( ix, :, 17)
     1 + 1081920 * g_quadpt ( ix, :, 18)
     1 + 537600 * g_quadpt ( ix, :, 19)
     1 + 537600 * g_quadpt ( ix, :, 20)

     1 - 1735650 * h_quadpt ( ix, :, 01)
     1 + 574560 * h_quadpt ( ix, :, 02)
     1 + 715680 * h_quadpt ( ix, :, 03)
     1 + 84000 * h_quadpt ( ix, :, 04)
     1 + 953820 * h_quadpt ( ix, :, 05)
     1 + 477120 * h_quadpt ( ix, :, 06)
     1 + 916020 * h_quadpt ( ix, :, 07)
     1 + 671790 * h_quadpt ( ix, :, 08)
     1 - 853440 * h_quadpt ( ix, :, 09)
     1 + 211680 * h_quadpt ( ix, :, 10)
     1 - 641760 * h_quadpt ( ix, :, 11)
     1 + 30240 * h_quadpt ( ix, :, 12)
     1 - 611520 * h_quadpt ( ix, :, 13)
     1 - 26460 * h_quadpt ( ix, :, 14)
     1 - 268800 * h_quadpt ( ix, :, 15)
     1 - 268800 * h_quadpt ( ix, :, 16)
     1 - 147840 * h_quadpt ( ix, :, 17)
     1 - 147840 * h_quadpt ( ix, :, 18)
     1 + 608160 * h_quadpt ( ix, :, 19)
     1 - 540960 * h_quadpt ( ix, :, 20)

     1 + 474069 * s_quadpt ( ix, :, 01)
     1 - 114912 * s_quadpt ( ix, :, 02)
     1 - 145131 * s_quadpt ( ix, :, 03)
     1 - 14539 * s_quadpt ( ix, :, 04)
     1 - 411810 * s_quadpt ( ix, :, 05)
     1 - 192514 * s_quadpt ( ix, :, 06)
     1 - 120561 * s_quadpt ( ix, :, 07)
     1 + 6755 * s_quadpt ( ix, :, 08)
     1 + 159782 * s_quadpt ( ix, :, 09)
     1 - 1638 * s_quadpt ( ix, :, 10)
     1 + 122234 * s_quadpt ( ix, :, 11)
     1 + 13104 * s_quadpt ( ix, :, 12)
     1 + 171248 * s_quadpt ( ix, :, 13)
     1 - 11466 * s_quadpt ( ix, :, 14)
     1 + 29288 * s_quadpt ( ix, :, 15)
     1 + 29288 * s_quadpt ( ix, :, 16)
     1 - 119392 * s_quadpt ( ix, :, 17)
     1 + 110432 * s_quadpt ( ix, :, 18)
     1 + 35840 * s_quadpt ( ix, :, 19)
     1 + 35840 * s_quadpt ( ix, :, 20) ) / 646380.0

      tu ( ix, :, 25) = u_quadpt ( ix, :, 06)

     1 + ( - 2118690 * f_quadpt ( ix, :, 01)
     1 + 574560 * f_quadpt ( ix, :, 02)
     1 + 1353240 * f_quadpt ( ix, :, 03)
     1 + 725130 * f_quadpt ( ix, :, 04)
     1 + 880740 * f_quadpt ( ix, :, 05)
     1 + 88200 * f_quadpt ( ix, :, 06)
     1 + 1141560 * f_quadpt ( ix, :, 07)
     1 + 443520 * f_quadpt ( ix, :, 08)
     1 - 584640 * f_quadpt ( ix, :, 09)
     1 - 268380 * f_quadpt ( ix, :, 10)
     1 - 853020 * f_quadpt ( ix, :, 11)
     1 - 7560 * f_quadpt ( ix, :, 12)
     1 - 860580 * f_quadpt ( ix, :, 13)
     1 + 60480 * f_quadpt ( ix, :, 14)
     1 + 609840 * f_quadpt ( ix, :, 15)
     1 - 539280 * f_quadpt ( ix, :, 16)
     1 - 282240 * f_quadpt ( ix, :, 17)
     1 - 282240 * f_quadpt ( ix, :, 18)
     1 - 40320 * f_quadpt ( ix, :, 19)
     1 - 40320 * f_quadpt ( ix, :, 20)

     1 + 3363570 * g_quadpt ( ix, :, 01)
     1 - 1149120 * g_quadpt ( ix, :, 02)
     1 - 2338560 * g_quadpt ( ix, :, 03)
     1 - 954240 * g_quadpt ( ix, :, 04)
     1 - 1832040 * g_quadpt ( ix, :, 05)
     1 - 1235850 * g_quadpt ( ix, :, 06)
     1 - 1431360 * g_quadpt ( ix, :, 07)
     1 - 168000 * g_quadpt ( ix, :, 08)
     1 + 1283520 * g_quadpt ( ix, :, 09)
     1 - 60480 * g_quadpt ( ix, :, 10)
     1 + 1223040 * g_quadpt ( ix, :, 11)
     1 + 483840 * g_quadpt ( ix, :, 12)
     1 + 1706880 * g_quadpt ( ix, :, 13)
     1 - 423360 * g_quadpt ( ix, :, 14)
     1 + 295680 * g_quadpt ( ix, :, 15)
     1 + 295680 * g_quadpt ( ix, :, 16)
     1 - 1216320 * g_quadpt ( ix, :, 17)
     1 + 1081920 * g_quadpt ( ix, :, 18)
     1 + 537600 * g_quadpt ( ix, :, 19)
     1 + 537600 * g_quadpt ( ix, :, 20)

     1 - 1208970 * h_quadpt ( ix, :, 01)
     1 + 574560 * h_quadpt ( ix, :, 02)
     1 + 660240 * h_quadpt ( ix, :, 03)
     1 + 16800 * h_quadpt ( ix, :, 04)
     1 + 966420 * h_quadpt ( ix, :, 05)
     1 + 152880 * h_quadpt ( ix, :, 06)
     1 + 815220 * h_quadpt ( ix, :, 07)
     1 + 249270 * h_quadpt ( ix, :, 08)
     1 - 860160 * h_quadpt ( ix, :, 09)
     1 + 272160 * h_quadpt ( ix, :, 10)
     1 - 588000 * h_quadpt ( ix, :, 11)
     1 - 22680 * h_quadpt ( ix, :, 12)
     1 - 610680 * h_quadpt ( ix, :, 13)
     1 - 34020 * h_quadpt ( ix, :, 14)
     1 - 53760 * h_quadpt ( ix, :, 15)
     1 - 53760 * h_quadpt ( ix, :, 16)
     1 - 144480 * h_quadpt ( ix, :, 17)
     1 - 144480 * h_quadpt ( ix, :, 18)
     1 + 581280 * h_quadpt ( ix, :, 19)
     1 - 567840 * h_quadpt ( ix, :, 20)

     1 - 577695 * s_quadpt ( ix, :, 01)
     1 + 114912 * s_quadpt ( ix, :, 02)
     1 + 252861 * s_quadpt ( ix, :, 03)
     1 + 122269 * s_quadpt ( ix, :, 04)
     1 + 304080 * s_quadpt ( ix, :, 05)
     1 + 84784 * s_quadpt ( ix, :, 06)
     1 + 228291 * s_quadpt ( ix, :, 07)
     1 + 100975 * s_quadpt ( ix, :, 08)
     1 - 159782 * s_quadpt ( ix, :, 09)
     1 + 1638 * s_quadpt ( ix, :, 10)
     1 - 122234 * s_quadpt ( ix, :, 11)
     1 - 13104 * s_quadpt ( ix, :, 12)
     1 - 171248 * s_quadpt ( ix, :, 13)
     1 + 11466 * s_quadpt ( ix, :, 14)
     1 - 29288 * s_quadpt ( ix, :, 15)
     1 - 29288 * s_quadpt ( ix, :, 16)
     1 + 119392 * s_quadpt ( ix, :, 17)
     1 - 110432 * s_quadpt ( ix, :, 18)
     1 - 35840 * s_quadpt ( ix, :, 19)
     1 - 35840 * s_quadpt ( ix, :, 20) ) / 646380.0

      tu ( ix, :, 26) = u_quadpt ( ix, :, 07)
     1 + ( - 1735650 * f_quadpt ( ix, :, 01)
     1 + 574560 * f_quadpt ( ix, :, 02)
     1 + 916020 * f_quadpt ( ix, :, 03)
     1 + 671790 * f_quadpt ( ix, :, 04)
     1 + 715680 * f_quadpt ( ix, :, 05)
     1 + 84000 * f_quadpt ( ix, :, 06)
     1 + 953820 * f_quadpt ( ix, :, 07)
     1 + 477120 * f_quadpt ( ix, :, 08)
     1 - 611520 * f_quadpt ( ix, :, 09)
     1 - 26460 * f_quadpt ( ix, :, 10)
     1 - 853440 * f_quadpt ( ix, :, 11)
     1 + 211680 * f_quadpt ( ix, :, 12)
     1 - 641760 * f_quadpt ( ix, :, 13)
     1 + 30240 * f_quadpt ( ix, :, 14)
     1 + 608160 * f_quadpt ( ix, :, 15)
     1 - 540960 * f_quadpt ( ix, :, 16)
     1 - 268800 * f_quadpt ( ix, :, 17)
     1 - 268800 * f_quadpt ( ix, :, 18)
     1 - 147840 * f_quadpt ( ix, :, 19)
     1 - 147840 * f_quadpt ( ix, :, 20)

     1 - 2645370 * g_quadpt ( ix, :, 01)
     1 + 574560 * g_quadpt ( ix, :, 02)
     1 + 1197000 * g_quadpt ( ix, :, 03)
     1 + 510720 * g_quadpt ( ix, :, 04)
     1 + 1340640 * g_quadpt ( ix, :, 05)
     1 + 1049370 * g_quadpt ( ix, :, 06)
     1 + 981540 * g_quadpt ( ix, :, 07)
     1 + 510720 * g_quadpt ( ix, :, 08)
     1 - 853860 * g_quadpt ( ix, :, 09)
     1 - 638400 * g_quadpt ( ix, :, 11)
     1 - 215460 * g_quadpt ( ix, :, 12)
     1 - 853860 * g_quadpt ( ix, :, 13)
     1 - 255360 * g_quadpt ( ix, :, 15)
     1 - 255360 * g_quadpt ( ix, :, 16)
     1 + 606480 * g_quadpt ( ix, :, 17)
     1 - 542640 * g_quadpt ( ix, :, 18)
     1 - 255360 * g_quadpt ( ix, :, 19)
     1 - 255360 * g_quadpt ( ix, :, 20)

     1 + 4225410 * h_quadpt ( ix, :, 01)
     1 - 1149120 * h_quadpt ( ix, :, 02)
     1 - 1431360 * h_quadpt ( ix, :, 03)
     1 - 168000 * h_quadpt ( ix, :, 04)
     1 - 2338560 * h_quadpt ( ix, :, 05)
     1 - 954240 * h_quadpt ( ix, :, 06)
     1 - 2262960 * h_quadpt ( ix, :, 07)
     1 - 1666770 * h_quadpt ( ix, :, 08)
     1 + 1706880 * h_quadpt ( ix, :, 09)
     1 - 423360 * h_quadpt ( ix, :, 10)
     1 + 1283520 * h_quadpt ( ix, :, 11)
     1 - 60480 * h_quadpt ( ix, :, 12)
     1 + 1223040 * h_quadpt ( ix, :, 13)
     1 + 483840 * h_quadpt ( ix, :, 14)
     1 + 537600 * h_quadpt ( ix, :, 15)
     1 + 537600 * h_quadpt ( ix, :, 16)
     1 + 295680 * h_quadpt ( ix, :, 17)
     1 + 295680 * h_quadpt ( ix, :, 18)
     1 - 1216320 * h_quadpt ( ix, :, 19)
     1 + 1081920 * h_quadpt ( ix, :, 20)

     1 + 474069 * s_quadpt ( ix, :, 01)
     1 - 114912 * s_quadpt ( ix, :, 02)
     1 - 120561 * s_quadpt ( ix, :, 03)
     1 + 6755 * s_quadpt ( ix, :, 04)
     1 - 145131 * s_quadpt ( ix, :, 05)
     1 - 14539 * s_quadpt ( ix, :, 06)
     1 - 411810 * s_quadpt ( ix, :, 07)
     1 - 192514 * s_quadpt ( ix, :, 08)
     1 + 171248 * s_quadpt ( ix, :, 09)
     1 - 11466 * s_quadpt ( ix, :, 10)
     1 + 159782 * s_quadpt ( ix, :, 11)
     1 - 1638 * s_quadpt ( ix, :, 12)
     1 + 122234 * s_quadpt ( ix, :, 13)
     1 + 13104 * s_quadpt ( ix, :, 14)
     1 + 35840 * s_quadpt ( ix, :, 15)
     1 + 35840 * s_quadpt ( ix, :, 16)
     1 + 29288 * s_quadpt ( ix, :, 17)
     1 + 29288 * s_quadpt ( ix, :, 18)
     1 - 119392 * s_quadpt ( ix, :, 19)
     1 + 110432 * s_quadpt ( ix, :, 20) ) / 646380.0

      tu ( ix, :, 27) = u_quadpt ( ix, :, 08)

     1 + ( - 1208970 * f_quadpt ( ix, :, 01)
     1 + 574560 * f_quadpt ( ix, :, 02)
     1 + 815220 * f_quadpt ( ix, :, 03)
     1 + 249270 * f_quadpt ( ix, :, 04)
     1 + 660240 * f_quadpt ( ix, :, 05)
     1 + 16800 * f_quadpt ( ix, :, 06)
     1 + 966420 * f_quadpt ( ix, :, 07)
     1 + 152880 * f_quadpt ( ix, :, 08)
     1 - 610680 * f_quadpt ( ix, :, 09)
     1 - 34020 * f_quadpt ( ix, :, 10)
     1 - 860160 * f_quadpt ( ix, :, 11)
     1 + 272160 * f_quadpt ( ix, :, 12)
     1 - 588000 * f_quadpt ( ix, :, 13)
     1 - 22680 * f_quadpt ( ix, :, 14)
     1 + 581280 * f_quadpt ( ix, :, 15)
     1 - 567840 * f_quadpt ( ix, :, 16)
     1 - 53760 * f_quadpt ( ix, :, 17)
     1 - 53760 * f_quadpt ( ix, :, 18)
     1 - 144480 * f_quadpt ( ix, :, 19)
     1 - 144480 * f_quadpt ( ix, :, 20)

     1 - 2118690 * g_quadpt ( ix, :, 01)
     1 + 574560 * g_quadpt ( ix, :, 02)
     1 + 1141560 * g_quadpt ( ix, :, 03)
     1 + 443520 * g_quadpt ( ix, :, 04)
     1 + 1353240 * g_quadpt ( ix, :, 05)
     1 + 725130 * g_quadpt ( ix, :, 06)
     1 + 880740 * g_quadpt ( ix, :, 07)
     1 + 88200 * g_quadpt ( ix, :, 08)
     1 - 860580 * g_quadpt ( ix, :, 09)
     1 + 60480 * g_quadpt ( ix, :, 10)
     1 - 584640 * g_quadpt ( ix, :, 11)
     1 - 268380 * g_quadpt ( ix, :, 12)
     1 - 853020 * g_quadpt ( ix, :, 13)
     1 - 7560 * g_quadpt ( ix, :, 14)
     1 - 40320 * g_quadpt ( ix, :, 15)
     1 - 40320 * g_quadpt ( ix, :, 16)
     1 + 609840 * g_quadpt ( ix, :, 17)
     1 - 539280 * g_quadpt ( ix, :, 18)
     1 - 282240 * g_quadpt ( ix, :, 19)
     1 - 282240 * g_quadpt ( ix, :, 20)

     1 + 3363570 * h_quadpt ( ix, :, 01)
     1 - 1149120 * h_quadpt ( ix, :, 02)
     1 - 1431360 * h_quadpt ( ix, :, 03)
     1 - 168000 * h_quadpt ( ix, :, 04)
     1 - 2338560 * h_quadpt ( ix, :, 05)
     1 - 954240 * h_quadpt ( ix, :, 06)
     1 - 1832040 * h_quadpt ( ix, :, 07)
     1 - 1235850 * h_quadpt ( ix, :, 08)
     1 + 1706880 * h_quadpt ( ix, :, 09)
     1 - 423360 * h_quadpt ( ix, :, 10)
     1 + 1283520 * h_quadpt ( ix, :, 11)
     1 - 60480 * h_quadpt ( ix, :, 12)
     1 + 1223040 * h_quadpt ( ix, :, 13)
     1 + 483840 * h_quadpt ( ix, :, 14)
     1 + 537600 * h_quadpt ( ix, :, 15)
     1 + 537600 * h_quadpt ( ix, :, 16)
     1 + 295680 * h_quadpt ( ix, :, 17)
     1 + 295680 * h_quadpt ( ix, :, 18)
     1 - 1216320 * h_quadpt ( ix, :, 19)
     1 + 1081920 * h_quadpt ( ix, :, 20)

     1 - 577695 * s_quadpt ( ix, :, 01)
     1 + 114912 * s_quadpt ( ix, :, 02)
     1 + 228291 * s_quadpt ( ix, :, 03)
     1 + 100975 * s_quadpt ( ix, :, 04)
     1 + 252861 * s_quadpt ( ix, :, 05)
     1 + 122269 * s_quadpt ( ix, :, 06)
     1 + 304080 * s_quadpt ( ix, :, 07)
     1 + 84784 * s_quadpt ( ix, :, 08)
     1 - 171248 * s_quadpt ( ix, :, 09)
     1 + 11466 * s_quadpt ( ix, :, 10)
     1 - 159782 * s_quadpt ( ix, :, 11)
     1 + 1638 * s_quadpt ( ix, :, 12)
     1 - 122234 * s_quadpt ( ix, :, 13)
     1 - 13104 * s_quadpt ( ix, :, 14)
     1 - 35840 * s_quadpt ( ix, :, 15)
     1 - 35840 * s_quadpt ( ix, :, 16)
     1 - 29288 * s_quadpt ( ix, :, 17)
     1 - 29288 * s_quadpt ( ix, :, 18)
     1 + 119392 * s_quadpt ( ix, :, 19)
     1 - 110432 * s_quadpt ( ix, :, 20) ) / 646380.0

      tu ( ix, :, 28) = u_quadpt ( ix, :, 09)

     1 + ( - 1998990 * f_quadpt ( ix, :, 01)
     1 + 143640 * f_quadpt ( ix, :, 02)
     1 + 909720 * f_quadpt ( ix, :, 03)
     1 + 833910 * f_quadpt ( ix, :, 04)
     1 + 550620 * f_quadpt ( ix, :, 05)
     1 + 510720 * f_quadpt ( ix, :, 06)
     1 + 550620 * f_quadpt ( ix, :, 07)
     1 + 510720 * f_quadpt ( ix, :, 08)
     1 - 207480 * f_quadpt ( ix, :, 09)
     1 - 422940 * f_quadpt ( ix, :, 11)
     1 - 422940 * f_quadpt ( ix, :, 13)
     1 + 606480 * f_quadpt ( ix, :, 15)
     1 - 542640 * f_quadpt ( ix, :, 16)
     1 - 255360 * f_quadpt ( ix, :, 17)
     1 - 255360 * f_quadpt ( ix, :, 18)
     1 - 255360 * f_quadpt ( ix, :, 19)
     1 - 255360 * f_quadpt ( ix, :, 20)

     1 + 4488750 * g_quadpt ( ix, :, 01)
     1 - 1149120 * g_quadpt ( ix, :, 02)
     1 - 2366280 * g_quadpt ( ix, :, 03)
     1 - 987840 * g_quadpt ( ix, :, 04)
     1 - 2256660 * g_quadpt ( ix, :, 05)
     1 - 1828890 * g_quadpt ( ix, :, 06)
     1 - 1266300 * g_quadpt ( ix, :, 07)
     1 - 594720 * g_quadpt ( ix, :, 08)
     1 + 1064700 * g_quadpt ( ix, :, 09)
     1 - 30240 * g_quadpt ( ix, :, 10)
     1 + 1249920 * g_quadpt ( ix, :, 11)
     1 + 457380 * g_quadpt ( ix, :, 12)
     1 + 1707300 * g_quadpt ( ix, :, 13)
     1 - 211680 * g_quadpt ( ix, :, 14)
     1 + 403200 * g_quadpt ( ix, :, 15)
     1 + 403200 * g_quadpt ( ix, :, 16)
     1 - 1214640 * g_quadpt ( ix, :, 17)
     1 + 1083600 * g_quadpt ( ix, :, 18)
     1 + 524160 * g_quadpt ( ix, :, 19)
     1 + 524160 * g_quadpt ( ix, :, 20)

     1 + 3579030 * h_quadpt ( ix, :, 01)
     1 - 1149120 * h_quadpt ( ix, :, 02)
     1 - 1431360 * h_quadpt ( ix, :, 03)
     1 - 168000 * h_quadpt ( ix, :, 04)
     1 - 1692180 * h_quadpt ( ix, :, 05)
     1 - 954240 * h_quadpt ( ix, :, 06)
     1 - 1832040 * h_quadpt ( ix, :, 07)
     1 - 1451310 * h_quadpt ( ix, :, 08)
     1 + 1275960 * h_quadpt ( ix, :, 09)
     1 - 423360 * h_quadpt ( ix, :, 10)
     1 + 1283520 * h_quadpt ( ix, :, 11)
     1 - 60480 * h_quadpt ( ix, :, 12)
     1 + 1223040 * h_quadpt ( ix, :, 13)
     1 + 268380 * h_quadpt ( ix, :, 14)
     1 + 537600 * h_quadpt ( ix, :, 15)
     1 + 537600 * h_quadpt ( ix, :, 16)
     1 + 295680 * h_quadpt ( ix, :, 17)
     1 + 295680 * h_quadpt ( ix, :, 18)
     1 - 1216320 * h_quadpt ( ix, :, 19)
     1 + 1081920 * h_quadpt ( ix, :, 20)

     1 + 1023891 * s_quadpt ( ix, :, 01)
     1 - 229824 * s_quadpt ( ix, :, 02)
     1 - 333417 * s_quadpt ( ix, :, 03)
     1 - 78449 * s_quadpt ( ix, :, 04)
     1 - 284466 * s_quadpt ( ix, :, 05)
     1 - 341978 * s_quadpt ( ix, :, 06)
     1 - 288246 * s_quadpt ( ix, :, 07)
     1 - 345254 * s_quadpt ( ix, :, 08)
     1 + 6160 * s_quadpt ( ix, :, 09)
     1 + 2016 * s_quadpt ( ix, :, 10)
     1 + 295456 * s_quadpt ( ix, :, 11)
     1 - 1764 * s_quadpt ( ix, :, 12)
     1 + 293692 * s_quadpt ( ix, :, 13)
     1 - 252 * s_quadpt ( ix, :, 14)
     1 + 118888 * s_quadpt ( ix, :, 15)
     1 + 118888 * s_quadpt ( ix, :, 16)
     1 - 89264 * s_quadpt ( ix, :, 17)
     1 + 140560 * s_quadpt ( ix, :, 18)
     1 - 90272 * s_quadpt ( ix, :, 19)
     1 + 139552 * s_quadpt ( ix, :, 20)  ) / 646380.0

      tu ( ix, :, 29) = u_quadpt ( ix, :, 11)

     1 + ( 3579030 * f_quadpt ( ix, :, 01)
     1 - 1149120 * f_quadpt ( ix, :, 02)
     1 - 1832040 * f_quadpt ( ix, :, 03)
     1 - 1451310 * f_quadpt ( ix, :, 04)
     1 - 1431360 * f_quadpt ( ix, :, 05)
     1 - 168000 * f_quadpt ( ix, :, 06)
     1 - 1692180 * f_quadpt ( ix, :, 07)
     1 - 954240 * f_quadpt ( ix, :, 08)
     1 + 1223040 * f_quadpt ( ix, :, 09)
     1 + 268380 * f_quadpt ( ix, :, 10)
     1 + 1275960 * f_quadpt ( ix, :, 11)
     1 - 423360 * f_quadpt ( ix, :, 12)
     1 + 1283520 * f_quadpt ( ix, :, 13)
     1 - 60480 * f_quadpt ( ix, :, 14)
     1 - 1216320 * f_quadpt ( ix, :, 15)
     1 + 1081920 * f_quadpt ( ix, :, 16)
     1 + 537600 * f_quadpt ( ix, :, 17)
     1 + 537600 * f_quadpt ( ix, :, 18)
     1 + 295680 * f_quadpt ( ix, :, 19)
     1 + 295680 * f_quadpt ( ix, :, 20)

     1 - 1998990 * g_quadpt ( ix, :, 01)
     1 + 143640 * g_quadpt ( ix, :, 02)
     1 + 550620 * g_quadpt ( ix, :, 03)
     1 + 510720 * g_quadpt ( ix, :, 04)
     1 + 909720 * g_quadpt ( ix, :, 05)
     1 + 833910 * g_quadpt ( ix, :, 06)
     1 + 550620 * g_quadpt ( ix, :, 07)
     1 + 510720 * g_quadpt ( ix, :, 08)
     1 - 422940 * g_quadpt ( ix, :, 09)
     1 - 207480 * g_quadpt ( ix, :, 11)
     1 - 422940 * g_quadpt ( ix, :, 13)
     1 - 255360 * g_quadpt ( ix, :, 15)
     1 - 255360 * g_quadpt ( ix, :, 16)
     1 + 606480 * g_quadpt ( ix, :, 17)
     1 - 542640 * g_quadpt ( ix, :, 18)
     1 - 255360 * g_quadpt ( ix, :, 19)
     1 - 255360 * g_quadpt ( ix, :, 20)

     1 + 4488750 * h_quadpt ( ix, :, 01)
     1 - 1149120 * h_quadpt ( ix, :, 02)
     1 - 1266300 * h_quadpt ( ix, :, 03)
     1 - 594720 * h_quadpt ( ix, :, 04)
     1 - 2366280 * h_quadpt ( ix, :, 05)
     1 - 987840 * h_quadpt ( ix, :, 06)
     1 - 2256660 * h_quadpt ( ix, :, 07)
     1 - 1828890 * h_quadpt ( ix, :, 08)
     1 + 1707300 * h_quadpt ( ix, :, 09)
     1 - 211680 * h_quadpt ( ix, :, 10)
     1 + 1064700 * h_quadpt ( ix, :, 11)
     1 - 30240 * h_quadpt ( ix, :, 12)
     1 + 1249920 * h_quadpt ( ix, :, 13)
     1 + 457380 * h_quadpt ( ix, :, 14)
     1 + 524160 * h_quadpt ( ix, :, 15)
     1 + 524160 * h_quadpt ( ix, :, 16)
     1 + 403200 * h_quadpt ( ix, :, 17)
     1 + 403200 * h_quadpt ( ix, :, 18)
     1 - 1214640 * h_quadpt ( ix, :, 19)
     1 + 1083600 * h_quadpt ( ix, :, 20)

     1 + 1023891 * s_quadpt ( ix, :, 01)
     1 - 229824 * s_quadpt ( ix, :, 02)
     1 - 288246 * s_quadpt ( ix, :, 03)
     1 - 345254 * s_quadpt ( ix, :, 04)
     1 - 333417 * s_quadpt ( ix, :, 05)
     1 - 78449 * s_quadpt ( ix, :, 06)
     1 - 284466 * s_quadpt ( ix, :, 07)
     1 - 341978 * s_quadpt ( ix, :, 08)
     1 + 293692 * s_quadpt ( ix, :, 09)
     1 - 252 * s_quadpt ( ix, :, 10)
     1 + 6160 * s_quadpt ( ix, :, 11)
     1 + 2016 * s_quadpt ( ix, :, 12)
     1 + 295456 * s_quadpt ( ix, :, 13)
     1 - 1764 * s_quadpt ( ix, :, 14)
     1 - 90272 * s_quadpt ( ix, :, 15)
     1 + 139552 * s_quadpt ( ix, :, 16)
     1 + 118888 * s_quadpt ( ix, :, 17)
     1 + 118888 * s_quadpt ( ix, :, 18)
     1 - 89264 * s_quadpt ( ix, :, 19)
     1 + 140560 * s_quadpt ( ix, :, 20)  ) / 646380.0

      tu ( ix, :, 30) = u_quadpt ( ix, :, 13)

     1 + ( 4488750 * f_quadpt ( ix, :, 01)
     1 - 1149120 * f_quadpt ( ix, :, 02)
     1 - 2256660 * f_quadpt ( ix, :, 03)
     1 - 1828890 * f_quadpt ( ix, :, 04)
     1 - 1266300 * f_quadpt ( ix, :, 05)
     1 - 594720 * f_quadpt ( ix, :, 06)
     1 - 2366280 * f_quadpt ( ix, :, 07)
     1 - 987840 * f_quadpt ( ix, :, 08)
     1 + 1249920 * f_quadpt ( ix, :, 09)
     1 + 457380 * f_quadpt ( ix, :, 10)
     1 + 1707300 * f_quadpt ( ix, :, 11)
     1 - 211680 * f_quadpt ( ix, :, 12)
     1 + 1064700 * f_quadpt ( ix, :, 13)
     1 - 30240 * f_quadpt ( ix, :, 14)
     1 - 1214640 * f_quadpt ( ix, :, 15)
     1 + 1083600 * f_quadpt ( ix, :, 16)
     1 + 524160 * f_quadpt ( ix, :, 17)
     1 + 524160 * f_quadpt ( ix, :, 18)
     1 + 403200 * f_quadpt ( ix, :, 19)
     1 + 403200 * f_quadpt ( ix, :, 20)

     1 + 3579030 * g_quadpt ( ix, :, 01)
     1 - 1149120 * g_quadpt ( ix, :, 02)
     1 - 1692180 * g_quadpt ( ix, :, 03)
     1 - 954240 * g_quadpt ( ix, :, 04)
     1 - 1832040 * g_quadpt ( ix, :, 05)
     1 - 1451310 * g_quadpt ( ix, :, 06)
     1 - 1431360 * g_quadpt ( ix, :, 07)
     1 - 168000 * g_quadpt ( ix, :, 08)
     1 + 1283520 * g_quadpt ( ix, :, 09)
     1 - 60480 * g_quadpt ( ix, :, 10)
     1 + 1223040 * g_quadpt ( ix, :, 11)
     1 + 268380 * g_quadpt ( ix, :, 12)
     1 + 1275960 * g_quadpt ( ix, :, 13)
     1 - 423360 * g_quadpt ( ix, :, 14)
     1 + 295680 * g_quadpt ( ix, :, 15)
     1 + 295680 * g_quadpt ( ix, :, 16)
     1 - 1216320 * g_quadpt ( ix, :, 17)
     1 + 1081920 * g_quadpt ( ix, :, 18)
     1 + 537600 * g_quadpt ( ix, :, 19)
     1 + 537600 * g_quadpt ( ix, :, 20)

     1 - 1998990 * h_quadpt ( ix, :, 01)
     1 + 143640 * h_quadpt ( ix, :, 02)
     1 + 550620 * h_quadpt ( ix, :, 03)
     1 + 510720 * h_quadpt ( ix, :, 04)
     1 + 550620 * h_quadpt ( ix, :, 05)
     1 + 510720 * h_quadpt ( ix, :, 06)
     1 + 909720 * h_quadpt ( ix, :, 07)
     1 + 833910 * h_quadpt ( ix, :, 08)
     1 - 422940 * h_quadpt ( ix, :, 09)
     1 - 422940 * h_quadpt ( ix, :, 11)
     1 - 207480 * h_quadpt ( ix, :, 13)
     1 - 255360 * h_quadpt ( ix, :, 15)
     1 - 255360 * h_quadpt ( ix, :, 16)
     1 - 255360 * h_quadpt ( ix, :, 17)
     1 - 255360 * h_quadpt ( ix, :, 18)
     1 + 606480 * h_quadpt ( ix, :, 19)
     1 - 542640 * h_quadpt ( ix, :, 20)

     1 + 1023891 * s_quadpt ( ix, :, 01)
     1 - 229824 * s_quadpt ( ix, :, 02)
     1 - 284466 * s_quadpt ( ix, :, 03)
     1 - 341978 * s_quadpt ( ix, :, 04)
     1 - 288246 * s_quadpt ( ix, :, 05)
     1 - 345254 * s_quadpt ( ix, :, 06)
     1 - 333417 * s_quadpt ( ix, :, 07)
     1 - 78449 * s_quadpt ( ix, :, 08)
     1 + 295456 * s_quadpt ( ix, :, 09)
     1 - 1764 * s_quadpt ( ix, :, 10)
     1 + 293692 * s_quadpt ( ix, :, 11)
     1 - 252 * s_quadpt ( ix, :, 12)
     1 + 6160 * s_quadpt ( ix, :, 13)
     1 + 2016 * s_quadpt ( ix, :, 14)
     1 - 89264 * s_quadpt ( ix, :, 15)
     1 + 140560 * s_quadpt ( ix, :, 16)
     1 - 90272 * s_quadpt ( ix, :, 17)
     1 + 139552 * s_quadpt ( ix, :, 18)
     1 + 118888 * s_quadpt ( ix, :, 19)
     1 + 118888 * s_quadpt ( ix, :, 20)  ) / 646380.0



      tu ( ix, :, 31) = u_quadpt ( ix, :, 01)

     1 + ( - 307230 * f_quadpt ( ix, :, 01)
     1 + 95760 * f_quadpt ( ix, :, 02)
     1 + 152670 * f_quadpt ( ix, :, 03)
     1 + 129920 * f_quadpt ( ix, :, 04)
     1 + 119280 * f_quadpt ( ix, :, 05)
     1 + 14000 * f_quadpt ( ix, :, 06)
     1 + 194880 * f_quadpt ( ix, :, 07)
     1 + 79520 * f_quadpt ( ix, :, 08)
     1 - 101920 * f_quadpt ( ix, :, 09)
     1 - 40320 * f_quadpt ( ix, :, 10)
     1 - 142240 * f_quadpt ( ix, :, 11)
     1 + 35280 * f_quadpt ( ix, :, 12)
     1 - 106960 * f_quadpt ( ix, :, 13)
     1 + 5040 * f_quadpt ( ix, :, 14)
     1 + 101360 * f_quadpt ( ix, :, 15)
     1 - 90160 * f_quadpt ( ix, :, 16)
     1 - 44800 * f_quadpt ( ix, :, 17)
     1 - 44800 * f_quadpt ( ix, :, 18)
     1 - 24640 * f_quadpt ( ix, :, 19)
     1 - 24640 * f_quadpt ( ix, :, 20)

     1 - 307230 * g_quadpt ( ix, :, 01)
     1 + 95760 * g_quadpt ( ix, :, 02)
     1 + 194880 * g_quadpt ( ix, :, 03)
     1 + 79520 * g_quadpt ( ix, :, 04)
     1 + 152670 * g_quadpt ( ix, :, 05)
     1 + 129920 * g_quadpt ( ix, :, 06)
     1 + 119280 * g_quadpt ( ix, :, 07)
     1 + 14000 * g_quadpt ( ix, :, 08)
     1 - 106960 * g_quadpt ( ix, :, 09)
     1 + 5040 * g_quadpt ( ix, :, 10)
     1 - 101920 * g_quadpt ( ix, :, 11)
     1 - 40320 * g_quadpt ( ix, :, 12)
     1 - 142240 * g_quadpt ( ix, :, 13)
     1 + 35280 * g_quadpt ( ix, :, 14)
     1 - 24640 * g_quadpt ( ix, :, 15)
     1 - 24640 * g_quadpt ( ix, :, 16)
     1 + 101360 * g_quadpt ( ix, :, 17)
     1 - 90160 * g_quadpt ( ix, :, 18)
     1 - 44800 * g_quadpt ( ix, :, 19)
     1 - 44800 * g_quadpt ( ix, :, 20)

     1 - 307230 * h_quadpt ( ix, :, 01)
     1 + 95760 * h_quadpt ( ix, :, 02)
     1 + 119280 * h_quadpt ( ix, :, 03)
     1 + 14000 * h_quadpt ( ix, :, 04)
     1 + 194880 * h_quadpt ( ix, :, 05)
     1 + 79520 * h_quadpt ( ix, :, 06)
     1 + 152670 * h_quadpt ( ix, :, 07)
     1 + 129920 * h_quadpt ( ix, :, 08)
     1 - 142240 * h_quadpt ( ix, :, 09)
     1 + 35280 * h_quadpt ( ix, :, 10)
     1 - 106960 * h_quadpt ( ix, :, 11)
     1 + 5040 * h_quadpt ( ix, :, 12)
     1 - 101920 * h_quadpt ( ix, :, 13)
     1 - 40320 * h_quadpt ( ix, :, 14)
     1 - 44800 * h_quadpt ( ix, :, 15)
     1 - 44800 * h_quadpt ( ix, :, 16)
     1 - 24640 * h_quadpt ( ix, :, 17)
     1 - 24640 * h_quadpt ( ix, :, 18)
     1 + 101360 * h_quadpt ( ix, :, 19)
     1 - 90160 * h_quadpt ( ix, :, 20)

     1 - 31293 * s_quadpt ( ix, :, 01)
     1 + 5985 * s_quadpt ( ix, :, 03)
     1 + 5985 * s_quadpt ( ix, :, 04)
     1 + 5985 * s_quadpt ( ix, :, 05)
     1 + 5985 * s_quadpt ( ix, :, 06)
     1 + 5985 * s_quadpt ( ix, :, 07)
     1 + 5985 * s_quadpt ( ix, :, 08)  ) / 53865.0

      tu ( ix, :, 32) = u_quadpt ( ix, :, 03)

     1 + ( 3902220 * f_quadpt ( ix, :, 01)
     1 - 1149120 * f_quadpt ( ix, :, 02)
     1 - 2047500 * f_quadpt ( ix, :, 03)
     1 - 1559040 * f_quadpt ( ix, :, 04)
     1 - 1431360 * f_quadpt ( ix, :, 05)
     1 - 168000 * f_quadpt ( ix, :, 06)
     1 - 2338560 * f_quadpt ( ix, :, 07)
     1 - 954240 * f_quadpt ( ix, :, 08)
     1 + 1223040 * f_quadpt ( ix, :, 09)
     1 + 483840 * f_quadpt ( ix, :, 10)
     1 + 1706880 * f_quadpt ( ix, :, 11)
     1 - 423360 * f_quadpt ( ix, :, 12)
     1 + 1283520 * f_quadpt ( ix, :, 13)
     1 - 60480 * f_quadpt ( ix, :, 14)
     1 - 1216320 * f_quadpt ( ix, :, 15)
     1 + 1081920 * f_quadpt ( ix, :, 16)
     1 + 537600 * f_quadpt ( ix, :, 17)
     1 + 537600 * f_quadpt ( ix, :, 18)
     1 + 295680 * f_quadpt ( ix, :, 19)
     1 + 295680 * f_quadpt ( ix, :, 20)

     1 - 1412460 * g_quadpt ( ix, :, 01)
     1 + 574560 * g_quadpt ( ix, :, 02)
     1 + 738360 * g_quadpt ( ix, :, 03)
     1 + 477120 * g_quadpt ( ix, :, 04)
     1 + 700560 * g_quadpt ( ix, :, 05)
     1 + 564060 * g_quadpt ( ix, :, 06)
     1 + 715680 * g_quadpt ( ix, :, 07)
     1 + 84000 * g_quadpt ( ix, :, 08)
     1 - 641760 * g_quadpt ( ix, :, 09)
     1 + 30240 * g_quadpt ( ix, :, 10)
     1 - 611520 * g_quadpt ( ix, :, 11)
     1 - 26460 * g_quadpt ( ix, :, 12)
     1 - 637980 * g_quadpt ( ix, :, 13)
     1 + 211680 * g_quadpt ( ix, :, 14)
     1 - 147840 * g_quadpt ( ix, :, 15)
     1 - 147840 * g_quadpt ( ix, :, 16)
     1 + 608160 * g_quadpt ( ix, :, 17)
     1 - 540960 * g_quadpt ( ix, :, 18)
     1 - 268800 * g_quadpt ( ix, :, 19)
     1 - 268800 * g_quadpt ( ix, :, 20)

     1 - 2322180 * h_quadpt ( ix, :, 01)
     1 + 574560 * h_quadpt ( ix, :, 02)
     1 + 766080 * h_quadpt ( ix, :, 03)
     1 + 510720 * h_quadpt ( ix, :, 04)
     1 + 1197000 * h_quadpt ( ix, :, 05)
     1 + 510720 * h_quadpt ( ix, :, 06)
     1 + 1125180 * h_quadpt ( ix, :, 07)
     1 + 941640 * h_quadpt ( ix, :, 08)
     1 - 853860 * h_quadpt ( ix, :, 09)
     1 - 638400 * h_quadpt ( ix, :, 11)
     1 - 638400 * h_quadpt ( ix, :, 13)
     1 - 215460 * h_quadpt ( ix, :, 14)
     1 - 255360 * h_quadpt ( ix, :, 15)
     1 - 255360 * h_quadpt ( ix, :, 16)
     1 - 255360 * h_quadpt ( ix, :, 17)
     1 - 255360 * h_quadpt ( ix, :, 18)
     1 + 606480 * h_quadpt ( ix, :, 19)
     1 - 542640 * h_quadpt ( ix, :, 20)

     1 + 413535 * s_quadpt ( ix, :, 01)
     1 - 114912 * s_quadpt ( ix, :, 02)
     1 - 343581 * s_quadpt ( ix, :, 03)
     1 - 156604 * s_quadpt ( ix, :, 04)
     1 - 138516 * s_quadpt ( ix, :, 05)
     1 - 11200 * s_quadpt ( ix, :, 06)
     1 - 163086 * s_quadpt ( ix, :, 07)
     1 - 32494 * s_quadpt ( ix, :, 08)
     1 + 122234 * s_quadpt ( ix, :, 09)
     1 + 13104 * s_quadpt ( ix, :, 10)
     1 + 171248 * s_quadpt ( ix, :, 11)
     1 - 11466 * s_quadpt ( ix, :, 12)
     1 + 159782 * s_quadpt ( ix, :, 13)
     1 - 1638 * s_quadpt ( ix, :, 14)
     1 - 119392 * s_quadpt ( ix, :, 15)
     1 + 110432 * s_quadpt ( ix, :, 16)
     1 + 35840 * s_quadpt ( ix, :, 17)
     1 + 35840 * s_quadpt ( ix, :, 18)
     1 + 29288 * s_quadpt ( ix, :, 19)
     1 + 29288 * s_quadpt ( ix, :, 20)  ) / 323190.0

      tu ( ix, :, 33) = u_quadpt ( ix, :, 05)

     1 + ( - 2322180 * f_quadpt ( ix, :, 01)
     1 + 574560 * f_quadpt ( ix, :, 02)
     1 + 1125180 * f_quadpt ( ix, :, 03)
     1 + 941640 * f_quadpt ( ix, :, 04)
     1 + 766080 * f_quadpt ( ix, :, 05)
     1 + 510720 * f_quadpt ( ix, :, 06)
     1 + 1197000 * f_quadpt ( ix, :, 07)
     1 + 510720 * f_quadpt ( ix, :, 08)
     1 - 638400 * f_quadpt ( ix, :, 09)
     1 - 215460 * f_quadpt ( ix, :, 10)
     1 - 853860 * f_quadpt ( ix, :, 11)
     1 - 638400 * f_quadpt ( ix, :, 13)
     1 + 606480 * f_quadpt ( ix, :, 15)
     1 - 542640 * f_quadpt ( ix, :, 16)
     1 - 255360 * f_quadpt ( ix, :, 17)
     1 - 255360 * f_quadpt ( ix, :, 18)
     1 - 255360 * f_quadpt ( ix, :, 19)
     1 - 255360 * f_quadpt ( ix, :, 20)

     1 + 3902220 * g_quadpt ( ix, :, 01)
     1 - 1149120 * g_quadpt ( ix, :, 02)
     1 - 2338560 * g_quadpt ( ix, :, 03)
     1 - 954240 * g_quadpt ( ix, :, 04)
     1 - 2047500 * g_quadpt ( ix, :, 05)
     1 - 1559040 * g_quadpt ( ix, :, 06)
     1 - 1431360 * g_quadpt ( ix, :, 07)
     1 - 168000 * g_quadpt ( ix, :, 08)
     1 + 1283520 * g_quadpt ( ix, :, 09)
     1 - 60480 * g_quadpt ( ix, :, 10)
     1 + 1223040 * g_quadpt ( ix, :, 11)
     1 + 483840 * g_quadpt ( ix, :, 12)
     1 + 1706880 * g_quadpt ( ix, :, 13)
     1 - 423360 * g_quadpt ( ix, :, 14)
     1 + 295680 * g_quadpt ( ix, :, 15)
     1 + 295680 * g_quadpt ( ix, :, 16)
     1 - 1216320 * g_quadpt ( ix, :, 17)
     1 + 1081920 * g_quadpt ( ix, :, 18)
     1 + 537600 * g_quadpt ( ix, :, 19)
     1 + 537600 * g_quadpt ( ix, :, 20)

     1 - 1412460 * h_quadpt ( ix, :, 01)
     1 + 574560 * h_quadpt ( ix, :, 02)
     1 + 715680 * h_quadpt ( ix, :, 03)
     1 + 84000 * h_quadpt ( ix, :, 04)
     1 + 738360 * h_quadpt ( ix, :, 05)
     1 + 477120 * h_quadpt ( ix, :, 06)
     1 + 700560 * h_quadpt ( ix, :, 07)
     1 + 564060 * h_quadpt ( ix, :, 08)
     1 - 637980 * h_quadpt ( ix, :, 09)
     1 + 211680 * h_quadpt ( ix, :, 10)
     1 - 641760 * h_quadpt ( ix, :, 11)
     1 + 30240 * h_quadpt ( ix, :, 12)
     1 - 611520 * h_quadpt ( ix, :, 13)
     1 - 26460 * h_quadpt ( ix, :, 14)
     1 - 268800 * h_quadpt ( ix, :, 15)
     1 - 268800 * h_quadpt ( ix, :, 16)
     1 - 147840 * h_quadpt ( ix, :, 17)
     1 - 147840 * h_quadpt ( ix, :, 18)
     1 + 608160 * h_quadpt ( ix, :, 19)
     1 - 540960 * h_quadpt ( ix, :, 20)

     1 + 413535 * s_quadpt ( ix, :, 01)
     1 - 114912 * s_quadpt ( ix, :, 02)
     1 - 163086 * s_quadpt ( ix, :, 03)
     1 - 32494 * s_quadpt ( ix, :, 04)
     1 - 343581 * s_quadpt ( ix, :, 05)
     1 - 156604 * s_quadpt ( ix, :, 06)
     1 - 138516 * s_quadpt ( ix, :, 07)
     1 - 11200 * s_quadpt ( ix, :, 08)
     1 + 159782 * s_quadpt ( ix, :, 09)
     1 - 1638 * s_quadpt ( ix, :, 10)
     1 + 122234 * s_quadpt ( ix, :, 11)
     1 + 13104 * s_quadpt ( ix, :, 12)
     1 + 171248 * s_quadpt ( ix, :, 13)
     1 - 11466 * s_quadpt ( ix, :, 14)
     1 + 29288 * s_quadpt ( ix, :, 15)
     1 + 29288 * s_quadpt ( ix, :, 16)
     1 - 119392 * s_quadpt ( ix, :, 17)
     1 + 110432 * s_quadpt ( ix, :, 18)
     1 + 35840 * s_quadpt ( ix, :, 19)
     1 + 35840 * s_quadpt ( ix, :, 20)  ) / 323190.0

      tu ( ix, :, 34) = u_quadpt ( ix, :, 07)

     1 + ( - 1412460 * f_quadpt ( ix, :, 01)
     1 + 574560 * f_quadpt ( ix, :, 02)
     1 + 700560 * f_quadpt ( ix, :, 03)
     1 + 564060 * f_quadpt ( ix, :, 04)
     1 + 715680 * f_quadpt ( ix, :, 05)
     1 + 84000 * f_quadpt ( ix, :, 06)
     1 + 738360 * f_quadpt ( ix, :, 07)
     1 + 477120 * f_quadpt ( ix, :, 08)
     1 - 611520 * f_quadpt ( ix, :, 09)
     1 - 26460 * f_quadpt ( ix, :, 10)
     1 - 637980 * f_quadpt ( ix, :, 11)
     1 + 211680 * f_quadpt ( ix, :, 12)
     1 - 641760 * f_quadpt ( ix, :, 13)
     1 + 30240 * f_quadpt ( ix, :, 14)
     1 + 608160 * f_quadpt ( ix, :, 15)
     1 - 540960 * f_quadpt ( ix, :, 16)
     1 - 268800 * f_quadpt ( ix, :, 17)
     1 - 268800 * f_quadpt ( ix, :, 18)
     1 - 147840 * f_quadpt ( ix, :, 19)
     1 - 147840 * f_quadpt ( ix, :, 20)

     1 - 2322180 * g_quadpt ( ix, :, 01)
     1 + 574560 * g_quadpt ( ix, :, 02)
     1 + 1197000 * g_quadpt ( ix, :, 03)
     1 + 510720 * g_quadpt ( ix, :, 04)
     1 + 1125180 * g_quadpt ( ix, :, 05)
     1 + 941640 * g_quadpt ( ix, :, 06)
     1 + 766080 * g_quadpt ( ix, :, 07)
     1 + 510720 * g_quadpt ( ix, :, 08)
     1 - 638400 * g_quadpt ( ix, :, 09)
     1 - 638400 * g_quadpt ( ix, :, 11)
     1 - 215460 * g_quadpt ( ix, :, 12)
     1 - 853860 * g_quadpt ( ix, :, 13)
     1 - 255360 * g_quadpt ( ix, :, 15)
     1 - 255360 * g_quadpt ( ix, :, 16)
     1 + 606480 * g_quadpt ( ix, :, 17)
     1 - 542640 * g_quadpt ( ix, :, 18)
     1 - 255360 * g_quadpt ( ix, :, 19)
     1 - 255360 * g_quadpt ( ix, :, 20)

     1 + 3902220 * h_quadpt ( ix, :, 01)
     1 - 1149120 * h_quadpt ( ix, :, 02)
     1 - 1431360 * h_quadpt ( ix, :, 03)
     1 - 168000 * h_quadpt ( ix, :, 04)
     1 - 2338560 * h_quadpt ( ix, :, 05)
     1 - 954240 * h_quadpt ( ix, :, 06)
     1 - 2047500 * h_quadpt ( ix, :, 07)
     1 - 1559040 * h_quadpt ( ix, :, 08)
     1 + 1706880 * h_quadpt ( ix, :, 09)
     1 - 423360 * h_quadpt ( ix, :, 10)
     1 + 1283520 * h_quadpt ( ix, :, 11)
     1 - 60480 * h_quadpt ( ix, :, 12)
     1 + 1223040 * h_quadpt ( ix, :, 13)
     1 + 483840 * h_quadpt ( ix, :, 14)
     1 + 537600 * h_quadpt ( ix, :, 15)
     1 + 537600 * h_quadpt ( ix, :, 16)
     1 + 295680 * h_quadpt ( ix, :, 17)
     1 + 295680 * h_quadpt ( ix, :, 18)
     1 - 1216320 * h_quadpt ( ix, :, 19)
     1 + 1081920 * h_quadpt ( ix, :, 20)

     1 + 413535 * s_quadpt ( ix, :, 01)
     1 - 114912 * s_quadpt ( ix, :, 02)
     1 - 138516 * s_quadpt ( ix, :, 03)
     1 - 11200 * s_quadpt ( ix, :, 04)
     1 - 163086 * s_quadpt ( ix, :, 05)
     1 - 32494 * s_quadpt ( ix, :, 06)
     1 - 343581 * s_quadpt ( ix, :, 07)
     1 - 156604 * s_quadpt ( ix, :, 08)
     1 + 171248 * s_quadpt ( ix, :, 09)
     1 - 11466 * s_quadpt ( ix, :, 10)
     1 + 159782 * s_quadpt ( ix, :, 11)
     1 - 1638 * s_quadpt ( ix, :, 12)
     1 + 122234 * s_quadpt ( ix, :, 13)
     1 + 13104 * s_quadpt ( ix, :, 14)
     1 + 35840 * s_quadpt ( ix, :, 15)
     1 + 35840 * s_quadpt ( ix, :, 16)
     1 + 29288 * s_quadpt ( ix, :, 17)
     1 + 29288 * s_quadpt ( ix, :, 18)
     1 - 119392 * s_quadpt ( ix, :, 19)
     1 + 110432 * s_quadpt ( ix, :, 20)  ) / 323190.0

      tu ( ix, :, 35) = u_quadpt ( ix, :, 01)

     1 + ( - 530670 * f_quadpt ( ix, :, 01)
     1 + 191520 * f_quadpt ( ix, :, 02)
     1 + 257460 * f_quadpt ( ix, :, 03)
     1 + 223930 * f_quadpt ( ix, :, 04)
     1 + 238560 * f_quadpt ( ix, :, 05)
     1 + 28000 * f_quadpt ( ix, :, 06)
     1 + 389760 * f_quadpt ( ix, :, 07)
     1 + 159040 * f_quadpt ( ix, :, 08)
     1 - 203840 * f_quadpt ( ix, :, 09)
     1 - 80640 * f_quadpt ( ix, :, 10)
     1 - 284480 * f_quadpt ( ix, :, 11)
     1 + 70560 * f_quadpt ( ix, :, 12)
     1 - 213920 * f_quadpt ( ix, :, 13)
     1 + 10080 * f_quadpt ( ix, :, 14)
     1 + 202720 * f_quadpt ( ix, :, 15)
     1 - 180320 * f_quadpt ( ix, :, 16)
     1 - 89600 * f_quadpt ( ix, :, 17)
     1 - 89600 * f_quadpt ( ix, :, 18)
     1 - 49280 * f_quadpt ( ix, :, 19)
     1 - 49280 * f_quadpt ( ix, :, 20)

     1 - 530670 * g_quadpt ( ix, :, 01)
     1 + 191520 * g_quadpt ( ix, :, 02)
     1 + 389760 * g_quadpt ( ix, :, 03)
     1 + 159040 * g_quadpt ( ix, :, 04)
     1 + 257460 * g_quadpt ( ix, :, 05)
     1 + 223930 * g_quadpt ( ix, :, 06)
     1 + 238560 * g_quadpt ( ix, :, 07)
     1 + 28000 * g_quadpt ( ix, :, 08)
     1 - 213920 * g_quadpt ( ix, :, 09)
     1 + 10080 * g_quadpt ( ix, :, 10)
     1 - 203840 * g_quadpt ( ix, :, 11)
     1 - 80640 * g_quadpt ( ix, :, 12)
     1 - 284480 * g_quadpt ( ix, :, 13)
     1 + 70560 * g_quadpt ( ix, :, 14)
     1 - 49280 * g_quadpt ( ix, :, 15)
     1 - 49280 * g_quadpt ( ix, :, 16)
     1 + 202720 * g_quadpt ( ix, :, 17)
     1 - 180320 * g_quadpt ( ix, :, 18)
     1 - 89600 * g_quadpt ( ix, :, 19)
     1 - 89600 * g_quadpt ( ix, :, 20)

     1 - 530670 * h_quadpt ( ix, :, 01)
     1 + 191520 * h_quadpt ( ix, :, 02)
     1 + 238560 * h_quadpt ( ix, :, 03)
     1 + 28000 * h_quadpt ( ix, :, 04)
     1 + 389760 * h_quadpt ( ix, :, 05)
     1 + 159040 * h_quadpt ( ix, :, 06)
     1 + 257460 * h_quadpt ( ix, :, 07)
     1 + 223930 * h_quadpt ( ix, :, 08)
     1 - 284480 * h_quadpt ( ix, :, 09)
     1 + 70560 * h_quadpt ( ix, :, 10)
     1 - 213920 * h_quadpt ( ix, :, 11)
     1 + 10080 * h_quadpt ( ix, :, 12)
     1 - 203840 * h_quadpt ( ix, :, 13)
     1 - 80640 * h_quadpt ( ix, :, 14)
     1 - 89600 * h_quadpt ( ix, :, 15)
     1 - 89600 * h_quadpt ( ix, :, 16)
     1 - 49280 * h_quadpt ( ix, :, 17)
     1 - 49280 * h_quadpt ( ix, :, 18)
     1 + 202720 * h_quadpt ( ix, :, 19)
     1 - 180320 * h_quadpt ( ix, :, 20)

     1 - 29241 * s_quadpt ( ix, :, 01)
     1 + 5985 * s_quadpt ( ix, :, 03)
     1 + 5985 * s_quadpt ( ix, :, 04)
     1 + 5985 * s_quadpt ( ix, :, 05)
     1 + 5985 * s_quadpt ( ix, :, 06)
     1 + 5985 * s_quadpt ( ix, :, 07)
     1 + 5985 * s_quadpt ( ix, :, 08)  ) / 71820.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Start the outer ADER iteration here.

      DO iter_ader = 1, max_iter_ader_cg

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 7) Get the primitive variables from the conserved variables at
! the spatial nodes at t>0. The density and pressure variables are
! given positive floor values to ensure that they can be used for physically
! meaningful source terms.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, 1, iquadpt) =
     1                   AMAX1 ( u_quadpt ( ix, 1, iquadpt), rhofloor)

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempaa = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempaa
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempaa
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempaa

      bfx ( ix, iquadpt) = u_quadpt ( ix, 6, iquadpt)
      bfy ( ix, iquadpt) = u_quadpt ( ix, 7, iquadpt)
      bfz ( ix, iquadpt) = u_quadpt ( ix, 8, iquadpt)

      tempab = 0.5 * rho ( ix, iquadpt) * ( vlx ( ix, iquadpt)**2
     1         + vly ( ix, iquadpt)**2 + vlz ( ix, iquadpt)**2)

      tempac = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      eng_t ( ix, iquadpt) = u_quadpt ( ix, 5, iquadpt)
     1                     - tempab - tempac

      prs ( ix, iquadpt) = eng_t ( ix, iquadpt)
     1                   * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

#if ( ENFORCE_POSITIVITY == 1)
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = entropy ( ix, iquadpt)
     1* rho ( ix, iquadpt)**( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)
      prs ( ix, iquadpt) = AMAX1 ( prs ( ix, iquadpt), prsfloor)

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#else
      IF ( prs ( ix, iquadpt) .LE. prsfloor) THEN

      prs ( ix, iquadpt) = prsfloor

      eng_t ( ix, iquadpt) = prs ( ix, iquadpt)
     1                 / ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      u_quadpt ( ix, 5, iquadpt) = eng_t ( ix, iquadpt)
     1                           + tempab + tempac

      END IF
#endif

#if ( NFLUID >= 1)
      specie ( ix, 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( ix, 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( ix, 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( ix, 1: NFLUID, iquadpt) * rho ( ix, iquadpt)
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, ivar + 1, iquadpt) =
     1         AMAX1 ( u_quadpt ( ix, ivar + 1, iquadpt), rhofloor)

      rho_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 1, iquadpt)

      tempaa = 1.0 / rho_euler ( ix, iquadpt, i_euler)

      vlx_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 2, iquadpt) * tempaa
      vly_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 3, iquadpt) * tempaa
      vlz_euler ( ix, iquadpt, i_euler) =
     1                       u_quadpt ( ix, ivar + 4, iquadpt) * tempaa

      tempab = 0.5 * rho_euler ( ix, iquadpt, i_euler)
     1           * ( vlx_euler ( ix, iquadpt, i_euler)**2
     1             + vly_euler ( ix, iquadpt, i_euler)**2
     1             + vlz_euler ( ix, iquadpt, i_euler)**2)

      prs_euler ( ix, iquadpt, i_euler) =
     1         ( u_quadpt ( ix, ivar + 5, iquadpt) - tempab)
     1       * ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      prs_euler ( ix, iquadpt, i_euler) =
     1    AMAX1 ( prs_euler ( ix, iquadpt, i_euler), prsfloor)

#if ( ISOTHERMAL_EQNS == 1)
      prs_euler ( ix, iquadpt, i_euler) =
     1     isospeed**2 * rho_euler ( ix, iquadpt, i_euler)
     1   + barotropic_stiffness
     1   * rho_euler ( ix, iquadpt, i_euler)**gamma
#endif

      eng_t_euler ( ix, iquadpt, i_euler) =
     1         prs_euler ( ix, iquadpt, i_euler)
     1       / ( gamma_eos_quadpt ( ix, i_euler + 1, iquadpt) - 1.0)

      IF ( prs_euler ( ix, iquadpt, i_euler) .LT. prsfloor * 1.00000001)
     1   u_quadpt ( ix, ivar + 5, iquadpt) =
     1                   eng_t_euler ( ix, iquadpt, i_euler) + tempab

#if ( NFLUID_EULER >= 1)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) = tempaa
     1       * u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1                    ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)

      specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler) =
     1    AMIN1 ( AMAX1 (
     1    specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler),
     1            0.0), 1.0 )

      u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1       specie_euler ( ix, 1: NFLUID_EULER, iquadpt, i_euler)
     1     * rho_euler ( ix, iquadpt, i_euler)

#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 8.a) Build fluxes ( f_quadpt, g_quadpt, h_quadpt) from the primitive
! variables that have now become available at spatial nodes at t=0.

! Multiply "f_quadpt" by "dtmove/dx" and so on.

! Do this for non-relativistic and relativistic MHD.

! ------------------------------

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! -----

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

! -----

      f_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 2, iquadpt)

      f_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfx ( ix, iquadpt)**2 * ob4pi

      f_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vly ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi

      f_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 2, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi

      f_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab

      f_quadpt ( ix, 6, iquadpt) = 0.0

      f_quadpt ( ix, 7, iquadpt) =
     1           vlx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         - vly ( ix, iquadpt) * bfx ( ix, iquadpt)

      f_quadpt ( ix, 8, iquadpt) =
     1           vlx ( ix, iquadpt) * bfz ( ix, iquadpt)
     1         - vlz ( ix, iquadpt) * bfx ( ix, iquadpt)

#if ( NFLUID >= 1)
      f_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlx ( ix, iquadpt)
#endif

! -----

      g_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 3, iquadpt)

      g_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         - bfy ( ix, iquadpt) * bfx ( ix, iquadpt) * ob4pi

      g_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vly ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfy ( ix, iquadpt)**2 * ob4pi

      g_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 3, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         - bfy ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi

      g_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vly ( ix, iquadpt)
     1                       - bfy ( ix, iquadpt) * tempab

      g_quadpt ( ix, 6, iquadpt) =
     1           vly ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         - vlx ( ix, iquadpt) * bfy ( ix, iquadpt)

      g_quadpt ( ix, 7, iquadpt) = 0.0

      g_quadpt ( ix, 8, iquadpt) =
     1           vly ( ix, iquadpt) * bfz ( ix, iquadpt)
     1         - vlz ( ix, iquadpt) * bfy ( ix, iquadpt)

#if ( NFLUID >= 1)
      g_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vly ( ix, iquadpt)
#endif

! -----

      h_quadpt ( ix, 1, iquadpt) = u_quadpt ( ix, 4, iquadpt)

      h_quadpt ( ix, 2, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vlx ( ix, iquadpt)
     1         - bfz ( ix, iquadpt) * bfx ( ix, iquadpt) * ob4pi

      h_quadpt ( ix, 3, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vly ( ix, iquadpt)
     1         - bfz ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi

      h_quadpt ( ix, 4, iquadpt) = u_quadpt ( ix, 4, iquadpt)
     1         * vlz ( ix, iquadpt)
     1         + prs ( ix, iquadpt)
     1         + tempaa - bfz ( ix, iquadpt)**2 * ob4pi

      h_quadpt ( ix, 5, iquadpt) = ( u_quadpt ( ix, 5, iquadpt)
     1         + prs ( ix, iquadpt) + tempaa) * vlz ( ix, iquadpt)
     1                       - bfz ( ix, iquadpt) * tempab

      h_quadpt ( ix, 6, iquadpt) =
     1           vlz ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         - vlx ( ix, iquadpt) * bfz ( ix, iquadpt)

      h_quadpt ( ix, 7, iquadpt) =
     1           vlz ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         - vly ( ix, iquadpt) * bfz ( ix, iquadpt)

      h_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      h_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1u_quadpt ( ix,  NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
     1         * vlz ( ix, iquadpt)
#endif

! -----

#if ( ISOTHERMAL_EQNS == 1)
      f_quadpt ( ix, 5, iquadpt) = 0.0
      g_quadpt ( ix, 5, iquadpt) = 0.0
      h_quadpt ( ix, 5, iquadpt) = 0.0
#endif

! -----

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! -----

      f_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)

      f_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 2, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)

      f_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vlx_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      f_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vlx_euler ( ix, iquadpt, i_euler)
#endif

! -----

      g_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)

      g_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 3, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)

      g_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vly_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      g_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vly_euler ( ix, iquadpt, i_euler)
#endif

! -----

      h_quadpt ( ix, ivar + 1, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)

      h_quadpt ( ix, ivar + 2, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vlx_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 3, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vly_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 4, iquadpt) =
     1                          u_quadpt ( ix, ivar + 4, iquadpt)
     1                       * vlz_euler ( ix, iquadpt, i_euler)
     1                       + prs_euler ( ix, iquadpt, i_euler)

      h_quadpt ( ix, ivar + 5, iquadpt) =
     1                        ( u_quadpt ( ix, ivar + 5, iquadpt)
     1                       + prs_euler ( ix, iquadpt, i_euler) )
     1                       * vlz_euler ( ix, iquadpt, i_euler)

#if ( NFLUID_EULER >= 1)
      h_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) =
     1u_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt)
     1         * vlz_euler ( ix, iquadpt, i_euler)
#endif

! -----

#if ( ISOTHERMAL_EQNS == 1)
      f_quadpt ( ix, ivar + 5, iquadpt) = 0.0
      g_quadpt ( ix, ivar + 5, iquadpt) = 0.0
      h_quadpt ( ix, ivar + 5, iquadpt) = 0.0
#endif

! -----

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ------------------------------

! To make cylindrical and spherical geometries accessible to non-rel.
! and relativistic flow, we put this part after either flux construction.

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      IF ( igeom == 2 .or. igeom == 3) THEN

         g_quadpt ( ix, :, iquadpt) = g_quadpt ( ix, :, iquadpt) * obx1

      END IF

      IF ( igeom == 3) THEN

         sin_x2 = 1.0 / ( x1 * sin ( x2))

         h_quadpt ( ix, :, iquadpt) = h_quadpt ( ix, :, iquadpt)
     1                              * sin_x2

      END IF

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 8.b) Optional step needed only for viscous, resistive or conductive MHD.

! ASSUME that "du_quadpt_dx, du_quadpt_dy and du_quadpt_dz" have been built at
! all orders (NOT DONE YET).
! We set about building the derivatives of the primitive
! variables and various auxiliary variables at the nodal points. These are
! needed only in the presence of viscous, resistive or conduction terms.

! These derivatives are then used to build non-ideal fluxes at the quadrature
! points given by "fni_quadpt, gni_quadpt and hni_quadpt".

! Lastly, the non-ideal fluxes are added into the ideal fluxes

! ----------------------------------------------------------------------

      IF ( ( ( viscosity_coef .GT. 0.0) .OR.
     1       ( resistivity_coef .GT. 0.0)) .AND. ( NUMORDER >= 3)) THEN

! ------------------------------

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = 1.0 / rho ( ix, iquadpt)

      drho_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 1, iquadpt)
      drho_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 1, iquadpt)
      drho_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 1, iquadpt)


      dvlx_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvlx_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvlx_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 2, iquadpt)
     1                   - vlx ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dvly_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvly_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvly_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 3, iquadpt)
     1                   - vly ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dvlz_dx ( ix, iquadpt) = tempaa * ( du_quadpt_dx ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dx ( ix, iquadpt) )

      dvlz_dy ( ix, iquadpt) = tempaa * ( du_quadpt_dy ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dy ( ix, iquadpt) )

      dvlz_dz ( ix, iquadpt) = tempaa * ( du_quadpt_dz ( ix, 4, iquadpt)
     1                   - vlz ( ix, iquadpt) * drho_dz ( ix, iquadpt) )


      dbfx_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 6, iquadpt)
      dbfx_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 6, iquadpt)
      dbfx_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 6, iquadpt)

      dbfy_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 7, iquadpt)
      dbfy_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 7, iquadpt)
      dbfy_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 7, iquadpt)

      dbfz_dx ( ix, iquadpt) = du_quadpt_dx ( ix, 8, iquadpt)
      dbfz_dy ( ix, iquadpt) = du_quadpt_dy ( ix, 8, iquadpt)
      dbfz_dz ( ix, iquadpt) = du_quadpt_dz ( ix, 8, iquadpt)


      dv2_dx ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dx ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dx ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dx ( ix, iquadpt) )

      dv2_dy ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dy ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dy ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dy ( ix, iquadpt) )

      dv2_dz ( ix, iquadpt) =
     1            2.0 * ( vlx ( ix, iquadpt) * dvlx_dz ( ix, iquadpt)
     1                  + vly ( ix, iquadpt) * dvly_dz ( ix, iquadpt)
     1                  + vlz ( ix, iquadpt) * dvlz_dz ( ix, iquadpt) )

      db2_dx ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dx ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dx ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dx ( ix, iquadpt) )

      db2_dy ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dy ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dy ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dy ( ix, iquadpt) )

      db2_dz ( ix, iquadpt) =
     1            2.0 * ( bfx ( ix, iquadpt) * dbfx_dz ( ix, iquadpt)
     1                  + bfy ( ix, iquadpt) * dbfy_dz ( ix, iquadpt)
     1                  + bfz ( ix, iquadpt) * dbfz_dz ( ix, iquadpt) )


      tempac = vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1       + vlz ( ix, iquadpt)**2

      tempab = molewt_quadpt ( ix, 1, iquadpt)
     1       / ( gasconst * rho ( ix, iquadpt) )

      temperature ( ix, iquadpt) = prs ( ix, iquadpt) * tempab

      tempab = tempab * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0)

      dtemperature_dx ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dx ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dx ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dx ( ix, iquadpt)
     1    - ob8pi * db2_dx ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dx ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

      dtemperature_dy ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dy ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dy ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dy ( ix, iquadpt)
     1    - ob8pi * db2_dy ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dy ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

      dtemperature_dz ( ix, iquadpt) =
     1      tempab * ( du_quadpt_dz ( ix, 5, iquadpt)
     1    - 0.5 * rho ( ix, iquadpt) * dv2_dz ( ix, iquadpt)
     1    - 0.5 * tempac * drho_dz ( ix, iquadpt)
     1    - ob8pi * db2_dz ( ix, iquadpt) )
     1 - temperature ( ix, iquadpt) * drho_dz ( ix, iquadpt)
     1 / rho ( ix, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

! The derivatives are within the reference element's local coordinates.
! To get true physical fluxes, we need to divide by the size of the zone here.

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = 1.0 / dx ( ix)
      tempab = 1.0 / dy ( ix)
      tempac = 1.0 / dz ( ix)

      drho_dx ( ix, iquadpt) = tempaa * drho_dx ( ix, iquadpt)
      dtemperature_dx ( ix, iquadpt) =
     1                 tempaa * dtemperature_dx ( ix, iquadpt)
      dvlx_dx ( ix, iquadpt) = tempaa * dvlx_dx ( ix, iquadpt)
      dvly_dx ( ix, iquadpt) = tempaa * dvly_dx ( ix, iquadpt)
      dvlz_dx ( ix, iquadpt) = tempaa * dvlz_dx ( ix, iquadpt)
      dbfx_dx ( ix, iquadpt) = tempaa * dbfx_dx ( ix, iquadpt)
      dbfy_dx ( ix, iquadpt) = tempaa * dbfy_dx ( ix, iquadpt)
      dbfz_dx ( ix, iquadpt) = tempaa * dbfz_dx ( ix, iquadpt)
      dv2_dx ( ix, iquadpt) = tempaa * dv2_dx ( ix, iquadpt)
      db2_dx ( ix, iquadpt) = tempaa * db2_dx ( ix, iquadpt)

      drho_dy ( ix, iquadpt) = tempab * drho_dy ( ix, iquadpt)
      dtemperature_dy ( ix, iquadpt) =
     1                 tempab * dtemperature_dy ( ix, iquadpt)
      dvlx_dy ( ix, iquadpt) = tempab * dvlx_dy ( ix, iquadpt)
      dvly_dy ( ix, iquadpt) = tempab * dvly_dy ( ix, iquadpt)
      dvlz_dy ( ix, iquadpt) = tempab * dvlz_dy ( ix, iquadpt)
      dbfx_dy ( ix, iquadpt) = tempab * dbfx_dy ( ix, iquadpt)
      dbfy_dy ( ix, iquadpt) = tempab * dbfy_dy ( ix, iquadpt)
      dbfz_dy ( ix, iquadpt) = tempab * dbfz_dy ( ix, iquadpt)
      dv2_dy ( ix, iquadpt) = tempab * dv2_dy ( ix, iquadpt)
      db2_dy ( ix, iquadpt) = tempab * db2_dy ( ix, iquadpt)

      drho_dz ( ix, iquadpt) = tempac * drho_dz ( ix, iquadpt)
      dtemperature_dz ( ix, iquadpt) =
     1                 tempac * dtemperature_dz ( ix, iquadpt)
      dvlx_dz ( ix, iquadpt) = tempac * dvlx_dz ( ix, iquadpt)
      dvly_dz ( ix, iquadpt) = tempac * dvly_dz ( ix, iquadpt)
      dvlz_dz ( ix, iquadpt) = tempac * dvlz_dz ( ix, iquadpt)
      dbfx_dz ( ix, iquadpt) = tempac * dbfx_dz ( ix, iquadpt)
      dbfy_dz ( ix, iquadpt) = tempac * dbfy_dz ( ix, iquadpt)
      dbfz_dz ( ix, iquadpt) = tempac * dbfz_dz ( ix, iquadpt)
      dv2_dz ( ix, iquadpt) = tempac * dv2_dz ( ix, iquadpt)
      db2_dz ( ix, iquadpt) = tempac * db2_dz ( ix, iquadpt)

      div_vl ( ix, iquadpt) = dvlx_dx ( ix, iquadpt)
     1                      + dvly_dy ( ix, iquadpt)
     1                      + dvlz_dz ( ix, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = rho ( ix, iquadpt) * viscosity_coef

! ----------

      fni_quadpt ( ix, 1, iquadpt) = 0.0

      fni_quadpt ( ix, 2, iquadpt) =
     1               2.0 * tempaa * ( dvlx_dx ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      fni_quadpt ( ix, 3, iquadpt) =
     1                     tempaa * ( dvlx_dy ( ix, iquadpt)
     1                              + dvly_dx ( ix, iquadpt) )

      fni_quadpt ( ix, 4, iquadpt) =
     1                     tempaa * ( dvlx_dz ( ix, iquadpt)
     1                              + dvlz_dx ( ix, iquadpt) )

      fni_quadpt ( ix, 5, iquadpt) =
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vlx ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvlx_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvlx_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvlx_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dx ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dx ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfx_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfx_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfx_dz ( ix, iquadpt) )

      fni_quadpt ( ix, 6, iquadpt) = 0.0

      fni_quadpt ( ix, 7, iquadpt) =
     1   resistivity_coef * ( dbfy_dx ( ix, iquadpt)
     1                      - dbfx_dy ( ix, iquadpt) )

      fni_quadpt ( ix, 8, iquadpt) =
     1   resistivity_coef * ( dbfz_dx ( ix, iquadpt)
     1                      - dbfx_dz ( ix, iquadpt) )

#if ( NVAR > 8)
      fni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      gni_quadpt ( ix, 1, iquadpt) = 0.0

      gni_quadpt ( ix, 2, iquadpt) =
     1                     tempaa * ( dvly_dx ( ix, iquadpt)
     1                              + dvlx_dy ( ix, iquadpt) )

      gni_quadpt ( ix, 3, iquadpt) =
     1               2.0 * tempaa * ( dvly_dy ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      gni_quadpt ( ix, 4, iquadpt) =
     1                     tempaa * ( dvly_dz ( ix, iquadpt)
     1                              + dvlz_dy ( ix, iquadpt) )

      gni_quadpt ( ix, 5, iquadpt) = 
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vly ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvly_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvly_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvly_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dy ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dy ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfy_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfy_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfy_dz ( ix, iquadpt) )

      gni_quadpt ( ix, 6, iquadpt) =
     1   resistivity_coef * ( dbfx_dy ( ix, iquadpt)
     1                      - dbfy_dx ( ix, iquadpt) )

      gni_quadpt ( ix, 7, iquadpt) = 0.0

      gni_quadpt ( ix, 8, iquadpt) =
     1   resistivity_coef * ( dbfz_dy ( ix, iquadpt)
     1                      - dbfy_dz ( ix, iquadpt) )

#if ( NVAR > 8)
      gni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      hni_quadpt ( ix, 1, iquadpt) = 0.0

      hni_quadpt ( ix, 2, iquadpt) =
     1                     tempaa * ( dvlz_dx ( ix, iquadpt)
     1                              + dvlx_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 3, iquadpt) =
     1                     tempaa * ( dvlz_dy ( ix, iquadpt)
     1                              + dvly_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 4, iquadpt) =
     1               2.0 * tempaa * ( dvlz_dz ( ix, iquadpt)
     1                         - ob3 * div_vl ( ix, iquadpt) )

      hni_quadpt ( ix, 5, iquadpt) =
     1   tempaa * ( - 2.0 * ob3 * div_vl ( ix, iquadpt)
     1                             * vlz ( ix, iquadpt)
     1            + vlx ( ix, iquadpt) * dvlz_dx ( ix, iquadpt)
     1            + vly ( ix, iquadpt) * dvlz_dy ( ix, iquadpt)
     1            + vlz ( ix, iquadpt) * dvlz_dz ( ix, iquadpt)
     1            + 0.5 * dv2_dz ( ix, iquadpt) )
     1 + resistivity_coef * ( 0.5 * db2_dz ( ix, iquadpt)
     1            - bfx ( ix, iquadpt) * dbfz_dx ( ix, iquadpt)
     1            - bfy ( ix, iquadpt) * dbfz_dy ( ix, iquadpt)
     1            - bfz ( ix, iquadpt) * dbfz_dz ( ix, iquadpt) )

      hni_quadpt ( ix, 6, iquadpt) =
     1   resistivity_coef * ( dbfx_dz ( ix, iquadpt)
     1                      - dbfz_dx ( ix, iquadpt) )

      hni_quadpt ( ix, 7, iquadpt) =
     1   resistivity_coef * ( dbfy_dz ( ix, iquadpt)
     1                      - dbfz_dy ( ix, iquadpt) )

      hni_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NVAR > 8)
      hni_quadpt ( ix, 9:NVAR , iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      f_quadpt ( ix, :, iquadpt) = f_quadpt ( ix, :, iquadpt)
     1                       + fni_quadpt ( ix, :, iquadpt)

      g_quadpt ( ix, :, iquadpt) = g_quadpt ( ix, :, iquadpt)
     1                       + gni_quadpt ( ix, :, iquadpt)

      h_quadpt ( ix, :, iquadpt) = h_quadpt ( ix, :, iquadpt)
     1                       + hni_quadpt ( ix, :, iquadpt)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

! ------------------------------

      END IF
! End of "( viscosity_coef .GT. 0.0) .OR. ( resistivity_coef .GT. 0.0)".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tempaa = dtmove / dx ( ix)
      f_quadpt ( ix, :, NUM_ZONECEN_BASES+1: NUM_SPACETIME_BASES) =
     1   tempaa
     1 * f_quadpt ( ix, :, NUM_ZONECEN_BASES+1: NUM_SPACETIME_BASES)

      tempab = dtmove / dy ( ix)
      g_quadpt ( ix, :, NUM_ZONECEN_BASES+1: NUM_SPACETIME_BASES) =
     1   tempab
     1 * g_quadpt ( ix, :, NUM_ZONECEN_BASES+1: NUM_SPACETIME_BASES)

      tempac = dtmove / dz ( ix)
      h_quadpt ( ix, :, NUM_ZONECEN_BASES+1: NUM_SPACETIME_BASES) =
     1   tempac
     1 * h_quadpt ( ix, :, NUM_ZONECEN_BASES+1: NUM_SPACETIME_BASES)

! ----------

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 9) Build source terms "s_quadpt" at the t>0 quadrature points.
! Remember that the primitive variables are still available from the previous
! loop at each of the t=0 quadrature points for constructing the source terms.

! Default is to build source terms at t>0 . In some ADER-CG updates
! (on Cartesian meshes), we may not require any source terms.
#if ( BUILD_SOURCE_TERMS == 1)
! ----------------------------------------------------------------------

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      s_quadpt ( oned_begin - oned_off: oned_end + oned_off,
     1                                           :, iquadpt) = 0.0

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

#if ( 0 == 1)

! For single-fluid adiabatic MHD turbulence, we may want to provide a gradual
! (explicit) heating/cooling term which ensures that pressure always
! tends to a mean pressure "prs_bar" in a time "tau_bar" (comparable
! to a turnover time). Such a source term can be implemented explicitly
! and its implementation is given below.
! Remember to declare "prs_bar, tau_bar".

      prs_bar = 1.0
      tau_bar = 1.0

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 5, iquadpt) = - ( prs ( ix, iquadpt) - prs_bar)
     1     / ( tau_bar * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0) )

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

#endif

! ------------------------------

! Build two-fluid source terms. In this part, the source terms are
! hard-coded for one family of ionized fluid and one family of neutrals.

#if ( N_EULER == 1) 

      IF ( ( igeom .EQ. 1) .AND. ( N_EULER .EQ. 1) .AND.
     1     ( ambipolar_diffusion_drag_coef .GT. 0.0)) THEN

      ivar = NFIELD + 1 + NFLUID

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 1, iquadpt) = 0.0

      tempaa = ambipolar_diffusion_drag_coef
     1       * rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

      s_quadpt ( ix, 2, iquadpt) = tempaa
     1           * ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt))

      s_quadpt ( ix, 3, iquadpt) = tempaa
     1           * ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt))

      s_quadpt ( ix, 4, iquadpt) = tempaa
     1           * ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt))

      s_quadpt ( ix, 5, iquadpt) = ( ambipolar_diffusion_drag_coef
     1   / ( molewt_ion + molewt_neut) )

     1 * ( rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

     1   * ( ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt) )
     1     * ( molewt_neut * vlx_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlx ( ix, iquadpt) )
     1     + ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt) )
     1     * ( molewt_neut * vly_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vly ( ix, iquadpt) )
     1     + ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt) )
     1     * ( molewt_neut * vlz_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlz ( ix, iquadpt) )
     1     )

     1   + 3.0 * ( molewt_neut * prs_euler ( ix, iquadpt, 1)
     1           * rho ( ix, iquadpt)
     1           - molewt_ion * prs ( ix, iquadpt)
     1           * rho_euler ( ix, iquadpt, 1) )
     1   )

      s_quadpt ( ix, 6, iquadpt) = 0.0
      s_quadpt ( ix, 7, iquadpt) = 0.0
      s_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      s_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) = 0.0
#endif


      s_quadpt ( ix, ivar + 1, iquadpt) = 0.0
      s_quadpt ( ix, ivar + 2, iquadpt) = - s_quadpt ( ix, 2, iquadpt)
      s_quadpt ( ix, ivar + 3, iquadpt) = - s_quadpt ( ix, 3, iquadpt)
      s_quadpt ( ix, ivar + 4, iquadpt) = - s_quadpt ( ix, 4, iquadpt)
      s_quadpt ( ix, ivar + 5, iquadpt) = - s_quadpt ( ix, 5, iquadpt)

#if ( NFLUID_EULER >= 1)
      s_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

      END IF
! End of "ambipolar_diffusion_drag_coef .GT. 0.0" part

#endif
! End of "N_EULER == 1"

! ------------------------------

! Cylindrical Geometry

      IF ( igeom .EQ. 2) THEN

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1                           * vlx ( ix, iquadpt) * obx1

      s_quadpt ( ix, 2, iquadpt) =
     1 ( rho ( ix, iquadpt) * vly ( ix, iquadpt)**2
     1 - rho ( ix, iquadpt) * vlx ( ix, iquadpt)**2
     1 + bfx ( ix, iquadpt)**2 * ob4pi
     1 - bfy ( ix, iquadpt)**2 * ob4pi) * obx1

      s_quadpt ( ix, 3, iquadpt) =
     1 - 2.0 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                              * vly ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                        * vlz ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 5, iquadpt) = - ( tempac * vlx ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = 0.0

      s_quadpt ( ix, 7, iquadpt) = 0.0

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Spherical Geometry

      IF ( igeom .EQ. 3) THEN

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      cos_x2 = cos ( x2)
      sin_x2 = sin ( x2)
      cot_x2 = cos_x2 / sin_x2

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 2.0 * vlx ( ix, iquadpt)
     1 + cot_x2 * vly ( ix, iquadpt) ) * obx1

      s_quadpt ( ix, 2, iquadpt) = rho ( ix, iquadpt)
     1 * ( - 2.0 * vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1   + vlz ( ix, iquadpt)**2) * obx1
     1 - cot_x2
     1 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1   - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1
     1 - ( - 2.0 * bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1   + bfz ( ix, iquadpt)**2) * ob4pi * obx1

      s_quadpt ( ix, 3, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1 + cot_x2
     1 * ( vly ( ix, iquadpt)**2 - vlz ( ix, iquadpt)**2)) * obx1
     1 + ( 3.0 * bfx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1   + cot_x2
     1   * ( bfy ( ix, iquadpt)**2 - bfz ( ix, iquadpt)**2))
     1 * ob4pi * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - rho ( ix, iquadpt) * vlz ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) + 2.0 * cot_x2
     1   * vly ( ix, iquadpt) ) * obx1
     1 + bfz ( ix, iquadpt)
     1 * ( 3.0 * bfx ( ix, iquadpt) + 2.0 * cot_x2
     1   * bfy ( ix, iquadpt) ) * ob4pi * obx1

      s_quadpt ( ix, 5, iquadpt) =
     1 - 2.0 * ( tempac * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab) * obx1
     1 - cot_x2 * ( tempac * vly ( ix, iquadpt)
     1 - bfy ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = - cot_x2
     1                           * g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 7, iquadpt) = g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Multiply "s_quadpt" by "dtmove".

      DO ix = oned_begin - oned_off, oned_end + oned_off
      s_quadpt ( ix, :, NUM_ZONECEN_BASES + 1: NUM_SPACETIME_BASES)=
     1 dtmove
     1*s_quadpt ( ix, :, NUM_ZONECEN_BASES + 1: NUM_SPACETIME_BASES)
      END DO

! ----------------------------------------------------------------------
#endif
! End of "BUILD_SOURCE_TERMS == 1"
! ----------------------------------------------------------------------

! STEP 10) By this point, the nodal flux coefficients
! ( f_quadpt, g_quadpt, h_quadpt) and the nodal source coefficients
! "s_quadpt" that enable us to define the fluxes and source terms
! any where in space and time within the zone have already been loaded.
! Now we can use the analytical inversion of the matrix to obtain
! the next iterate for the conserved variables "u_quadpt" . The nodal
! conserved variables "u_quadpt" also enable us to construct the conserved
! variables at any space-time point in the zone.

! ----------------------------------------------------------------------

#if ( NUMORDER == 2)
! cg, 2nd order : Build "u_quadpt" at t > 0

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 05) = tu ( ix, :, 01)
     1 + 2.0 * s_quadpt ( ix, :, 05) / 3.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2".

! ----------------------------------------------------------------------

#if ( NUMORDER == 3)
! cg, 3rd order : Build "u_quadpt" at t > 0

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      tu ( ix, :, 01) =
     1   0.5 * ( f_quadpt ( ix, :, 11) - f_quadpt ( ix, :, 12) )
     1 + 0.5 * ( g_quadpt ( ix, :, 11) - g_quadpt ( ix, :, 13) )
     1 + 0.5 * ( h_quadpt ( ix, :, 11) - h_quadpt ( ix, :, 14) )

      tu ( ix, :, 02) = tu ( ix, :, 01) + ( -4.0 * s_quadpt ( ix, :, 11)
     1 - 3.0 * s_quadpt ( ix, :, 15) ) / 24.0

      u_quadpt ( ix, :, 11) = tu ( ix, :, 11) + tu ( ix, :, 01)
     1 + 0.125 * ( 4.0 * s_quadpt ( ix, :, 11) - s_quadpt ( ix, :, 15) )

      u_quadpt ( ix, :, 12) = tu ( ix, :, 12) + tu ( ix, :, 02)
     1 + 2.0 * s_quadpt ( ix, :, 12) / 3.0

      u_quadpt ( ix, :, 13) = tu ( ix, :, 13) + tu ( ix, :, 02)
     1 + 2.0 * s_quadpt ( ix, :, 13) / 3.0

      u_quadpt ( ix, :, 14) = tu ( ix, :, 14) + tu ( ix, :, 02)
     1 + 2.0 * s_quadpt ( ix, :, 14) / 3.0

      u_quadpt ( ix, :, 15) = tu ( ix, :, 15) + 4.0 * tu ( ix, :, 01)
     1 + 0.1 * ( 8.0 * s_quadpt ( ix, :, 11) + s_quadpt ( ix, :, 15) )

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3".

! ----------------------------------------------------------------------

#if ( NUMORDER == 4)
! cg, 4th order : Build "u_quadpt" at t > 0

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 21) = tu ( ix, :, 21)

     1 + ( 23940 * f_quadpt ( ix, :, 21)
     1 - 59850 * f_quadpt ( ix, :, 22)
     1 + 35910 * f_quadpt ( ix, :, 23) 
     1 - 11970 * f_quadpt ( ix, :, 31)
     1 + 11970 * f_quadpt ( ix, :, 32)

     1 + 23940 * g_quadpt ( ix, :, 21)
     1 - 59850 * g_quadpt ( ix, :, 24)
     1 + 35910 * g_quadpt ( ix, :, 25)
     1 - 11970 * g_quadpt ( ix, :, 31)
     1 + 11970 * g_quadpt ( ix, :, 33)

     1 + 23940 * h_quadpt ( ix, :, 21)
     1 - 59850 * h_quadpt ( ix, :, 26)
     1 + 35910 * h_quadpt ( ix, :, 27)
     1 - 11970 * h_quadpt ( ix, :, 31)
     1 + 11970 * h_quadpt ( ix, :, 34)

     1 + 189468 * s_quadpt ( ix, :, 21)
     1 - 17955 * s_quadpt ( ix, :, 22)
     1 - 17955 * s_quadpt ( ix, :, 23)
     1 - 17955 * s_quadpt ( ix, :, 24)
     1 - 17955 * s_quadpt ( ix, :, 25)
     1 - 17955 * s_quadpt ( ix, :, 26)
     1 - 17955 * s_quadpt ( ix, :, 27)
     1 - 39843 * s_quadpt ( ix, :, 31)
     1 + 11286 * s_quadpt ( ix, :, 35) ) / 215460.0

      u_quadpt ( ix, :, 22) = tu ( ix, :, 22)

     1 + ( 502740 * f_quadpt ( ix, :, 21)
     1 - 395010 * f_quadpt ( ix, :, 22)
     1 - 107730 * f_quadpt ( ix, :, 23)
     1 - 35910 * f_quadpt ( ix, :, 31)
     1 + 35910 * f_quadpt ( ix, :, 32)

     1 - 143640 * g_quadpt ( ix, :, 21)
     1 + 215460 * g_quadpt ( ix, :, 22)
     1 + 35910 * g_quadpt ( ix, :, 24)
     1 + 107730 * g_quadpt ( ix, :, 25)
     1 - 215460 * g_quadpt ( ix, :, 30)
     1 - 35910 * g_quadpt ( ix, :, 31)
     1 + 35910 * g_quadpt ( ix, :, 33)

     1 - 143640 * h_quadpt ( ix, :, 21)
     1 + 215460 * h_quadpt ( ix, :, 22)
     1 + 35910 * h_quadpt ( ix, :, 26)
     1 + 107730 * h_quadpt ( ix, :, 27)
     1 - 215460 * h_quadpt ( ix, :, 29)
     1 - 35910 * h_quadpt ( ix, :, 31)
     1 + 35910 * h_quadpt ( ix, :, 34)

     1 - 164160 * s_quadpt ( ix, :, 21)
     1 + 517104 * s_quadpt ( ix, :, 22)
     1 + 107730 * s_quadpt ( ix, :, 23)
     1 - 53865 * s_quadpt ( ix, :, 24)
     1 - 53865 * s_quadpt ( ix, :, 25)
     1 - 53865 * s_quadpt ( ix, :, 26)
     1 - 53865 * s_quadpt ( ix, :, 27)
     1 + 31293 * s_quadpt ( ix, :, 31)
     1 - 150822 * s_quadpt ( ix, :, 32)
     1 + 33858 * s_quadpt ( ix, :, 35) ) / 646380.0

      u_quadpt ( ix, :, 23) = tu ( ix, :, 23)

     1 + ( - 359100 * f_quadpt ( ix, :, 21)
     1 + 35910 * f_quadpt ( ix, :, 22)
     1 + 323190 * f_quadpt ( ix, :, 23)
     1 - 35910 * f_quadpt ( ix, :, 31)
     1 + 35910 * f_quadpt ( ix, :, 32)

     1 + 287280 * g_quadpt ( ix, :, 21)
     1 - 215460 * g_quadpt ( ix, :, 22)
     1 - 395010 * g_quadpt ( ix, :, 24)
     1 + 107730 * g_quadpt ( ix, :, 25)
     1 + 215460 * g_quadpt ( ix, :, 30)
     1 - 35910 * g_quadpt ( ix, :, 31)
     1 + 35910 * g_quadpt ( ix, :, 33)

     1 + 287280 * h_quadpt ( ix, :, 21)
     1 - 215460 * h_quadpt ( ix, :, 22)
     1 - 395010 * h_quadpt ( ix, :, 26)
     1 + 107730 * h_quadpt ( ix, :, 27)
     1 + 215460 * h_quadpt ( ix, :, 29)
     1 - 35910 * h_quadpt ( ix, :, 31)
     1 + 35910 * h_quadpt ( ix, :, 34)

     1 + 439128 * s_quadpt ( ix, :, 21)
     1 - 193914 * s_quadpt ( ix, :, 22)
     1 + 215460 * s_quadpt ( ix, :, 23)
     1 - 53865 * s_quadpt ( ix, :, 24)
     1 - 53865 * s_quadpt ( ix, :, 25)
     1 - 53865 * s_quadpt ( ix, :, 26)
     1 - 53865 * s_quadpt ( ix, :, 27)
     1 - 270351 * s_quadpt ( ix, :, 31)
     1 + 150822 * s_quadpt ( ix, :, 32)
     1 + 33858 * s_quadpt ( ix, :, 35) ) / 646380.0

      u_quadpt ( ix, :, 24) = tu ( ix, :, 24)

     1 + ( - 143640 * f_quadpt ( ix, :, 21)
     1 + 35910 * f_quadpt ( ix, :, 22)
     1 + 107730 * f_quadpt ( ix, :, 23)
     1 + 215460 * f_quadpt ( ix, :, 24)
     1 - 215460 * f_quadpt ( ix, :, 30)
     1 - 35910 * f_quadpt ( ix, :, 31)
     1 + 35910 * f_quadpt ( ix, :, 32)

     1 + 502740 * g_quadpt ( ix, :, 21)
     1 - 395010 * g_quadpt ( ix, :, 24)
     1 - 107730 * g_quadpt ( ix, :, 25)
     1 - 35910 * g_quadpt ( ix, :, 31)
     1 + 35910 * g_quadpt ( ix, :, 33)

     1 - 143640 * h_quadpt ( ix, :, 21)
     1 + 215460 * h_quadpt ( ix, :, 24)
     1 + 35910 * h_quadpt ( ix, :, 26)
     1 + 107730 * h_quadpt ( ix, :, 27)
     1 - 215460 * h_quadpt ( ix, :, 28)
     1 - 35910 * h_quadpt ( ix, :, 31)
     1 + 35910 * h_quadpt ( ix, :, 34)

     1 - 164160 * s_quadpt ( ix, :, 21)
     1 - 53865 * s_quadpt ( ix, :, 22)
     1 - 53865 * s_quadpt ( ix, :, 23)
     1 + 517104 * s_quadpt ( ix, :, 24)
     1 + 107730 * s_quadpt ( ix, :, 25)
     1 - 53865 * s_quadpt ( ix, :, 26)
     1 - 53865 * s_quadpt ( ix, :, 27)
     1 + 31293 * s_quadpt ( ix, :, 31)
     1 - 150822 * s_quadpt ( ix, :, 33)
     1 + 33858 * s_quadpt ( ix, :, 35) ) / 646380.0

      u_quadpt ( ix, :, 25) = tu ( ix, :, 25)

     1 + ( 287280 * f_quadpt ( ix, :, 21)
     1 - 395010 * f_quadpt ( ix, :, 22)
     1 + 107730 * f_quadpt ( ix, :, 23)
     1 - 215460 * f_quadpt ( ix, :, 24)
     1 + 215460 * f_quadpt ( ix, :, 30)
     1 - 35910 * f_quadpt ( ix, :, 31)
     1 + 35910 * f_quadpt ( ix, :, 32)

     1 - 359100 * g_quadpt ( ix, :, 21)
     1 + 35910 * g_quadpt ( ix, :, 24)
     1 + 323190 * g_quadpt ( ix, :, 25)
     1 - 35910 * g_quadpt ( ix, :, 31)
     1 + 35910 * g_quadpt ( ix, :, 33)

     1 + 287280 * h_quadpt ( ix, :, 21)
     1 - 215460 * h_quadpt ( ix, :, 24)
     1 - 395010 * h_quadpt ( ix, :, 26)
     1 + 107730 * h_quadpt ( ix, :, 27)
     1 + 215460 * h_quadpt ( ix, :, 28)
     1 - 35910 * h_quadpt ( ix, :, 31)
     1 + 35910 * h_quadpt ( ix, :, 34)

     1 + 439128 * s_quadpt ( ix, :, 21)
     1 - 53865 * s_quadpt ( ix, :, 22)
     1 - 53865 * s_quadpt ( ix, :, 23)
     1 - 193914 * s_quadpt ( ix, :, 24)
     1 + 215460 * s_quadpt ( ix, :, 25)
     1 - 53865 * s_quadpt ( ix, :, 26)
     1 - 53865 * s_quadpt ( ix, :, 27)
     1 - 270351 * s_quadpt ( ix, :, 31)
     1 + 150822 * s_quadpt ( ix, :, 33)
     1 + 33858 * s_quadpt ( ix, :, 35) ) / 646380.0

      u_quadpt ( ix, :, 26) = tu ( ix, :, 26)

     1 + ( - 143640 * f_quadpt ( ix, :, 21)
     1 + 35910 * f_quadpt ( ix, :, 22)
     1 + 107730 * f_quadpt ( ix, :, 23)
     1 + 215460 * f_quadpt ( ix, :, 26)
     1 - 215460 * f_quadpt ( ix, :, 29)
     1 - 35910 * f_quadpt ( ix, :, 31)
     1 + 35910 * f_quadpt ( ix, :, 32)

     1 - 143640 * g_quadpt ( ix, :, 21)
     1 + 35910 * g_quadpt ( ix, :, 24)
     1 + 107730 * g_quadpt ( ix, :, 25)
     1 + 215460 * g_quadpt ( ix, :, 26)
     1 - 215460 * g_quadpt ( ix, :, 28)
     1 - 35910 * g_quadpt ( ix, :, 31)
     1 + 35910 * g_quadpt ( ix, :, 33)

     1 + 502740 * h_quadpt ( ix, :, 21)
     1 - 395010 * h_quadpt ( ix, :, 26)
     1 - 107730 * h_quadpt ( ix, :, 27)
     1 - 35910 * h_quadpt ( ix, :, 31)
     1 + 35910 * h_quadpt ( ix, :, 34)

     1 - 164160 * s_quadpt ( ix, :, 21)
     1 - 53865 * s_quadpt ( ix, :, 22)
     1 - 53865 * s_quadpt ( ix, :, 23)
     1 - 53865 * s_quadpt ( ix, :, 24)
     1 - 53865 * s_quadpt ( ix, :, 25)
     1 + 517104 * s_quadpt ( ix, :, 26)
     1 + 107730 * s_quadpt ( ix, :, 27)
     1 + 31293 * s_quadpt ( ix, :, 31)
     1 - 150822 * s_quadpt ( ix, :, 34)
     1 + 33858 * s_quadpt ( ix, :, 35) ) / 646380.0

      u_quadpt ( ix, :, 27) = tu ( ix, :, 27)

     1 + ( 287280 * f_quadpt ( ix, :, 21)
     1 - 395010 * f_quadpt ( ix, :, 22)
     1 + 107730 * f_quadpt ( ix, :, 23)
     1 - 215460 * f_quadpt ( ix, :, 26)
     1 + 215460 * f_quadpt ( ix, :, 29)
     1 - 35910 * f_quadpt ( ix, :, 31)
     1 + 35910 * f_quadpt ( ix, :, 32)

     1 + 287280 * g_quadpt ( ix, :, 21)
     1 - 395010 * g_quadpt ( ix, :, 24)
     1 + 107730 * g_quadpt ( ix, :, 25)
     1 - 215460 * g_quadpt ( ix, :, 26)
     1 + 215460 * g_quadpt ( ix, :, 28)
     1 - 35910 * g_quadpt ( ix, :, 31)
     1 + 35910 * g_quadpt ( ix, :, 33)

     1 - 359100 * h_quadpt ( ix, :, 21)
     1 + 35910 * h_quadpt ( ix, :, 26)
     1 + 323190 * h_quadpt ( ix, :, 27)
     1 - 35910 * h_quadpt ( ix, :, 31)
     1 + 35910 * h_quadpt ( ix, :, 34)

     1 + 439128 * s_quadpt ( ix, :, 21)
     1 - 53865 * s_quadpt ( ix, :, 22)
     1 - 53865 * s_quadpt ( ix, :, 23)
     1 - 53865 * s_quadpt ( ix, :, 24)
     1 - 53865 * s_quadpt ( ix, :, 25)
     1 - 193914 * s_quadpt ( ix, :, 26)
     1 + 215460 * s_quadpt ( ix, :, 27)
     1 - 270351 * s_quadpt ( ix, :, 31)
     1 + 150822 * s_quadpt ( ix, :, 34)
     1 + 33858 * s_quadpt ( ix, :, 35) ) / 646380.0

      u_quadpt ( ix, :, 28) = tu ( ix, :, 28)

     1 + ( - 359100 * f_quadpt ( ix, :, 21)
     1 + 251370 * f_quadpt ( ix, :, 22)
     1 + 107730 * f_quadpt ( ix, :, 23)
     1 + 215460 * f_quadpt ( ix, :, 24)
     1 + 215460 * f_quadpt ( ix, :, 26)
     1 - 215460 * f_quadpt ( ix, :, 29)
     1 - 215460 * f_quadpt ( ix, :, 30)
     1 - 35910 * f_quadpt ( ix, :, 31)
     1 + 35910 * f_quadpt ( ix, :, 32)

     1 + 287280 * g_quadpt ( ix, :, 21)
     1 - 179550 * g_quadpt ( ix, :, 24)
     1 - 107730 * g_quadpt ( ix, :, 25)
     1 + 215460 * g_quadpt ( ix, :, 26)
     1 - 215460 * g_quadpt ( ix, :, 28)
     1 - 35910 * g_quadpt ( ix, :, 31)
     1 + 35910 * g_quadpt ( ix, :, 33)

     1 + 287280 * h_quadpt ( ix, :, 21)
     1 + 215460 * h_quadpt ( ix, :, 24)
     1 - 179550 * h_quadpt ( ix, :, 26)
     1 - 107730 * h_quadpt ( ix, :, 27)
     1 - 215460 * h_quadpt ( ix, :, 28)
     1 - 35910 * h_quadpt ( ix, :, 31)
     1 + 35910 * h_quadpt ( ix, :, 34)

     1 - 465804 * s_quadpt ( ix, :, 21)
     1 - 53865 * s_quadpt ( ix, :, 22)
     1 - 53865 * s_quadpt ( ix, :, 23)
     1 + 86184 * s_quadpt ( ix, :, 24)
     1 + 107730 * s_quadpt ( ix, :, 25)
     1 + 86184 * s_quadpt ( ix, :, 26)
     1 + 107730 * s_quadpt ( ix, :, 27)
     1 + 430920 * s_quadpt ( ix, :, 28)
     1 + 182115 * s_quadpt ( ix, :, 31)
     1 - 150822 * s_quadpt ( ix, :, 33)
     1 - 150822 * s_quadpt ( ix, :, 34)
     1 + 33858 * s_quadpt ( ix, :, 35) ) / 646380.0

      u_quadpt ( ix, :, 29) = tu ( ix, :, 29)

     1 + ( 287280 * f_quadpt ( ix, :, 21)
     1 - 179550 * f_quadpt ( ix, :, 22)
     1 - 107730 * f_quadpt ( ix, :, 23)
     1 + 215460 * f_quadpt ( ix, :, 26)
     1 - 215460 * f_quadpt ( ix, :, 29)
     1 - 35910 * f_quadpt ( ix, :, 31)
     1 + 35910 * f_quadpt ( ix, :, 32)

     1 - 359100 * g_quadpt ( ix, :, 21)
     1 + 215460 * g_quadpt ( ix, :, 22)
     1 + 251370 * g_quadpt ( ix, :, 24)
     1 + 107730 * g_quadpt ( ix, :, 25)
     1 + 215460 * g_quadpt ( ix, :, 26)
     1 - 215460 * g_quadpt ( ix, :, 28)
     1 - 215460 * g_quadpt ( ix, :, 30)
     1 - 35910 * g_quadpt ( ix, :, 31)
     1 + 35910 * g_quadpt ( ix, :, 33)

     1 + 287280 * h_quadpt ( ix, :, 21)
     1 + 215460 * h_quadpt ( ix, :, 22)
     1 - 179550 * h_quadpt ( ix, :, 26)
     1 - 107730 * h_quadpt ( ix, :, 27)
     1 - 215460 * h_quadpt ( ix, :, 29)
     1 - 35910 * h_quadpt ( ix, :, 31)
     1 + 35910 * h_quadpt ( ix, :, 34)

     1 - 465804 * s_quadpt ( ix, :, 21)
     1 + 86184 * s_quadpt ( ix, :, 22)
     1 + 107730 * s_quadpt ( ix, :, 23)
     1 - 53865 * s_quadpt ( ix, :, 24)
     1 - 53865 * s_quadpt ( ix, :, 25)
     1 + 86184 * s_quadpt ( ix, :, 26)
     1 + 107730 * s_quadpt ( ix, :, 27)
     1 + 430920 * s_quadpt ( ix, :, 29)
     1 + 182115 * s_quadpt ( ix, :, 31)
     1 - 150822 * s_quadpt ( ix, :, 32)
     1 - 150822 * s_quadpt ( ix, :, 34)
     1 + 33858 * s_quadpt ( ix, :, 35) ) / 646380.0

      u_quadpt ( ix, :, 30) = tu ( ix, :, 30)

     1 + ( 287280 * f_quadpt ( ix, :, 21)
     1 - 179550 * f_quadpt ( ix, :, 22)
     1 - 107730 * f_quadpt ( ix, :, 23)
     1 + 215460 * f_quadpt ( ix, :, 24)
     1 - 215460 * f_quadpt ( ix, :, 30)
     1 - 35910 * f_quadpt ( ix, :, 31)
     1 + 35910 * f_quadpt ( ix, :, 32)

     1 + 287280 * g_quadpt ( ix, :, 21)
     1 + 215460 * g_quadpt ( ix, :, 22)
     1 - 179550 * g_quadpt ( ix, :, 24)
     1 - 107730 * g_quadpt ( ix, :, 25)
     1 - 215460 * g_quadpt ( ix, :, 30)
     1 - 35910 * g_quadpt ( ix, :, 31)
     1 + 35910 * g_quadpt ( ix, :, 33)

     1 - 359100 * h_quadpt ( ix, :, 21)
     1 + 215460 * h_quadpt ( ix, :, 22)
     1 + 215460 * h_quadpt ( ix, :, 24)
     1 + 251370 * h_quadpt ( ix, :, 26)
     1 + 107730 * h_quadpt ( ix, :, 27)
     1 - 215460 * h_quadpt ( ix, :, 28)
     1 - 215460 * h_quadpt ( ix, :, 29)
     1 - 35910 * h_quadpt ( ix, :, 31)
     1 + 35910 * h_quadpt ( ix, :, 34)

     1 - 465804 * s_quadpt ( ix, :, 21)
     1 + 86184 * s_quadpt ( ix, :, 22)
     1 + 107730 * s_quadpt ( ix, :, 23)
     1 + 86184 * s_quadpt ( ix, :, 24)
     1 + 107730 * s_quadpt ( ix, :, 25)
     1 - 53865 * s_quadpt ( ix, :, 26)
     1 - 53865 * s_quadpt ( ix, :, 27)
     1 + 430920 * s_quadpt ( ix, :, 30)
     1 + 182115 * s_quadpt ( ix, :, 31)
     1 - 150822 * s_quadpt ( ix, :, 32)
     1 - 150822 * s_quadpt ( ix, :, 33)
     1 + 33858 * s_quadpt ( ix, :, 35)  ) / 646380.0

      u_quadpt ( ix, :, 31) = tu ( ix, :, 31)

     1 + ( - 23940 * f_quadpt ( ix, :, 21)
     1 - 11970 * f_quadpt ( ix, :, 22)
     1 + 35910 * f_quadpt ( ix, :, 23)
     1 + 11970 * f_quadpt ( ix, :, 31)
     1 - 11970 * f_quadpt ( ix, :, 32)

     1 - 23940 * g_quadpt ( ix, :, 21)
     1 - 11970 * g_quadpt ( ix, :, 24)
     1 + 35910 * g_quadpt ( ix, :, 25)
     1 + 11970 * g_quadpt ( ix, :, 31)
     1 - 11970 * g_quadpt ( ix, :, 33)

     1 - 23940 * h_quadpt ( ix, :, 21)
     1 - 11970 * h_quadpt ( ix, :, 26)
     1 + 35910 * h_quadpt ( ix, :, 27)
     1 + 11970 * h_quadpt ( ix, :, 31)
     1 - 11970 * h_quadpt ( ix, :, 34)

     1 + 63954 * s_quadpt ( ix, :, 21)
     1 - 5985 * s_quadpt ( ix, :, 22)
     1 - 5985 * s_quadpt ( ix, :, 23)
     1 - 5985 * s_quadpt ( ix, :, 24)
     1 - 5985 * s_quadpt ( ix, :, 25)
     1 - 5985 * s_quadpt ( ix, :, 26)
     1 - 5985 * s_quadpt ( ix, :, 27)
     1 + 1881 * s_quadpt ( ix, :, 31)
     1 + 1368 * s_quadpt ( ix, :, 35)  ) / 53865.0

      u_quadpt ( ix, :, 32) = tu ( ix, :, 32)

     1 + ( 718200 * f_quadpt ( ix, :, 21)
     1 - 502740 * f_quadpt ( ix, :, 22)
     1 - 215460 * f_quadpt ( ix, :, 23)
     1 + 71820 * f_quadpt ( ix, :, 31)
     1 - 71820 * f_quadpt ( ix, :, 32)

     1 - 574560 * g_quadpt ( ix, :, 21)
     1 + 430920 * g_quadpt ( ix, :, 22)
     1 + 359100 * g_quadpt ( ix, :, 24)
     1 + 215460 * g_quadpt ( ix, :, 25)
     1 - 430920 * g_quadpt ( ix, :, 30)
     1 + 71820 * g_quadpt ( ix, :, 31)
     1 - 71820 * g_quadpt ( ix, :, 33)

     1 - 574560 * h_quadpt ( ix, :, 21)
     1 + 430920 * h_quadpt ( ix, :, 22)
     1 + 359100 * h_quadpt ( ix, :, 26)
     1 + 215460 * h_quadpt ( ix, :, 27)
     1 - 430920 * h_quadpt ( ix, :, 29)
     1 + 71820 * h_quadpt ( ix, :, 31)
     1 - 71820 * h_quadpt ( ix, :, 34)

     1 - 68742 * s_quadpt ( ix, :, 21)
     1 + 308826 * s_quadpt ( ix, :, 22)
     1 + 71820 * s_quadpt ( ix, :, 23)
     1 - 35910 * s_quadpt ( ix, :, 24)
     1 - 35910 * s_quadpt ( ix, :, 25)
     1 - 35910 * s_quadpt ( ix, :, 26)
     1 - 35910 * s_quadpt ( ix, :, 27)
     1 + 22059 * s_quadpt ( ix, :, 31)
     1 - 10773 * s_quadpt ( ix, :, 32)
     1 + 8208 * s_quadpt ( ix, :, 35) ) / 323190.0

      u_quadpt ( ix, :, 33) = tu ( ix, :, 33)

     1 + ( - 574560 * f_quadpt ( ix, :, 21)
     1 + 359100 * f_quadpt ( ix, :, 22)
     1 + 215460 * f_quadpt ( ix, :, 23)
     1 + 430920 * f_quadpt ( ix, :, 24)
     1 - 430920 * f_quadpt ( ix, :, 30)
     1 + 71820 * f_quadpt ( ix, :, 31)
     1 - 71820 * f_quadpt ( ix, :, 32)

     1 + 718200 * g_quadpt ( ix, :, 21)
     1 - 502740 * g_quadpt ( ix, :, 24)
     1 - 215460 * g_quadpt ( ix, :, 25)
     1 + 71820 * g_quadpt ( ix, :, 31)
     1 - 71820 * g_quadpt ( ix, :, 33)

     1 - 574560 * h_quadpt ( ix, :, 21)
     1 + 430920 * h_quadpt ( ix, :, 24)
     1 + 359100 * h_quadpt ( ix, :, 26)
     1 + 215460 * h_quadpt ( ix, :, 27)
     1 - 430920 * h_quadpt ( ix, :, 28)
     1 + 71820 * h_quadpt ( ix, :, 31)
     1 - 71820 * h_quadpt ( ix, :, 34)

     1 - 68742 * s_quadpt ( ix, :, 21)
     1 - 35910 * s_quadpt ( ix, :, 22)
     1 - 35910 * s_quadpt ( ix, :, 23)
     1 + 308826 * s_quadpt ( ix, :, 24)
     1 + 71820 * s_quadpt ( ix, :, 25)
     1 - 35910 * s_quadpt ( ix, :, 26)
     1 - 35910 * s_quadpt ( ix, :, 27)
     1 + 22059 * s_quadpt ( ix, :, 31)
     1 - 10773 * s_quadpt ( ix, :, 33)
     1 + 8208 * s_quadpt ( ix, :, 35) ) / 323190.0

      u_quadpt ( ix, :, 34) = tu ( ix, :, 34)

     1 + ( - 574560 * f_quadpt ( ix, :, 21)
     1 + 359100 * f_quadpt ( ix, :, 22)
     1 + 215460 * f_quadpt ( ix, :, 23)
     1 + 430920 * f_quadpt ( ix, :, 26)
     1 - 430920 * f_quadpt ( ix, :, 29)
     1 + 71820 * f_quadpt ( ix, :, 31)
     1 - 71820 * f_quadpt ( ix, :, 32)

     1 - 574560 * g_quadpt ( ix, :, 21)
     1 + 359100 * g_quadpt ( ix, :, 24)
     1 + 215460 * g_quadpt ( ix, :, 25)
     1 + 430920 * g_quadpt ( ix, :, 26)
     1 - 430920 * g_quadpt ( ix, :, 28)
     1 + 71820 * g_quadpt ( ix, :, 31)
     1 - 71820 * g_quadpt ( ix, :, 33)

     1 + 718200 * h_quadpt ( ix, :, 21)
     1 - 502740 * h_quadpt ( ix, :, 26)
     1 - 215460 * h_quadpt ( ix, :, 27)
     1 + 71820 * h_quadpt ( ix, :, 31)
     1 - 71820 * h_quadpt ( ix, :, 34)

     1 - 68742 * s_quadpt ( ix, :, 21)
     1 - 35910 * s_quadpt ( ix, :, 22)
     1 - 35910 * s_quadpt ( ix, :, 23)
     1 - 35910 * s_quadpt ( ix, :, 24)
     1 - 35910 * s_quadpt ( ix, :, 25)
     1 + 308826 * s_quadpt ( ix, :, 26)
     1 + 71820 * s_quadpt ( ix, :, 27)
     1 + 22059 * s_quadpt ( ix, :, 31)
     1 - 10773 * s_quadpt ( ix, :, 34)
     1 + 8208 * s_quadpt ( ix, :, 35) ) / 323190.0

      u_quadpt ( ix, :, 35) = tu ( ix, :, 35)

     1 + ( - 215460 * f_quadpt ( ix, :, 21)
     1 + 107730 * f_quadpt ( ix, :, 22)
     1 + 107730 * f_quadpt ( ix, :, 23)
     1 + 107730 * f_quadpt ( ix, :, 31)
     1 - 107730 * f_quadpt ( ix, :, 32)

     1 - 215460 * g_quadpt ( ix, :, 21)
     1 + 107730 * g_quadpt ( ix, :, 24)
     1 + 107730 * g_quadpt ( ix, :, 25)
     1 + 107730 * g_quadpt ( ix, :, 31)
     1 - 107730 * g_quadpt ( ix, :, 33)

     1 - 215460 * h_quadpt ( ix, :, 21)
     1 + 107730 * h_quadpt ( ix, :, 26)
     1 + 107730 * h_quadpt ( ix, :, 27)
     1 + 107730 * h_quadpt ( ix, :, 31)
     1 - 107730 * h_quadpt ( ix, :, 34)

     1 + 69768 * s_quadpt ( ix, :, 21)
     1 - 5985 * s_quadpt ( ix, :, 22)
     1 - 5985 * s_quadpt ( ix, :, 23)
     1 - 5985 * s_quadpt ( ix, :, 24)
     1 - 5985 * s_quadpt ( ix, :, 25)
     1 - 5985 * s_quadpt ( ix, :, 26)
     1 - 5985 * s_quadpt ( ix, :, 27)
     1 + 20007 * s_quadpt ( ix, :, 31)
     1 + 11286 * s_quadpt ( ix, :, 35) ) / 71820.0

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO
! End of "DO iter_ader =" loop.

! End the outer ADER iteration here.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 11) In order to build a more accurate version of the source terms,
! the t > 0 source terms, i.e. "s_quadpt", have to be rebuilt here.

#if ( BUILD_SOURCE_TERMS == 1)
! ----------------------------------------------------------------------

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      s_quadpt ( oned_begin - oned_off: oned_end + oned_off,
     1                                           :, iquadpt) = 0.0

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

#if ( 0 == 1)

! For single-fluid adiabatic MHD turbulence, we may want to provide a gradual
! (explicit) heating/cooling term which ensures that pressure always
! tends to a mean pressure "prs_bar" in a time "tau_bar" (comparable
! to a turnover time). Such a source term can be implemented explicitly
! and its implementation is given below.
! Remember to declare "prs_bar, tau_bar".

      prs_bar = 1.0
      tau_bar = 1.0

      IF ( igeom .EQ. 1) THEN

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 5, iquadpt) = - ( prs ( ix, iquadpt) - prs_bar)
     1     / ( tau_bar * ( gamma_eos_quadpt ( ix, 1, iquadpt) - 1.0) )

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

#endif

! ------------------------------

! Build two-fluid source terms. In this part, the source terms are
! hard-coded for one family of ionized fluid and one family of neutrals.

#if ( N_EULER == 1) 

      IF ( ( igeom .EQ. 1) .AND. ( N_EULER .EQ. 1) .AND.
     1     ( ambipolar_diffusion_drag_coef .GT. 0.0)) THEN

      ivar = NFIELD + 1 + NFLUID

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      s_quadpt ( ix, 1, iquadpt) = 0.0

      tempaa = ambipolar_diffusion_drag_coef
     1       * rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

      s_quadpt ( ix, 2, iquadpt) = tempaa
     1           * ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt))

      s_quadpt ( ix, 3, iquadpt) = tempaa
     1           * ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt))

      s_quadpt ( ix, 4, iquadpt) = tempaa
     1           * ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt))

      s_quadpt ( ix, 5, iquadpt) = ( ambipolar_diffusion_drag_coef
     1   / ( molewt_ion + molewt_neut) )

     1 * ( rho ( ix, iquadpt) * rho_euler ( ix, iquadpt, 1)

     1   * ( ( vlx_euler ( ix, iquadpt, 1) - vlx ( ix, iquadpt) )
     1     * ( molewt_neut * vlx_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlx ( ix, iquadpt) )
     1     + ( vly_euler ( ix, iquadpt, 1) - vly ( ix, iquadpt) )
     1     * ( molewt_neut * vly_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vly ( ix, iquadpt) )
     1     + ( vlz_euler ( ix, iquadpt, 1) - vlz ( ix, iquadpt) )
     1     * ( molewt_neut * vlz_euler ( ix, iquadpt, 1)
     1       + molewt_ion * vlz ( ix, iquadpt) )
     1     )

     1   + 3.0 * ( molewt_neut * prs_euler ( ix, iquadpt, 1)
     1           * rho ( ix, iquadpt)
     1           - molewt_ion * prs ( ix, iquadpt)
     1           * rho_euler ( ix, iquadpt, 1) )
     1   )

      s_quadpt ( ix, 6, iquadpt) = 0.0
      s_quadpt ( ix, 7, iquadpt) = 0.0
      s_quadpt ( ix, 8, iquadpt) = 0.0

#if ( NFLUID >= 1)
      s_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) = 0.0
#endif


      s_quadpt ( ix, ivar + 1, iquadpt) = 0.0
      s_quadpt ( ix, ivar + 2, iquadpt) = - s_quadpt ( ix, 2, iquadpt)
      s_quadpt ( ix, ivar + 3, iquadpt) = - s_quadpt ( ix, 3, iquadpt)
      s_quadpt ( ix, ivar + 4, iquadpt) = - s_quadpt ( ix, 4, iquadpt)
      s_quadpt ( ix, ivar + 5, iquadpt) = - s_quadpt ( ix, 5, iquadpt)

#if ( NFLUID_EULER >= 1)
      s_quadpt ( ix, ivar + NFIELD_EULER + 1:
     1           ivar + NFIELD_EULER + NFLUID_EULER, iquadpt) = 0.0
#endif

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO iquadpt ="

      END IF
! End of "ambipolar_diffusion_drag_coef .GT. 0.0" part

#endif
! End of "N_EULER == 1"

! ------------------------------

! Cylindrical Geometry

      IF ( igeom .EQ. 2) THEN

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1                           * vlx ( ix, iquadpt) * obx1

      s_quadpt ( ix, 2, iquadpt) =
     1 ( rho ( ix, iquadpt) * vly ( ix, iquadpt)**2
     1 - rho ( ix, iquadpt) * vlx ( ix, iquadpt)**2
     1 + bfx ( ix, iquadpt)**2 * ob4pi
     1 - bfy ( ix, iquadpt)**2 * ob4pi) * obx1

      s_quadpt ( ix, 3, iquadpt) =
     1 - 2.0 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                              * vly ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - ( rho ( ix, iquadpt) * vlx ( ix, iquadpt)
     1                        * vlz ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * bfz ( ix, iquadpt) * ob4pi) * obx1

      s_quadpt ( ix, 5, iquadpt) = - ( tempac * vlx ( ix, iquadpt)
     1 - bfx ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = 0.0

      s_quadpt ( ix, 7, iquadpt) = 0.0

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Spherical Geometry

      IF ( igeom .EQ. 3) THEN

      DO iquadpt = NUM_ZONECEN_BASES + 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      x1 = nodal_points ( ix, 1, iquadpt)
      x2 = nodal_points ( ix, 2, iquadpt)
      x3 = nodal_points ( ix, 3, iquadpt)

      obx1 = 1.0 / x1

      cos_x2 = cos ( x2)
      sin_x2 = sin ( x2)
      cot_x2 = cos_x2 / sin_x2

      tempaa = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      tempab = ( vlx ( ix, iquadpt) * bfx ( ix, iquadpt)
     1         + vly ( ix, iquadpt) * bfy ( ix, iquadpt)
     1         + vlz ( ix, iquadpt) * bfz ( ix, iquadpt)) * ob4pi

      tempac = u_quadpt ( ix, 5, iquadpt) + prs ( ix, iquadpt) + tempaa

! ----------

      s_quadpt ( ix, 1, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 2.0 * vlx ( ix, iquadpt)
     1 + cot_x2 * vly ( ix, iquadpt) ) * obx1

      s_quadpt ( ix, 2, iquadpt) = rho ( ix, iquadpt)
     1 * ( - 2.0 * vlx ( ix, iquadpt)**2 + vly ( ix, iquadpt)**2
     1   + vlz ( ix, iquadpt)**2) * obx1
     1 - cot_x2
     1 * ( rho ( ix, iquadpt) * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1   - bfx ( ix, iquadpt) * bfy ( ix, iquadpt) * ob4pi) * obx1
     1 - ( - 2.0 * bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1   + bfz ( ix, iquadpt)**2) * ob4pi * obx1

      s_quadpt ( ix, 3, iquadpt) = - rho ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) * vly ( ix, iquadpt)
     1 + cot_x2
     1 * ( vly ( ix, iquadpt)**2 - vlz ( ix, iquadpt)**2)) * obx1
     1 + ( 3.0 * bfx ( ix, iquadpt) * bfy ( ix, iquadpt)
     1   + cot_x2
     1   * ( bfy ( ix, iquadpt)**2 - bfz ( ix, iquadpt)**2))
     1 * ob4pi * obx1

      s_quadpt ( ix, 4, iquadpt) =
     1 - rho ( ix, iquadpt) * vlz ( ix, iquadpt)
     1 * ( 3.0 * vlx ( ix, iquadpt) + 2.0 * cot_x2
     1   * vly ( ix, iquadpt) ) * obx1
     1 + bfz ( ix, iquadpt)
     1 * ( 3.0 * bfx ( ix, iquadpt) + 2.0 * cot_x2
     1   * bfy ( ix, iquadpt) ) * ob4pi * obx1

      s_quadpt ( ix, 5, iquadpt) =
     1 - 2.0 * ( tempac * vlx ( ix, iquadpt)
     1         - bfx ( ix, iquadpt) * tempab) * obx1
     1 - cot_x2 * ( tempac * vly ( ix, iquadpt)
     1 - bfy ( ix, iquadpt) * tempab) * obx1

      s_quadpt ( ix, 6, iquadpt) = - cot_x2
     1                           * g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 7, iquadpt) = g_quadpt ( ix, 6, iquadpt) * obx1

      s_quadpt ( ix, 8, iquadpt) = - f_quadpt ( ix, 8, iquadpt) * obx1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END IF

! ------------------------------

! Multiply "s_quadpt" by "dtmove".

      DO ix = oned_begin - oned_off, oned_end + oned_off
      s_quadpt ( ix, :, NUM_ZONECEN_BASES + 1: NUM_SPACETIME_BASES)=
     1 dtmove
     1*s_quadpt ( ix, :, NUM_ZONECEN_BASES + 1: NUM_SPACETIME_BASES)
      END DO

! ----------------------------------------------------------------------
#endif
! End of "BUILD_SOURCE_TERMS == 1"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 12) The t > 0 nodal values of "u_quadpt" have converged by this
! point. Use them to build the t > 0 modes of "u".

! ----------------------------------------------------------------------

#if ( NUMORDER == 2)
! cg, 2nd order : Build "u" modes at t > 0

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u ( ix, :, 05) = u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 01)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2".

! ----------------------------------------------------------------------

#if ( NUMORDER == 3)
! cg, 3rd order : Build "u" modes at t > 0

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u ( ix, :, 12) = 2.0 * ( u_quadpt ( ix, :, 15)
     1 - 2.0 * u_quadpt ( ix, :, 11) + u_quadpt ( ix, :, 01) )

      u ( ix, :, 11) = 2.0 * ( u_quadpt ( ix, :, 11)
     1 - u_quadpt ( ix, :, 01) - 0.25 * u ( ix, :, 12) )

      u ( ix, :, 13) = 4.0 * ( u_quadpt ( ix, :, 12)
     1 - u_quadpt ( ix, :, 11) - u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 01) )

      u ( ix, :, 14) = 4.0 * ( u_quadpt ( ix, :, 13)
     1 - u_quadpt ( ix, :, 11) - u_quadpt ( ix, :, 04)
     1 + u_quadpt ( ix, :, 01) )

      u ( ix, :, 15) = 4.0 * ( u_quadpt ( ix, :, 14)
     1 - u_quadpt ( ix, :, 11) - u_quadpt ( ix, :, 06)
     1 + u_quadpt ( ix, :, 01) )

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3".

! ----------------------------------------------------------------------

#if ( NUMORDER == 4)
! cg, 4th order : Build "u" modes at t > 0

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u ( ix, :, 21) = 0.5 * ( - 5.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 03) - u_quadpt ( ix, :, 04)
     1 - u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 06)
     1 - u_quadpt ( ix, :, 07) - u_quadpt ( ix, :, 08)
     1 + 12.0 * u_quadpt ( ix, :, 21) + u_quadpt ( ix, :, 22)
     1 + u_quadpt ( ix, :, 23) + u_quadpt ( ix, :, 24)
     1 + u_quadpt ( ix, :, 25) + u_quadpt ( ix, :, 26)
     1 + u_quadpt ( ix, :, 27) - 9.0 * u_quadpt ( ix, :, 31)
     1 + 2 * u_quadpt ( ix, :, 35) )

      u ( ix, :, 22) = 4.5 * ( 2.0 * u_quadpt ( ix, :, 01)
     1 - 5.0 * u_quadpt ( ix, :, 21) + 4.0 * u_quadpt ( ix, :, 31)
     1 - u_quadpt ( ix, :, 35) )

      u ( ix, :, 23) = - 4.5 * ( u_quadpt ( ix, :, 01)
     1 - 3.0 * u_quadpt ( ix, :, 21) + 3.0 * u_quadpt ( ix, :, 31)
     1 - u_quadpt ( ix, :, 35) )

      u ( ix, :, 24) = 3.0 * ( u_quadpt ( ix, :, 01)
     1 - 2.0 * u_quadpt ( ix, :, 03) + u_quadpt ( ix, :, 04)
     1 - 2.0 * u_quadpt ( ix, :, 21) + 3.0 * u_quadpt ( ix, :, 22)
     1 - u_quadpt ( ix, :, 23) + u_quadpt ( ix, :, 31)
     1 - u_quadpt ( ix, :, 32) )

      u ( ix, :, 25) = 3.0 * ( u_quadpt ( ix, :, 01)
     1 - 2.0 * u_quadpt ( ix, :, 05) + u_quadpt ( ix, :, 06)
     1 - 2.0 * u_quadpt ( ix, :, 21) + 3.0 * u_quadpt ( ix, :, 24)
     1 - u_quadpt ( ix, :, 25) + u_quadpt ( ix, :, 31)
     1 - u_quadpt ( ix, :, 33) )

      u ( ix, :, 26) = 3.0 * ( u_quadpt ( ix, :, 01)
     1 - 2.0 * u_quadpt ( ix, :, 07) + u_quadpt ( ix, :, 08)
     1 - 2.0 * u_quadpt ( ix, :, 21) + 3.0 * u_quadpt ( ix, :, 26)
     1 - u_quadpt ( ix, :, 27) + u_quadpt ( ix, :, 31)
     1 - u_quadpt ( ix, :, 34) )

      u ( ix, :, 27) = - 9.0 * ( u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 03) - 2.0 * u_quadpt ( ix, :, 21)
     1 + 2.0 * u_quadpt ( ix, :, 22) + u_quadpt ( ix, :, 31)
     1 - u_quadpt ( ix, :, 32) )

      u ( ix, :, 28) = - 9.0 * ( u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 05) - 2.0 * u_quadpt ( ix, :, 21)
     1 + 2.0 * u_quadpt ( ix, :, 24) + u_quadpt ( ix, :, 31)
     1 - u_quadpt ( ix, :, 33) )

      u ( ix, :, 29) = - 9.0 * ( u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 07) - 2.0 * u_quadpt ( ix, :, 21)
     1 + 2.0 * u_quadpt ( ix, :, 26) + u_quadpt ( ix, :, 31)
     1 - u_quadpt ( ix, :, 34) )

      u ( ix, :, 30) = 6.0 * ( 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 03) - u_quadpt ( ix, :, 04)
     1 - 2.0 * u_quadpt ( ix, :, 21) + u_quadpt ( ix, :, 22)
     1 + u_quadpt ( ix, :, 23) )

      u ( ix, :, 31) = 6.0 * ( 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 06)
     1 - 2.0 * u_quadpt ( ix, :, 21) + u_quadpt ( ix, :, 24)
     1 + u_quadpt ( ix, :, 25) )

      u ( ix, :, 32) = 6.0 * ( 2.0 * u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 07) - u_quadpt ( ix, :, 08)
     1 - 2.0 * u_quadpt ( ix, :, 21) + u_quadpt ( ix, :, 26)
     1 + u_quadpt ( ix, :, 27) )

      u ( ix, :, 33) = - 12.0 * ( u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 03) - u_quadpt ( ix, :, 05)
     1 + u_quadpt ( ix, :, 13) - u_quadpt ( ix, :, 21)
     1 + u_quadpt ( ix, :, 22) + u_quadpt ( ix, :, 24)
     1 - u_quadpt ( ix, :, 30) )

      u ( ix, :, 34) = - 12.0 * ( u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 07)
     1 + u_quadpt ( ix, :, 09) - u_quadpt ( ix, :, 21)
     1 + u_quadpt ( ix, :, 24) + u_quadpt ( ix, :, 26)
     1 - u_quadpt ( ix, :, 28) )

      u ( ix, :, 35) = - 12.0 * ( u_quadpt ( ix, :, 01)
     1 - u_quadpt ( ix, :, 03) - u_quadpt ( ix, :, 07)
     1 + u_quadpt ( ix, :, 11) - u_quadpt ( ix, :, 21)
     1 + u_quadpt ( ix, :, 22) + u_quadpt ( ix, :, 26)
     1 - u_quadpt ( ix, :, 29) )

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! STEP 13) Returning the space-time integral of source terms that are treated
! time-explicitly is easy and is done here.

! Because the "dudt" terms in the subroutine that calls this one already
! carry the source terms, we leave the "dtmove" multiplied into the
! "s_int" terms.

! IF "BUILD_SOURCE_TERMS == 0"then "s_quadpt" should be 0, so "s_int" will be 0.

! In principle, the factor of "dtmove", that was previously multiplied into
! the source terms, is divided out of "s_int" here. However,
! because the "dudt" terms in the subroutine that calls this one already
! carry the source terms, we leave the "dtmove" multiplied into the
! "s_int" terms.

! ------------------------------
#if ( NUMORDER == 2)
      DO ix = oned_begin - oned_off, oned_end + oned_off

      s_int ( ix, :) = 0.5 * ( s_quadpt ( ix, :, 01)
     1 + s_quadpt ( ix, :, 05) )

      END DO
#endif
! End of "NUMORDER == 2"
! ------------------------------
#if ( NUMORDER == 3)
      DO ix = oned_begin - oned_off, oned_end + oned_off

      s_int ( ix, :) = ( -5.0 * s_quadpt ( ix, :, 01)
     1 + 4.0 * s_quadpt ( ix, :, 11) + s_quadpt ( ix, :, 15)
     1 + s_quadpt ( ix, :, 02) + s_quadpt ( ix, :, 03)
     1 + s_quadpt ( ix, :, 04) + s_quadpt ( ix, :, 05)
     1 + s_quadpt ( ix, :, 06) + s_quadpt ( ix, :, 07) ) / 6.0

      END DO
#endif
! End of "NUMORDER == 3"
! ------------------------------
#if ( NUMORDER == 4)
      DO ix = oned_begin - oned_off, oned_end + oned_off

      s_int ( ix, :) = ( 15.0 * s_quadpt ( ix, :, 01)
     1 - 2.0 * s_quadpt ( ix, :, 03)
     1 - 2.0 * s_quadpt ( ix, :, 04)
     1 - 2.0 * s_quadpt ( ix, :, 05)
     1 - 2.0 * s_quadpt ( ix, :, 06)
     1 - 2.0 * s_quadpt ( ix, :, 07)
     1 - 2.0 * s_quadpt ( ix, :, 08)
     1 - 27.0 * s_quadpt ( ix, :, 21)
     1 + 6.0 * s_quadpt ( ix, :, 22)
     1 + 6.0 * s_quadpt ( ix, :, 23)
     1 + 6.0 * s_quadpt ( ix, :, 24)
     1 + 6.0 * s_quadpt ( ix, :, 25)
     1 + 6.0 * s_quadpt ( ix, :, 26)
     1 + 6.0 * s_quadpt ( ix, :, 27)
     1 + 9.0 * s_quadpt ( ix, :, 31)
     1 + 3.0 * s_quadpt ( ix, :, 35) ) / 24.0

      END DO
#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------

! STEP 14) Here we give ourselves the option of trapping a zone that has
! not converged. This may consist of zeroing out higher order moments, i.e.
! we fall back to first order here if we absolutely must.

! This is also where we can use Dumbser's trick of zeroing out the
! higher moments in zones that produce negative pressures or densities.

      corrupted_zone ( oned_begin - oned_off: oned_end + oned_off) = 0

      DO iquadpt = 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

        IF ( ( rho ( ix, iquadpt) .LT. rhofloor) .OR.
     1       ( prs ( ix, iquadpt) .LT. prsfloor) )
     1                                      corrupted_zone ( ix) = 1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO

! ------------------------------

#if ( N_EULER >= 1)

      ivar_lo = NFIELD + 1 + NFLUID

      DO i_euler = 1, N_EULER

      ivar = ivar_lo
     1     + ( i_euler - 1) * ( NFIELD_EULER + NFLUID_EULER)

      DO iquadpt = 1, NUM_SPACETIME_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

        IF ( ( rho_euler ( ix, iquadpt, i_euler) .LT. rhofloor) .OR.
     1       ( prs_euler ( ix, iquadpt, i_euler) .LT. prsfloor) )
     1                                      corrupted_zone ( ix) = 1

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO i_euler ="

#endif
! End of "N_EULER >= 1"

! ------------------------------

! For ADER_CG ( but not for ADER_DG) "f ( ix, :, 1)" is just the
! zone-centered flux at t=0. For ADER_DG, this will have to change.

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( corrupted_zone ( ix) .EQ. 1) THEN

! -----

      u ( ix, :, 1) = w ( ix, :, 1)

      s_int ( ix, :) = 0.0

! -----

      w ( ix, :, 2: NUM_ZONECEN_BASES) = 0.0

      u ( ix, :, 2: NUM_SPACETIME_BASES) = 0.0

! -----

      END IF

! ----------

      END DO
! End of "DO ix =" loop.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ADER_CG_NODAL_3D_TIME_EVOL_MHD

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

