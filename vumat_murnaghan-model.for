     ! Created by: Sunia Tanweer, National University of Sciences
     ! and Technology (NUST), Islamabad, Pakistan
     ! Supervised by: Dr. Xiang Zhang, University of Wyoming,
     ! Laramie WY, USA
     ! Computations for Advanced Materials and Manufacturing Lab, 
     ! University of Wyoming, USA
     ! VUMAT to be used with ABAQUS/Explicit 

     ! Five inputs required, refer to: 
     ! Numerical implementation of the Murnaghan material model in ABAQUS/Standard
     ! Stanisław  Jemioło, Aleksander  Franus
     ! MATEC Web Conf. 196 01042 (2018)
     ! DOI: 10.1051/matecconf/201819601042


      subroutine vumat(                                                         &
! Read only (unmodifiable)variables -
        nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,                  &
        stepTime, totalTime, dt, cmname, coordMp, charLength,                   &
        props, density, strainInc, relSpinInc,                                  &
        tempOld, stretchOld, defgradOld, fieldOld,                              &
        stressOld, stateOld, enerInternOld, enerInelasOld,                      &
        tempNew, stretchNew, defgradNew, fieldNew,                              &
! Write only (modifiable) variables -
        stressNew, stateNew, enerInternNew, enerInelasNew )
!
        implicit none
        integer j_sys_Dimension, n_vec_Length, maxblk
        !---additioinal parameters from vaba_param.inc
        parameter (j_sys_Dimension = 2)
        parameter( n_vec_Length = 136  )
        parameter( maxblk = n_vec_Length  )

!------Start my explicit variable type
        character*80 cmname
        integer nprops, nblock, nshr, nstatev, nfieldv, ndir,lanneal
        real*8 stepTime, totalTime, dt
!
        real*8 props(nprops), density(nblock), coordMp(nblock,*),              &
        strainInc(nblock,ndir+nshr),                        &
        relSpinInc(nblock,nshr), tempOld(nblock),                               &
        stretchOld(nblock,ndir+nshr),                                           &
        defgradOld(nblock,ndir+nshr+nshr),                                      &
        fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),                  &
        stateOld(nblock,nstatev), enerInternOld(nblock),                        &
        enerInelasOld(nblock), tempNew(nblock),                                 &
        stretchNew(nblock,ndir+nshr),                                           &
        defgradNew(nblock,ndir+nshr+nshr),                                      &
        fieldNew(nblock,nfieldv),                                               &
        stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),                  &
        enerInternNew(nblock), enerInelasNew(nblock), charLength(nblock) 

!

!------------Define parameters---------------------------------------------------
       integer i,j,k
      real*8 a1,a2,a3,a4,a5,a6,a7,lambda,mu,v1,v2,v3
      real*8 l,m,n,cofu2,cofu4,cofidem,lam
      
      real*8 evol,bxx,byy,bzz,bxy,byz,bxz,detu, detF, bi1,bi2,trace_bsq,xpow

      real*8 u(3,3), u2(3,3), u4(3,3), corot(3,3), F(3,3)

      real*8 zero,half,two,four,fourthird,three,twothird,nine, &
                one,onethird,six,twoninth,fourninth,eightthird, &
                twelve,eightninth,seventhird,fivethird, &
                eighttwentyseventh,eighth,fourth,twentyfourth, Ident2nd(3,3)

      parameter(zero=0.d0,half=1.d0/2.d0,two=2.d0,four=4.d0,&
                        twothird=2.d0/3.d0,nine=9.d0,fourthird=4.d0/3.d0,&
                        three=3.d0,one=1.d0,onethird=1.d0/3.d0,six=6.d0,&
                        twoninth=2.d0/9.d0,fourninth=4.d0/9.d0,&
                        eightthird=8.d0/3.d0,twelve=12.d0,&
                        eightninth=8.d0/9.d0,seventhird=7.d0/3.d0,&
                        fivethird=5.d0/3.d0,eighttwentyseventh=8.d0/27.d0,&
                  eighth=1.d0/8.d0,fourth=1.d0/4.d0,twentyfourth=1.d0/24.d0)

      parameter (Ident2nd =reshape([1.d0, 0.d0, 0.d0,           & 
                          0.d0, 1.d0, 0.d0,         &
                          0.d0, 0.d0, 1.d0],[3,3]))

!------------Input material properties---------------------------------------------------

      lambda = props(1)
      mu = props(2)
      v1 = props(3)
      v2 = props(4)
      v3 = props(5)
      k = 1
      i = 1
      j = 1

!------------Calculate l, m and n---------------------------------------------------

      l = (v1*half) + v2
      m = v2 + two*v3 
      n = four*v3

!------------Calculate a_i for the stored energy potential function--------------

      a1 = (- six*lambda - four*mu + n + nine*l)*eighth
      a2 = (lambda + two*mu - two*m - three*l)*eighth
      a3 = (- four*mu + six*m - n)*eighth
      a4 = - m*fourth
      a5 = (l + two*m)*twentyfourth
      a6 = n*eighth
      a7 = (nine*lambda + six*mu - nine*l)*eighth

!    calculate b matrix, stresses, etc. in block form
!
     do k = 1, nblock

!---Obtain U tensor (F=RU)

    u(1,1) = stretchNew(k,1)
    u(2,2) = stretchNew(k,2)
    u(3,3) = stretchNew(k,3)
    u(1,2) = stretchNew(k,4)

!--Obtain deformation gradient tensor F, see http://ivt-abaqusdoc.ivt.ntnu.no:2080/v6.14/books/sub/default.htm for the ordering of F from 3x3 to 9x1

    F(1,1) = defgradNew(k,1)
    F(2,2) = defgradNew(k,2)
    F(3,3) = defgradNew(k,3)
    F(1,2) = defgradNew(k,4)
    F(2,1) = defgradNew(k,5)

!--Calculate B tensor (B should be U^2)

     bxx = u(1,1)*u(1,1) + u(1,2)*u(1,2)
     byy = u(1,2)*u(1,2) + u(2,2)*u(2,2)
     bzz = u(3,3)*u(3,3)
     bxy = u(1,1)*u(1,2) + u(1,2)*u(2,2)
     bxz = zero
     byz = zero

!-- determint of U
    detu = u(3,3)*(u(1,1)*u(2,2) - u(1,2)*u(1,2))

!-- determint of F
   detF= F(3,3)*(F(1,1)*F(2,2) - F(1,2)*F(2,1))

!----for 3D case----------------
    if (nshr > one ) then

       u(1,1) = stretchNew(k,1)
       u(2,2) = stretchNew(k,2)
       u(3,3) = stretchNew(k,3)
       u(1,2) = stretchNew(k,4)
       u(2,3) = stretchNew(k,5)
       u(3,1) = stretchNew(k,6)
       u(1,3) = stretchNew(k,6)
       u(3,2) = stretchNew(k,5)
       u(2,1) = stretchNew(k,4)

!--Obtain deformation gradient tensor F, see http://ivt-abaqusdoc.ivt.ntnu.no:2080/v6.14/books/sub/default.htm for the ordering of F from 3x3 to 9x1
       F(2,3) = defgradNew(k,5)
       F(3,1) = defgradNew(k,6)
       F(2,1) = defgradNew(k,7)
       F(3,2) = defgradNew(k,8)
       F(1,3) = defgradNew(k,9)

! ------ find u^2 by matmul ---------------
       u2=matmul(u,u)

! -------- equate cauchy green tensor with u^2 -------------
       bxx = u2(1,1)
       byy = u2(2,2)
       bzz = u2(3,3)
       bxy = u2(1,2)
       bxz = u2(1,3)
       byz = u2(2,3)

!-- determint of U
       detu = detu &
            + u(3,1)*(u(1,2)*u(2,3) - u(2,2)*u(1,3)) &  
            - u(3,2)*(u(2,3)*u(1,1) - u(2,1)*u(1,3))

!-- determint of F
       detF = detF &
            + F(3,1)*(F(1,2)*F(2,3) - F(2,2)*F(1,3)) &  
            - F(3,2)*(F(2,3)*F(1,1) - F(2,1)*F(1,3))
    end if

    xpow = (one/(detF**two)**onethird)
    bxx = bxx * xpow 
    byy = byy * xpow 
    bzz = bzz * xpow
    bxy = bxy * xpow
    bxz = bxz * xpow
    byz = byz * xpow

    u(1,1) = u(1,1)*(xpow**half)
    u(1,2) = u(1,2)*(xpow**half)
    u(1,3) = u(1,3)*(xpow**half)
    u(2,1) = u(2,1)*(xpow**half)
    u(2,2) = u(2,2)*(xpow**half)
    u(2,3) = u(2,3)*(xpow**half)
    u(3,1) = u(3,1)*(xpow**half)
    u(3,2) = u(3,2)*(xpow**half)
    u(3,3) = u(3,3)*(xpow**half)

    u2=matmul(u,u)
    u4=matmul(u2, u2)

!   BI1 ( first invariant of B )
!
    bi1 = bxx + byy + bzz

    trace_bsq = u4(1,1) + u4(2,2) + u4(3,3)

!   BI2 ( second invariant of B )

    bi2 = half*((bi1**two) - trace_bsq)

!--- coefficient of the B (or U^2) term
    cofu2 = (two*a1/(detF**onethird)) + four*a2*bi1*(detF**onethird) &
                 + two*a4*detF*bi2 + six*a5*detF*(bi1**two) + two*a3*bi1*(detF**onethird) &
                 + two*a4*detF*(bi1**2)

!--- coefficient of the B^2 (or U^4) term
    cofu4 = two*a3*(detF**onethird) + two*a4*detF*bi1

!--- coefficient of the identity tensor term
    cofidem = two*a6*detF

!---- corotational stress for murnaghan model ---------- 
    corot=u2*cofu2- u4*cofu4 + cofidem*Ident2nd

! update stress for the element ---------------

    stressNew(k,1) = corot(1,1)
    stressNew(k,2) = corot(2,2) 
    stressNew(k,3) = corot(3,3)
    stressNew(k,4) = corot(1,2) 
    if ( nshr > 1 ) then
        stressNew(k,5) =corot(2,3)
        stressNew(k,6) =corot(3,1)
    end if

! internal energy updated is equal to stored energy potential function

    enerInternNew(k) = a1*bi1*((detu**two)**onethird) &
            + a2*(bi1**two)*(detu**fourthird) + &
            a3*bi2*(detu**fourthird) + a4*bi1*bi2*(detu**two) + &
            a5*(detu**two)*(bi1**three) - a6 + a6*(detu**two) + a7

    enerInternNew(k) = enerInternNew(k)/density(k) 

     end do

     return
     end

