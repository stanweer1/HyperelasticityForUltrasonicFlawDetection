     ! Created by: Sunia Tanweer, National University of Sciences
     ! and Technology (NUST), Islamabad, Pakistan
     ! Supervised by: Dr. Xiang Zhang, University of Wyoming,
     ! Laramie WY, USA
     ! Computations for Advanced Materials and Manufacturing Lab, 
     ! University of Wyoming, USA
     ! VUANISOHYPER to be used with ABAQUS/EXPLICIT 

     ! Five inputs required, refer to: 
     ! Numerical implementation of the Murnaghan material model in ABAQUS/Standard
     ! Stanisław  Jemioło, Aleksander  Franus
     ! MATEC Web Conf. 196 01042 (2018)
     ! DOI: 10.1051/matecconf/201819601042

      subroutine vuanisohyper_inv (nblock, nFiber, nInv, &
        jElem, kIntPt, kLayer, kSecPt, cmname, nstatev, &
        nfieldv, nprops, props, tempOld, tempNew, &
        fieldOld, fieldNew, stateOld, sInvariant, zeta, &  
        uDev, duDi, d2uDiDi, stateNew)
!
      include 'vaba_param.inc'
!
      dimension props(nprops),  tempOld(nblock), &
       fieldOld(nblock,nfieldv), stateOld(nblock,nstatev),& 
       tempNew(nblock), fieldNew(nblock,nfieldv),&
       sInvariant(nblock,nInv), zeta(nblock,nFiber*(nFiber-1)/2),&
       uDev(nblock), duDi(nblock,nInv), &
       d2uDiDi(nblock,nInv*(nInv+1)/2), stateNew(nblock,nstatev)
!
      character*80 cmname
!

!    Declare variables
!
      integer i,j,k
      real*8 a1,a2,a3,a4,a5,a6,a7,lambda,mu,v1,v2,v3
      real*8 AJ,BI1,BI2,l,m,n,U

      real*8 zero,half,two,four,fourthird,three,twothird,nine, &
                one,onethird,six,twoninth,fourninth,eightthird, &
                twelve,eightninth,seventhird,fivethird, &
                eighttwentyseventh,eighth,fourth,twentyfourth

      parameter(zero=0.d0,half=1.d0/2.d0,two=2.d0,four=4.d0,&
                        twothird=2.d0/3.d0,nine=9.d0,fourthird=4.d0/3.d0,&
                        three=3.d0,one=1.d0,onethird=1.d0/3.d0,six=6.d0,&
                        twoninth=2.d0/9.d0,fourninth=4.d0/9.d0,&
                        eightthird=8.d0/3.d0,twelve=12.d0,&
                        eightninth=8.d0/9.d0,seventhird=7.d0/3.d0,&
                        fivethird=5.d0/3.d0,eighttwentyseventh=8.d0/27.d0,&
                  eighth=1.d0/8.d0,fourth=1.d0/4.d0,twentyfourth=1.d0/24.d0)

      ! Obtain material properties
      !

      lambda = props(1)
      mu = props(2)
      v1 = props(3)
      v2 = props(4)
      v3 = props(5)
      k = 1
      i = 1
      j = 1

     ! Calculate l, m and n 
     !
      l = (v1*half) + v2
      m = v2 + two*v3 
      n = four*v3

      ! Calculate a_i for the stored energy potential function
      !
      a1 = (- six*lambda - four*mu + n + nine*l)*eighth
      a2 = (lambda + two*mu - two*m - three*l)*eighth
      a3 = (- four*mu + six*m - n)*eighth
      a4 = - m*fourth
      a5 = (l + two*m)*twentyfourth
      a6 = n*eighth
      a7 = (nine*lambda + six*mu - nine*l)*eighth

     ! Enter do loop -- vectorized code
     ! 
     do k = 1, nblock

     ! Equate invariants with given variables in vuanisohyper
     !
     BI1 = sInvariant(k,1)
     BI2 = sInvariant(k,2)
     AJ = sInvariant(k,3)

     ! Stored energy potential function for Murnaghan model of hyperelasticity
     !
     U = a1*BI1*((AJ**two)**onethird) + a2*(BI1**two)*(AJ**fourthird) + &
            a3*BI2*(AJ**fourthird) + a4*BI1*BI2*(AJ**two) + &
            a5*(AJ**two)*(BI1**three) + a6*(-one + (AJ**two)) + a7

      ! first derivatives
      !
      duDi(k,1) =  a1*((AJ**two)**onethird) + two*a2*BI1*(AJ**fourthird) + &
                    a4*BI2*(AJ**two) + three*a5*(BI1**two)*(AJ**two)

      duDi(k,2) = a3*(AJ**fourthird) + a4*BI1*(AJ**two)

      duDi(k,3) = twothird*a1*BI1/(AJ**onethird) + &
                    (fourthird)*a2*(BI1**two)*(AJ**onethird) + &
                    (fourthird)*a3*BI2*(AJ**onethird) + two*a4*BI1*BI2*AJ &
                    + two*a5*AJ*(BI1**three) + two*a6*AJ
      !
      ! second derivatives
      !
      d2uDiDi(k,1) = two*a2*(AJ**fourthird) + six*a5*BI1*(AJ**two)

      d2uDiDi(k,3) =  zero

      d2uDiDi(k,6) =  (-twoninth)*a1*BI1*(AJ**(-fourthird)) + &
                (fourninth)*a2*(BI1**two)*((AJ**(-two))**onethird) + &
                (fourninth)*a3*BI2*((AJ**(-two))**onethird) + two*a4*BI1*BI2 &
                + two*a5*(BI1**three) + two*a6

      d2uDiDi(k,2) =  a4*(AJ**two)

      d2uDiDi(k,4) =  (twothird)*a1*(AJ**(-onethird)) + &
                (eightthird)*a2*BI1*(AJ**onethird) &
                + two*a4*BI2*AJ + six*a5*(BI1**two)*AJ

      d2uDiDi(k,5) = (fourthird)*a3*(AJ**onethird) + two*a4*BI1*AJ

     ! Update state 
     !
      stateNew(k, 1) = BI1
      stateNew(k, 2) = BI2
      stateNew(k, 3) = AJ

     end do

      return
      end