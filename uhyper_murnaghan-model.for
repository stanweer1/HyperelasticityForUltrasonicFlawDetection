     ! Created by: Sunia Tanweer, National University of Sciences
     ! and Technology (NUST), Islamabad, Pakistan
     ! Supervised by: Dr. Xiang Zhang, University of Wyoming,
     ! Laramie WY, USA
     ! Computations for Advanced Materials and Manufacturing Lab, 
     ! University of Wyoming, USA
     ! UHYPER to be used with ABAQUS/Standard 

     ! Five inputs required, refer to: 
     ! Numerical implementation of the Murnaghan material model in ABAQUS/Standard
     ! Stanisław  Jemioło, Aleksander  Franus
     ! MATEC Web Conf. 196 01042 (2018)
     ! DOI: 10.1051/matecconf/201819601042

      subroutine uhyper(BI1,BI2,AJ,U,UI1,UI2,UI3,TEMP, &
      NOEL,CMNAME,INCMPFLAG,NUMSTATEV, &
      STATEV,NUMFIELDV,FIELDV,FIELDVINC, &
      NUMPROPS,PROPS)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION UI1(3),UI2(6),UI3(6),STATEV(*), &
      FIELDV(*),FIELDVINC(*),PROPS(*)

      ! declare variables
      !
      real*8 a1,a2,a3,a4,a5,a6,a7,lambda,mu,v1,v2,v3,l,m,n 

     ! declare variables to initialize as constants    
     !
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

      ! initalize
      !
      U = zero
      UI1 = zero
      UI2 = zero
      UI3 = zero

     ! Stored energy potential function for Murnaghan model of hyperelasticity
     !
     U = a1*BI1*((AJ**two)**onethird) + a2*(BI1**two)*(AJ**fourthird) + &
            a3*BI2*(AJ**fourthird) + a4*BI1*BI2*(AJ**two) + &
            a5*(AJ**two)*(BI1**three) + a6*(-one + (AJ**two)) + a7

      ! first derivatives
      !
      UI1(1) =  a1*((AJ**two)**onethird) + two*a2*BI1*(AJ**fourthird) + &
                    a4*BI2*(AJ**two) + three*a5*(BI1**two)*(AJ**two)

      UI1(2) = a3*(AJ**fourthird) + a4*BI1*(AJ**two)

      UI1(3) = (twothird)*a1*BI1*AJ**(-onethird) + &
                    (fourthird)*a2*(BI1**two)*(AJ**onethird) + &
                    (fourthird)*a3*BI2*(AJ**onethird) + two*a4*BI1*BI2*AJ &
                    + two*a5*AJ*(BI1**three) + two*a6*AJ
      !
      ! second derivatives
      !
      UI2(1) = two*a2*(AJ**fourthird) + six*a5*BI1*(AJ**two)

      UI2(2) =  zero

      UI2(3) =  (-twoninth)*a1*BI1*(AJ**(-fourthird)) + &
                (fourninth)*a2*(BI1**two)*((AJ**(-two))**onethird) + &
                (fourninth)*a3*BI2*((AJ**(-two))**onethird) + two*a4*BI1*BI2 &
                + two*a5*(BI1**three) + two*a6

      UI2(4) =  a4*(AJ**two)

      UI2(5) =  (twothird)*a1*(AJ**(-onethird)) + &
                (eightthird)*a2*BI1*(AJ**onethird) &
                + two*a4*BI2*AJ + six*a5*(BI1**two)*AJ

      UI2(6) = (fourthird)*a3*(AJ**onethird) + two*a4*BI1*AJ

      !
      ! third derivatives
      !
 
      UI3(1) = (eightthird)*a2*(AJ**onethird) + twelve*a5*BI1*AJ

      UI3(2) = zero

      UI3(3) = two*a4*AJ

      UI3(4) =  (-twoninth)*a1*(AJ**(-fourthird)) + &
                (eightninth)*a2*BI1*((AJ**(-two))**onethird) &
                + two*a4*BI2 + six*a5*(BI1**two)

      UI3(5) =  (fourninth)*a3*((AJ**(-two))**onethird) + two*a4*BI1

      UI3(6) =  (eighttwentyseventh)*a1*BI1*(AJ**(-seventhird)) &
                + (-eighttwentyseventh)*a2*(BI1**two)*(AJ**(-fivethird)) &
                + (-eighttwentyseventh)*a3*BI2*(AJ**(-fivethird))

     ! Compute the effective stretch (not necessary)
     !
      effStr = dsqrt(dabs(BI1*onethird))


      ! Update the state variables
      !
      statev(1) = BI1
      statev(2) = BI2
      statev(3) = AJ
      statev(4) = effStr

      return
      end subroutine uhyper
