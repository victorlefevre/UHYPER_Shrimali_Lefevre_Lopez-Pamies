!**********************************************************************
! Legal notice: UHYPER_Shrimali_Lefevre_Lopez-Pamies.for (Windows) 
!
! Copyright (C) 2018 Bhavesh Shrimali (bshrima2@illinois.edu)
!                    Victor Lefèvre (victor.lefevre@northwestern.edu)
!                    Oscar Lopez-Pamies (pamies@illinois.edu)
!
! This ABAQUS UHYPER subroutine implements the hyperelastic energy  
! density derived in [1] for the macroscopic elastic response of
! non-Gaussian elastomers weakened by an isotropic and non-percolative 
! distribution of equiaxed pores. This result is valid for any choice 
! of I1-based incompressible energy density characterizing 
! the non-Gaussian isotropic elastic response of the underlying   
! elastomer. The present subroutine is implemented for the  
! choice of strain energy density proposed in [2].
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see https://www.gnu.org/licenses/
!
!**********************************************************************
! Usage:
!
! The subroutine is to be used as an compressible USER hyperelastic 
! model with 5 material properties, e.g.,
! *HYPERELASTIC, USER, TYPE=COMPRESSIBLE, PROPERTIES=5
! in the input (.inp) file.
!
! The 5 materials properties for the model to be provided as input to
! the subroutine via the PROPS array are listed in the table below:

!  AMU1    = PROPS(1)  ! PARAMETER #1 OF THE ELASTOMER
!  ALPHA1  = PROPS(2)  ! EXPONENT #1 OF THE ELASTOMER
!  AMU2    = PROPS(3)  ! PARAMETER #2 OF THE ELASTOMER
!  ALPHA2  = PROPS(4)  ! EXPONENT #2 OF THE ELASTOMER 
!  AF0      = PROPS(5)  ! INITIAL POROSITY 
!
! The two material parameters AMU1, AMU2 characterizing the elastic    
! behavior of the underlying elastomer are non-negative real numbers 
! (AMU1 >= 0, AMU2 >= 0) with strictly positive sum (AMU1 + AMU2 > 0). 
! The two exponents ALPHA1, ALPHA2 are non-zero real numbers 
! (ALPHA1 ≠ 0, ALPHA2 ≠ 0) leading to a strongly elliptic strain 
! energy (see eq. (22) in [2]). This is left to the user to check.
!
! The initial porosity (AF0) must satisfy 0 <= AF0 <= 1. 
!
! As expected from physical considerations, this macroscopic energy
! remains finite so long the determinant of the deformation 
! gradient (AJ) satisfies the condition AJ - 1 + AF0 > 0. This
! inequality constraint is enforced through a MOREAU-YOSIDA
! regularization (see, e.g. [3]). The underlying weight (ANU) is set
! here by default to ANU = 1.0e15. The subroutine issues two kinds
! of messages regarding the constraint:
!  -- a WARNING message when AJ - 1 + AF0 < 1e-9 which allows the job
!      to carry on
!  -- an ERROR message when AJ - 1 + AF0 < -0.01 and TERMINATES the job
! In both cases, please treat the results with caution and check that 
! the current local porosity given by (AJ - 1 + AF0) / AJ 
! remains positive; see below on how to request it.
!
! The porosity in the deformed configuration (AJ - 1 + AF0) / AJ
! is required to be output (to check the results for instance, 
! see above) as a solution-dependent state variable (SDV), e.g.,
! using the following lines in the input (.inp) file:
! *DEPVAR
! 1
! 1, Porosity, Current local porosity
! The solution-dependent state variable may be initialized using the 
! following lines in the input (.inp) file:
! *INITIAL CONDITIONS, TYPE=SOLUTION
! <some element set>, AF0
!
!**********************************************************************
! Additional information:
!
! This subroutine creates a solution-dependent state variable for
! for the current porosity (see above) but does not create predefined 
! field variables. 
!
! Examples can be found in the article posted in the SIMULIA Learning 
! Community: https://r1132100503382-eu1-3dswym.3dexperience.3ds.com/
! #community:39/post:mRdxC3xkRzajJ6LVk0SgwA
!
! Please consult the ABAQUS Documentation for additional references 
! regarding the use of compressible USER hyperelastic models with
! the UHYPER subroutine and the use of solution-dependent 
! state variables.
!
! Due the near-incompressible nature of this model at low porosities,
! the use of hybrid elements is recommended.
!
!**********************************************************************
! References:
!
! [1] Shrimali, B., Lefèvre, V., Lopez-Pamies, O. 2019. A simple
!     explicit homogenization solution for the macroscopic elastic 
!     response of isotropic porous elastomers. J. Mech. Phys. Solids 
!     122, 364--380.
! [2] Lopez-Pamies, O., 2010. A new I1-based hyperelastic model for 
!     rubber elastic materials. C. R. Mec. 338, 3--11.
! [3] Parikh, N., Boyd, S., 2013. Proximal algorithms. Found. Trends
!     Optim. 1, 123--231.
!
!**********************************************************************
!      
      SUBROUTINE UHYPER(BI1,BI2,AJ,U,UI1,UI2,UI3,TEMP,NOEL,
     1 CMNAME,INCMPFLAG,NUMSTATEV,STATEV,NUMFIELDV,FIELDV,
     2 FIELDVINC,NUMPROPS,PROPS)
!
      INCLUDE 'ABA_PARAM.INC'
#INCLUDE <SMAASPUSERSUBROUTINES.HDR>
!
      CHARACTER*80 CMNAME
      DIMENSION U(2),UI1(3),UI2(6),UI3(6),STATEV(*),FIELDV(*),
     1 FIELDVINC(*),PROPS(*)     
!
!     STDB_ABQERR AND GET_THREAD_ID INITIALIZATION
!
      DIMENSION INTV(1),REALV(3)
      CHARACTER*8 CHARV(1)
      CHARACTER*100 STRING1, STRING2, STRING3
      CHARACTER*300 STRING
!
      DATA LWRITE /1/
!      
      INTEGER MYTHREADID  
!
      INTV(1)=0
      REALV(1)=0.
      REALV(2)=0.
      REALV(3)=0.
      CHARV(1)=''
!
      MYTHREADID = GET_THREAD_ID()
!      
!     INPUT CHECKS
!
      IF (MYTHREADID.EQ.0) THEN
        IF (INCMPFLAG.EQ.1) THEN
          STRING1='INCOMPRESSIBILITY FLAG IS 1. THE MODEL IS COMPRES'
          STRING2='SIBLE. SET USER TYPE=COMPRESSIBLE.'
          STRING = TRIM(STRING1) // TRIM(STRING2)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE IF (NUMSTATEV.NE.1) THEN  
          INTV(1)=NUMSTATEV
          STRING1='RECEIVED REQUEST FOR %I SOLUTION-DEPENDENT STATE'
          STRING2=' VARIABLES. THE SUBROUTINE CREATES 1 SOLUTION'
          STRING3='-DEPENDENT STATE VARIABLE.'
          STRING = TRIM(STRING1) // TRIM(STRING2) // TRIM(STRING3)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE IF (NUMFIELDV.NE.0) THEN 
          INTV(1)=NUMFIELDV 
          STRING1='RECEIVED REQUEST FOR %I PREDEFINED FIELD  VARI'
          STRING2='ABLES. THE SUBROUTINE DOES NOT CREATE PREDEFINED'
          STRING3=' FIELD VARIABLES.'
          STRING = TRIM(STRING1) // TRIM(STRING2) // TRIM(STRING3)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE IF (NUMPROPS.NE.5) THEN    
          INTV(1)=NUMPROPS  
          STRING1='RECEIVED %I MATERIAL PROPERTIES. THE SUBROUTINE'
          STRING2=' REQUIRES 5 MATERIAL PROPERTIES.'
          STRING = TRIM(STRING1) // TRIM(STRING2)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF
      END IF     
!
!     MATERIAL PARAMETERS
!  
      AMU1    = PROPS(1)  ! PARAMETER #1 OF THE ELASTOMER
      ALPHA1  = PROPS(2)  ! I1 EXPONENT #1 OF THE ELASTOMER
      AMU2    = PROPS(3)  ! PARAMETER #2 OF THE ELASTOMER
      ALPHA2  = PROPS(4)  ! I1 EXPONENT #2 OF THE ELASTOMER 
      AF0      = PROPS(5)  ! INITIAL POROSITY
!
!     PARTIAL MATERIAL PARAMETERS CHECKS
!     
      IF (((AMU1.LT.0.).OR.(AMU2.LT.0.).OR.(AMU1+AMU2.LE.0.))
     1 .AND.(MYTHREADID.EQ.0)) THEN
        REALV(1)=AMU1      
        REALV(2)=AMU2      
        REALV(3)=AMU1+AMU2      
        STRING1='RECEIVED AMU1 = %R AND AMU2 = %R, AMU1 + AMU2 = %R.'
        STRING2=' THE PARAMETERS AMU1 AND AMU2 MUST BE NON-NEGATIVE'
        STRING3=' AND AMU1 + AMU2, MUST BE GREATER THAT ZERO.'
        STRING = TRIM(STRING1) // TRIM(STRING2) // TRIM(STRING3)
        CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
      END IF            
!      
      IF (((ALPHA1.EQ.0.).OR.(ALPHA2.EQ.0.)).AND.
     1 (MYTHREADID.EQ.0)) THEN
        REALV(1)=ALPHA1      
        REALV(2)=ALPHA2      
        STRING1='RECEIVED ALPHA1 = %R AND ALPHA2 = %R.'
        STRING2=' THE EXPONENTS ALPHA1 AND ALPHA2 MUST BE NON-ZERO.'
        STRING = TRIM(STRING1) // TRIM(STRING2)
        CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
      END IF             
!      
      IF (AF0.LT.0.) THEN 
        IF (MYTHREADID.EQ.0) THEN
          REALV(1)=AF0  
          STRING1='RECEIVED AF0 = %R. THE INITIAL POROSITY IS NEGATIVE.'
          STRING = TRIM(STRING1)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE
          CALL XIT
        END IF
      ELSE IF (AF0.GT.1.) THEN 
        IF (MYTHREADID.EQ.0) THEN
          REALV(1)=AF0  
          STRING1='RECEIVED AF0 = %R. THE INITIAL POROSITY' 
          STRING2=' IS GREATER THAN 1.'
          STRING = TRIM(STRING1) // TRIM(STRING2)   
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE
          CALL XIT
        END IF
      END IF  
!
!     CONSTRAINT AJ + AF0 - 1 > 0 CHECK: WARNING MESSAGE
!
      IF (AJ+AF0.LT.1.0 .AND. LWRITE.EQ.1 ) THEN 
        LWRITE = 0
        INTV(1)=NOEL  
        STRING1='THE POROSITY IS NEGATIVE IN ELEMENT #%I.'    
        STRING2=' TREAT THE RESULTS WITH CAUTION.'    
        STRING = TRIM(STRING1) // TRIM(STRING2)   
        CALL STDB_ABQERR(-1,STRING,INTV,REALV,CHARV)
      END IF 
!
!     CONSTRAINT AJ + AF0 - 1 > 0 CHECK: JOB TERMINATION
! 
      IF (AJ+AF0.LE.0.99) THEN  
        INTV(1)=NOEL    
        STRING1='THE POROSITY IS NEGATIVE IN ELEMENT #%I.'    
        STRING2=' JOB TERMINATED.'    
        STRING = TRIM(STRING1) // TRIM(STRING2) 
        CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
      END IF
!
!     RECURRING RATIOS AND FACTORS
!      
      AOT = 1./3.
      ATT = 2./3.
      AFT = 4./3.
      AST = 7./3.
      ATET = 10./3. 
      AFN = 4./9.
      AOMAC = 1.-AF0
      ACT = AJ-AOMAC ! AJ + AF0 - 1
!
!     ENFORCE LOWER LIMIT TO ACT FOR THE EVALUATION OF ASI1
!
      ACT_MIN = 1.0e-9
      ACT0 = MAX(ACT, ACT_MIN)
      DACT0DJ = 1.0
      IF (ACT.LT.ACT_MIN) DACT0DJ = 0.0
!
!     FIRST INVARIANT AI1 = F.F AND PARTIAL DERIVATIVES
!      
      AI1=BI1*AJ**ATT
      DAI1DBI1=AJ**ATT
      DAI1DAJ=ATT*BI1*AJ**(-AOT)
      DAI1DBI1AJ=ATT*AJ**(-AOT)
      DAI1DAJAJ=-AOT*ATT*BI1*AJ**(-AFT)
      DAI1DBI1AJAJ=-AOT*ATT*AJ**(-AFT)
      DAI1DAJAJAJ=AFT*AOT*ATT*BI1*AJ**(-AST)
!
!     EQ (19) IN [1] AND PARTIAL DERIVATIVES
!      
      ASI1=3.*AOMAC*(AI1-3.)/(3.+2.*AF0) + (3.*(2.*AJ-1. -
     1 AOMAC*(2.*AF0+3.*AJ**ATT)*AJ**AOT/(3.+2.*AF0) - 
     2 AF0**AOT*AJ**AOT*(2.*ACT0-AF0)/ACT0**AOT))/AJ**AOT
!      
      DASI1DAI1=3.*AOMAC/(3.+2.*AF0)
!      
      DASI1DAJ=AJ**(-AFT)+(2.*(3.+7.*AF0))/((3.+2.*AF0)*AJ**AOT) - 
     1 AF0**AOT*(4.*ACT0+AF0)*DACT0DJ/ACT0**AFT  
!      
      DASI1DAJAJ=(2.*AF0**AOT*(AF0+ACT0)*DACT0DJ**2.0/ACT0**AST 
     1 -2./AJ**AST - (3.+7.*AF0)/((3.+2.*AF0)*AJ**AFT))*ATT
!      
      DASI1DAJAJAJ=(7./AJ**ATET + 2.*(3.+7.*AF0)/((3.+2.*AF0)*AJ**AST) - 
     1 AF0**AOT*(4.*ACT0+7.*AF0)*DACT0DJ**3.0/ACT0**ATET)*AFN
!
!     ARGUMENT IN EQ (18) IN [1] AND PARTIAL DERIVATIVES
!      
      AII1=ABS(ASI1)/AOMAC+3.
!      
      DAII1DASI1=1./AOMAC      
!      
      DAII1DBI1=DAII1DASI1*DASI1DAI1*DAI1DBI1
!      
      DAII1DAJ=DAII1DASI1*(DASI1DAI1*DAI1DAJ+DASI1DAJ)
!      
      DAII1DAJAJ=DAII1DASI1*(DASI1DAI1*DAI1DAJAJ+DASI1DAJAJ)
!      
      DAII1DBI1AJ=DAII1DASI1*DASI1DAI1*DAI1DBI1AJ
!      
      DAII1DBI1AJAJ=DAII1DASI1*DASI1DAI1*DAI1DBI1AJAJ
!      
      DAII1DAJAJAJ=DAII1DASI1*(DASI1DAI1*DAI1DAJAJAJ+DASI1DAJAJAJ)
!      
!     NON-GAUSSIAN HYPERELASTIC MODEL FOR THE ELASTOMER [2] 
!     AND DERIVATIVES
!      
      AP1 = AMU1*3.**(1.-ALPHA1)*0.5
      AP2 = AMU2*3.**(1.-ALPHA2)*0.5   
!      
      PSI = AP1/ALPHA1*(AII1**ALPHA1-3.**ALPHA1)+
     1      AP2/ALPHA2*(AII1**ALPHA2-3.**ALPHA2)
!     
      DPSI = AP1*AII1**(ALPHA1-1.)+AP2*AII1**(ALPHA2-1.)
!     
      DDPSI = AP1*(ALPHA1-1.)*AII1**(ALPHA1-2.)+
     1        AP2*(ALPHA2-1.)*AII1**(ALPHA2-2.)  
!     
      DDDPSI = AP1*(ALPHA1-1.)*(ALPHA1-2.)*AII1**(ALPHA1-3.)+
     1         AP2*(ALPHA2-1.)*(ALPHA2-2.)*AII1**(ALPHA2-3.)  
!      
!     MOREAU-YOSIDA REGULARIZATION FOR THE CONSTRAINT ACT>0
!      
      ANU=1.0e15      
!      
      AMYR=ANU/2*(ABS(ACT)-ACT)**2.
!      
      DAMYRDAJ=ANU*2*(ACT-ABS(ACT)) 
!           
      DAMYRDAJAJ = 0
      IF (ACT.LE.0) DAMYRDAJAJ = ANU*4.
!            
      DAMYRDAJAJAJ=0.
!
!     MACROSCOPIC STRAIN ENERGY DENSITY FUNCTION FOR 
!     POROUS ELASTOMER [1]
!            
      U(1) = AOMAC*PSI + AMYR
      U(2) = 0.
!
!     FIRST PARTIAL DERIVATIVES
!
      UI1(1) = AOMAC*DPSI*DAII1DBI1
      UI1(2) = 0.
      UI1(3) = AOMAC*DPSI*DAII1DAJ + DAMYRDAJ
!
!     SECOND PARTIAL DERIVATIVES
!
      UI2(1) = AOMAC*DDPSI*DAII1DBI1**2.
      UI2(2) = 0.
      UI2(3) = AOMAC*DDPSI*DAII1DAJ**2.+
     1  AOMAC*DPSI*DAII1DAJAJ + DAMYRDAJAJ
      UI2(4) = 0.
      UI2(5) = AOMAC*DDPSI*DAII1DBI1*DAII1DAJ+
     1 AOMAC*DPSI*DAII1DBI1AJ
      UI2(6) = 0.
!
!     THIRD PARTIAL DERIVATIVES
!
      UI3(1) = AOMAC*DDDPSI*DAII1DBI1**2.*DAII1DAJ+
     1     2.*AOMAC*DDPSI*DAII1DBI1*DAII1DBI1AJ
      UI3(2) = 0.
      UI3(3) = 0.
      UI3(4) = AOMAC*DDDPSI*DAII1DBI1*DAII1DAJ**2.+
     1     2.*AOMAC*DDPSI*DAII1DBI1AJ*DAII1DAJ+
     2     AOMAC*DDPSI*DAII1DBI1*DAII1DAJAJ+
     3     AOMAC*DPSI*DAII1DBI1AJAJ
      UI3(5) = 0.
      UI3(6) = AOMAC*DDDPSI*DAII1DAJ**3.+
     1     3.*AOMAC*DDPSI*DAII1DAJ*DAII1DAJAJ+
     2    AOMAC*DPSI*DAII1DAJAJAJ + DAMYRDAJAJAJ
!
!     SOLUTION-DEPENDENT STATE VARIABLE
!     CURRENT POROSITY 
!      
      STATEV(1)=ACT/AJ 
!      
      RETURN
      END
!
!**********************************************************************
