!**********************************************************************
! Legal notice: UHYPER_Shrimali_Lefevre_Lopez-Pamies.for (Windows) and
!               UHYPER_Shrimali_Lefevre_Lopez-Pamies.for (Linux)
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
! along with this program.  If not, see https://www.gnu.org/licenses/.
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
!  AC      = PROPS(5)  ! INITIAL POROSITY 
!
! The two material parameters AMU1, AMU2 characterizing the elastic    
! behavior of the underlying elastomer are non-negative real numbers 
! (AMU1 >= 0, AMU2 >= 0). The two exponents ALPHA1, ALPHA2 are real 
! numbers leading to a strongly elliptic strain energy (see eq. (22)
! in [2]). This is left to the user to check.
!
! The initial porosity (AC) must satisfy 0 <= AC <= 1. 
!
! As expected from physical considerations, this macroscopic energy
! remains finite so long the determinant of the deformation 
! gradient (AJ) satisfies the condition AJ - 1 + AC > 0. This
! inequality constraint is enforced through a MOREAU-YOSIDA
! regularization (see, e.g. [3]). The underlying weight (ANU) is set
! here by default to ANU = 1,000,000.
!
! The porosity in the deformed configuration (AJ - 1 + AC) / AJ
! may be output (to be plotted for instance) as a solution-dependent 
! state variable (SDV), e.g., using the lines
! *DEPVAR
! 1
! 1, Porosity, Current local porosity
! in the input (.inp) file. If needed, this solution-dependent state 
! variable can be initialized using the lines
! *INITIAL CONDITIONS, TYPE=SOLUTION
! <some element set>, AC
! in the input (.inp) file.
!
!**********************************************************************
! Additional information:
!
! This subroutine can create a solution-dependent state variable for
! for the current porosity (see above) but does not create predefined 
! field variables. 
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