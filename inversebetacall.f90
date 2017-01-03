! ==============================================================================
! 		    Inverse Beta Example Calling Code version 1.0	       !
! ==============================================================================
!
! AUTHOR: David Kipping
!         Harvard-Smithsonian Center for Astrophysics
!         Please report any problems to: dkipping@cfa.harvard.edu
!
! CITATION: If using this code, please cite:
!           [1] Kipping, D. M.,  2013, 'Parametrizing the exoplanet eccentricity 
!               distribution with the Beta distribution', MNRAS, accepted
!           [2] Allan Macleod, Algorithm AS 245, 'A Robust and Reliable 
!               Algorithm for the Logarithm of the Gamma Function', Applied 
!               Statistics, Volume 38, Number 2, 1989, pages 397-402.
!
! DESCRIPTION: This is an example routine which calls the inversebeta.f module,
!              which in turn computes the inverse cumulative density function 
!              (CDF) of a Beta distribution. Evaluating the inverse Beta CDF for 
!              a range of uniform random variates will reproduce a Beta
!              distribution and hence allows for direct sampling from a Beta
!              prior (see APPLICATIONS).
!
! COMPILING: gfortran -c inversebeta.f
!            gfortran -o inversebetacall inversebetacall.f90 inversebeta.o
!
! EXECUTING: ./inversebetacall
!
! INPUTS: a, b, and n should be changed as required
!
! OUTPUTS: A file "beta_variates.dat" will be outputted, containing n random
!          variates drawn from a Beta distribution with shape parameters a & b.
!          If you have Mathematica, you can plot the resulting distribution of
!          Beta variates using exampleplot.nb
!
! APPLICATIONS: One may enforce a Beta prior on a fitted parameter, such as
!               an exoplanet's orbital eccentricity, by fitting for parameter z 
!               in a uniformly between 0<z<1. At each realization, this z value 
!               is converted to x (where x=eccentricity in the aforementioned 
!               example), for a given choice of a & b. These two latter terms 
!               define the shape of the Beta prior. a = 0.867 and b = 3.030 is 
!               recommended as a general prior for exoplanet orbital 
!               eccentricity.
!
! CHANGES: v1.0 Initial version released

PROGRAM inversebetacall

use inversebetamod

  implicit none

 REAL(8), PARAMETER :: a =  0.867D0 ! 1st shape parameter of Beta distribution
 REAL(8), PARAMETER :: b =  3.030D0 ! 2nd shape parameter of Beta distribution
 REAL(8) :: beta_log
 INTEGER :: i, ifault
 INTEGER, PARAMETER :: n = 1D5 ! Number of draws to compute
 REAL(8), DIMENSION(n) :: z ! Inputted z value 
 REAL(8), DIMENSION(n) :: x ! Outputted Beta variate

 ! Generate some z values, in this case simply a random uniform variate
 call random_flat(n,z) 

 OPEN(UNIT=11,FILE='beta_variates.dat',STATUS='UNKNOWN')
 ! Compute inverse beta
 beta_log = alngam(a,ifault)+alngam(b,ifault)-alngam(a+b,ifault)
 DO i=1,n
   x(i) = xinbta ( a, b, beta_log, z(i), ifault )
   write(11,*) x(i) ! Output
 END DO
 CLOSE(11)

CONTAINS

! =======================================================
 SUBROUTINE random_flat(n,r)
 INTEGER i, n
 DOUBLE PRECISION seeda, r(n)

 call random_seed()
 DO i=1,n
   call random_number(seeda)
   r(i)=seeda
 END DO

 END SUBROUTINE random_flat
! =======================================================

END PROGRAM inversebetacall
