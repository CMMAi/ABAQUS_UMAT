      SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd,
     1 rpl, ddsddt, drplde, drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c  WRITE (6,*) '
c  NOTE:  MODIFICATIONS TO *UMAT FOR ABAQUS VERSION 5.3 (14 APR '94)
c 
c  (1)  The list of variables above defining the *UMAT subroutine, 
c  and the first (standard) block of variables dimensioned below, 
c  have variable names added compared to earlier ABAQUS versions. 
c
c  (2)  The statement: include 'aba_param.inc' must be added as below.
c
c  (3)  As of version 5.3, ABAQUS files use double precision only.
c  The file aba_param.inc has a line "implicit real*8" and, since 
c  it is included in the main subroutine, it will define the variables
c  there as double precision.  But other subroutines still need the
c  definition "implicit real*8" since there may be variables that are
c  not passed to them through the list or common block.
c
c  (4)  This is current as of version 5.6 of ABAQUS.
c
c  (5)  Note added by J. W. Kysar (4 November 1997).  This UMAT has been
c   modified to keep track of the cumulative shear strain in each 
c   individual slip system.  This information is needed to correct an
c   error in the implementation of the Bassani and Wu hardening law.
c   Any line of code which has been added or modified is preceded
c   immediately by a line beginning CFIXA and succeeded by a line 
c   beginning CFIXB.  Any comment line added or modified will begin 
c   with CFIX.
c
c   The hardening law by Bassani and Wu was implemented incorrectly.
c   This law is a function of both hyperbolic secant squared and hyperbolic
c   tangent.  However, the arguments of sech and tanh are related to the *total*
c   slip on individual slip systems.  Formerly, the UMAT implemented this
c   hardening law by using the *current* slip on each slip system.  Therein 
c   lay the problem. The UMAT did not restrict the current slip to be a 
c   positive value.  So when a slip with a negative sign was encountered, the 
c   term containing tanh led to a negative hardening rate (since tanh is an 
c   odd function).

c   The UMAT has been fixed by adding state variables to keep track of the 
c   *total* slip on each slip system by integrating up the absolute value 
c   of slip rates for each individual slip system.  These "solution dependent 
c   variables" are available for postprocessing.  The only required change 
c   in the input file is that the DEPVAR command must be changed.
c 
C-----  Use single precision on Cray by
C     (1) deleting the statement "IMPLICIT*8 (A-H,O-Z)";
C     (2) changing "REAL*8 FUNCTION" to "FUNCTION";
C     (3) changing double precision functions DSIGN to SIGN.
C
C-----  Subroutines:
C
C       ROTATION     -- forming rotation matrix, i.e. the direction 
C                       cosines of cubic crystal [100], [010] and [001]
C                       directions in global system at the initial 
C                       state
C
C       SLIPSYS      -- calculating number of slip systems, unit 
C                       vectors in slip directions and unit normals to 
C                       slip planes in a cubic crystal at the initial 
C                       state
C
C       GSLPINIT     -- calculating initial value of current strengths 
C                       at initial state
C
C       STRAINRATE   -- based on current values of resolved shear 
C                       stresses and current strength, calculating 
C                       shear strain-rates in slip systems
C
C       LATENTHARDEN -- forming self- and latent-hardening matrix 
C
C       ITERATION    -- generating arrays for the Newton-Rhapson 
C                       iteration
C
C       LUDCMP       -- LU decomposition
C
C       LUBKSB       -- linear equation solver based on LU 
C                       decomposition method (must call LUDCMP first)


C-----  Function subprogram:

C       F -- shear strain-rates in slip systems


C-----  Variables:
C
C       STRESS -- stresses (INPUT & OUTPUT)
C                 Cauchy stresses for finite deformation
C       STATEV -- solution dependent state variables (INPUT & OUTPUT)
C       DDSDDE -- Jacobian matrix (OUTPUT)

C-----  Variables passed in for information:
C
C       STRAN  -- strains
C                 logarithmic strain for finite deformation 
C                 (actually, integral of the symmetric part of velocity
C                  gradient with respect to time)
C       DSTRAN -- increments of strains
C       CMNAME -- name given in the *MATERIAL option
C       NDI    -- number of direct stress components
C       NSHR   -- number of engineering shear stress components
C       NTENS  -- NDI+NSHR
C       NSTATV -- number of solution dependent state variables (as 
C                 defined in the *DEPVAR option)
C       PROPS  -- material constants entered in the *USER MATERIAL 
C                 option
C       NPROPS -- number of material constants
C

C-----  This subroutine provides the plastic constitutive relation of 
C     single crystals for finite element code ABAQUS.  The plastic slip
C     of single crystal obeys the Schmid law.  The program gives the 
C     choice of small deformation theory and theory of finite rotation 
C     and finite strain.
C       The strain increment is composed of elastic part and plastic 
C     part.  The elastic strain increment corresponds to lattice 
C     stretching, the plastic part is the sum over all slip systems of 
C     plastic slip.  The shear strain increment for each slip system is
C     assumed a function of the ratio of corresponding resolved shear 
C     stress over current strength, and of the time step.  The resolved
C     shear stress is the double product of stress tensor with the slip
C     deformation tensor (Schmid factor), and the increment of current 
C     strength is related to shear strain increments over all slip 
C     systems through self- and latent-hardening functions.

C-----  The implicit integration method proposed by Peirce, Shih and 
C     Needleman (1984) is used here.  The subroutine provides an option
C     of iteration to solve stresses and solution dependent state 
C     variables within each increment.

C-----  The present program is for a single CUBIC crystal.  However, 
C     this code can be generalized for other crystals (e.g. HCP, 
C     Tetragonal, Orthotropic, etc.).  Only subroutines ROTATION and 
C     SLIPSYS need to be modified to include the effect of crystal 
C     aspect ratio.
C

C-----  Important notice:
C
C     (1) The number of state variables NSTATV must be larger than (or 
CFIX      equal to) TEN (10) times the total number of slip systems in
C         all sets, NSLPTL, plus FIVE (5)
CFIX           NSTATV >= 10 * NSLPTL + 5
C         Denote s as a slip direction and m as normal to a slip plane.
C         Here (s,-m), (-s,m) and (-s,-m) are NOT considered 
C         independent of (s,m).  The number of slip systems in each set
C         could be either 6, 12, 24 or 48 for a cubic crystal, e.g. 12 
C         for {110}<111>.
C
C         Users who need more parameters to characterize the 
C         constitutive law of single crystal, e.g. the framework 
C         proposed by Zarka, should make NSTATV larger than (or equal 
C         to) the number of those parameters NPARMT plus nine times 
C         the total number of slip systems, NSLPTL, plus five
CFIX           NSTATV >= NPARMT + 10 * NSLPTL + 5
C
C     (2) The tangent stiffness matrix in general is not symmetric if 
C         latent hardening is considered.  Users must declare "UNSYMM" 
C         in the input file, at the *USER MATERIAL card.
C

      PARAMETER (ND=150)
C-----  The parameter ND determines the dimensions of the arrays in 
C     this subroutine.  The current choice 150 is a upper bound for a 
C     cubic crystal with up to three sets of slip systems activated.  
C     Users may reduce the parameter ND to any number as long as larger
C     than or equal to the total number of slip systems in all sets.  
C     For example, if {110}<111> is the only set of slip system 
C     potentially activated, ND could be taken as twelve (12).  
c
      include 'aba_param.inc'
c
      CHARACTER*8 CMNAME
      EXTERNAL F

      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      DIMENSION ISPDIR(3), ISPNOR(3), NSLIP(3), 
     2          SLPDIR(3,ND), SLPNOR(3,ND), SLPDEF(6,ND), 
     3          SLPSPN(3,ND), DSPDIR(3,ND), DSPNOR(3,ND), 
     4          DLOCAL(6,6), D(6,6), ROTD(6,6), ROTATE(3,3), 
     5          FSLIP(ND), DFDXSP(ND), DDEMSD(6,ND), 
     6          H(ND,ND), DDGDDE(ND,6), 
     7          DSTRES(6), DELATS(6), DSPIN(3), DVGRAD(3,3),
     8          DGAMMA(ND), DTAUSP(ND), DGSLIP(ND), 
     9          WORKST(ND,ND), INDX(ND), TERM(3,3), TRM0(3,3), ITRM(3)

      DIMENSION FSLIP1(ND), STRES1(6), GAMMA1(ND), TAUSP1(ND), 
     2          GSLP1(ND), SPNOR1(3,ND), SPDIR1(3,ND), DDSDE1(6,6),
     3          DSOLD(6), DGAMOD(ND), DTAUOD(ND), DGSPOD(ND), 
     4          DSPNRO(3,ND), DSPDRO(3,ND), 
     5          DHDGDG(ND,ND)


C-----  NSLIP  -- number of slip systems in each set
C-----  SLPDIR -- slip directions (unit vectors in the initial state)
C-----  SLPNOR -- normals to slip planes (unit normals in the initial 
C                 state)
C-----  SLPDEF -- slip deformation tensors (Schmid factors)
C                 SLPDEF(1,i) -- SLPDIR(1,i)*SLPNOR(1,i)
C                 SLPDEF(2,i) -- SLPDIR(2,i)*SLPNOR(2,i)
C                 SLPDEF(3,i) -- SLPDIR(3,i)*SLPNOR(3,i)
C                 SLPDEF(4,i) -- SLPDIR(1,i)*SLPNOR(2,i)+
C                                SLPDIR(2,i)*SLPNOR(1,i)
C                 SLPDEF(5,i) -- SLPDIR(1,i)*SLPNOR(3,i)+
C                                SLPDIR(3,i)*SLPNOR(1,i)
C                 SLPDEF(6,i) -- SLPDIR(2,i)*SLPNOR(3,i)+
C                                SLPDIR(3,i)*SLPNOR(2,i)
C                 where index i corresponds to the ith slip system
C-----  SLPSPN -- slip spin tensors (only needed for finite rotation)
C                 SLPSPN(1,i) -- [SLPDIR(1,i)*SLPNOR(2,i)-
C                                 SLPDIR(2,i)*SLPNOR(1,i)]/2
C                 SLPSPN(2,i) -- [SLPDIR(3,i)*SLPNOR(1,i)-
C                                 SLPDIR(1,i)*SLPNOR(3,i)]/2
C                 SLPSPN(3,i) -- [SLPDIR(2,i)*SLPNOR(3,i)-
C                                 SLPDIR(3,i)*SLPNOR(2,i)]/2
C                 where index i corresponds to the ith slip system
C-----  DSPDIR -- increments of slip directions
C-----  DSPNOR -- increments of normals to slip planes
C
C-----  DLOCAL -- elastic matrix in local cubic crystal system
C-----  D      -- elastic matrix in global system
C-----  ROTD   -- rotation matrix transforming DLOCAL to D
C
C-----  ROTATE -- rotation matrix, direction cosines of [100], [010] 
C                 and [001] of cubic crystal in global system
C
C-----  FSLIP  -- shear strain-rates in slip systems
C-----  DFDXSP -- derivatives of FSLIP w.r.t x=TAUSLP/GSLIP, where 
C                 TAUSLP is the resolved shear stress and GSLIP is the 
C                 current strength
C
C-----  DDEMSD -- double dot product of the elastic moduli tensor with 
C                 the slip deformation tensor plus, only for finite 
C                 rotation, the dot product of slip spin tensor with 
C                 the stress
C
C-----  H      -- self- and latent-hardening matrix
C                 H(i,i) -- self hardening modulus of the ith slip 
C                           system (no sum over i)
C                 H(i,j) -- latent hardening molulus of the ith slip 
C                           system due to a slip in the jth slip system
C                           (i not equal j)
C
C-----  DDGDDE -- derivatice of the shear strain increments in slip 
C                 systems w.r.t. the increment of strains
C
C-----  DSTRES -- Jaumann increments of stresses, i.e. corotational 
C                 stress-increments formed on axes spinning with the 
C                 material
C-----  DELATS -- strain-increments associated with lattice stretching
C                 DELATS(1) - DELATS(3) -- normal strain increments
C                 DELATS(4) - DELATS(6) -- engineering shear strain 
C                                          increments
C-----  DSPIN  -- spin-increments associated with the material element
C                 DSPIN(1) -- component 12 of the spin tensor
C                 DSPIN(2) -- component 31 of the spin tensor
C                 DSPIN(3) -- component 23 of the spin tensor
C
C-----  DVGRAD -- increments of deformation gradient in the current 
C                 state, i.e. velocity gradient times the increment of 
C                 time
C
C-----  DGAMMA -- increment of shear strains in slip systems
C-----  DTAUSP -- increment of resolved shear stresses in slip systems 
C-----  DGSLIP -- increment of current strengths in slip systems
C
C
C-----  Arrays for iteration:
C
C            FSLIP1, STRES1, GAMMA1, TAUSP1, GSLP1 , SPNOR1, SPDIR1, 
C            DDSDE1, DSOLD , DGAMOD, DTAUOD, DGSPOD, DSPNRO, DSPDRO,
C            DHDGDG
C
C
C-----  Solution dependent state variable STATEV:
C            Denote the number of total slip systems by NSLPTL, which 
C            will be calculated in this code.
C
C       Array STATEV:
C       1          - NSLPTL    :  current strength in slip systems
C       NSLPTL+1   - 2*NSLPTL  :  shear strain in slip systems
C       2*NSLPTL+1 - 3*NSLPTL  :  resolved shear stress in slip systems
C
C       3*NSLPTL+1 - 6*NSLPTL  :  current components of normals to slip
C                                 slip planes
C       6*NSLPTL+1 - 9*NSLPTL  :  current components of slip directions
C
CFIX    9*NSLPTL+1 - 10*NSLPTL :  total cumulative shear strain on each 
CFIX                              slip system (sum of the absolute 
CFIX                              values of shear strains in each slip 
CFIX                              system individually)
CFIX
CFIX    10*NSLPTL+1             : total cumulative shear strain on all 
C                                 slip systems (sum of the absolute 
C                                 values of shear strains in all slip 
C                                 systems)
C
CFIX    10*NSLPTL+2 - NSTATV-4  : additional parameters users may need 
C                                 to characterize the constitutive law 
C                                 of a single crystal (if there are 
C                                 any).
C
C       NSTATV-3               :  number of slip systems in the 1st set
C       NSTATV-2               :  number of slip systems in the 2nd set
C       NSTATV-1               :  number of slip systems in the 3rd set
C       NSTATV                 :  total number of slip systems in all 
C                                 sets
C
C
C-----  Material constants PROPS:
C
C       PROPS(1) - PROPS(21) -- elastic constants for a general elastic
C                               anisotropic material
C
C            isotropic   : PROPS(i)=0  for  i>2
C                          PROPS(1) -- Young's modulus
C                          PROPS(2) -- Poisson's ratio
C
C            cubic       : PROPS(i)=0  for i>3
C                          PROPS(1) -- c11
C                          PROPS(2) -- c12
C                          PROPS(3) -- c44
C
C            orthotropic : PORPS(i)=0  for  i>9
C                          PROPS(1) - PROPS(9) are D1111, D1122, D2222,
C                          D1133, D2233, D3333, D1212, D1313, D2323, 
C                          respectively, which has the same definition 
C                          as ABAQUS for orthotropic materials
C                          (see *ELASTIC card)
C
C            anisotropic : PROPS(1) - PROPS(21) are D1111, D1122, 
C                          D2222, D1133, D2233, D3333, D1112, D2212, 
C                          D3312, D1212, D1113, D2213, D3313, D1213, 
C                          D1313, D1123, D2223, D3323, D1223, D1323, 
C                          D2323, respectively, which has the same 
C                          definition as ABAQUS for anisotropic 
C                          materials (see *ELASTIC card)
C
C
C       PROPS(25) - PROPS(56) -- parameters characterizing all slip 
C                                systems to be activated in a cubic 
C                                crystal
C
C            PROPS(25) -- number of sets of slip systems (maximum 3), 
C                         e.g. (110)[1-11] and (101)[11-1] are in the 
C                         same set of slip systems, (110)[1-11] and 
C                         (121)[1-11] belong to different sets of slip 
C                         systems
C                         (It must be a real number, e.g. 3., not 3 !)
C
C            PROPS(33) - PROPS(35) -- normal to a typical slip plane in
C                                     the first set of slip systems, 
C                                     e.g. (1 1 0)
C                                     (They must be real numbers, e.g. 
C                                      1. 1. 0., not 1 1 0 !)
C            PROPS(36) - PROPS(38) -- a typical slip direction in the 
C                                     first set of slip systems, e.g. 
C                                     [1 1 1]
C                                     (They must be real numbers, e.g. 
C                                      1. 1. 1., not 1 1 1 !)
C
C            PROPS(41) - PROPS(43) -- normal to a typical slip plane in
C                                     the second set of slip systems
C                                     (real numbers)
C            PROPS(44) - PROPS(46) -- a typical slip direction in the 
C                                     second set of slip systems
C                                     (real numbers)
C
C            PROPS(49) - PROPS(51) -- normal to a typical slip plane in
C                                     the third set of slip systems
C                                     (real numbers)
C            PROPS(52) - PROPS(54) -- a typical slip direction in the 
C                                     third set of slip systems
C                                     (real numbers)
C
C
C       PROPS(57) - PROPS(72) -- parameters characterizing the initial 
C                                orientation of a single crystal in 
C                                global system
C            The directions in global system and directions in local 
C            cubic crystal system of two nonparallel vectors are needed
C            to determine the crystal orientation.
C
C            PROPS(57) - PROPS(59) -- [p1 p2 p3], direction of first 
C                                     vector in local cubic crystal 
C                                     system, e.g. [1 1 0]
C                                     (They must be real numbers, e.g. 
C                                      1. 1. 0., not 1 1 0 !)
C            PROPS(60) - PROPS(62) -- [P1 P2 P3], direction of first 
C                                     vector in global system, e.g. 
C                                     [2. 1. 0.]
C                                     (It does not have to be a unit 
C                                      vector)
C
C            PROPS(65) - PROPS(67) -- direction of second vector in 
C                                     local cubic crystal system (real 
C                                     numbers)
C            PROPS(68) - PROPS(70) -- direction of second vector in 
C                                     global system
C
C
C       PROPS(73) - PROPS(96) -- parameters characterizing the visco-
C                                plastic constitutive law (shear 
C                                strain-rate vs. resolved shear 
C                                stress), e.g. a power-law relation
C
C            PROPS(73) - PROPS(80) -- parameters for the first set of 
C                                     slip systems
C            PROPS(81) - PROPS(88) -- parameters for the second set of 
C                                     slip systems
C            PROPS(89) - PROPS(96) -- parameters for the third set of 
C                                     slip systems
C
C
C       PROPS(97) - PROPS(144)-- parameters characterizing the self-
C                                and latent-hardening laws of slip 
C                                systems
C
C            PROPS(97) - PROPS(104)-- self-hardening parameters for the
C                                     first set of slip systems
C            PROPS(105)- PROPS(112)-- latent-hardening parameters for 
C                                     the first set of slip systems and
C                                     interaction with other sets of 
C                                     slip systems
C
C            PROPS(113)- PROPS(120)-- self-hardening parameters for the
C                                     second set of slip systems
C            PROPS(121)- PROPS(128)-- latent-hardening parameters for 
C                                     the second set of slip systems 
C                                     and interaction with other sets 
C                                     of slip systems
C
C            PROPS(129)- PROPS(136)-- self-hardening parameters for the
C                                     third set of slip systems
C            PROPS(137)- PROPS(144)-- latent-hardening parameters for 
C                                     the third set of slip systems and
C                                     interaction with other sets of
C                                     slip systems
C
C
C       PROPS(145)- PROPS(152)-- parameters characterizing forward time
C                                integration scheme and finite 
C                                deformation
C
C            PROPS(145) -- parameter theta controlling the implicit 
C                          integration, which is between 0 and 1
C                          0.  : explicit integration
C                          0.5 : recommended value
C                          1.  : fully implicit integration
C
C            PROPS(146) -- parameter NLGEOM controlling whether the 
C                          effect of finite rotation and finite strain 
C                          of crystal is considered,
C                          0.        : small deformation theory
C                          otherwise : theory of finite rotation and 
C                                      finite strain
C
C
C       PROPS(153)- PROPS(160)-- parameters characterizing iteration 
C                                method
C
C            PROPS(153) -- parameter ITRATN controlling whether the 
C                          iteration method is used, 
C                          0.        : no iteration
C                          otherwise : iteration
C
C            PROPS(154) -- maximum number of iteration ITRMAX 
C
C            PROPS(155) -- absolute error of shear strains in slip 
C                          systems GAMERR
C


C-----  Elastic matrix in local cubic crystal system: DLOCAL
      DO J=1,6
         DO I=1,6
            DLOCAL(I,J)=0.
         END DO
      END DO

      CHECK=0.
      DO J=10,21
         CHECK=CHECK+ABS(PROPS(J))
      END DO

      IF (CHECK.EQ.0.) THEN
         DO J=4,9
            CHECK=CHECK+ABS(PROPS(J))
         END DO

         IF (CHECK.EQ.0.) THEN

            IF (PROPS(3).EQ.0.) THEN

C-----  Isotropic material
               GSHEAR=PROPS(1)/2./(1.+PROPS(2))
               E11=2.*GSHEAR*(1.-PROPS(2))/(1.-2.*PROPS(2))
               E12=2.*GSHEAR*PROPS(2)/(1.-2.*PROPS(2))

               DO J=1,3
                  DLOCAL(J,J)=E11

                  DO I=1,3
                     IF (I.NE.J) DLOCAL(I,J)=E12
                  END DO

                  DLOCAL(J+3,J+3)=GSHEAR
               END DO

            ELSE

C-----  Cubic material
               DO J=1,3
                  DLOCAL(J,J)=PROPS(1)

                  DO I=1,3
                     IF (I.NE.J) DLOCAL(I,J)=PROPS(2)
                  END DO

                  DLOCAL(J+3,J+3)=PROPS(3)
               END DO

            END IF

         ELSE

C-----  Orthotropic metarial
            DLOCAL(1,1)=PROPS(1)
            DLOCAL(1,2)=PROPS(2)
            DLOCAL(2,1)=PROPS(2)
            DLOCAL(2,2)=PROPS(3)

            DLOCAL(1,3)=PROPS(4)
            DLOCAL(3,1)=PROPS(4)
            DLOCAL(2,3)=PROPS(5)
            DLOCAL(3,2)=PROPS(5)
            DLOCAL(3,3)=PROPS(6)

            DLOCAL(4,4)=PROPS(7)
            DLOCAL(5,5)=PROPS(8)
            DLOCAL(6,6)=PROPS(9)

         END IF

      ELSE

C-----  General anisotropic material
         ID=0
         DO J=1,6
            DO I=1,J
               ID=ID+1
               DLOCAL(I,J)=PROPS(ID)
               DLOCAL(J,I)=DLOCAL(I,J)
            END DO
         END DO
      END IF

C-----  Rotation matrix: ROTATE, i.e. direction cosines of [100], [010]
C     and [001] of a cubic crystal in global system
C
      CALL ROTATION (PROPS(57), ROTATE)

C-----  Rotation matrix: ROTD to transform local elastic matrix DLOCAL 
C     to global elastic matrix D
C
      DO J=1,3
         J1=1+J/3
         J2=2+J/2

         DO I=1,3
            I1=1+I/3
            I2=2+I/2

            ROTD(I,J)=ROTATE(I,J)**2
            ROTD(I,J+3)=2.*ROTATE(I,J1)*ROTATE(I,J2)
            ROTD(I+3,J)=ROTATE(I1,J)*ROTATE(I2,J)
            ROTD(I+3,J+3)=ROTATE(I1,J1)*ROTATE(I2,J2)+
     2                    ROTATE(I1,J2)*ROTATE(I2,J1)

         END DO
      END DO

C-----  Elastic matrix in global system: D
C     {D} = {ROTD} * {DLOCAL} * {ROTD}transpose
C
      DO J=1,6
         DO I=1,6
            D(I,J)=0.
         END DO
      END DO

      DO J=1,6
         DO I=1,J

            DO K=1,6
               DO L=1,6
                  D(I,J)=D(I,J)+DLOCAL(K,L)*ROTD(I,K)*ROTD(J,L)
               END DO
            END DO

            D(J,I)=D(I,J)

         END DO
      END DO

C-----  Total number of sets of slip systems: NSET
      NSET=NINT(PROPS(25))
      IF (NSET.LT.1) THEN
         WRITE (6,*) '***ERROR - zero sets of slip systems'
         STOP
      ELSE IF (NSET.GT.3) THEN
         WRITE (6,*) 
     2     '***ERROR - more than three sets of slip systems'
         STOP
      END IF

C-----  Implicit integration parameter: THETA
      THETA=PROPS(145)

C-----  Finite deformation ?
C-----  NLGEOM = 0,   small deformation theory
C       otherwise, theory of finite rotation and finite strain, Users 
C     must declare "NLGEOM" in the input file, at the *STEP card
C
      IF (PROPS(146).EQ.0.) THEN
         NLGEOM=0
      ELSE
         NLGEOM=1
      END IF

C-----  Iteration?
C-----  ITRATN = 0, no iteration
C       otherwise, iteration (solving increments of stresses and 
C     solution dependent state variables)
C
      IF (PROPS(153).EQ.0.) THEN
         ITRATN=0
      ELSE
         ITRATN=1
      END IF

      ITRMAX=NINT(PROPS(154))
      GAMERR=PROPS(155)

      NITRTN=-1

      DO I=1,NTENS
         DSOLD(I)=0.
      END DO

      DO J=1,ND
         DGAMOD(J)=0.
         DTAUOD(J)=0.
         DGSPOD(J)=0.
         DO I=1,3
            DSPNRO(I,J)=0.
            DSPDRO(I,J)=0.
         END DO
      END DO

C-----  Increment of spin associated with the material element: DSPIN
C     (only needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               TERM(I,J)=DROT(J,I)
               TRM0(I,J)=DROT(J,I)
            END DO

            TERM(J,J)=TERM(J,J)+1.D0
            TRM0(J,J)=TRM0(J,J)-1.D0
         END DO

         CALL LUDCMP (TERM, 3, 3, ITRM, DDCMP)

         DO J=1,3
            CALL LUBKSB (TERM, 3, 3, ITRM, TRM0(1,J))
         END DO

         DSPIN(1)=TRM0(2,1)-TRM0(1,2)
         DSPIN(2)=TRM0(1,3)-TRM0(3,1)
         DSPIN(3)=TRM0(3,2)-TRM0(2,3)

      END IF

C-----  Increment of dilatational strain: DEV
      DEV=0.D0
      DO I=1,NDI
         DEV=DEV+DSTRAN(I)
      END DO

C-----  Iteration starts (only when iteration method is used)
1000  CONTINUE

C-----  Parameter NITRTN: number of iterations
C       NITRTN = 0 --- no-iteration solution
C
      NITRTN=NITRTN+1

C-----  Check whether the current stress state is the initial state
      IF (STATEV(1).EQ.0.) THEN

C-----  Initial state
C
C-----  Generating the following parameters and variables at initial 
C     state:
C          Total number of slip systems in all the sets NSLPTL
C          Number of slip systems in each set NSLIP
C          Unit vectors in initial slip directions SLPDIR
C          Unit normals to initial slip planes SLPNOR
C
         NSLPTL=0
         DO I=1,NSET
            ISPNOR(1)=NINT(PROPS(25+8*I))
            ISPNOR(2)=NINT(PROPS(26+8*I))
            ISPNOR(3)=NINT(PROPS(27+8*I))

            ISPDIR(1)=NINT(PROPS(28+8*I))
            ISPDIR(2)=NINT(PROPS(29+8*I))
            ISPDIR(3)=NINT(PROPS(30+8*I))

            CALL SLIPSYS (ISPDIR, ISPNOR, NSLIP(I), SLPDIR(1,NSLPTL+1), 
     2                    SLPNOR(1,NSLPTL+1), ROTATE)

            NSLPTL=NSLPTL+NSLIP(I)
         END DO

         IF (ND.LT.NSLPTL) THEN
            WRITE (6,*) 
     2 '***ERROR - parameter ND chosen by the present user is less than
     3             the total number of slip systems NSLPTL'
            STOP
         END IF

C-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
         END DO

C-----  Initial value of state variables: unit normal to a slip plane 
C     and unit vector in a slip direction
C
         STATEV(NSTATV)=FLOAT(NSLPTL)
         DO I=1,NSET
            STATEV(NSTATV-4+I)=FLOAT(NSLIP(I))
         END DO

         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               STATEV(IDNOR)=SLPNOR(I,J)

               IDDIR=IDDIR+1
               STATEV(IDDIR)=SLPDIR(I,J)
            END DO
         END DO

C-----  Initial value of the current strength for all slip systems
C
         CALL GSLPINIT (STATEV(1), NSLIP, NSLPTL, NSET, PROPS(97))

C-----  Initial value of shear strain in slip systems
CFIX--  Initial value of cumulative shear strain in each slip systems

         DO I=1,NSLPTL
            STATEV(NSLPTL+I)=0.
CFIXA
            STATEV(9*NSLPTL+I)=0.
CFIXB
         END DO

CFIXA
         STATEV(10*NSLPTL+1)=0.
CFIXB

C-----  Initial value of the resolved shear stress in slip systems
         DO I=1,NSLPTL
            TERM1=0.

            DO J=1,NTENS
               IF (J.LE.NDI) THEN
                  TERM1=TERM1+SLPDEF(J,I)*STRESS(J)
               ELSE
                  TERM1=TERM1+SLPDEF(J-NDI+3,I)*STRESS(J)
               END IF
            END DO

            STATEV(2*NSLPTL+I)=TERM1
         END DO

      ELSE

C-----  Current stress state
C
C-----  Copying from the array of state variables STATVE the following
C          parameters and variables at current stress state:
C          Total number of slip systems in all the sets NSLPTL
C          Number of slip systems in each set NSLIP
C          Current slip directions SLPDIR
C          Normals to current slip planes SLPNOR
C
         NSLPTL=NINT(STATEV(NSTATV))
         DO I=1,NSET
            NSLIP(I)=NINT(STATEV(NSTATV-4+I))
         END DO

         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               SLPNOR(I,J)=STATEV(IDNOR)

               IDDIR=IDDIR+1
               SLPDIR(I,J)=STATEV(IDDIR)
            END DO
         END DO

C-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
         END DO

      END IF

C-----  Slip spin tensor: SLPSPN (only needed for finite rotation)
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            SLPSPN(1,J)=0.5*(SLPDIR(1,J)*SLPNOR(2,J)-
     2                       SLPDIR(2,J)*SLPNOR(1,J))
            SLPSPN(2,J)=0.5*(SLPDIR(3,J)*SLPNOR(1,J)-
     2                       SLPDIR(1,J)*SLPNOR(3,J))
            SLPSPN(3,J)=0.5*(SLPDIR(2,J)*SLPNOR(3,J)-
     2                       SLPDIR(3,J)*SLPNOR(2,J))
         END DO
      END IF

C-----  Double dot product of elastic moduli tensor with the slip 
C     deformation tensor (Schmid factors) plus, only for finite 
C     rotation, the dot product of slip spin tensor with the stress: 
C     DDEMSD
C
      DO J=1,NSLPTL
         DO I=1,6
            DDEMSD(I,J)=0.
            DO K=1,6
               DDEMSD(I,J)=DDEMSD(I,J)+D(K,I)*SLPDEF(K,J)
            END DO
         END DO
      END DO

      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL

            DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(1,J)*STRESS(1)
            DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(2,J)*STRESS(1)

            IF (NDI.GT.1) THEN
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(1,J)*STRESS(2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(3,J)*STRESS(2)
            END IF

            IF (NDI.GT.2) THEN
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(2,J)*STRESS(3)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(3,J)*STRESS(3)
            END IF

            IF (NSHR.GE.1) THEN
               DDEMSD(1,J)=DDEMSD(1,J)+SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(2,J)=DDEMSD(2,J)-SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(3,J)*STRESS(NDI+1)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(2,J)*STRESS(NDI+1)
            END IF

            IF (NSHR.GE.2) THEN
               DDEMSD(1,J)=DDEMSD(1,J)-SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(3,J)=DDEMSD(3,J)+SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(3,J)*STRESS(NDI+2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(1,J)*STRESS(NDI+2)
            END IF

            IF (NSHR.EQ.3) THEN
               DDEMSD(2,J)=DDEMSD(2,J)+SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(3,J)=DDEMSD(3,J)-SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(2,J)*STRESS(NDI+3)
               DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(1,J)*STRESS(NDI+3)
            END IF

         END DO
      END IF

C-----  Shear strain-rate in a slip system at the start of increment: 
C     FSLIP, and its derivative: DFDXSP
C
      ID=1
      DO I=1,NSET
         IF (I.GT.1) ID=ID+NSLIP(I-1)
         CALL STRAINRATE (STATEV(NSLPTL+ID), STATEV(2*NSLPTL+ID), 
     2                    STATEV(ID), NSLIP(I), FSLIP(ID), DFDXSP(ID), 
     3                    PROPS(65+8*I))
      END DO

C-----  Self- and latent-hardening laws
CFIXA  
       CALL LATENTHARDEN (STATEV(NSLPTL+1), STATEV(2*NSLPTL+1), 
     2                   STATEV(1), STATEV(9*NSLPTL+1),
     3                   STATEV(10*NSLPTL+1), NSLIP, NSLPTL, 
     4                   NSET, H(1,1), PROPS(97), ND)
CFIXB

C-----  LU decomposition to solve the increment of shear strain in a 
C     slip system
C
      TERM1=THETA*DTIME
      DO I=1,NSLPTL
         TAUSLP=STATEV(2*NSLPTL+I)
         GSLIP=STATEV(I)
         X=TAUSLP/GSLIP
         TERM2=TERM1*DFDXSP(I)/GSLIP
         TERM3=TERM1*X*DFDXSP(I)/GSLIP

         DO J=1,NSLPTL
            TERM4=0.
            DO K=1,6
               TERM4=TERM4+DDEMSD(K,I)*SLPDEF(K,J)
            END DO

            WORKST(I,J)=TERM2*TERM4+H(I,J)*TERM3*DSIGN(1.D0,FSLIP(J))

            IF (NITRTN.GT.0) WORKST(I,J)=WORKST(I,J)+TERM3*DHDGDG(I,J)

         END DO

         WORKST(I,I)=WORKST(I,I)+1.
      END DO

      CALL LUDCMP (WORKST, NSLPTL, ND, INDX, DDCMP)

C-----  Increment of shear strain in a slip system: DGAMMA
      TERM1=THETA*DTIME
      DO I=1,NSLPTL

         IF (NITRTN.EQ.0) THEN
            TAUSLP=STATEV(2*NSLPTL+I)
            GSLIP=STATEV(I)
            X=TAUSLP/GSLIP
            TERM2=TERM1*DFDXSP(I)/GSLIP

            DGAMMA(I)=0.
            DO J=1,NDI
               DGAMMA(I)=DGAMMA(I)+DDEMSD(J,I)*DSTRAN(J)
            END DO

            IF (NSHR.GT.0) THEN
               DO J=1,NSHR
                  DGAMMA(I)=DGAMMA(I)+DDEMSD(J+3,I)*DSTRAN(J+NDI)
               END DO
            END IF

            DGAMMA(I)=DGAMMA(I)*TERM2+FSLIP(I)*DTIME

         ELSE
            DGAMMA(I)=TERM1*(FSLIP(I)-FSLIP1(I))+FSLIP1(I)*DTIME
     2                -DGAMOD(I)

         END IF

      END DO

      CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DGAMMA)

      DO I=1,NSLPTL
         DGAMMA(I)=DGAMMA(I)+DGAMOD(I)
      END DO

C-----  Update the shear strain in a slip system: STATEV(NSLPTL+1) - 
C     STATEV(2*NSLPTL)
C
      DO I=1,NSLPTL
         STATEV(NSLPTL+I)=STATEV(NSLPTL+I)+DGAMMA(I)-DGAMOD(I)
      END DO

C-----  Increment of current strength in a slip system: DGSLIP
      DO I=1,NSLPTL
         DGSLIP(I)=0.
         DO J=1,NSLPTL
            DGSLIP(I)=DGSLIP(I)+H(I,J)*ABS(DGAMMA(J))
         END DO
      END DO

C-----  Update the current strength in a slip system: STATEV(1) - 
C     STATEV(NSLPTL)
C
      DO I=1,NSLPTL
         STATEV(I)=STATEV(I)+DGSLIP(I)-DGSPOD(I)
      END DO

C-----  Increment of strain associated with lattice stretching: DELATS
      DO J=1,6
         DELATS(J)=0.
      END DO

      DO J=1,3
         IF (J.LE.NDI) DELATS(J)=DSTRAN(J)
         DO I=1,NSLPTL
            DELATS(J)=DELATS(J)-SLPDEF(J,I)*DGAMMA(I)
         END DO
      END DO

      DO J=1,3
         IF (J.LE.NSHR) DELATS(J+3)=DSTRAN(J+NDI)
         DO I=1,NSLPTL
            DELATS(J+3)=DELATS(J+3)-SLPDEF(J+3,I)*DGAMMA(I)
         END DO
      END DO

C-----  Increment of deformation gradient associated with lattice 
C     stretching in the current state, i.e. the velocity gradient 
C     (associated with lattice stretching) times the increment of time:
C     DVGRAD (only needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               IF (I.EQ.J) THEN
                  DVGRAD(I,J)=DELATS(I)
               ELSE
                  DVGRAD(I,J)=DELATS(I+J+1)
               END IF
            END DO
         END DO

         DO J=1,3
            DO I=1,J
               IF (J.GT.I) THEN
                  IJ2=I+J-2
                  IF (MOD(IJ2,2).EQ.1) THEN
                     TERM1=1.
                  ELSE
                     TERM1=-1.
                  END IF

                  DVGRAD(I,J)=DVGRAD(I,J)+TERM1*DSPIN(IJ2)
                  DVGRAD(J,I)=DVGRAD(J,I)-TERM1*DSPIN(IJ2)

                  DO K=1,NSLPTL
                     DVGRAD(I,J)=DVGRAD(I,J)-TERM1*DGAMMA(K)*
     2                                       SLPSPN(IJ2,K)
                     DVGRAD(J,I)=DVGRAD(J,I)+TERM1*DGAMMA(K)*
     2                                       SLPSPN(IJ2,K)
                  END DO
               END IF

            END DO
         END DO

      END IF

C-----  Increment of resolved shear stress in a slip system: DTAUSP
      DO I=1,NSLPTL
         DTAUSP(I)=0.
         DO J=1,6
            DTAUSP(I)=DTAUSP(I)+DDEMSD(J,I)*DELATS(J)
         END DO
      END DO

C-----  Update the resolved shear stress in a slip system: 
C     STATEV(2*NSLPTL+1) - STATEV(3*NSLPTL)
C
      DO I=1,NSLPTL
         STATEV(2*NSLPTL+I)=STATEV(2*NSLPTL+I)+DTAUSP(I)-DTAUOD(I)
      END DO

C-----  Increment of stress: DSTRES
      IF (NLGEOM.EQ.0) THEN
         DO I=1,NTENS
            DSTRES(I)=0.
         END DO
      ELSE
         DO I=1,NTENS
            DSTRES(I)=-STRESS(I)*DEV
         END DO
      END IF

      DO I=1,NDI
         DO J=1,NDI
            DSTRES(I)=DSTRES(I)+D(I,J)*DSTRAN(J)
         END DO

         IF (NSHR.GT.0) THEN
            DO J=1,NSHR
               DSTRES(I)=DSTRES(I)+D(I,J+3)*DSTRAN(J+NDI)
            END DO
         END IF

         DO J=1,NSLPTL
            DSTRES(I)=DSTRES(I)-DDEMSD(I,J)*DGAMMA(J)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO I=1,NSHR

            DO J=1,NDI
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J)*DSTRAN(J)
            END DO

            DO J=1,NSHR
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J+3)*DSTRAN(J+NDI)
            END DO

            DO J=1,NSLPTL
               DSTRES(I+NDI)=DSTRES(I+NDI)-DDEMSD(I+3,J)*DGAMMA(J)
            END DO

         END DO
      END IF

C-----  Update the stress: STRESS
      DO I=1,NTENS
         STRESS(I)=STRESS(I)+DSTRES(I)-DSOLD(I)
      END DO

C-----  Increment of normal to a slip plane and a slip direction (only 
C     needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            DO I=1,3
               DSPNOR(I,J)=0.
               DSPDIR(I,J)=0.

               DO K=1,3
                  DSPNOR(I,J)=DSPNOR(I,J)-SLPNOR(K,J)*DVGRAD(K,I)
                  DSPDIR(I,J)=DSPDIR(I,J)+SLPDIR(K,J)*DVGRAD(I,K)
               END DO

            END DO
         END DO

C-----  Update the normal to a slip plane and a slip direction (only 
C     needed for finite rotation)
C
         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               STATEV(IDNOR)=STATEV(IDNOR)+DSPNOR(I,J)-DSPNRO(I,J)

               IDDIR=IDDIR+1
               STATEV(IDDIR)=STATEV(IDDIR)+DSPDIR(I,J)-DSPDRO(I,J)
            END DO
         END DO

      END IF

C-----  Derivative of shear strain increment in a slip system w.r.t. 
C     strain increment: DDGDDE
C
      TERM1=THETA*DTIME
      DO I=1,NTENS
         DO J=1,NSLPTL
            TAUSLP=STATEV(2*NSLPTL+J)
            GSLIP=STATEV(J)
            X=TAUSLP/GSLIP
            TERM2=TERM1*DFDXSP(J)/GSLIP
            IF (I.LE.NDI) THEN
               DDGDDE(J,I)=TERM2*DDEMSD(I,J)
            ELSE
               DDGDDE(J,I)=TERM2*DDEMSD(I-NDI+3,J)
            END IF
         END DO

         CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DDGDDE(1,I))

      END DO

C-----  Derivative of stress increment w.r.t. strain increment, i.e. 
C     Jacobian matrix
C
C-----  Jacobian matrix: elastic part
      DO J=1,NTENS
         DO I=1,NTENS
            DDSDDE(I,J)=0.
         END DO
      END DO

      DO J=1,NDI
         DO I=1,NDI
            DDSDDE(I,J)=D(I,J)
            IF (NLGEOM.NE.0) DDSDDE(I,J)=DDSDDE(I,J)-STRESS(I)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR
            DO I=1,NSHR
               DDSDDE(I+NDI,J+NDI)=D(I+3,J+3)
            END DO

            DO I=1,NDI
               DDSDDE(I,J+NDI)=D(I,J+3)
               DDSDDE(J+NDI,I)=D(J+3,I)
               IF (NLGEOM.NE.0)
     2            DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-STRESS(J+NDI)
            END DO
         END DO
      END IF

C-----  Jacobian matrix: plastic part (slip)
      DO J=1,NDI
         DO I=1,NDI
            DO K=1,NSLPTL
               DDSDDE(I,J)=DDSDDE(I,J)-DDEMSD(I,K)*DDGDDE(K,J)
            END DO
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR

            DO I=1,NSHR
               DO K=1,NSLPTL
                  DDSDDE(I+NDI,J+NDI)=DDSDDE(I+NDI,J+NDI)-
     2                                DDEMSD(I+3,K)*DDGDDE(K,J+NDI)
               END DO
            END DO

            DO I=1,NDI
               DO K=1,NSLPTL
                  DDSDDE(I,J+NDI)=DDSDDE(I,J+NDI)-
     2                            DDEMSD(I,K)*DDGDDE(K,J+NDI)
                  DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-
     2                            DDEMSD(J+3,K)*DDGDDE(K,I)
               END DO
            END DO

         END DO
      END IF

      IF (ITRATN.NE.0) THEN
         DO J=1,NTENS
            DO I=1,NTENS
               DDSDDE(I,J)=DDSDDE(I,J)/(1.+DEV)
            END DO
         END DO
      END IF

C-----  Iteration ?
      IF (ITRATN.NE.0) THEN

C-----  Save solutions (without iteration):
C            Shear strain-rate in a slip system FSLIP1
C            Current strength in a slip system GSLP1
C            Shear strain in a slip system GAMMA1
C            Resolved shear stress in a slip system TAUSP1
C            Normal to a slip plane SPNOR1
C            Slip direction SPDIR1
C            Stress STRES1
C            Jacobian matrix DDSDE1
C
         IF (NITRTN.EQ.0) THEN

            IDNOR=3*NSLPTL
            IDDIR=6*NSLPTL
            DO J=1,NSLPTL
               FSLIP1(J)=FSLIP(J)
               GSLP1(J)=STATEV(J)
               GAMMA1(J)=STATEV(NSLPTL+J)
               TAUSP1(J)=STATEV(2*NSLPTL+J)
               DO I=1,3
                  IDNOR=IDNOR+1
                  SPNOR1(I,J)=STATEV(IDNOR)

                  IDDIR=IDDIR+1
                  SPDIR1(I,J)=STATEV(IDDIR)
               END DO
            END DO

            DO J=1,NTENS
               STRES1(J)=STRESS(J)
               DO I=1,NTENS
                  DDSDE1(I,J)=DDSDDE(I,J)
               END DO
            END DO

         END IF

C-----  Increments of stress DSOLD, and solution dependent state 
C     variables DGAMOD, DTAUOD, DGSPOD, DSPNRO, DSPDRO (for the next 
C     iteration)
C
         DO I=1,NTENS
            DSOLD(I)=DSTRES(I)
         END DO

         DO J=1,NSLPTL
            DGAMOD(J)=DGAMMA(J)
            DTAUOD(J)=DTAUSP(J)
            DGSPOD(J)=DGSLIP(J)
            DO I=1,3
               DSPNRO(I,J)=DSPNOR(I,J)
               DSPDRO(I,J)=DSPDIR(I,J)
            END DO
         END DO

C-----  Check if the iteration solution converges
         IDBACK=0
         ID=0
         DO I=1,NSET
            DO J=1,NSLIP(I)
               ID=ID+1
               X=STATEV(2*NSLPTL+ID)/STATEV(ID)
               RESIDU=THETA*DTIME*F(X,PROPS(65+8*I))+DTIME*(1.0-THETA)*
     2                FSLIP1(ID)-DGAMMA(ID)
               IF (ABS(RESIDU).GT.GAMERR) IDBACK=1
            END DO
         END DO

         IF (IDBACK.NE.0.AND.NITRTN.LT.ITRMAX) THEN
C-----  Iteration: arrays for iteration
CFIXA
            CALL ITERATION (STATEV(NSLPTL+1), STATEV(2*NSLPTL+1), 
     2                      STATEV(1), STATEV(9*NSLPTL+1), 
     3                      STATEV(10*NSLPTL+1), NSLPTL, 
     4                      NSET, NSLIP, ND, PROPS(97), DGAMOD,
     5                      DHDGDG)
CFIXB

            GO TO 1000

         ELSE IF (NITRTN.GE.ITRMAX) THEN
C-----  Solution not converge within maximum number of iteration (the 
C     solution without iteration will be used)
C
            DO J=1,NTENS
               STRESS(J)=STRES1(J)
               DO I=1,NTENS
                  DDSDDE(I,J)=DDSDE1(I,J)
               END DO
            END DO

            IDNOR=3*NSLPTL
            IDDIR=6*NSLPTL
            DO J=1,NSLPTL
               STATEV(J)=GSLP1(J)
               STATEV(NSLPTL+J)=GAMMA1(J)
               STATEV(2*NSLPTL+J)=TAUSP1(J)

               DO I=1,3
                  IDNOR=IDNOR+1
                  STATEV(IDNOR)=SPNOR1(I,J)

                  IDDIR=IDDIR+1
                  STATEV(IDDIR)=SPDIR1(I,J)
               END DO
            END DO

         END IF

      END IF

C-----  Total cumulative shear strains on all slip systems (sum of the 
C       absolute values of shear strains in all slip systems)
CFIX--  Total cumulative shear strains on each slip system (sum of the 
CFIX    absolute values of shear strains in each individual slip system)
C
      DO I=1,NSLPTL
CFIXA
         STATEV(10*NSLPTL+1)=STATEV(10*NSLPTL+1)+ABS(DGAMMA(I))
         STATEV(9*NSLPTL+I)=STATEV(9*NSLPTL+I)+ABS(DGAMMA(I))
CFIXB
      END DO

      RETURN
      END


C----------------------------------------------------------------------


      SUBROUTINE ROTATION (PROP, ROTATE)

C-----  This subroutine calculates the rotation matrix, i.e. the 
C     direction cosines of cubic crystal [100], [010] and [001] 
C     directions in global system

C-----  The rotation matrix is stored in the array ROTATE.

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PROP(16), ROTATE(3,3), TERM1(3,3), TERM2(3,3), INDX(3) 

C-----  Subroutines:
C
C       CROSS  -- cross product of two vectors
C
C       LUDCMP -- LU decomposition
C
C       LUBKSB -- linear equation solver based on LU decomposition 
C                 method (must call LUDCMP first)


C-----  PROP -- constants characterizing the crystal orientation 
C               (INPUT)
C
C            PROP(1) - PROP(3) -- direction of the first vector in 
C                                 local cubic crystal system
C            PROP(4) - PROP(6) -- direction of the first vector in 
C                                 global system
C
C            PROP(9) - PROP(11)-- direction of the second vector in 
C                                 local cubic crystal system
C            PROP(12)- PROP(14)-- direction of the second vector in 
C                                 global system
C
C-----  ROTATE -- rotation matrix (OUTPUT):
C
C            ROTATE(i,1) -- direction cosines of direction [1 0 0] in 
C                           local cubic crystal system
C            ROTATE(i,2) -- direction cosines of direction [0 1 0] in 
C                           local cubic crystal system
C            ROTATE(i,3) -- direction cosines of direction [0 0 1] in 
C                           local cubic crystal system

C-----  local matrix: TERM1
      CALL CROSS (PROP(1), PROP(9), TERM1, ANGLE1)

C-----  LU decomposition of TERM1
      CALL LUDCMP (TERM1, 3, 3, INDX, DCMP)

C-----  inverse matrix of TERM1: TERM2
      DO J=1,3
         DO I=1,3
            IF (I.EQ.J) THEN
               TERM2(I,J)=1.
            ELSE
               TERM2(I,J)=0.
            END IF
         END DO
      END DO

      DO J=1,3
         CALL LUBKSB (TERM1, 3, 3, INDX, TERM2(1,J))
      END DO

C-----  global matrix: TERM1
      CALL CROSS (PROP(4), PROP(12), TERM1, ANGLE2)

C-----  Check: the angle between first and second vector in local and 
C     global systems must be the same.  The relative difference must be
C     less than 0.1%.
C
      IF (ABS(ANGLE1/ANGLE2-1.).GT.0.001) THEN 
         WRITE (6,*) 
     2      '***ERROR - angles between two vectors are not the same'
         STOP
      END IF

C-----  rotation matrix: ROTATE
      DO J=1,3
         DO I=1,3
            ROTATE(I,J)=0.
            DO K=1,3
               ROTATE(I,J)=ROTATE(I,J)+TERM1(I,K)*TERM2(K,J)
            END DO
         END DO
      END DO

      RETURN
      END


C-----------------------------------


           SUBROUTINE CROSS (A, B, C, ANGLE)

C-----  (1) normalize vectors A and B to unit vectors
C       (2) store A, B and A*B (cross product) in C

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION A(3), B(3), C(3,3)

           SUM1=SQRT(A(1)**2+A(2)**2+A(3)**2)
           SUM2=SQRT(B(1)**2+B(2)**2+B(3)**2)

           IF (SUM1.EQ.0.) THEN
              WRITE (6,*) '***ERROR - first vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,1)=A(I)/SUM1
              END DO
           END IF

           IF (SUM2.EQ.0.) THEN
              WRITE (6,*) '***ERROR - second vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,2)=B(I)/SUM2
              END DO
           END IF

           ANGLE=0.
           DO I=1,3
              ANGLE=ANGLE+C(I,1)*C(I,2)
           END DO
           ANGLE=ACOS(ANGLE)

           C(1,3)=C(2,1)*C(3,2)-C(3,1)*C(2,2)
           C(2,3)=C(3,1)*C(1,2)-C(1,1)*C(3,2)
           C(3,3)=C(1,1)*C(2,2)-C(2,1)*C(1,2)
           SUM3=SQRT(C(1,3)**2+C(2,3)**2+C(3,3)**2)
           IF (SUM3.LT.1.E-8) THEN
              WRITE (6,*) 
     2           '***ERROR - first and second vectors are parallel'
               STOP
            END IF

           RETURN
           END


C----------------------------------------------------------------------


      SUBROUTINE SLIPSYS (ISPDIR, ISPNOR, NSLIP, SLPDIR, SLPNOR, 
     2                    ROTATE)

C-----  This subroutine generates all slip systems in the same set for 
C     a CUBIC crystal.  For other crystals (e.g., HCP, Tetragonal, 
C     Orthotropic, ...), it has to be modified to include the effect of
C     crystal aspect ratio.

C-----  Denote s as a slip direction and m as normal to a slip plane.  
C     In a cubic crystal, (s,-m), (-s,m) and (-s,-m) are NOT considered
C     independent of (s,m).

C-----  Subroutines:  LINE1 and LINE

C-----  Variables:
C
C     ISPDIR -- a typical slip direction in this set of slip systems 
C               (integer)  (INPUT)
C     ISPNOR -- a typical normal to slip plane in this set of slip 
C               systems (integer)  (INPUT)
C     NSLIP  -- number of independent slip systems in this set 
C               (OUTPUT)
C     SLPDIR -- unit vectors of all slip directions  (OUTPUT)
C     SLPNOR -- unit normals to all slip planes  (OUTPUT)
C     ROTATE -- rotation matrix (INPUT)
C          ROTATE(i,1) -- direction cosines of [100] in global system
C          ROTATE(i,2) -- direction cosines of [010] in global system
C          ROTATE(i,3) -- direction cosines of [001] in global system
C
C     NSPDIR -- number of all possible slip directions in this set
C     NSPNOR -- number of all possible slip planes in this set
C     IWKDIR -- all possible slip directions (integer)
C     IWKNOR -- all possible slip planes (integer)


C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ISPDIR(3), ISPNOR(3), SLPDIR(3,50), SLPNOR(3,50), 
     *          ROTATE(3,3), IWKDIR(3,24), IWKNOR(3,24), TERM(3)

      NSLIP=0
      NSPDIR=0
      NSPNOR=0

C-----  Generating all possible slip directions in this set
C
C       Denote the slip direction by [lmn].  I1 is the minimum of the 
C     absolute value of l, m and n, I3 is the maximum and I2 is the 
C     mode, e.g. (1 -3 2), I1=1, I2=2 and I3=3.  I1<=I2<=I3.

      I1=MIN(IABS(ISPDIR(1)),IABS(ISPDIR(2)),IABS(ISPDIR(3)))
      I3=MAX(IABS(ISPDIR(1)),IABS(ISPDIR(2)),IABS(ISPDIR(3)))
      I2=IABS(ISPDIR(1))+IABS(ISPDIR(2))+IABS(ISPDIR(3))-I1-I3

      RMODIR=SQRT(FLOAT(I1*I1+I2*I2+I3*I3))

C     I1=I2=I3=0
      IF (I3.EQ.0) THEN 
         WRITE (6,*) '***ERROR - slip direction is [000]'
         STOP

C     I1=I2=0, I3>0   ---   [001] type
      ELSE IF (I2.EQ.0) THEN
         NSPDIR=3
         DO J=1,3
            DO I=1,3
               IWKDIR(I,J)=0
               IF (I.EQ.J) IWKDIR(I,J)=I3
            END DO
         END DO

C     I1=0, I3>=I2>0
      ELSE IF (I1.EQ.0) THEN

C        I1=0, I3=I2>0   ---   [011] type
         IF (I2.EQ.I3) THEN
            NSPDIR=6
            DO J=1,6
               DO I=1,3
                  IWKDIR(I,J)=I2
                  IF (I.EQ.J.OR.J-I.EQ.3) IWKDIR(I,J)=0
                  IWKDIR(1,6)=-I2
                  IWKDIR(2,4)=-I2
                  IWKDIR(3,5)=-I2
               END DO
            END DO

C        I1=0, I3>I2>0   ---   [012] type
         ELSE
            NSPDIR=12
            CALL LINE1 (I2, I3, IWKDIR(1,1), 1)
            CALL LINE1 (I3, I2, IWKDIR(1,3), 1)
            CALL LINE1 (I2, I3, IWKDIR(1,5), 2)
            CALL LINE1 (I3, I2, IWKDIR(1,7), 2)
            CALL LINE1 (I2, I3, IWKDIR(1,9), 3)
            CALL LINE1 (I3, I2, IWKDIR(1,11), 3)

         END IF

C     I1=I2=I3>0   ---   [111] type
      ELSE IF (I1.EQ.I3) THEN
         NSPDIR=4
         CALL LINE (I1, I1, I1, IWKDIR)

C     I3>I2=I1>0   ---   [112] type
      ELSE IF (I1.EQ.I2) THEN
         NSPDIR=12
         CALL LINE (I1, I1, I3, IWKDIR(1,1))
         CALL LINE (I1, I3, I1, IWKDIR(1,5))
         CALL LINE (I3, I1, I1, IWKDIR(1,9))

C     I3=I2>I1>0   ---   [122] type
      ELSE IF (I2.EQ.I3) THEN
         NSPDIR=12
         CALL LINE (I1, I2, I2, IWKDIR(1,1))
         CALL LINE (I2, I1, I2, IWKDIR(1,5))
         CALL LINE (I2, I2, I1, IWKDIR(1,9))

C     I3>I2>I1>0   ---   [123] type
      ELSE
         NSPDIR=24
         CALL LINE (I1, I2, I3, IWKDIR(1,1))
         CALL LINE (I3, I1, I2, IWKDIR(1,5))
         CALL LINE (I2, I3, I1, IWKDIR(1,9))
         CALL LINE (I1, I3, I2, IWKDIR(1,13))
         CALL LINE (I2, I1, I3, IWKDIR(1,17))
         CALL LINE (I3, I2, I1, IWKDIR(1,21))

      END IF

C-----  Generating all possible slip planes in this set
C
C       Denote the normal to slip plane by (pqr).  J1 is the minimum of
C     the absolute value of p, q and r, J3 is the maximum and J2 is the
C     mode, e.g. (1 -2 1), J1=1, J2=1 and J3=2.  J1<=J2<=J3.

      J1=MIN(IABS(ISPNOR(1)),IABS(ISPNOR(2)),IABS(ISPNOR(3)))
      J3=MAX(IABS(ISPNOR(1)),IABS(ISPNOR(2)),IABS(ISPNOR(3)))
      J2=IABS(ISPNOR(1))+IABS(ISPNOR(2))+IABS(ISPNOR(3))-J1-J3

      RMONOR=SQRT(FLOAT(J1*J1+J2*J2+J3*J3))

      IF (J3.EQ.0) THEN 
         WRITE (6,*) '***ERROR - slip plane is [000]'
         STOP

C     (001) type
      ELSE IF (J2.EQ.0) THEN
         NSPNOR=3
         DO J=1,3
            DO I=1,3
               IWKNOR(I,J)=0
               IF (I.EQ.J) IWKNOR(I,J)=J3
            END DO
         END DO

      ELSE IF (J1.EQ.0) THEN

C     (011) type
         IF (J2.EQ.J3) THEN
            NSPNOR=6
            DO J=1,6
               DO I=1,3
                  IWKNOR(I,J)=J2
                  IF (I.EQ.J.OR.J-I.EQ.3) IWKNOR(I,J)=0
                  IWKNOR(1,6)=-J2
                  IWKNOR(2,4)=-J2
                  IWKNOR(3,5)=-J2
               END DO
            END DO

C     (012) type
         ELSE
            NSPNOR=12
            CALL LINE1 (J2, J3, IWKNOR(1,1), 1)
            CALL LINE1 (J3, J2, IWKNOR(1,3), 1)
            CALL LINE1 (J2, J3, IWKNOR(1,5), 2)
            CALL LINE1 (J3, J2, IWKNOR(1,7), 2)
            CALL LINE1 (J2, J3, IWKNOR(1,9), 3)
            CALL LINE1 (J3, J2, IWKNOR(1,11), 3)

         END IF

C     (111) type
      ELSE IF (J1.EQ.J3) THEN
         NSPNOR=4
         CALL LINE (J1, J1, J1, IWKNOR)

C     (112) type
      ELSE IF (J1.EQ.J2) THEN
         NSPNOR=12
         CALL LINE (J1, J1, J3, IWKNOR(1,1))
         CALL LINE (J1, J3, J1, IWKNOR(1,5))
         CALL LINE (J3, J1, J1, IWKNOR(1,9))

C     (122) type
      ELSE IF (J2.EQ.J3) THEN
         NSPNOR=12
         CALL LINE (J1, J2, J2, IWKNOR(1,1))
         CALL LINE (J2, J1, J2, IWKNOR(1,5))
         CALL LINE (J2, J2, J1, IWKNOR(1,9))

C     (123) type
      ELSE
         NSPNOR=24
         CALL LINE (J1, J2, J3, IWKNOR(1,1))
         CALL LINE (J3, J1, J2, IWKNOR(1,5))
         CALL LINE (J2, J3, J1, IWKNOR(1,9))
         CALL LINE (J1, J3, J2, IWKNOR(1,13))
         CALL LINE (J2, J1, J3, IWKNOR(1,17))
         CALL LINE (J3, J2, J1, IWKNOR(1,21))

      END IF

C-----  Generating all slip systems in this set
C
C-----  Unit vectors in slip directions: SLPDIR, and unit normals to 
C     slip planes: SLPNOR in local cubic crystal system
C
      WRITE (6,*) '          '
      WRITE (6,*) ' #          Slip plane          Slip direction'

      DO J=1,NSPNOR
         DO I=1,NSPDIR

            IDOT=0
            DO K=1,3
               IDOT=IDOT+IWKDIR(K,I)*IWKNOR(K,J)
            END DO

            IF (IDOT.EQ.0) THEN
               NSLIP=NSLIP+1
               DO K=1,3
                  SLPDIR(K,NSLIP)=IWKDIR(K,I)/RMODIR
                  SLPNOR(K,NSLIP)=IWKNOR(K,J)/RMONOR
               END DO

               WRITE (6,10) NSLIP, 
     2                      (IWKNOR(K,J),K=1,3), (IWKDIR(K,I),K=1,3)

            END IF

         END DO
      END DO
10    FORMAT(1X,I2,9X,'(',3(1X,I2),1X,')',10X,'[',3(1X,I2),1X,']')

      WRITE (6,*) 'Number of slip systems in this set = ',NSLIP
      WRITE (6,*) '          '

      IF (NSLIP.EQ.0) THEN
         WRITE (6,*) 
     *      'There is no slip direction normal to the slip planes!'
         STOP

      ELSE

C-----  Unit vectors in slip directions: SLPDIR, and unit normals to 
C     slip planes: SLPNOR in global system
C
         DO J=1,NSLIP
            DO I=1,3
               TERM(I)=0.
               DO K=1,3
                  TERM(I)=TERM(I)+ROTATE(I,K)*SLPDIR(K,J)
               END DO
            END DO
            DO I=1,3
               SLPDIR(I,J)=TERM(I)
            END DO

            DO I=1,3
               TERM(I)=0.
               DO K=1,3
                  TERM(I)=TERM(I)+ROTATE(I,K)*SLPNOR(K,J)
               END DO
            END DO
            DO I=1,3
               SLPNOR(I,J)=TERM(I)
            END DO
         END DO

      END IF

      RETURN
      END


C----------------------------------


           SUBROUTINE LINE (I1, I2, I3, IARRAY)

C-----  Generating all possible slip directions <lmn> (or slip planes 
C     {lmn}) for a cubic crystal, where l,m,n are not zeros.

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION IARRAY(3,4)

           DO J=1,4
              IARRAY(1,J)=I1
              IARRAY(2,J)=I2
              IARRAY(3,J)=I3
           END DO

           DO I=1,3
              DO J=1,4
                 IF (J.EQ.I+1) IARRAY(I,J)=-IARRAY(I,J)
              END DO
           END DO

           RETURN
           END


C-----------------------------------


           SUBROUTINE LINE1 (J1, J2, IARRAY, ID)

C-----  Generating all possible slip directions <0mn> (or slip planes 
C     {0mn}) for a cubic crystal, where m,n are not zeros and m does 
C     not equal n.

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION IARRAY(3,2)

           IARRAY(ID,1)=0
           IARRAY(ID,2)=0

           ID1=ID+1
           IF (ID1.GT.3) ID1=ID1-3
           IARRAY(ID1,1)=J1
           IARRAY(ID1,2)=J1

           ID2=ID+2
           IF (ID2.GT.3) ID2=ID2-3
           IARRAY(ID2,1)=J2
           IARRAY(ID2,2)=-J2
  
           RETURN
           END


C----------------------------------------------------------------------


      SUBROUTINE GSLPINIT (GSLIP0, NSLIP, NSLPTL, NSET, PROP)

C-----  This subroutine calculates the initial value of current 
C     strength for each slip system in a rate-dependent single crystal.
C     Two sets of initial values, proposed by Asaro, Pierce et al, and 
C     by Bassani, respectively, are used here.  Both sets assume that 
C     the initial values for all slip systems are the same (initially 
C     isotropic).

C-----  These initial values are assumed the same for all slip systems 
C     in each set, though they could be different from set to set, e.g.
C     <110>{111} and <110>{100}.

C-----  Users who want to use their own initial values may change the 
C     function subprogram GSLP0.  The parameters characterizing these 
C     initial values are passed into GSLP0 through array PROP.

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL GSLP0
      DIMENSION GSLIP0(NSLPTL), NSLIP(NSET), PROP(16,NSET)

C-----  Function subprograms:
C
C       GSLP0 -- User-supplied function subprogram given the initial 
C                value of current strength at initial state

C-----  Variables:
C
C     GSLIP0 -- initial value of current strength (OUTPUT)
C
C     NSLIP  -- number of slip systems in each set (INPUT)
C     NSLPTL -- total number of slip systems in all the sets (INPUT)
C     NSET   -- number of sets of slip systems (INPUT)
C
C     PROP   -- material constants characterizing the initial value of 
C               current strength (INPUT)
C
C               For Asaro, Pierce et al's law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- saturation stress TAUs in the ith set of  
C                            slip systems
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C
C               For Bassani's law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- stage I stress TAUI in the ith set of  
C                            slip systems (or the breakthrough stress 
C                            where large plastic flow initiates)
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C

      ID=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ID=ID+1
            GSLIP0(ID)=GSLP0(NSLPTL,NSET,NSLIP,PROP(1,I),ID,ISET)
         END DO
      END DO

      RETURN
      END


C----------------------------------


C-----  Use single precision on cray
C
           REAL*8 FUNCTION GSLP0(NSLPTL,NSET,NSLIP,PROP,ISLIP,ISET)

C-----     User-supplied function subprogram given the initial value of
C        current strength at initial state

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION NSLIP(NSET), PROP(16)

           GSLP0=PROP(3)

           RETURN
           END


C----------------------------------------------------------------------


      SUBROUTINE STRAINRATE (GAMMA, TAUSLP, GSLIP, NSLIP, FSLIP, 
     2                       DFDXSP, PROP)

C-----  This subroutine calculates the shear strain-rate in each slip 
C     system for a rate-dependent single crystal.  The POWER LAW 
C     relation between shear strain-rate and resolved shear stress 
C     proposed by Hutchinson, Pan and Rice, is used here.

C-----  The power law exponents are assumed the same for all slip 
C     systems in each set, though they could be different from set to 
C     set, e.g. <110>{111} and <110>{100}.  The strain-rate coefficient
C     in front of the power law form are also assumed the same for all 
C     slip systems in each set. 

C-----  Users who want to use their own constitutive relation may 
C     change the function subprograms F and its derivative DFDX, 
C     where F is the strain hardening law, dGAMMA/dt = F(X), 
C     X=TAUSLP/GSLIP.  The parameters characterizing F are passed into 
C     F and DFDX through array PROP.

C-----  Function subprograms:
C
C       F    -- User-supplied function subprogram which gives shear 
C               strain-rate for each slip system based on current 
C               values of resolved shear stress and current strength
C
C       DFDX -- User-supplied function subprogram dF/dX, where x is the
C               ratio of resolved shear stress over current strength

C-----  Variables:
C
C     GAMMA  -- shear strain in each slip system at the start of time 
C               step  (INPUT)
C     TAUSLP -- resolved shear stress in each slip system (INPUT)
C     GSLIP  -- current strength (INPUT)
C     NSLIP  -- number of slip systems in this set (INPUT)
C
C     FSLIP  -- current value of F for each slip system (OUTPUT)
C     DFDXSP -- current value of DFDX for each slip system (OUTPUT)
C
C     PROP   -- material constants characterizing the strain hardening 
C               law (INPUT)
C
C               For the current power law strain hardening law 
C               PROP(1) -- power law hardening exponent
C               PROP(1) = infinity corresponds to a rate-independent 
C               material
C               PROP(2) -- coefficient in front of power law hardening


C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F, DFDX
      DIMENSION GAMMA(NSLIP), TAUSLP(NSLIP), GSLIP(NSLIP), 
     2          FSLIP(NSLIP), DFDXSP(NSLIP), PROP(8)

      DO I=1,NSLIP
         X=TAUSLP(I)/GSLIP(I)
         FSLIP(I)=F(X,PROP)
         DFDXSP(I)=DFDX(X,PROP)
      END DO

      RETURN
      END


C-----------------------------------


C-----  Use single precision on cray
C
           REAL*8 FUNCTION F(X,PROP)

C-----     User-supplied function subprogram which gives shear 
C        strain-rate for each slip system based on current values of 
C        resolved shear stress and current strength
C
C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION PROP(8)

           F=PROP(2)*(ABS(X))**PROP(1)*DSIGN(1.D0,X)

           RETURN
           END


C-----------------------------------


C-----  Use single precision on cray
C
           REAL*8 FUNCTION DFDX(X,PROP)

C-----     User-supplied function subprogram dF/dX, where x is the 
C        ratio of resolved shear stress over current strength

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION PROP(8)

           DFDX=PROP(1)*PROP(2)*(ABS(X))**(PROP(1)-1.)

           RETURN
           END


C----------------------------------------------------------------------

CFIXA
      SUBROUTINE LATENTHARDEN (GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL, 
     2                         NSLIP, NSLPTL, NSET, H, PROP, ND)
CFIXB

C-----  This subroutine calculates the current self- and latent-
C     hardening moduli for all slip systems in a rate-dependent single 
C     crystal.  Two kinds of hardening law are used here.  The first 
C     law, proposed by Asaro, and Pierce et al, assumes a HYPER SECANT 
C     relation between self- and latent-hardening moduli and overall 
C     shear strain.  The Bauschinger effect has been neglected.  The 
C     second is Bassani's hardening law, which gives an explicit 
C     expression of slip interactions between slip systems.  The 
C     classical three stage hardening for FCC single crystal could be 
C     simulated.

C-----  The hardening coefficients are assumed the same for all slip 
C     systems in each set, though they could be different from set to 
C     set, e.g. <110>{111} and <110>{100}.

C-----  Users who want to use their own self- and latent-hardening law 
C     may change the function subprograms HSELF (self hardening) and 
C     HLATNT (latent hardening).  The parameters characterizing these 
C     hardening laws are passed into HSELF and HLATNT through array 
C     PROP.


C-----  Function subprograms:
C
C       HSELF  -- User-supplied self-hardening function in a slip 
C                 system
C
C       HLATNT -- User-supplied latent-hardening function

C-----  Variables:
C
C     GAMMA  -- shear strain in all slip systems at the start of time 
C               step  (INPUT)
C     TAUSLP -- resolved shear stress in all slip systems (INPUT)
C     GSLIP  -- current strength (INPUT)
CFIX  GMSLTL -- total cumulative shear strains on each individual slip system 
CFIX            (INPUT)
C     GAMTOL -- total cumulative shear strains over all slip systems 
C               (INPUT)
C     NSLIP  -- number of slip systems in each set (INPUT)
C     NSLPTL -- total number of slip systems in all the sets (INPUT)
C     NSET   -- number of sets of slip systems (INPUT)
C
C     H      -- current value of self- and latent-hardening moduli 
C               (OUTPUT)
C               H(i,i) -- self-hardening modulus of the ith slip system
C                         (no sum over i)
C               H(i,j) -- latent-hardening molulus of the ith slip 
C                         system due to a slip in the jth slip system 
C                         (i not equal j)
C
C     PROP   -- material constants characterizing the self- and latent-
C               hardening law (INPUT)
C
C               For the HYPER SECANT hardening law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- saturation stress TAUs in the ith set of  
C                            slip systems
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets 
C                            of slip systems to self-hardening in the 
C                            ith set of slip systems Q1
C
C               For Bassani's hardening law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- stage I stress TAUI in the ith set of  
C                            slip systems (or the breakthrough stress 
C                            where large plastic flow initiates)
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C               PROP(4,i) -- hardening modulus during easy glide Hs in 
C                            the ith set of slip systems
C               PROP(5,i) -- amount of slip Gamma0 after which a given 
C                            interaction between slip systems in the 
C                            ith set reaches peak strength
C               PROP(6,i) -- amount of slip Gamma0 after which a given 
C                            interaction between slip systems in the 
C                            ith set and jth set (i not equal j) 
C                            reaches peak strength
C               PROP(7,i) -- representing the magnitude of the strength
C                            of interaction in the ith set of slip 
C                            system
C               PROP(8,i) -- representing the magnitude of the strength
C                            of interaction between the ith set and jth
C                            set of system
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets 
C                            of slip systems to self-hardening in the 
C                            ith set of slip systems Q1
C
C     ND     -- leading dimension of arrays defined in subroutine UMAT 
C               (INPUT) 


C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL HSELF, HLATNT
CFIXA
      DIMENSION GAMMA(NSLPTL), TAUSLP(NSLPTL), GMSLTL(NSLPTL),
     2          GSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET), 
     3          H(ND,NSLPTL)
CFIXB

      CHECK=0.
      DO I=1,NSET
         DO J=4,8
            CHECK=CHECK+ABS(PROP(J,I))
         END DO
      END DO

C-----  CHECK=0   --  HYPER SECANT hardening law
C       otherwise --  Bassani's hardening law

      ISELF=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ISELF=ISELF+1

            DO LATENT=1,NSLPTL
               IF (LATENT.EQ.ISELF) THEN
CFIXA
                  H(LATENT,ISELF)=HSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,
     2                                  NSET,NSLIP,PROP(1,I),CHECK,
     3                                  ISELF,ISET)
CFIXB
               ELSE
CFIXA
                  H(LATENT,ISELF)=HLATNT(GAMMA,GMSLTL,GAMTOL,NSLPTL,
     2                                   NSET,NSLIP,PROP(1,I),CHECK,
     3                                   ISELF,ISET,LATENT)
CFIXB

               END IF
            END DO

         END DO
      END DO

      RETURN
      END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION HSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP,CHECK,ISELF,ISET)
CFIXB

C-----     User-supplied self-hardening function in a slip system

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), NSLIP(NSET), PROP(16),
     2               GMSLTL(NSLPTL)
CFIXB

           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              HSELF=PROP(1)*TERM2**2

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)

              ID=0
              G=1.
              DO I=1,NSET
                 IF (I.EQ.ISET) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

                 DO J=1,NSLIP(I)
                    ID=ID+1
                    IF (ID.NE.ISELF) THEN
CFIXA
		       G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
		    END IF

                 END DO
              END DO

              HSELF=F*G

           END IF

           RETURN
           END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION HLATNT(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,LATENT)
CFIXB

C-----     User-supplied latent-hardening function

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), NSLIP(NSET), PROP(16),
     2               GMSLTL(NSLPTL)
CFIXB

           ILOWER=0
           IUPPER=NSLIP(1)
           IF (ISET.GT.1) THEN
              DO K=2,ISET
                 ILOWER=ILOWER+NSLIP(K-1)
                 IUPPER=IUPPER+NSLIP(K)
              END DO
           END IF

           IF (LATENT.GT.ILOWER.AND.LATENT.LE.IUPPER) THEN
              Q=PROP(9)
           ELSE
              Q=PROP(10)
           END IF

           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              HLATNT=PROP(1)*TERM2**2*Q

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)

              ID=0
              G=1.
              DO I=1,NSET
                 IF (I.EQ.ISET) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

                 DO J=1,NSLIP(I)
                    ID=ID+1
                    IF (ID.NE.ISELF) THEN
CFIXA
		       G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
		    END IF

                 END DO
              END DO

              HLATNT=F*G*Q

           END IF

           RETURN
           END


C----------------------------------------------------------------------

CFIXA
      SUBROUTINE ITERATION (GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL, 
     2                      NSLPTL, NSET, NSLIP, ND, PROP, DGAMOD, 
     3                      DHDGDG)
CFIXB

C-----  This subroutine generates arrays for the Newton-Rhapson 
C     iteration method.

C-----  Users who want to use their own self- and latent-hardening law 
C     may change the function subprograms DHSELF (self hardening) and 
C     DHLATN (latent hardening).  The parameters characterizing these 
C     hardening laws are passed into DHSELF and DHLATN through array 
C     PROP.


C-----  Function subprograms:
C
C       DHSELF -- User-supplied function of the derivative of self-
C                 hardening moduli
C
C       DHLATN -- User-supplied function of the derivative of latent-
C                 hardening moduli

C-----  Variables:
C
C     GAMMA  -- shear strain in all slip systems at the start of time 
C               step  (INPUT)
C     TAUSLP -- resolved shear stress in all slip systems (INPUT)
C     GSLIP  -- current strength (INPUT)
CFIX  GMSLTL -- total cumulative shear strains on each individual slip system 
CFIX            (INPUT)
C     GAMTOL -- total cumulative shear strains over all slip systems 
C               (INPUT)
C     NSLPTL -- total number of slip systems in all the sets (INPUT)
C     NSET   -- number of sets of slip systems (INPUT)
C     NSLIP  -- number of slip systems in each set (INPUT)
C     ND     -- leading dimension of arrays defined in subroutine UMAT 
C               (INPUT) 
C
C     PROP   -- material constants characterizing the self- and latent-
C               hardening law (INPUT)
C
C               For the HYPER SECANT hardening law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- saturation stress TAUs in the ith set of  
C                            slip systems
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets 
C                            of slip systems to self-hardening in the 
C                            ith set of slip systems Q1
C
C               For Bassani's hardening law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- stage I stress TAUI in the ith set of  
C                            slip systems (or the breakthrough stress 
C                            where large plastic flow initiates)
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C               PROP(4,i) -- hardening modulus during easy glide Hs in 
C                            the ith set of slip systems
C               PROP(5,i) -- amount of slip Gamma0 after which a given 
C                            interaction between slip systems in the 
C                            ith set reaches peak strength
C               PROP(6,i) -- amount of slip Gamma0 after which a given 
C                            interaction between slip systems in the 
C                            ith set and jth set (i not equal j) 
C                            reaches peak strength
C               PROP(7,i) -- representing the magnitude of the strength
C                            of interaction in the ith set of slip 
C                            system
C               PROP(8,i) -- representing the magnitude of the strength
C                            of interaction between the ith set and jth
C                            set of system
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets 
C                            of slip systems to self-hardening in the 
C                            ith set of slip systems Q1
C
C-----  Arrays for iteration:
C
C       DGAMOD (INPUT)
C
C       DHDGDG (OUTPUT)
C

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL DHSELF, DHLATN
CFIXA
      DIMENSION GAMMA(NSLPTL), TAUSLP(NSLPTL), GMSLTL(NSLPTL),
     2          GSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET), 
     3          DGAMOD(NSLPTL), DHDGDG(ND,NSLPTL)
CFIXB

      CHECK=0.
      DO I=1,NSET
         DO J=4,8
            CHECK=CHECK+ABS(PROP(J,I))
         END DO
      END DO

C-----  CHECK=0   --  HYPER SECANT hardening law
C       otherwise --  Bassani's hardening law

      ISELF=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ISELF=ISELF+1

            DO KDERIV=1,NSLPTL
               DHDGDG(ISELF,KDERIV)=0.

               DO LATENT=1,NSLPTL
                  IF (LATENT.EQ.ISELF) THEN
CFIXA
                     DHDG=DHSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP(1,I),CHECK,ISELF,ISET,
     3                           KDERIV)
CFIXB
                  ELSE
CFIXA
                     DHDG=DHLATN(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP(1,I),CHECK,ISELF,ISET,
     3                           LATENT,KDERIV)
CFIXB
                  END IF

                  DHDGDG(ISELF,KDERIV)=DHDGDG(ISELF,KDERIV)+
     2                                 DHDG*ABS(DGAMOD(LATENT))
               END DO

            END DO
         END DO
      END DO

      RETURN
      END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION DHSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,
     3                            KDERIV)
CFIXB

C-----  User-supplied function of the derivative of self-hardening
C     moduli

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), GMSLTL(NSLPTL), 
     2               NSLIP(NSET), PROP(16)
CFIXB

           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=PROP(1)/(PROP(2)-PROP(3))*DSIGN(1.D0,GAMMA(KDERIV))
              DHSELF=-2.*PROP(1)*TERM2**2*TANH(TERM1)*TERM3

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=(PROP(1)-PROP(4))/(PROP(2)-PROP(3))

              IF (KDERIV.EQ.ISELF) THEN
                 F=-2.*(PROP(1)-PROP(4))*TERM2**2*TANH(TERM1)*TERM3
                 ID=0
                 G=1.
                 DO I=1,NSET
                    IF (I.EQ.ISET) THEN
                       GAMMA0=PROP(5)
                       FAB=PROP(7)
                    ELSE
                       GAMMA0=PROP(6)
                       FAB=PROP(8)
                    END IF

                    DO J=1,NSLIP(I)
                       ID=ID+1
CFIXA
                       IF (ID.NE.ISELF) G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
                    END DO
                 END DO

              ELSE
                 F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)
                 ILOWER=0
                 IUPPER=NSLIP(1)
                 IF (ISET.GT.1) THEN
                    DO K=2,ISET
                       ILOWER=ILOWER+NSLIP(K-1)
                       IUPPER=IUPPER+NSLIP(K)
                    END DO
                 END IF

                 IF (KDERIV.GT.ILOWER.AND.KDERIV.LE.IUPPER) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

CFIXA
                 TERM4=GMSLTL(KDERIV)/GAMMA0
CFIXB
                 TERM5=2.*EXP(-TERM4)/(1.+EXP(-2.*TERM4))
                 G=FAB/GAMMA0*TERM5**2

              END IF

              DHSELF=F*G

           END IF

           RETURN
           END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION DHLATN(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,LATENT,
     3                            KDERIV)
CFIXB

C-----  User-supplied function of the derivative of latent-hardening 
C     moduli

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), GMSLTL(NSLPTL), NSLIP(NSET), 
     2               PROP(16)
CFIXB

           ILOWER=0
           IUPPER=NSLIP(1)
           IF (ISET.GT.1) THEN
              DO K=2,ISET
                 ILOWER=ILOWER+NSLIP(K-1)
                 IUPPER=IUPPER+NSLIP(K)
              END DO
           END IF

           IF (LATENT.GT.ILOWER.AND.LATENT.LE.IUPPER) THEN
              Q=PROP(9)
           ELSE
              Q=PROP(10)
           END IF

           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=PROP(1)/(PROP(2)-PROP(3))*DSIGN(1.D0,GAMMA(KDERIV))
              DHLATN=-2.*PROP(1)*TERM2**2*TANH(TERM1)*TERM3*Q

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=(PROP(1)-PROP(4))/(PROP(2)-PROP(3))

              IF (KDERIV.EQ.ISELF) THEN
                 F=-2.*(PROP(1)-PROP(4))*TERM2**2*TANH(TERM1)*TERM3
                 ID=0
                 G=1.
                 DO I=1,NSET
                    IF (I.EQ.ISET) THEN
                       GAMMA0=PROP(5)
                       FAB=PROP(7)
                    ELSE
                       GAMMA0=PROP(6)
                       FAB=PROP(8)
                    END IF

                    DO J=1,NSLIP(I)
                       ID=ID+1
CFIXA
                       IF (ID.NE.ISELF) G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
                    END DO
                 END DO

              ELSE
                 F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)
                 ILOWER=0
                 IUPPER=NSLIP(1)
                 IF (ISET.GT.1) THEN
                    DO K=2,ISET
                       ILOWER=ILOWER+NSLIP(K-1)
                       IUPPER=IUPPER+NSLIP(K)
                    END DO
                 END IF

                 IF (KDERIV.GT.ILOWER.AND.KDERIV.LE.IUPPER) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF
CFIXA
                 TERM4=GMSLTL(KDERIV)/GAMMA0
CFIXB
                 TERM5=2.*EXP(-TERM4)/(1.+EXP(-2.*TERM4))
                 G=FAB/GAMMA0*TERM5**2

              END IF

              DHLATN=F*G*Q

           END IF

           RETURN
           END


C----------------------------------------------------------------------


      SUBROUTINE LUDCMP (A, N, NP, INDX, D)

C-----  LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=200, TINY=1.0E-20)
      DIMENSION A(NP,NP), INDX(N), VV(NMAX)

      D=1.
      DO I=1,N
         AAMAX=0.

         DO J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
         END DO

         IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
         VV(I)=1./AAMAX
      END DO

      DO J=1,N
         DO I=1,J-1
            SUM=A(I,J)

            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
         END DO
         AAMAX=0.

         DO I=J,N
            SUM=A(I,J)

            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            END IF
         END DO

         IF (J.NE.IMAX) THEN
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            END DO

            D=-D
            VV(IMAX)=VV(J)
         END IF

         INDX(J)=IMAX
         IF (A(J,J).EQ.0.) A(J,J)=TINY
         IF (J.NE.N) THEN
            DUM=1./A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*DUM
            END DO
         END IF

      END DO

      RETURN
      END


C----------------------------------------------------------------------


      SUBROUTINE LUBKSB (A, N, NP, INDX, B)

C-----  Linear equation solver based on LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), INDX(N), B(N)

      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)

         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            END DO
         ELSE IF (SUM.NE.0.) THEN
            II=I
         END IF

         B(I)=SUM
      END DO

      DO I=N,1,-1
         SUM=B(I)

         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            END DO
         END IF

         B(I)=SUM/A(I,I)
      END DO

      RETURN
      END
