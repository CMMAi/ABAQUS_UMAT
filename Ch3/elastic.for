      subroutine umat(
C *** �ѼƦC ********************************************************* C
     1  STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     2  RPL,DDSDDT,DRPLDE,DRPLDT,
     3  STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     4  NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     5  CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C **** �ŧi�Ѽ� ****************************************************** C
      include 'ABA_PARAM.INC'
C -------------------------------------------------------------------- C
      character*80 CMNAME
      dimension STRESS(NTENS),STATEV(NSTATV),
     1  DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2  STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3  PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4  JSTEP(4)
C -------------------------------------------------------------------- C
      real*8 E, nu, mu, lam, mu2
      integer*4 i, j
C *** ���Ƽҫ��D�{�� ************************************************* C
C --- ���ƨϥο��~�A������R
      if (NTENS .EQ. 1) then
        write(7, *) '���~�G���u�u�ʼҫ����䴩�@�������C'
        call XIT
      endif
C --- ���ƫl�ׯx�}�k�s
      do i= 1, NTENS
        do j=1,NTENS
          DDSDDE(i,j)=0.d0
        enddo
      enddo
C --- ���ưѼ�Ū���P�p��
      E   = PROPS(1)
      nu  = PROPS(2)
      mu  = E/(2.d0+2.d0*nu)
      lam = mu*nu/(0.5d0-nu)
C --- �����O�M�l�ׯx�}
      do j= 1,NSHR
        i = j+NDI
        DDSDDE(i,i) = mu
        STRESS(i)   = STRESS(i)+mu*DSTRAN(i)
      enddo
C --- ���V���O�M�l�ׯx�}�]�������ܩM�T�����ơ^
      mu2 = 2.d0*mu
      do i= 1,NDI
        do j= 1,NDI
          DDSDDE(i,j) = lam
          STRESS(i)   = STRESS(i)+ lam*DSTRAN(j)
        enddo
        DDSDDE(i,i) = lam + mu2
        STRESS(i)   = STRESS(i)+ mu2*DSTRAN(i)
      enddo
C ******************************************************************** C
      return
      end subroutine umat
