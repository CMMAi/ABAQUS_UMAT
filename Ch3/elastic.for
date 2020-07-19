      subroutine umat(
C *** 參數列 ********************************************************* C
     1  STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     2  RPL,DDSDDT,DRPLDE,DRPLDT,
     3  STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     4  NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     5  CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C **** 宣告參數 ****************************************************** C
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
C *** 材料模型主程式 ************************************************* C
C --- 材料使用錯誤，停止分析
      if (NTENS .EQ. 1) then
        write(7, *) '錯誤：本線彈性模型不支援一維元素。'
        call XIT
      endif
C --- 材料勁度矩陣歸零
      do i= 1, NTENS
        do j=1,NTENS
          DDSDDE(i,j)=0.d0
        enddo
      enddo
C --- 材料參數讀取與計算
      E   = PROPS(1)
      nu  = PROPS(2)
      mu  = E/(2.d0+2.d0*nu)
      lam = mu*nu/(0.5d0-nu)
C --- 剪應力和勁度矩陣
      do j= 1,NSHR
        i = j+NDI
        DDSDDE(i,i) = mu
        STRESS(i)   = STRESS(i)+mu*DSTRAN(i)
      enddo
C --- 正向應力和勁度矩陣（平面應變和三維材料）
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
