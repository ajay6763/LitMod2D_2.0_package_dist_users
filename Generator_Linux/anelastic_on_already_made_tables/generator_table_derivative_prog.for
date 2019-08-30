CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Last change Mar-2014: Program converted to subroutine called by 
C     LitMod_v4.0 (MF)
C
      program anelastic_atten_deri_T_P
c**********************************************************************
C     THIS SUBROUTINE CALCULATES P- AND S-WAVE VELOCITIES INCLUDING 
C     ANELASTIC (ATTENUATION) EFFECTS. IT NEEDS AS INPUT THE T-P 
C     CONDITIONS AND THE VELOCITIES CALCULATED WITH THE ANHARMONIC
C     APPROXIMATION. THIS INFO IS READ FROM AN EXTERNAL FILE CREATED
C     IN A PREVIOUS STEP IN LITMOD.
C        
C     REFERENCES: Afonso et al. (2008) G3.
C**********************************************************************
C     DSIZE is the mantle grain size, typically 5 mm 
C     IOSPE is the oscillation period which can take the value of 
C     50, 75 or 100 sec.
C
C**********************************************************************
      IMPLICIT double PRECISION (a-h,o-z) 
      integer :: N, i, lines_in_file
      real, allocatable :: d(:,:),Vp_mat(:,:),Vs_mat(:,:),P(:,:),T(:,:)
      real, allocatable :: dVp_dT(:,:), dVs_dT(:,:), dVp_dP(:,:)
      real, allocatable :: dVs_dP(:,:)
      allocate(d(120000,10))
      allocate(Vp_mat(300,400))
      allocate(Vs_mat(300,400))
      allocate(P(300,400))
      allocate(T(300,400))
      allocate(dVp_dP(300,400))
      allocate(dVs_dP(300,400))
      allocate(dVp_dT(300,400))
      allocate(dVs_dT(300,400))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      open(1,file='TABLE_LITMOD_atten_corr',status='unknown')
      open(2,file='TABLE_LITMOD3',status='unknown')
c reading into arrays
      do i=1,120000
        read(1,1010,iostat=ios)d(i,1),d(i,2),d(i,3),d(i,4),d(i,5),
     *  d(i,6),d(i,7),d(i,8),d(i,9),d(i,10)
      enddo
1010  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6,1F12.3)
c @ ajay
c  reshaping Vp and Vs into T,P mesh to calculate gradients at P and T
      do i=1,400
        do j=1,300
            index=j+((i-1)*300)
            Vp_mat(j,i)=d(index,4)
            Vs_mat(j,i)=d(index,5)
c Here Pressure is not converted to pascal because in the the subroutine
c where effect of these derivatives are include pressure difference is
c is calculated in bars so I guess anharmonic derivatives from generated
c are bars, so just to be consistace I calculate derivatives using bars            
            P(j,i)=d(index,2)
            T(j,i)=d(index,1)
        enddo
       enddo 
c @ ajay
      write(*,*)"calculating temperature derivative"
c    dT
      do i=1,400
        do j=1,300
            
            dVp_dT(j,i)=(Vp_mat(j+1,i)-2.E0*Vp_mat(j,i)+Vp_mat(j-1,i))
     *  /(T(j,i)-T(j+1,i))**2.E0 
c            if (dVp_dT(j,i)==0.0) then
c                dVp_dT(j,i)=dVp_dT(j-1,i)
c            endif
            dVs_dT(j,i)=(Vs_mat(j+1,i)-2.E0*Vs_mat(j,i)+Vs_mat(j-1,i))
     *  /(T(j,i)-T(j+1,i))**2.E0 
c            if (dVs_dT(j,i)==0.0) then
c                dVs_dT(j,i) = dVs_dT(j-1,i)
c           endif
        write(*,*)dVs_dT(j,i),P(j,i),T(j,i)
        enddo
       enddo    
      write(*,*)"calculating pressure derivative"
c    dP
      do i=1,300
        do j=1,400
            dVp_dP(i,j)=(Vp_mat(i,j)-Vp_mat(i,j+1))/(P(i,j)-
     *  P(i,j+1))
       if (dVp_dP(i,j)==0.0) then
               dVp_dP(i,j)=dVp_dP(i,j-1)
        endif
            dVs_dP(i,j)=(Vs_mat(i,j)-Vs_mat(i,j+1))/(P(i,j)-
     *  P(i,j+1))
        if (dVs_dP(i,j)==0.0) then
                dVs_dP(i,j)=dVs_dP(i,j-1)
        endif
        write(*,*)dVs_dP(i,j),P(i,j),T(i,j)
        enddo
       enddo
       write(*,*)"size of vp,vs,P,T"
       write(*,*)size(Vp_mat),size(Vs_mat),size(P),size(T)
       write(*,*)"size of derivatie"
       write(*,*)size(dVp_dP),size(dVp_dT)
c writing into files
      do i=1,400
       do j=1,300
            index=j+((i-1)*300)
            write(2,1010)d(index,1),d(index,2),d(index,3),d(index,4),
c     *  d(index,5),d(index,6),d(index,7),d(index,8),d(index,9),
c     *  dVp_dT(i,j),dVs_dT(i,j),dVp_dP(i,j),dVs_dP(i,j)
     *  d(index,5),
     *  dVp_dT(i,j)+d(index,6),dVs_dT(i,j)+d(index,7),dVp_dP(i,j)
     *  +d(index,8),dVs_dP(i,j)+d(index,9),d(index,10)  
c             write(2,*)T(i,j),P(i,j),
c     *  dVp_dT(i,j),dVs_dT(i,j),dVp_dP(i,j),dVs_dP(i,j)
        enddo
       enddo 
       
c1010  format(F10.2,F10.2,F10.2,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6)
c 1003 format(2F8.2,2F6.3,2F8.4)   
c 1911  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6,1F8.2)


c    	endif
      
c  20    continue

30    close(1)
      close(2)

      Print *, 'Attenuation calculations are finished'
              return
      end
