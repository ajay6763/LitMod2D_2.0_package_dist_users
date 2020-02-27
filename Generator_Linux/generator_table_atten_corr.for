CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Last change Mar-2014: Program converted to subroutine called by 
C     LitMod_v4.0 (MF)
C
      SUBROUTINE generator_table_atten_corr(DSIZE,IOSPE)
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
      INTEGER cond 
c     previous data
c      Parameter(pi=3.1415926,AA1=750,alfa3=0.26,energi=424.0E3,
c     *          volexp=1.25E-05,rgas=8.314472)
c  updated from Jackshon et.al 2010
      Parameter(pi=3.1415926,AA1=816,alfa3=0.36,energi=293.0E3,
     *          volexp=1.20E-05,rgas=8.314472)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*)"Now including anelastic attenuation effects."
c      write(*,*)"Enter the grain size (5,10) in mm :"
c      DSIZE=5
c      read(*,*)DSIZE
c      write(*,*)"Enter oscillation period for anelastic effects(s):"
c      IOSPE=100
c      read(*,*)IOSPE
      open(1,file='TABLE_LITMOD_no_atten',status='unknown')
      open(2,file='TABLE_LITMOD_atten_corr',status='unknown')
      open(3,file='Attenuation_parameters_info',status='unknown')
      write(3,*)"Grain size(mm) =" ,DSIZE
      write(3,*)"Time period (secods) =", IOSPE
      
        do 20
      read(1,1010,iostat=ios)aT,aP,aden,avp,avs,avpdt,avsdt,avpdp,avsdp
c1010  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6)
c1911  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6,1F8.2)
c1912  format(F9.3,1F17.4,1F15.3,F9.3)
1010  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6)
      if(IOS.lt.0.0)then
        goto 30
      else 
c    ...exponential term...
          
       parexp=DEXP((-(energi+(volexp*aP*1.E+05)))/(rgas*(aT)))
       
       selectcase (iospe)
        case (50)
         sqatt50=AA1*(((50.0E0*(1.0E0/(dsize*1000.0E0)))*parexp))**alfa3
         cots50= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqatt50)*0.5E0
         cotp50= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqatt50)*(2.0E0/9.0E0)
c        cotp50=0.001
         vsf50=avs*(1.0E0-cots50)
         vpf50=avp*(1.0E0-cotp50)
         vpat=vpf50
         vsat=vsf50
        case (75)
         sqatt75=AA1*(((75.0E0*(1.0E0/(dsize*1000.0E0)))*parexp))**alfa3
         cots75= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqatt75)*0.5E0
         cotp75= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqatt75)*(2.0E0/9.0E0)
c         cotp75=0.001
         vsf75=avs*(1.0E0-cots75)
         vpf75=avp*(1.0E0-cotp75)
         vpat=vpf75
         vsat=vsf75
        case (100)
         sqa100=AA1*(((100.0E0*(1.0E0/(dsize*1000.0E0)))*parexp))**alfa3
         cots100= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqa100)*0.5E0
         cotp100= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqa100)*(2.0E0/9.0E0)
c         cotp100=0.001
         vsf100=avs*(1.0E0-cots100)
         vpf100=avp*(1.0E0-cotp100)
         vpat=vpf100
         vsat=vsf100
        case default
         print*,'Oscillation period IOSPE must be 50, 75 or 100 s'
       endselect
c      difp=(vpat/avp)-1
c      difs=(vsat/avs)-1    
c       write(2,*)aT,aP,aden,avp,avs,avpdt,difqp,avsdt,difs
c 1010 format(F10.2,F10.2,F10.2,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6)
c 1003 format(2F8.2,2F6.3,2F8.4)   
c 1911  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6,1F8.2)
c      write(2,1003)aT,aP,vpat,vsat,difp,difs,avp,avs
       write(2,1010)aT,aP,aden,vpat,vsat,avpdt,avsdt,avpdp,avsdp
c 1003 format(F9.3,1F17.4,2F4.4,2F4.4,2F4.4,2F4.4,2F4.4,2F4.4,2F4.4)
c1010  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6)
      endif
  20    continue
30    close(1)
      close(2)
      Print *, 'Attenuation calculations are finished'
              return
      end
