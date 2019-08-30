CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Last change Mar-2014: Program converted to subroutine called by 
C     LitMod_v4.0 (MF)
C
      program generator_table_atten_corr
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
      REAL k
c     previous data
c      Parameter(pi=3.1415926,AA1=750,alfa3=0.26,energi=424.0E3,
c     *          volexp=1.25E-05,rgas=8.314472)
c  updated from Jackshon et.al 2010
      Parameter(pi=3.1415926,AA1=68,alfa3=0.34,energi=293.0E03,
     *          volexp=1.25E-05,rgas=8.314472)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*)"Now including anelastic attenuation effects."
      write(*,*)"Enter the grain size (5,10) in mm :"
c      DSIZE=5
      read(*,*)DSIZE
      write(*,*)"Enter the oscillation period for anelastic effects(s):"
c      IOSPE=100
      read(*,*)IOSPE
      open(1,file='TABLE_LITMOD2',status='unknown')
      open(2,file='TABLE_LITMOD_anelastic_deri',status='unknown')
        do 20
      
c      read(3,1912,iostat=ios)aT,aP,cond
      read(1,1010,iostat=ios)aT,aP,ade,avp,avs,avpdt,avsdt,avpdp,avsdp,k
c1010  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6)
c1911  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6,1F8.2)
c1912  format(F9.3,1F17.4,I8.2)
c1010  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6)
1010  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6,1F12.3)
1011  format(F9.3,1F17.4,1F15.3,2F12.3,6E17.6,1F12.3)
      if(IOS.lt.0.0)then
        goto 30
      else 
c    ...exponential term...
          
c       parexp=DEXP(((-energi+(volexp*(aP*1.0E05))))/(rgas*(aT-273.15)))
       parexp=DEXP((-(energi+(volexp*aP*1.0E05)))/(rgas*(aT)))
c       p_=((-energi+(volexp*aP*1.0E05))*alfa3*2.0E0)/
c     * (rgas*(aT-273.15)**2.0E0)
c       p_=(1.0E0/aT)*(((-energi+(volexp*(aP*1.0E05)))/(rgas*(aT)))**
c#     Goes 2000
      p_=(((-energi+(volexp*(aP)))*alfa3)/
     * (rgas*(aT**2.0E0))) 
c#    derivative calulated by me 
c      p_=(((energi+(volexp*aP*1.0E05))*alfa3)/(rgas*(aT**2.E0)))/parexp
c      write(*,*)p_
c       parexp=DEXP((-(energi+(volexp*aP)))/(rgas*(aT-273.15)))
       selectcase (iospe)
        case (50)
         sqatt50=AA1*(((50.0E0*(1.0E0/(dsize*1000.0E0)))*parexp))**alfa3
         cots50= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqatt50)*0.5E0
         cotp50= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqatt50)*(2.0E0/9.0E0)
c        cotp50=0.001
c         der_f_s=(p_/DTAN((pi*alfa3)/2.0E0))*0.5E0
c         der_f_p=(p_/DTAN((pi*alfa3)/2.0E0))*(2.0E0/9.0E0)
         dVsdT=(1.0E0/parexp)*cots50*p_
         dVpdT=cotp50*p_
        case (75)
         sqatt75=AA1*(((75.0E0*(1.0E0/(dsize*1000.0E0)))*parexp))**alfa3
         cots75= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqatt75)*0.5E0
         cotp75= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqatt75)*(2.0E0/9.0E0)
c         cotp75=0.001
c         der_f=p_/DTAN((pi*alfa3)/2.0E0)
         dVsdT=cots75*p_
         dVpdT=cotp75*p_
        case (100)
         sqa100=AA1*(((100.0E0*(1.0E0/(dsize*1000.0E0)))*parexp))**alfa3
         cots100= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqa100)*0.5E0
         cotp100= ((1.0E0/DTAN((pi*alfa3)/2.0E0))*sqa100)*(2.0E0/9.0E0)
c         cotp100=0.001
c         der_f=p_/DTAN((pi*alfa3)/2.0E0)
         dVsdT=cots100*p_
         dVpdT=cotp100*p_
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
       write(2,1011)aT,aP,ade,vpat,vsat,avpdt,avsdt,avpdp,avsdp,dVpdT,
     *  dVsdT,k
c 1003 format(F9.3,1F17.4,2F4.4,2F4.4,2F4.4,2F4.4,2F4.4,2F4.4,2F4.4)
c1010  format(F9.3,1F17.4,1F15.3,2F12.3,4E17.6)
      endif
  20    continue
30    close(1)
      close(2)
      Print *, 'Attenuation calculations are finished'
              return
      end
