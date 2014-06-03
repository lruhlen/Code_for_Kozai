      PROGRAM STELLAR
      include 'parm.h'
      include 'var.h'
      DATA GRAV/6.6704D-8/
      save
c

      call start
      imod=0
      konv=0
      mode=0

      do i=1,100000
        imod=imod+1
        TIME=TIME+DTIME
        MODEL=MODEL+1
C
C  Store old values of time dependent variables in arrays with names
C  beginning with a "v".
C
        do j=1,N
          do k=1,NG
            vx(j,k)=x(j,k)
          enddo
          vhydrogen(j)=hydrogen(j)
          vhelium4(j)=helium4(j)
          vhelium3(j)=helium3(j)
          vdeuterium(j)=deuterium(j)
          vcarbon(j)=carbon(j)
          vtrogen(j)=trogen(j)
          voxygen(j)=oxygen(j)
        enddo
C
C  End of storing old values of time dependent variables
C
        konv=max(0,konv-1)
        call massflux(TIME,Zflux)
        if(konv.gt.0. or. MODEL.eq.1) then
          mode=999
        else
          mode=0
          call gridmov
        endif
        dTAUmx= taufactor*taujk
c       call elatim('NOPRINT')

cc  Calculate the efactor value (tidal heating rate) at this time step
cc  Call atmos to find Rplanet...?  Here...?
        rplanet = X(N,2)
        Efacmax = (63.0/4.0) * aval**(-7.5) * (grav * zstar)**1.5
        Efacmax = Efacmax * zstar * (rplanet**5.0) / qval 

c The following line is a kludge        
c        Efacmax = 1e26

        write(6,*) "Efacmax = ",Efacmax


c        eccNow = ecc * ( 1.0 + SIN(frequency * TIME) ) * 0.5 
        eccNow = ecc * SIN(frequency * TIME) 


        efac = Efacmax * (eccNow**2.0)



c        write(6,*) "In main, freq = ",frequency
c        write(6,*) "In main, time = ",time
c        write(6,*) "In main, freq * time = ",frequency*TIME
c        write(6,*) "In main, sin() = ",SIN(frequency*TIME)
c        write(6,*) "In main, ecc = ",ecc
        write(6,*) "In main, eccNow = ",eccNow
c
        write(6,*) "In main, efac = ",efac


cc  Call the Henyey subroutine 
        call henyey(mode)
        if(mode.ne.0) then
          konv=konv+3
c          goto 99
C
C  NO CONVERGENCE IN HENYEY.  TRY AGAIN
C
          TIME=TIME-DTIME
          MODEL=MODEL-1
          imod=imod-1
          if(abs(DTIME/DTMIN-1.).lt.0.01) stop 'NO CONVERGENCE'
          DTIME=DTIME/float(konv)
          DTIME=max(DTIME,DTMIN)
          write(6,*) "NEW DTIME value = ",DTIME
C
C  Reset time dependent variables with old values.
C
          do j=1,N
            do k=1,NG
              x(j,k)=vx(j,k)
            enddo
            hydrogen(j)=vhydrogen(j)
            helium4(j)=vhelium4(j)
            helium3(j)=vhelium3(j)
            deuterium(j)=vdeuterium(j)
            carbon(j)=vcarbon(j)
            trogen(j)=vtrogen(j)
            oxygen(j)=voxygen(j)
          enddo
C
C  End of variable reset.
C
        else
C
C  CONVERGENCE IN HENYEY.  ADVANCE THE COMPOSITION
C

          call compchange  
          call addmass
c          if(konv.lt.2) call addsub
          call printm(imod)
C
C  OPTIMIZE THE TIME STEP, TAKING INTO ACCOUNT THE PREVIOUS CONVERGENCE
C  BEHAVIOR AND THE NET CHANGES.
C
          if(CHANGE.lt.CHGMIN .and. ITER.lt.9 .and. konv.eq.0)
     *                         DTIME=DTIME*1.1
          if(CHANGE.gt.CHGMAX) DTIME=DTIME*0.8
          if(CHANGE.gt.2.*CHGMAX) DTIME=DTIME*0.8
          DTIME=max(DTIME,DTMIN)
          DTIME=min(DTIME,DTMAX)
        endif
c.        if((nmod.gt.0 .and. imod.ge.nmod) .or. konv.ge.7) goto 99
c.        write(6,*) "konv = ",konv
        if((nmod.gt.0 .and. imod.ge.nmod) .or. konv.ge.18) goto 99
      enddo
 99   call printm(-1)
c. 99   continue
      stop
      end

      SUBROUTINE PRINTM(imod)
C+
C  THE SUBROUTINE PRINTM WRITES OUT MODEL INFORMATION
C
C  AUTHOR: H.W. YORKE   30-AUG-02               (JPL / CALTECH)
C-
      include 'parm.h'
      include 'var.h'
      CHARACTER*1 CS(2)
      DATA MODWRT/-1/
      DATA CS/'*',' '/
      save
C
      if(mod(IMOD,NRIT).ne.0 .and. mod(MODEL,JRIT).ne.0) return
      if(IMOD.gt.0) then
        te4 = X(N,3)/7.125e-4/X(N,2)/X(N,2) 
        tee=sqrt(sqrt(te4))
c       Znet=1.d0-hydrogen(1)-helium4(1)-helium3(1)-deuterium(1)        &
c    &      -carbon(1)-trogen(1)-oxygen(1)
        if(iter.lt.abs(itmax/2)) then
          write (6,200)  MODEL,TIME,DTIME,CHANGE,                       &
     &                 X(1,1),X(1,4),X(N,2),X(N,3),tee,Zmass,Zflux 
c    &        hydrogen(1),helium4(1),helium3(1),deuterium(1),           &
c    &        carbon(1),oxygen(1),Znet,                                 &
200       format(' MODEL:',i6,'  TIME:',1p,D11.4,' DTIME:',D11.4,       &
     &     ' NET CHANGE: ',D9.2,                                        &
     &     /,' Pc,Tc,R,L,TE,M:',6D10.3,D11.3                            &
     &     /,1x,83("*"))
c    &     /,' H,He4,He3,D,C,O:',0p,7f8.5,1p,                           &
        else
          write (6,300)  MODEL,TIME,DTIME,CHANGE,                       &
     &                 X(1,1),X(1,4),X(N,2),X(N,3),tee,Zmass,Zflux
300       format(' MODEL:',i6,'  TIME:',1p,D11.4,' DTIME:',D11.4,       &
     &     ' NET CHANGE: ',D9.2,                                        &
     &     /,' XXXXXXXXXXXXXX:',6D10.3,D11.3                            &
     &     /,1x,83("*"))
        endif
      endif 


      if((mod(IMOD,nrit).eq.0 .or. IMOD.lt.0) .and. MODEL.ne.MODWRT)    &
     &  then
        MODWRT=MODEL
        if(deuterium(1).lt.1.d-37) then
          if(helium3(1).lt.1.d-37) then
             write(6,*) "Start of model ",MODEL
            write (6,401) MODEL,TIME
          else
             write(6,*) "Start of model ",MODEL
            write (6,301) MODEL,TIME
          endif
        else
           write(6,*) "Start of model ",MODEL
          write (6,201) MODEL,TIME
        endif
401     format(/,' MODEL:',i6,'  TIME:',1p,D12.4,/                      &
     &         '    J       dM         M         P        R    ',       &
     &         '        L         T         RHO    1-BETA     HYDR  ',  &
     &         '   HE4       C        N        O')
301     format(/,' MODEL:',i6,'  TIME:',1p,D12.4,/                      &
     &         '    J       dM         M         P        R    ',       &
     &         '        L         T         RHO    1-BETA     HYDR  ',  &
     &         '   HE4      HE3       N        O')
201     format(/,' MODEL:',i6,'  TIME:',1p,G25.15,/                     &
     &         '    J       dM         M         P        R    ',       &
     &         '        L         T         RHO    1-BETA     HYDR  ',  &
     &         '   HE4       D        N        O')
        do j=1,N
          call invstate(91,x(j,1),x(j,4),hydrogen(j),helium4(j),Crad,   &
     &         RHO,cP,alpha,beta,delta)
          if(deuterium(1).lt.1.d-37) then
            if(helium3(1).lt.1.d-37) then
            write (6,202) J,CS(1+ICONVECT(j)),dM(j),zM(j),              &
     &                    (X(j,k),k=1,NG),rho,beta,                     &
     &                    hydrogen(j),helium4(j),carbon(j),             &
     &                    trogen(j),oxygen(j)
            else
            write (6,202) J,CS(1+ICONVECT(j)),dM(J),zM(J),              &
     &                    (X(j,k),k=1,NG),rho,beta,                     &
     &                    hydrogen(j),helium4(j),helium3(j),            &
     &                    trogen(j),oxygen(j)
            endif
          else
            write (6,202) J,CS(1+ICONVECT(j)),dM(J),zM(J),              &
     &                    (X(j,k),k=1,NG),rho,beta,                     &
     &                    hydrogen(j),helium4(j),deuterium(j),          &
     &                    trogen(j),oxygen(j)
          endif
C202      format(I5,1x,a1,0p,F9.6,1p,12D9.2)
c202       format(I5,1x,a1,1p,D9.2,D13.6,D10.3,D9.2,D13.5,2D10.3,6D9.2)
 202      format(I5,1x,a1,1p,13G25.15)
        enddo
        write(6,*) "End of model ",MODEL
        write(6,*) "Start of atmos ",MODEL
        if(NATM.gt.0 .and. IMOD.ne.0) then
          call atmos(Tatm,RHOatm,Ratm,Patm,X(N,2),X(N,3),NATM)
        else
          call atmos(Tatm,RHOatm,Ratm,Patm,X(N,2),X(N,3),0)
        endif
        write(6,*) "End of atmos ",MODEL
        if(MODEL.ne.0) then
          write (JUNIT) MODEL,N,TIME,dTIME,taujk,                       &
     &     ((x(j,k),k=1,4),dM(j),                                       &
     &     hydrogen(j),helium3(j),helium4(j),deuterium(j),              &
     &     carbon(j),trogen(j),oxygen(j),                               &
     &     j=1,N)
c         backspace JUNIT
c         read (JUNIT)
          NREC=NREC+1
          write (6,204) MODEL,TIME,NREC,JUNIT
204       format(' MODEL:',i6,'  TIME:',1p,D12.4,' STORED AS RECORD',   &
     &     i5,' ON UNIT',i3)
        endif

      endif

      return
      end
      SUBROUTINE START
C+
C  THE SUBROUTINE START PREPARES FOR THE FIRST RUNNUNG OF HYDRO.
C
C  AUTHOR: H.W. YORKE   30-JAN-02               (JPL / CALTECH)
C-
      include 'parm.h'
      include 'var.h'
      character*50 SFILE,OFILE
      character*1  why(3)
      dimension ZH(14)
      dimension dens(MJ),temp(MJ), xl(MJ), r(MJ) ,pres(MJ)
      equivalence (X(1,1),pres) , (X(1,4),temp)
      equivalence (X(1,2),r) , (X(1,3),xl)
      COMMON/ABUND/XBA(14),H1(14),AH(14)
      COMMON/TIDALINTS/sigma,xmass
c      COMMON/TIDALDBLS/frequency
c      COMMON/ORBVALS/ecc,aval,orbn,qval,zstar,efac
      
      DATA ZH     /1.0081451,4.003874,12.0038156,23.,24.32,             &
     & 26.97,28.06,32.07,39.102,40.08,55.85,14.0075257,16.,19.99/
      DATA RSUN/6.96d10/
      DATA PI  /3.14159265359d0/ 

      save
c
      do i=1,14
        ah(i)=zh(i)
      enddo
c
      call comrd
      READ (5,200) IUNIT,SFILE
      WRITE(6,200) IUNIT,SFILE
      call comrd
      READ (5,200) JUNIT,OFILE
      WRITE(6,200) JUNIT,OFILE
c
      if(IUNIT.eq.1) then
        open(unit=IUNIT,file=SFILE, form='formatted' , status='old' ) 
      else
         write(6,*) SFILE
        open(unit=IUNIT,file=SFILE, form='unformatted' , status='old' )
      endif
      if(IUNIT.ne.JUNIT) then              
         open(unit=JUNIT,file=OFILE, form='unformatted' ,               &
     &  status='unknown' )
      endif
C
 200  FORMAT(20X,i3,1X,a)
 201  FORMAT(5(6X,I6))
 202  FORMAT(1P,4(6X,D9.2))
 203  FORMAT(4(6X,D9.2))

C
      call comrd
      READ (5,201) NREC,NMOD,NRIT,ITMIN,ITMAX
c     WRITE(6,201) NREC,NMOD,NRIT,ITMIN,ITMAX
      call comrd
      READ (5,201) JADD,JSUB,NATM
c     WRITE(6,201) JADD,JSUB,NATM
      call comrd
      READ (5,203) Atmx,Atmn,dTAUmx,RLH
c     WRITE(6,202) Atmx,Atmn,dTAUmx,RLH
      taufactor=dTAUmx
      call comrd
      READ (5,203) dLmx,dLmn,dXmx,dXmn
c     WRITE(6,202) dLmx,dLmn,dXmx,dXmn
      call comrd
      READ (5,203) dPmx,dPmn,Crad,Cwrk
c     WRITE(6,202) dPmx,dPmn,Crad,Cwrk
      call comrd
      READ (5,203) dZmax,dZmin,dZdt
      dZmax=dZmax*1.00001
c     WRITE(6,202) dZmax,dZmin,dZdt
      call comrd
      READ (5,203) EPS
c     WRITE(6,202) EPS
      call comrd
      READ (5,203) SMIN
c     WRITE(6,202) SMIN
      call comrd
      READ (5,203) SMAX
c     WRITE(6,202) SMAX
C
      call comrd
      READ (5,203) DTIME,FACTIM,DTMIN,DTMAX
c     WRITE(6,202) DTIME,FACTIM,DTMIN,DTMAX
      call comrd
      READ (5,203) CHGMIN,CHGMAX
c     WRITE(6,202) CHGMIN,CHGMAX

      call comrd

      READ (5,203) x1,x2
      xmass = x1
      sigma = x2
      call comrd
      READ(5,203) xperiod, zstar
      xperiod = 3.1536E+7 * xperiod
      frequency = 2.0*PI / xperiod

      call comrd
      READ(5,203) x1,x2,x3,x4
      ecc = x1
      aval = x2
      orbn = x3
      qval = x4

      call comrd
      READ (5,*)   (H1(J), J = 1,4) 
      call comrd
      READ (5,*)   (H1(J), J = 5,8) 
      call comrd
      READ (5,*)   (H1(J), J = 9,12) 
      call comrd
      READ (5,*)   (H1(J), J = 13,14) 
      call comrd
      NG=MH
      JRIT = max(1,NRIT/10)
C
      XH2=H1(2)
      H1(2)=abs(H1(2))
      SUM = 1.
      DO J = 1,14
        SUM = SUM - H1(J)
      END DO
      XX=H1(1)
      YY=SUM
C H1(2) is temporarily used as the primordial deuterium abundance and
C then replaced by the initial helium abundance.
      H1(2)=YY
      DO J = 1,14
        XBA(J) = H1(J) / AH(J)
      END DO
      ZZ= 1. - XX - YY 
      IF(IUNIT.EQ.1) THEN
        NREC=0
        MODEL=0
C SETUP STARTING MODEL
        do j=1,MJ
          read(IUNIT,*,end=34) zM(j),r(j),temp(j),xl(j),dens(j)
          if(r(j).le.0.d0) goto 35
        enddo
        J=MJ+1 
        GOTO 35
  34    backspace(IUNIT)
  35    N=min(J-1,MJ)
        close(IUNIT)
        Zmass=zM(N+1)
C**************
        amu = 2.*XX + 0.75*YY + 0.57*ZZ
        Rglog=log10(83145100.d0*amu)
        do i = 1,N
          r(i) = RSUN*r(i) 
          tlog = log10(temp(i))
          hydrogen(i) = XX
          helium4(i) = YY
          helium3(i) = 1.d-5 
          deuterium(i)=abs(XH2)
          carbon(i)=H1(3)
          trogen(i)=H1(12)
          oxygen(i)=H1(13)
          denlog = log10(dens(i))
          xl(i) = xl(i)*3.85d33 
          Plog = Rglog + denlog + Tlog
          pres(i)=10.d0**Plog
          zM(i)=zM(i)*zMass
        end do 
        dM (1)=zM(1)
        do i=2,N
          dM (i)=zM(i)-zM(i-1)
        enddo
        TIME=0.
        taujk=1.D11
        dTAUmx=taujk*taufactor
        call printm(MODEL)
      ELSE
        if(NREC.lt.0) NREC=99999
        nn=max(1,NREC)
        do imod=1,nn
          read(IUNIT,end=301) MODEL,N,TIME,dTIM,taujk,                  &
     &     ((x(j,k),k=1,4),dM(j),                                       &
     &     hydrogen(j),helium3(j),helium4(j),deuterium(j),              &
     &     carbon(j),trogen(j),oxygen(j),                               &
     &     j=1,N)
        enddo


        IMOD=NREC+1
        goto 302
 301    backspace(IUNIT)
 302    NREC=IMOD-1


       zM(1)=dM(1)
        do i=2,N
          zM(i)=zM(i-1)+dM(i)
        enddo
        if(NREC.eq.0) then
          TIME=0.d0
          MODEL=0
        endif
        write(6,209)  NREC,MODEL,N,TIME,dTIM
 209    format(' RECORD',i4,' MODEL:',i6,' N:',i4,' TIME:',             &
     &           1p,D11.4,' DTIME:',D11.4)
        if(JUNIT.ne.IUNIT) NREC=0
        if(DTIME.lt.0.0d0) DTIME=dTIM
        write(6,*) '... continuing calculations with MODEL=',MODEL
        Zmass=zM(N)
        write(*,*) 'Zmass=', Zmass
        if(Jadd.gt.0) then
          why(1)='I'
          call add(Jadd,why)
        endif
        if(Jsub.gt.0) call sub(Jsub)
        if(FACTIM.lt.0.d0) TIME=0.
      ENDIF

c. Following lines were added for the kozai energy input modifications to the code:
      TIME=0.0
c      DTIME = min(dtim, (xperiod / 20.0))
      DTMAX = xperiod/20.0
      dTIME = dtmax/2.0
      write(6,*) "Initial dtmax = ",dtmax
      write(6,*) "Initial DTIME = ",dtime
      xmass= xmass*zMass
      sigma = sigma*zMass
      RETURN
      END

      subroutine comrd
      character*1 line(1)
      character*80 lin
      equivalence (lin,line)
      save
c
1     read (5,203,end=999) lin
      if(line(1).eq.' ') then
      write(6,203)         lin
        backspace 5
        return
      else
        if(line(1).ne.'c' .and. line(1).ne.'C')                         &
     &  write(6,203)  lin
c202    format(1x,a)
      endif
      goto 1
 999  write(6,203) 'COMRD: END OF INPUT DATA'
 203  format(a)
      return
      end

      subroutine opacity(j,tlog,rholog,bkap)
      include 'parm.h' 
      include 'var.h' 
      parameter(NTVAL=60,NRVAL=14)
      DIMENSION TVAL(NTVAL),RVAL(NRVAL)
      DATA TVAL/
     +2.000,2.097,2.176,2.243,2.301,2.352,2.398,2.439,
     +2.477,2.512,2.544,2.574,2.602,2.628,2.653,2.677,
     +2.699,2.740,2.778,2.813,2.845,2.875,2.903,2.929,
     +2.954,2.978,3.000,3.021,3.041,3.061,3.079,3.097,
     +3.114,3.130,3.146,3.161,3.176,3.204,3.230,3.255,
     +3.279,3.301,3.350,3.400,3.450,3.500,3.550,3.600,
     +3.650,3.700,3.800,3.900,4.000,4.079,4.176,4.301,
     +4.477,4.699,4.845,5.000/
      DATA RVAL/-12.,-11.,-10.,-9.,-8.,-7.,-6.,-5.,-4.,-3.,-2.,-1.,0.,1.
     1/
      data ifirst/0/
      save

c     if(j.lt.100) write(6,*) 'opacity: j,tlog,rhol=',j,tlog,rhol
      if(ifirst.eq.0) then
        ifirst=1

C....read in the low temperature opacity 
        open(unit=52, file='opac.cool', status='old')
        do itemp=1,60
          READ(52,46) (ZKAP(idense,itemp), idense=1,7)
c         WRITE(6,46) (ZKAP(idense,itemp), idense=1,7)
          READ(52,46) (ZKAP(idense,itemp),idense=8,14)
c         WRITE(6,46) (ZKAP(idense,itemp),idense=8,14)
46        FORMAT(7F7.3)
        end do
        close(52)
        open(unit=50,file='hyd.cond',status='old')
        open(unit=51,file='hel.cond',status='old')


c.....read in hydrogen conductive opacities
        read(50,*) (rhyr(i),i=1,23)
        read(50,*) (thyr(i),i=1,28)
        do i=1,28
          read(50,*) (auxchyop(k,i),k=1,23)
        end do
c.....read in helium conductive opacities
        read(51,*) (rhel(i),i=1,25)
        read(51,*) (thel(i),i=1,32)
        do i=1,28
          read(51,*) (auxcheop(k,i),k=1,25)
        end do

        close(50)
        close(51) 
      endif

      YY=helium4(j)+helium3(j)

c.....Start here with logT,logRho, Rference logTs logRhos, Table,
c..    log Opacity is Bkap.
      RHOL=rholog
      IF(TLOG .GE. TVAL(NTVAL)) GO TO 206
      call bracket(TVAL,TLOG,NTVAL,mode,it1,it,iterx)
      RHOL=max( RVAL(1) , RHOL )
      RHOL=min( RVAL(NRVAL)-1.E-9 , RHOL )

      WS1 = TVAL(IT1) - TLOG                                                1900
      WS = TVAL(IT1) - TVAL(IT)                                             2000
      WS1 = WS1 / WS                                                        2100
      WS = 1. - WS1                                                         2200

      call bracket(RVAL,RHOL,NRVAL,mode,iw1,iw,iterx)
      WS2 = RHOL - RVAL(IW1)                                                2800
      WS3 = RVAL(IW) - RVAL(IW1)                                            2900
      WS2 = WS2 / WS3                                                       3000
      WS3 = 1. - WS2                                                        3100
      Z00 = ZKAP(IW,IT)                                                     3200
      Z10 = ZKAP(IW1,IT)                                                    3300
      Z11 = ZKAP(IW1,IT1)                                                   3400
      Z01 = ZKAP(IW,IT1)                                                    3500
      IF(Z00.EQ.0..OR.Z01.EQ.0..OR.Z11.EQ.0..OR.Z10.EQ.0.) GO TO 200        3600
      WS4 = WS2*(WS1*Z00 + WS*Z01)                                          3700
      WS5 = WS3*(WS1*Z10 + WS*Z11)                                          3800
      BKAP= (WS4+WS5)  

C  Slowly transition into OPAL opacities
      if(Tlog.gt.TVAL(NTVAL)-.1) then
        W1=(TVAL(NTVAL)-Tlog)*10.
        W2=1.d0-W1
        call opaltab(Tlog,RHOlog,hydrogen(j),YY,bkap2,ierror)
        if(ierror.eq.0)  BKAP=W1*BKAP+W2*bkap2
      endif
      GO TO 205

 206  continue
      call opaltab(Tlog,RHOlog,hydrogen(j),YY,bkap,ierror)
c     if(ierror.ne.0) Bkap= 4.7

c.....add in conductive opacity:
      if(rhol.gt.0. .and. tlog.gt.5.0)then
        if(hydrogen(j) .le. 1.e-5) then
          hyckap=1.
        else
          call chyopacity(tlog,rhol,hyckap)
        end if
        call cheopacity(tlog,rhol,heckap)
        ckap=hydrogen(j)*hyckap+(1.-hydrogen(j))*heckap
        ckap=10**ckap
        realbkap = 10.**bkap
        realbkap=1./(1./realbkap+1./Ckap)
        bkap=log10(realbkap)
      end if

205   continue                 
      RETURN                   
200   write(6,201) TLOG,RHOL    
201   FORMAT('  OUTSIDE TABLE, T = ', 1PE12.3, '  RHO= ',E12.3)
      STOP 'opacity'
      END

           Subroutine opacityhy(T,rho,realk)
      include 'parm.h' 
           COMMON/AUXhyOPACITY/RHOREF(18),TREF(15),rkapparef(18,15)
c..........This subroutine computes the opacity for a specified temperature
c..........and pressure by linearly interpolating between the proper points
c..........of a pure hydrogen opacity table.

           data its,irhs /15,18/
      save

           if ((T.le.Tref(1)).or.(T.ge.Tref(its)).or.(Rho.le.Rhoref(1))
     +     .or.(RhO.ge.Rhoref(irhs))) then
                  write(6,*) 'Temperature: ',T,' Density: ',Rho
                  write(6,*) 'Out of Bounds.'
                 realk=1.
C                goto 2104
                 stop 'opacityhy'
           endif

c..........Determine which of the reference points bracket the temperature.

           call bracket(Tref,T,its,mode,icool,ihot,iterx)
           call bracket(Rhoref,Rho,irhs,mode,jrare,jdense,iterx)

c..........Determine Krare and Kdense by a linear interpolation in temperature.

           rare = rkapparef(jrare,icool)
     +             +((rkapparef(jrare,ihot)-rkapparef(jrare,icool))/
     +             (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

           dense = rkapparef (jdense,icool)
     +              +((rkapparef(jdense,ihot)-rkapparef(jdense,icool))/
     +              (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

c.........Determine kappa by a linear interpolation in the density

          realk=rare+((dense-rare)/(rhoref(jdense)-rhoref(jrare)))
     +              *(rho - rhoref(jrare))

c.........We have thus determined k.
2104   continue
c       write(*,*)  ihot,icool,jdense,jrare,krare,kdense,k


c       write(*,*) 'rkapp(jr,ic),rk(jr,ih),rk(jd,ic),rk(jd,ih)'
c       write(*,*)  rkapparef(jrare,icool),rkapparef(jrare,ihot)
c       write(*,*)  rkapparef(jdense,icool),rkapparef(jdense,ihot)
c       write(*,*) 'tref(ihot),tref(icool),rhoref(jdense),rhoref(jrare)'

c       write(*,*) tref(ihot),tref(icool),rhoref(jdense),rhoref(jrare)


        return
        end
 

           Subroutine opacityhe(T,rho,realk)

c..........This subroutine computes the opacity for a specified temperature
c..........and pressure by linearly interpolating between the proper points
c..........of a pure helium opacity table.
      include 'parm.h' 
           COMMON/AUXheOPACITY/RHOREF(18),TREF(15),rkapparef(18,15)
           data its,irhs /15,18/
      save

           if ((T.le.Tref(1)).or.(T.ge.Tref(its)).or.(Rho.le.Rhoref(1))
     +     .or.(RhO.ge.Rhoref(irhs))) then
                  write(6,*) 'Temperature: ',T,' Density: ',Rho
                  write(6,*) 'Out of Bounds.'
                 realk=1.
C                goto 2104
                 stop 'opacityhe'
           endif

c..........Determine which of the reference points bracket the temperature.

           call bracket(Tref,T,its,mode,icool,ihot,iterx)
           call bracket(Rhoref,Rho,irhs,mode,jrare,jdense,iterx)

c..........Determine Krare and Kdense by a linear interpolation in temperature.

           rare = rkapparef(jrare,icool)
     +             +((rkapparef(jrare,ihot)-rkapparef(jrare,icool))/
     +             (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

           dense = rkapparef (jdense,icool)
     +              +((rkapparef(jdense,ihot)-rkapparef(jdense,icool))/
     +              (Tref(ihot)-Tref(icool)))*(T-Tref(icool))



c.........Determine kappa by a linear interpolation in the density

          realk=rare+((dense-rare)/(rhoref(jdense)-rhoref(jrare)))
     +              *(rho - rhoref(jrare))

c.........We have thus determined k.
2104   continue
c       write(*,*)  ihot,icool,jdense,jrare,krare,kdense,k
c       write(*,*) 'rkapp(jr,ic),rk(jr,ih),rk(jd,ic),rk(jd,ih)'
c       write(*,*)  rkapparef(jrare,icool),rkapparef(jrare,ihot)
c       write(*,*)  rkapparef(jdense,icool),rkapparef(jdense,ihot)
c       write(*,*) 'tref(ihot),tref(icool),rhoref(jdense),rhoref(jrare)'

c       write(*,*) tref(ihot),tref(icool),rhoref(jdense),rhoref(jrare)


        return
        end

           Subroutine chyopacity(T,rho,realk)

c..........This subroutine computes the opacity for a specified temperature
c..........and pressure by linearlly interpolating between the proper points
c..........of the Hubbard-Lampe (gL-padded) hydrogen conductive opacity table

      include 'parm.h' 
           common /ahycopacity/ Rhoref(23),Tref(28),rkapparef(23,28)
           data its,irhs /28,23/
      save

           condmin=rkapparef(1,its)
c          write(*,*) 'chyopacity: Temperature: ',T,' Density: ',Rho
           if (T.le.Tref(1).or.Rho.le.Rhoref(1)) then
             write(*,*) 'Temperature: ',T,' Density: ',Rho
             write(*,*) 'Out of Bounds in conductive hydrogen O-table.'
             stop 'chyopacity'
           endif

c..........Determine which of the reference points bracket the temperature.

           call bracket(Tref,T,its,mode,icool,ihot,iterx)
           call bracket(Rhoref,Rho,irhs,mode,jrare,jdense,iterx)

c..........Determine Krare and Kdense by a linear interpolation in temperature

           rare = rkapparef(jrare,icool)
     +             +((rkapparef(jrare,ihot)-rkapparef(jrare,icool))/
     +             (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

           dense = rkapparef (jdense,icool)
     +              +((rkapparef(jdense,ihot)-rkapparef(jdense,icool))/
     +              (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

c.........Determine kappa by a linear interpolation in the density

          realk=rare+((dense-rare)/(rhoref(jdense)-rhoref(jrare)))
     +              *(rho - rhoref(jrare))

          realk=min(realk,condmin)

c.........We have thus determined k.

          return
          end

           Subroutine cheopacity(T,rho,realk)

c..........This subroutine computes the opacity for a specified temperature
c..........and pressure by linearlly interpolating between the proper points
c..........of the Hubbard-Lampe (gL-padded) helium conductive opacity table

      include 'parm.h' 
           common /ahecopacity/ Rhoref(25),Tref(32),rkapparef(25,32)
           data its,irhs /32,25/
      save
c
             condmin=rkapparef(1,its)
           if(T.le.Tref(1).or.Rho.le.Rhoref(1)) then
             write(*,*) 'Temperature: ',T,' Density: ',Rho
             write(*,*) 'Out of Bounds in conductive helium O-table.'
             stop 'cheopacity'
           endif

c..........Determine which of the reference points bracket the temperature.

           call bracket(Tref,T,its,mode,icool,ihot,iterx)
           call bracket(Rhoref,Rho,irhs,mode,jrare,jdense,iterx)

c..........Determine Krare and Kdense by a linear interpolation in temperature

           rare = rkapparef(jrare,icool)
     +             +((rkapparef(jrare,ihot)-rkapparef(jrare,icool))/
     +             (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

           dense = rkapparef (jdense,icool)
     +              +((rkapparef(jdense,ihot)-rkapparef(jdense,icool))/
     +              (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

c.........Determine kappa by a linear interpolation in the density

          realk=rare+((dense-rare)/(rhoref(jdense)-rhoref(jrare)))
     +              *(rho - rhoref(jrare))

          realk=max(realk,condmin)

c.........We have thus determined k.

          return
          end
C+
C  NAME:    bracket (Version 1.0)
C  AUTHOR:  H.W. Yorke (JPL)
C  DATE:    24-Feb-05     V1.0
C  UPDATES:
C
C  The subroutine bracket finds the upper and lower indices of the
C  members of a monotonic increasing table which bracket a given number.
C
C  *** NOTE:    bracket does not check whether the table is monotonically
C               increasing.
C
C  USAGE: call bracket (Xarray , Xval, Narray , mode , I1 , I2 , ITER )
C
C  WHERE    Xarray   is the given monotonically increasing table (INPUT)
C           Xval     is the given number to be bracketed         (INPUT)
C           Narray   is the length of the table                  (INPUT)
C           mode     = -1, 0 , 1 depending on whether Xval is     OUTPUT
C                    less than the smallest value in the table
C                    (mode=-1), greater than the maximum value
C                    (mode=1) or lies within the table (mode=0)
C           I1       is the lower index of bracketing             OUTPUT
C           I2       is the upper index of bracketing             OUTPUT
C           ITER     is the number of iterations                  OUTPUT
C
C  PROGRAMS USED: none
C-
      subroutine bracket(Xarray,Xval,Narray,mode,I1,I2,iter)
      include 'parm.h' 
      dimension Xarray(Narray)
      save
c
      i=0
      mode=0
      if(Xval.lt.Xarray(1)) then
        I1=1
        I2=2
        mode=-1
        return
      endif
      if(Xval.gt.Xarray(Narray)) then
        I1=Narray-1
        I2=Narray
        mode=1
        return
      endif
      I1=1
      I2=Narray
      do i=1,Narray
        ID=(I1+I2)/2
        if(Xarray(ID).lt.Xval) then
          I1=ID
        else
          I2=ID
        endif
        if(I2-I1.eq.1) goto 90
      enddo
      stop 'bracket: Unable to find position in table'
 90   iter=i
      continue
      return
      end

      subroutine opaltab(Tlog,RHOlog,XH,Y,opac,ierror)
      include 'parm.h' 
      include 'var.h' 
      parameter (MX=14,MR=19,MT=70)
      character*5 TABLE
      dimension OPACTB(MX,MR,MT),XX(MX),RR(MR),TT(MT)
      COMMON/ABUND/XBA(14),H1(14),AH(14)
      data ifirst/0/
      data Z0,Z1/0.,1./
      save

      ierror=0
      if(ifirst.eq.0) then
        ifirst=1
        open(51,file='GN93hz',form='formatted',status='old',err=99)
        do i=1,239
          read(51,*)
        enddo
        do ix=1,MX
          read(51,*)
 201      format(a5,31x,f6.4,12x,f6.4,2(5x,f6.4))
          read(51,201,err=98) TABLE,XX(ix),ZZ,XC,XO
c         write(6,201)        TABLE,XX(ix),ZZ,XC,XO
          if(TABLE.ne.'TABLE') goto 98
C
C  XX(ix) is normally used to store the hydrogen mass content used to
C  generate the table.  However, for X=0 (no hydrogen) then -XX(ix) is
C  used to store the amount of additional oxygen and carbon.
C
          if(XC+XO.gt.0.0001) then
            XX(ix)=-(XC+XO)
          else
            ZZtab = ZZ
          endif
          read(51,*)
          read(51,*)
          read(51,*)
 202      format(4x,f6.1,18f7.1)
          read(51,202,err=98) RR
c         write(6,202) RR
          read(51,*)
          do it=1,MT
 203        format(f4.2,19f7.3)
            read(51,203,err=98) TT(it),(OPACTB(ix,ir,it),ir=1,MR)
c           write(6,203) TT(it),(OPACTB(ix,ir,it),ir=1,MR)
          enddo
c         write(6,203) TT
        enddo
        close(51)
        dr=(RR(MR)-RR(1))/float(MR-1)
        IXmax=1
        XXmax=XX(1)
        do ix=1,MX
          if(XX(ix).gt.XXmax) then
            IXmax=ix
            XXmax=XX(ix)
          endif
        enddo
        Znorm = Z1-H1(1)-H1(2)
      endif
 
      if(XH.lt.XX(1) .or. XH.gt.XXmax) ierror=1
      dZZ=((Z1-XH-Y)-Znorm) * (Z1-ZZtab)/(Z1-Znorm)
      if(dZZ.le.0.0014d0) then
        IX1=IXmax-1
        IX2=IXmax
        if(XH.lt.XX(IX1)) then
          IX1=1
          IX2=2
          if(XH.gt.XX(IX2)) then
            IUP=IXmax-1
            IDN=2
            do i=1,2
              IX1=(IUP+IDN)/2
              if(XH.gt.XX(IX1)) then
                IDN=IX1
              else
                IUP=IX1
              endif
           enddo
           IX1=IDN
           IX2=IUP
          endif
        endif
        DX2=(XH-XX(IX1))/(XX(IX2)-XX(IX1))
        DX1=Z1-DX2
      else
        if(XH.gt.1.d-20) then
          ierror=99
          write(6,*) 'You have no business being in this part of the',  &
     &     ' opacity table  ZZtab,dZZ=',ZZtab,dZZ
          write(6,*) 'X,Y,Znorm=',XH,Y,Znorm
          write(6,'(a,1p,14E11.3)') 'Abundances:',H1
          do j=1,N
            if(abs(XH-hydrogen(j)).lt.1.d-10) then
              write(6,'(a,i4,9F8.5)') 'J,lnT,X,Y4,Y3,D,C,N,O=',j,Tlog,  &
     &          XH,helium4(j),helium3(j),deuterium(j),carbon(j),        &
     &          trogen(j),oxygen(j),1.d0-XH-helium4(j)-helium3(j)-      &
     &          deuterium(j)-carbon(j)-trogen(j)-oxygen(j)
            endif
          enddo
          stop 'opaltab'
        endif
        IX1=MX-1
        IX2=MX
        if(dZZ.lt.-XX(IX1)) then
          IX1=1
          IX2=IXmax+1
          if(dZZ.gt.-XX(IX2)) then
            IUP=MX-1
            IDN=IXmax+1
            do i=1,2
              IX1=(IUP+IDN)/2
              if(dZZ.gt.-XX(IX1)) then
                IDN=IX1
              else
                IUP=IX1
              endif
           enddo
           IX1=IDN
           IX2=IUP
          endif
        endif
        DX2=(dZZ+XX(IX1))/(-XX(IX2)+XX(IX1))
        DX1=Z1-DX2
      endif
c     write(6,204) 'X:',XH,XX(ix1),XX(ix2),ix1,ix2,dx1,dx2
 204  format(a,1p,3E12.4,2i10,2E12.4,i3)

      Rlog=RHOlog-3.*(Tlog-6.)
      if(Rlog.lt.RR(1) .or. Rlog.gt.RR(MR)) ierror=ierror+2
      DR2=(Rlog-RR(1))/dr
      IR1=max(DR2,Z0)
      IR2=min(IR1+2,MR)
      IR1=IR2-1
      DR2=DR2-float(IR1-1)
      DR1=Z1-DR2
c     write(6,204) 'R:',Rlog,RR(ir1),RR(ir2),ir1,ir2,dr1,dr2

      if(Tlog.lt.TT(1) .or. Tlog.gt.TT(MT)) ierror=ierror+4
      call bracket(TT,Tlog,MT,mode,IT1,IT2,iterx)
      DT2=(Tlog-TT(IT1))/(TT(IT2)-TT(IT1))
      DT1=Z1-DT2
c     write(6,204) 'T:',Tlog,TT(it1),TT(it2),it1,it2,dt1,dt2,iterx

      opac=dx1*(dr1*(dt1*OPACTB(ix1,ir1,it1)+dt2*OPACTB(ix1,ir1,it2))
     &        + dr2*(dt1*OPACTB(ix1,ir2,it1)+dt2*OPACTB(ix1,ir2,it2)))
     &    +dx2*(dr1*(dt1*OPACTB(ix2,ir1,it1)+dt2*OPACTB(ix2,ir1,it2))
     &        + dr2*(dt1*OPACTB(ix2,ir2,it1)+dt2*OPACTB(ix2,ir2,it2)))

      return
 98   write(6,*) 'Error reading opacity tables'
      stop 'opaltab'
 99   write(6,*) 'Error opening opacity tables'
      stop 'opaltab'
      end

       subroutine nucrat(Rho,T,np,nd,nhe3,dtime,E,                      &
     &                  Rpp,Rhe3he3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,      &
     &                  Xhe4,Xc,Xn,Xo,xpp23,iflip,iwrite)
      include 'parm.h'

c......
       double precision k,np,nde,nd,N0,kb,Lambda0,mu,nhe3
       double precision n(14),nc,nn,no
       COMMON/ABUND/XCA(14),H1(14),A(14)

       dimension Xz(14),AA(14),XBA(14)
       data Xz/1.d0,2.d0,6.d0,11.d0,12.d0,13.d0,14.d0,16.d0,19.d0,20.d0,&
     &        26.d0,7.d0,8.d0,10.d0/
       DATA AA/1.0081451d0,4.003874d0,12.0038156d0,23.d0,24.32d0,       &
     &        26.97d0,28.06d0,32.07d0,39.102d0,40.08d0,55.85d0,         &
     &        14.0075257d0,16.d0,19.99d0/
       DATA Z0,Z1,Z13,HMASS/0.d0,1.d0,0.333333333333333333d0,1.6732D-24/
       DATA CV,CV2,CA,CA2/.934d0,.8724d0,0.5d0,0.25d0/
       DATA ifirst/0/
       DATA PI  /3.14159265359d0/ 
       COMMON/UsefulInfo/J2,Ntot,Niter
       COMMON/TimeInfo/TimeVal
       COMMON/TIDALINTS/sigma,xmass
       COMMON/TIDALDBLS/frequency
       COMMON/CURRENTVALS/currentMass
c       COMMON/tidalinfo/efac
      COMMON/ORBVALS/ecc,aval,orbn,qval,zstar,efac
      
c......recalculate helium abundance


c......Given the temperature, the density and an abundance set,
c......this routine returns the nuclear energy generation rate
c......in [ergs/sec/gm]. The reactions considered are:

c......proton:proton                  1H(p,e+ve)2H
c......deuterium:proton               2H(p,y)3He
c......He3:He3                        3He(3He,H+H)4He
c......Corrections for ppII & ppIII
c......Equilibrium CNO
c......He4:He4:He4                    4He(aa,y)12C
c......C:He4                          12C(a,y)16O
c......Various high T cooling processes

       save
c................................................................
c       write(6,*) "In nucrat, efac = ",efac
cc Linear energy input scaling
c       ratio = real(J2)/real(Ntot)
c       ratio = 1.0 - ratio
c       if (J2.GT.50) ratio = 0.0
c       E = 1.0 * ratio

cc Gaussian energy input scaling BY MASS
       temp = (currentMass - xmass)**2.0
       temp = temp/(2.0*sigma**2.0)

       temp = -1.0*temp
       temp = exp(temp)
       temp = temp / SQRT(2.0*PI*sigma**2.0)
c       amplitude = efactor * 0.5*(1.0+ SIN(TimeVal*frequency))
c       E = amplitude * temp
       E = efac * temp

cc Gaussian energy input scaling BY CELL NUMBER
c       jnought = int(xmass)
c       sigmasquared = sigma**2.0
c       temp = real(J2-jnought)**2.0 
c       temp = temp/(2.0*sigmasquared)
c       temp = -1.0*temp
c       temp = exp(temp)
c       temp = temp / SQRT(2.0*PI*sigmasquared)
c       amplitude = efactor * 0.5*(1.0+ SIN(TimeVal*frequency))
c       amplitude = efactor * 0.5
c       E = amplitude * temp
c       write(6,*) "E=",E
c. Following line is a kludge...
c       E = 0.5 * efactor



cc Set energy/heating rate to zero
c       E = 0.0
c       write(6,*) "In nucrat, returning E = ",E
       return
c................................................................

       j = abs(iwrite)
       if(ifirst.eq.0) then
         ifirst=1
         do i=1,14
           A(i)=AA(i)
           xba(i)=xca(i)
         enddo
c        xmetal=x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(14)
       endif
c
       Rpp = Z0
       Xpd = Z0
       Rhe3he3=Z0
       Rcno=Z0
       Xpc=Z0
       Xpn=Z0
       Xpo=Z0
       R3a=Z0
       R12a=Z0
       xpp23=Z1
       E=Z0
c. April 3, 2012:
c. Commented out the following if-statement and replaced it with a return statement of my own
c. in order to effectively shut off all the fusion effects within the code
c       if(T.lt.5.d5) return
c       return
c
c......recompute mass fractions
       xh=(np/rho)*HMASS
       xhe3=3.*(nhe3/rho)*HMASS
       xba(1) = xh/a(1) 
       xba(2) = xhe4/a(2) 
       xba(3) = xc/a(3) 
       xba(12)= xn/a(12) 
       xba(13)= xo/a(13) 
       fpd=Z1
       Zd=Z1
       Zp=Z1
       Zhe3 = 2.d0
       k=1.380658e-16

       zeta=Z0
       do i=1,14
         zeta = zeta+(Xz(i)**2+Xz(i))*xba(i)
       end do
C*******  
       if(T.ge.3.d7) then
          T8=T*1.d-8
C  f3a formula from Clayton (1968)
          f3a=exp(2.76d-3*sqrt(RHO/T8)/T8)
C  Triple alpha from Kippenhahn & Weigert (1991)
          E3a=5.1d11 * f3a * RHO**2 * (Xhe4/T8)**3 * exp(-44.027d0/T8)
          Q3a=1.164d-5
          R3a=E3a*RHO/Q3a
C  Carbon:alpha from Kippenhahn & Weigert(1991)
          T7=T8*10.d0
          f12a=exp(0.071d0*sqrt(zeta*RHO/T7**3))
          T823=T8**(2.d0/3.d0)
          g12a=((Z1+0.134d0*T823)/(Z1+0.01d0*T823))**2
          E12a=1.3d27 * f12a * Xc * Xhe4 * RHO / T8**2 * g12a *
     *         exp(-69.20d0/T8**Z13)
          Q12a=1.146d-5
          R12a=E12a * RHO / Q12a
          E  =E3a + E12a
C Neutrino loss rates Beaudet, Petrosian, and Salpeter  1967
          T9=T8/10.d0
          alam=T*1.683d-10
          axi=((rho/2.d0/1.d9)**(.33333333d0))/alam
          fpair1 = (6.002d19 + 2.084d20*axi + 1.872d21*axi**2)
     &     *exp(-5.5924*axi)
          fpair2=axi**3 + 9.383d-1/alam - 4.141d-1/alam/alam
     &     + 5.829d-2/alam**3
          fpair=fpair1/fpair2
          alam2=alam*alam
          alam4=alam2*alam2
          alam6=alam4*alam2
          alam8=alam6*alam2
          glam=1.d0 -13.04d0*alam2+133.5d0*alam4 
     &    + 1534.d0*alam6 + 918.6d0*alam8 
          epair=glam/rho*fpair*exp(-2.d0/alam)
C         Z2A=xhe4 + xc*3.d0 + xo*4.d0
C         ebrems=0.76d0*Z2A*T8**6
C         if(rho.lt. 1.d8) then
C         factor=2.5d0*(8.d0- log10(rho)) -1.d0
C    &      +0.25d0*log10(rho)
C         ebrems = ebrems/factor
C         end if 
          ebrems=0.
          fphot1 = (4.886d10 + 7.580d10*axi + 6.023d10*axi**2)
     &     *exp(-1.5654*axi)
          fphot2=axi**3 + 6.29d-3/alam + 7.483d-3/alam/alam
     &     + 3.061d-4/alam**3
          fphot=fphot1/fphot2
C    ********  mue=2  *****************
          ephot=0.5*alam**5*fphot
C         write(6,*) ' ephot = ', ephot, ' ebrems = ', ebrems
C         write(6,*) ' epair = ', epair
          fplas1 = (2.320d-7 + 8.449d-8*axi + 1.787d-8*axi**2)
     &     *exp(-.5646*axi)
          fplas2=axi**3 + 2.581d-2/alam + 1.734d-2/alam/alam
     &     + 6.990d-4/alam**3
          fplas=fplas1/fplas2
          eplas=rho*rho/8.d0*fplas 
C         write(6,*) ' eplas = ', eplas , ' eps = ', epsneu   
C     Munakata, Kohyama, and Itoh 1985 corrections
   
          eplas = eplas*CV2
          qphot=1.d0+rho/2.d0*(1.875d8*alam + 1.653d8*alam2
     & + 8.499d8*alam2*alam -1.604d8*alam4)**(-1.d0)
          qphot=0.666d0/qphot*(1.d0 + 2.045d0*alam)**(-2.066d0)
          fphot=fphot*0.893d9*alam4*alam4*(1.d0+143.8d0*alam**(3.555d0)
     &    )**(-0.3516d0)*axi**3*exp(0.556d0*axi**4.48d0/
     &    (150.d0+axi**3.3d0))
          fphot=fphot*0.5d0*((CV2+CA2))*(1.d0 - (CV2-CA2)/(CV2+CA2)
     &    *qphot)
          ephot=fphot/rho

C
          qpair=(10.7480d0*alam2 + 0.3967d0*sqrt(alam)
     &    +1.005d0)**(-1.d0)*(1.d0+rho/2.d0*(7.692d7*alam*alam2
     &    +9.715d6*sqrt(alam))**(-1.d0))**(-0.3d0)
          epair=epair*0.5d0*(CV2+CA2)*(1.d0+(CV2-CA2)/
     &    (CV2+CA2)*qpair) 
          
          

          epsneu=(eplas+ephot+ebrems+epair)

          E = E - epsneu
       endif

      if(xh.gt.1.d-37) then
c......1. Calculate screening factors fpp (weak&intermediate screening)
c.........a. Weak pp screening:

             kb=.5d0
             etahb=1.127d0
             b=1.d0
             zetab=2.d0
             mu=etahb**2/zeta
             Lambda0=1.88d+8*sqrt(Rho/(mu*T**3))
             H12=kb*etahb*zetab*(Lambda0**b)
             fppw=exp(H12)

c.........b. Weak he3he3 screening

             etahb=2.d0
             zetab=8.d0
             mu = etahb**2/zeta
             Lambda0=1.88d8*sqrt(Rho/(mu*T**3))
             h33= kb*etahb*zetab*(Lambda0**b)
             fhe3he3=exp(h33)
c             fhe3he3 = 1.1
c            if(H12.gt.(.1)) then
c..........c Intermediate pp screening:

                 b=0.860d0
                 kb=0.380d0

                 sumni=Z0
                 do i=1,14
                     n(i)=(xba(i)*Rho)/(HMASS)
                     sumni=sumni+n(i)
                 end do
                 avz3b=Z0
                 zbar=Z0
                 do i=1,14
                     avz3b=avz3b+Xz(i)**(3.d0*b-Z1)*n(i)
                     zbar=zbar+Xz(i)*n(i)
                 end do
                 avz3b=avz3b/sumni
                 zbar=zbar/sumni
                 etahb=avz3b/(1.127d0**(3.d0*b-2.d0)*zbar**(2.-2.*b))

                 zetab=1.63d0
                 H12=kb*etahb*zetab*(Lambda0**b)
                 fppi=exp(H12)
C
c            end if
          fpp = min(fppi,fppw)
C
c......2. Calculate effective cross section factors Spp and Spd, She3he3

       Ap=Z1
       Ad=2.d0
       Ahe3=3.d0
C
       N0=6.0225d+23
C
C       Np=Rho*(xh/Ap)*N0
C       nd=Rho*(xd/Ad)*N0

       App=.5d0
       Apd=(Ap*Ad)/(Ap+Ad)
       Ahe3he3 = Ahe3/2.d0
       S0pp=4.07d-22
       S0pd=2.5d-04
       S0he3he3 = 5.15d3

       dS0pp=4.52d-24
       dS0pd=7.9d-06
       dS0he3he3 = -0.9d0

       T6=T/1.d6
       taupp=42.487d0*(App/T6)**Z13
       taupd=42.487d0*(Apd/T6)**Z13
       tauhe3he3 =42.487d0*(16.d0*Ahe3he3/T6)**Z13

       E0pp=taupp*k*T*Z13
       E0pd=1.2204d0*(Apd*T6**2)**Z13
       E0he3he3 = tauhe3he3*k*T*Z13

       Spp=S0pp*(Z1+(5.d0/(12.d0*taupp)))+dS0pp*(E0pp+.97222d0*k*T)
       Spd=S0pd*(Z1+(5.d0/(12.d0*taupd)))+dS0pd*(E0pd+.97222d0*k*T)
       She3he3 = S0he3he3*(Z1+ (5.d0/(12.d0*tauhe3he3)))
     +   +dS0he3he3*(E0he3he3 + .97222d0*k*T)

c......Get the average products of cross section times velocity:

       sigvpp=1.3005d-15*(Z1/(App*T6**2))**Z13
       sigvpp=sigvpp*fpp*Spp*exp(-taupp)

       sigvpd=1.3005d-15*(Z1/(Apd*T6**2))**Z13
c      Xpd = Z0
       sigvpd=sigvpd*fpd*Spd*exp(-taupd)

       sigvhe3he3=1.3005d-15*(Zhe3**2/(Ahe3he3*T6**2))**Z13
       sigvhe3he3=sigvhe3he3*fhe3he3*She3he3*exp(-tauhe3he3)
c......Get the nuclear generation rates for the three reactions

       Rpp = np*np*sigvpp/2.d0
       Rhe3he3 = nhe3*nhe3*sigvhe3he3/2.d0
c......Correct nd for equilibrium deuterium burning
       Xpd=np*sigvpd
       if(Xpd*dtime.gt.100.d0) then
         nde=Z0
       else
         nde=nd*exp(-Xpd*dtime)
       endif
       if(Xpd*dtime.gt.0.01d0) then
         Rpd = (nd-nde)/dtime
       else
         Rpd = nde*Xpd
       endif

c......Get the Q values (adjusted for neutrino losses)

       Qpp = 2.310d-6 - 4.245d-7
       Qpd = 8.801d-6
       Qhe3he3 = 2.06d-5
c......Add on extra deuterium burning Q to the Qpp

       Qpp = Qpp + Qpd
c.......Add on he3 contribution to Qpp to form Qpp1
       Qpp1=Qhe3he3/2.d0 + Qpp

c......Get the energy generation rate

       if(iflip.eq.0) then
         E = E + Qpp*Rpp/Rho
         E = E + Qhe3he3*Rhe3he3/Rho
c......Add on contribution from primordial deuterium (Note that is is
c......implicitly assumed that deuterium is used up by the time Rpp is
c......significant.)
         E = E + Qpd*Rpd/Rho
         if(iwrite.gt.0)                                                &
     &     write(6,779) 'j2,nde,nd,Xpd,Rpd,E:',iwrite,nde,nd,Xpd,Rpd,E
 779       format(a,i5,1p,7E10.2)
       else

c....    He3 in statistical equilibrium.

c...     PPII and PPIII corrections due to Parker, Bahcall and Fowler, 1964.
C...     (Fit to curves for Y=0.1, 0.5, 0.9)
         xpp23=Z1+max(Z0,(T6-10.d0)/36.d0)
         xpp23=min(1.5d0,xpp23) + max(Z0,XHe4*(0.8d0-((T6-18.)/11.)**2))
         xpp23=min(1.95d0,xpp23)
c
         E = E + xpp23*Qpp1*Rpp/rho
c...     Compute CNO generation rate (from Bahcall 1989)
         Apc = 12./(12.+Z1)
         Apn = 14./(14.+Z1)
         Apo = 16./(16.+Z1)
         sigvpc=1.3005d-15*(6.d0/(Apc*T6**2))**Z13
         sigvpn=1.3005d-15*(7.d0/(Apn*T6**2))**Z13
         sigvpo=1.3005d-15*(8.d0/(Apo*T6**2))**Z13
         fpc=Z1
         fpn=Z1
         fpo=Z1
         taupc=136.93d0/T6**Z13
         taupn=152.31d0/T6**Z13
         taupo=166.96d0/T6**Z13
         S0pc=1.45
         S0pn=3.32
         S0po=9.4
         Spc=S0pc*(Z1+(5.d0/(12.d0*taupc)))
         Spn=S0pn*(Z1+(5.d0/(12.d0*taupn)))
         Spo=S0po*(Z1+(5.d0/(12.d0*taupo)))
         sigvpc=sigvpc*fpc*Spc*exp(-taupc)
         sigvpn=sigvpn*fpn*Spn*exp(-taupn)
         sigvpo=sigvpo*fpo*Spo*exp(-taupo)
         nc=xc*RHO/(HMASS*12.)
         nn=xn*RHO/(HMASS*14.)
         no=xo*RHO/(HMASS*16.)
         Qpc=1.7635d-5
         Qpn=2.2466d-5
         Qpo=5.6961d-6
         Xpc=np*sigvpc
         Xpn=np*sigvpn
         Xpo=np*sigvpo
         Rpc=nc*Xpc
         Rpn=nn*Xpn
         Rpo=no*Xpo
         Rcno=Rpc+Rpn+Rpo
         Ecno=(Qpc*Rpc+Qpn*Rpn+Qpo*Rpo)/RHO
c.  March 27, 2012:
c. Debugging effort, step 2.
c. UNCOMMENT the line below when you're
c. done testing out the effect of this modification
         E = E + Ecno

       end if
      endif
C       write(6,1414)E 
C       write(6,1414)T,Rhe3he3,fppw,fppi,nP,rho,rpP
1414   format(1x,1p,8E11.3)
       if(ifirst.eq.1 .and. T.gt.2.e7) then
         ifirst=2
         write(6,'(" xc,xn,xo=",1p,3E12.4," Ecno=",E12.4," E=",E12.4)')
     &     xc,xn,xo,Ecno,E
        write(6,1414) T,Rhe3he3,fppw,fppi,nP,rho,Rpp
       endif
c.  March 27, 2012:
c.  And also April 3, 2012
c. Debugging effort, step 1.
c. (comment out the line below after you're done
c.  testing the effect of this change.)
c.       E = 0.0
c.       E = 0.98 * E

c       write(6,*) "In nucrat, E = ",E
       return

       end

      subroutine compchange
      include 'parm.h'
      include 'var.h'
      dimension pp(5),AA(3,4)
c     dimension AA(3,4)
      COMMON/ABUND/XCA(14),H1(14),A(14)
      DATA Z0,Z1,Z13,HMASS/0.d0,1.d0,0.333333333333333333d0,1.6732D-24/
      data ifirst/0/
      save
c
C  Reset all abundances to original values when the deuterium
C  abundance is given as a negative value.
      if(XH2.lt.Z0) then
        H1(2)=abs(XH2)
        SUM = 1.
        DO J = 1,14
          SUM = SUM - H1(J)
        ENDDO
C H1(2) is temporarily used as the primordial deuterium abundance and
C then replaced by the initial helium abundance.
        H1(2)=SUM
        do i=1,N
          helium3(i) = 1.d-5
          hydrogen(i) = H1(1)+Z13*helium3(i)
          helium4(i) = H1(2)-4.d0/3.d0*helium3(i)
          deuterium(i)=abs(XH2)
          carbon(i)=H1(3)
          trogen(i)=H1(12)
          oxygen(i)=H1(13)
        enddo
        return
      endif
c
C  The following redistributes any inaccuracy of abundances to the
C  most abundant species.
c
      if(ifirst.eq.0) then
        ifirst=1
        XYZx = hydrogen(N) + helium4(N) + helium3(N) + deuterium(N)     &
     &       + oxygen(N) + trogen(N) + carbon(N)
        do i=1,N
          pp(1)=hydrogen(i)
          pp(2)=helium4(i)
          pp(3)=carbon(i)
          pp(4)=trogen(i)
          pp(5)=oxygen(i)
          XYZ= hydrogen(i) + helium4(i) + helium3(i) + deuterium(i)     &
     &       + oxygen(i) + trogen(i) + carbon(i)
          k=1
          do j=2,5
            if(pp(j).gt.pp(k)) k=j
          enddo
          pp(k)=pp(k) + XYZx-XYZ
          hydrogen(i) = pp(1)
          helium4(i)  = pp(2)
          carbon(i)   = pp(3)
          trogen(i)   = pp(4)
          oxygen(i)   = pp(5)
        enddo
      endif
      YYmax = hydrogen(N) + helium4(N) + helium3(N) + deuterium(N)
      do i=1,N 
        Tnuc = X(i,4)     
        olddeuterium=deuterium(i)
        oldhydrogen= hydrogen(i)
        oldhelium3 = helium3(i)
        oldhelium4 = helium4(i)
        oldcarbon  = carbon(i)
        oldtrogen  = trogen(i)
        oldoxygen  = oxygen(i)
        call flip(Tnuc,hydrogen(i),helium4(i),helium3(i),iflip)
        call invstate(41,x(i,1),x(i,4),hydrogen(i),helium4(i),Crad,RHO, &
     &               cP,alpha,beta,delta)
        QNh1 = hydrogen(i)*rho/HMASS
        QNh2 =deuterium(i)*rho/2.d0/HMASS
        QNhe3 = helium3(i)*rho/3.d0/HMASS

        call nucrat(Rho,Tnuc,QNh1,QNh2,QNhe3,dtime,ee,                  &
     &        Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                   &
     &        helium4(i),carbon(i),trogen(i),oxygen(i),                 &
     &        xpp23,iflip,0)

        if(Tnuc.gt.5.d7) then
          helium4(i)=helium4(i)+deuterium(i)+helium3(i)
          deuterium(i)=Z0
          helium3(i)=Z0
        else
          if(Xpd*dtime.gt.100.d0 .or. deuterium(i).lt.1.d-40) then
            deuterium(i)=Z0
          else
            deuterium(i)=deuterium(i)*exp(-Xpd*dtime)
          endif
          dxddt = (deuterium(i) - olddeuterium)/dtime
        endif
        dxhe4dt= -(3.d0*R3a+R12a)*4.d0*HMASS/RHO
        dxcdt  = (R3a-R12a)*12.d0*HMASS/RHO
        dxodt  = R12a*16.d0*HMASS/RHO
        if(iflip.eq.0) then
          dxhe3dt =(Rpp-2.d0*Rhe3)*3.d0*HMASS/rho-1.5d0*dxddt
          helium3(i) = dtime*dxhe3dt + helium3(i)
          if(helium3(i).lt.Z0) then
            helium3(i)=Z0
            dxhdt = (-3.d0*Rpp)*HMASS/rho
            dxhe4dt= dxhe4dt - dxhdt 
            helium4(i) = helium4(i) + oldhelium3
            dxh=min(oldhelium*Z13,hydrogen(i))
            hydrogen(i)= hydrogen(i) - dxh
            helium4 (i)= helium4 (i) + dxh
          else
            dxhdt  = (-3.d0*Rpp+2.d0*Rhe3)*HMASS/rho+.5d0*dxddt
            dxhe4dt= dxhe4dt + Rhe3 *4.d0*HMASS/rho
          endif
        else
          brem=Z1/dtime
          AA(1,1)=brem+Xpc
          AA(1,2)=-(Z1-4.d-4)*Xpn
          AA(1,3)=Z0
          AA(1,4)=vcarbon(i)*brem
          AA(2,1)=-Xpc
          AA(2,2)=brem+Xpn
          AA(2,3)=-Xpo
          AA(2,4)=vtrogen(i)*brem
          AA(3,1)=Z0
          AA(3,2)=-4.d-4*Xpn
          AA(3,3)=brem+Xpo
          AA(3,4)=voxygen(i)*brem
          CALL GIRL(AA,3,1)
          carbon(i)=AA(1,4)
          trogen(i)=AA(2,4)
          oxygen(i)=AA(3,4)
          dxhdt  = -(2.d0*xpp23*Rpp+4.d0*Rcno)*HMASS/rho
c         if(i.eq.1) then
c           write(6,*) xpp23,Rpp,HMASS,rho,Rcno,-dxhdt
c           stop
c         endif
          dxhe4dt= dxhe4dt - dxhdt 
        end if

        hydrogen(i)= dtime*dxhdt   + hydrogen(i)
        helium4(i) = dtime*dxhe4dt + helium4(i)
        if(oldhydrogen.gt.Z0 .and. hydrogen(i) .lt. 1.d-37) then
          hydrogen(i) = Z0
          helium3(i) = Z0
          helium4(i) = YYmax
        endif
        if(oldhelium4.gt.Z0 .and. helium4(i).lt.1.d-37) then
          helium4(i)= Z0
          carbon(i) = carbon(i) + oldhelium4
        else
          carbon(i) = carbon(i) + dxcdt*Dtime
        endif
        if(oldcarbon.gt.Z0 .and. carbon(i).lt.Z0) then
          carbon(i)= Z0
          oxygen(i) = oxygen(i) + carbon(i)-dxcdt*Dtime
        else
          oxygen(i) = oxygen(i) + dxodt*Dtime
        endif
c       if(i.eq.1) then
c         Znet=1.d0-oldhydrogen-oldhelium4-oldhelium3-olddeuterium      &
c    &      -oldcarbon-oldnigrogen-oldoxygen
c         write (6,200) i,oldhydrogen,oldhelium4,oldhelium3,            &
c    &        olddeuterium,oldcarbon,oldnigrogen,oldoxygen,Znet,iflip
c         Znet=1.d0-hydrogen(1)-helium4(1)-helium3(1)-deuterium(1)      &
c    &      -carbon(1)-trogen(1)-oxygen(1)
c         write (6,200) i,hydrogen(1),helium4(1),helium3(1),            &
c    &        deuterium(1),carbon(1),trogen(1),oxygen(1),Znet
c200      format(i4,' H,He4,He3,D,C,N:=',8F8.5,i2)
c       endif
      end do 
 
C..... MixinG

      iset=0
      ic1=0
      ic2=0
      do j=1,N
        if(iconvect(j) .eq. 0) then
          if(iset.eq.0) then
            ic1=j
            iset=1
          else
            ic2=j
          end if
        else
          if(iset.eq.1) then
            ic2=j-1
            iset=2
          end if
        end if
        if(iset.eq.2 .or. (iset.eq.1 .and. j.eq.N)) then
          if(ic1.lt.ic2) then
            runsum = Z0
            runsum1= Z0
            runsum2= Z0
            runsum3= Z0
            runsum4= Z0
            runsum5= Z0
            runsum6= Z0
            runsum7= Z0
            if(ic2 .eq. N .and. Zflux.gt.0.) then
              runsum  = Zflux*dtime
              runsum1 = H1(1)*runsum
              runsum2 = XH2 *runsum
              runsum3 = 1.d-5*runsum
              runsum4 = H1(2)*runsum
              runsum5 = H1(3)*runsum
              runsum6 = H1(12)*runsum
              runsum7 = H1(13)*runsum
            endif
            do i = ic1,ic2
              runsum1 = runsum1 + hydrogen(i)*dM(I) 
              runsum2 = runsum2 + deuterium(i)*dM(I) 
              runsum3 = runsum3 + helium3(i)*dM(I) 
              runsum4 = runsum4 + helium4(i)*dM(I)
              runsum5 = runsum5 + carbon(i)*dM(I)
              runsum6 = runsum6 + trogen(i)*dM(I)
              runsum7 = runsum7 + oxygen(i)*dM(I)
              runsum  = runsum  + dM(i)
            enddo
            Xh = runsum1/runsum
            xd = runsum2/runsum
            xhe3=runsum3/runsum
            xhe4=runsum4/runsum
            xc  =runsum5/runsum
            xn  =runsum6/runsum
            xo  =runsum7/runsum
C            write(6,*) xh,xd,xhe3,xhe4
            do i = ic1,ic2
              hydrogen(i) = xh
              deuterium(i)= xd
              helium3(i)  = xhe3
              helium4(i)  = xhe4
              carbon(i)   = xc
              trogen(i)   = xn
              oxygen(i)   = xo
            end do
          endif
          iset=0
          ic1=0
          ic2=0
        endif
      enddo
c     Znet=1.d0-hydrogen(1)-helium4(1)-helium3(1)-deuterium(1)          &
c    &      -carbon(1)-trogen(1)-oxygen(1)
c     write (6,200) hydrogen(1),helium4(1),helium3(1),deuterium(1),     &
c    &        carbon(1),trogen(1),oxygen(1),Znet
c
C     END MIXING 
      RETURN
      END 
      subroutine flip(T,hydrogen,helium4,helium3,iflip)
      include 'parm.h'
      data Z0/0.d0/
      save
c
      if(T.gt.1.0d+7 .or. helium3.eq.Z0) then
        iflip=1
        helium4=helium4+helium3
        dX=min(helium3/3.d0,hydrogen)
        hydrogen=hydrogen-dX
        helium4 =helium4 +dX
        helium3=Z0
      else
        iflip=0
      end if
      return
      end
      subroutine invstate(im,Pr,temp,X,Y,Crad,rho,cP,alpha,Rbeta1,delta)
C+
C  The subroutine invstate determines the density as a function of
C  pressure and temperature.  Other thermodynamic quantities are
C  also calculated.  This current version includes the effects of
C  degeneracy, including partial relativistic degeneracy; it does
C  not include partial ionization of the degenerate gas.
C
C  Author:  Harold W. Yorke      10-SEP-02               (JPL / CALTECH)
C-
      include 'parm.h'
      real*4 Plog,Tlog,XX,YY,Ymax,EOS(5)
      real*8 XCA,H1,A
      COMMON/ABUND/XCA(14),H1(14),A(14)
      data iorder/5/
      DATA Arad3/2.52197145D-15/
      DATA negative/0/
      DATA ifirst/0/
      save
C
      if(ifirst.eq.0) then
        ifirst=1
        Ymax=H1(1)+H1(2)
      endif
      Prad=Crad*Arad3*Temp**4
      Pgas=Pr-Prad
      if(Pgas.lt.0.d0) then
        write(6,'(a,i3,a)') 'INVSTATE',im,': Gas pressure is negative'
        write(6,'(i4,1x,a,1p,4e12.4)')                                  &
     &     negative,'Pr,temp,Prad,Pgas:',Pr,temp,Prad,Pgas
        negative=negative+1
        if(negative.gt.10) stop 'invstate'
        Pgas = 0.01*Pr
        Prad = 0.99*Pr
      endif
      Rbeta=Pgas/Pr
      Rbeta1=Prad/Pr
      i1=0
c
      Plog=log10(Pgas)
      Tlog=log10(Temp)
      XX=X
      YY=Y
      call eostab(Plog,Tlog,XX,YY,Ymax,EOS,iorder,ierror)
      if(ierror.ne.0) then
c       write(6,201) Plog,Tlog,EOS,iorder,ierror
 201    format(7f10.5,2i2)
c       stop 'invstate'
      endif
      rho     = 10.**eos(2)
      alphgas = eos(3)
      deltgas = eos(4)
      cVgas   = 10.**(eos(5)+6.)
      cPgas   = cVgas+deltgas**2*Pgas/(alphgas*temp*rho)
c     atg     = deltgas*Pgas/(rho*temp*cPgas)
c
c CHANGE FINAL RESULTS TO ACCOUNT FOR RADIATION PRESSURE
c
      alpha = alphgas/Rbeta
      delta = deltgas + 4.d0*alpha*Rbeta1
      cV    = cVgas + 12.d0*Prad/(temp*rho)
      cP    = cV + delta**2*Pr/(alpha*temp*rho)
c
      return
      end
      subroutine eostab(Plog,Tlog,X,Y,Ymax,EOS,iorder,ierror)
      implicit real*4(a-h,o-z) , integer (i-n)
      parameter (NZ=5,NX=5,MT=139,NTX=(NX+NZ)*MT,MPT=187165)
      dimension IPTR(NTX+1),EOSA(MPT,5)
      dimension EOS(5)
      dimension Pmin(4),DP1(4),DP2(4),XX(4),IP(4),IP1(4),IP2(4),NP(4)
      data ifirst/0/
      save

      ierror=0
      if(ifirst.eq.0) then
c         write(6,*) 'Reading in eostable for the first time'
        ifirst=1
        open(7,file='eospointer',form='formatted',status='old')
        read(7,*) IPTR,kx,kt,tmin,tmax,dP
        close(7)
        open(7,file='eostable',form='formatted',status='old')
c        write(6,*) 'Done reading in the eostable for the first time'
c
        if(kx.ne.NX .or. kt.ne.MT) then
          write(6,*) 'Mismatch of dimensions in eostab and eostable',
     &     ' NX=',NX,kx,' MT=',MT,kt
          stop 'EOSTAB'
        endif

        dx  = 1./float(NX)
        dz  = 1./float(NZ)
        Tlogmin=log10(Tmin)
        Tlogmax=log10(Tmax)
        dt  = (Tlogmax-Tlogmin)/float(MT-1)

        ITX = 0
        do ix=1,NX+NZ
           do it=1,MT
              ITX=ITX+1
              read (7,201) (EOSA(IPTR(ITX),i),i=1,5)
 201          format(4x,1p,E12.4,7x,f9.6,5x,f4.1,3x,f9.6,3x,f9.6)
              IPmin=IPTR(ITX)+1
              IPmax=IPTR(ITX+1)-1
              do j=IPmin,IPmax
                 read(7,203) (EOSA(j,i),i=1,5)
 203             format(f10.6,f11.6,3f9.6)
              enddo
           enddo
        enddo
      close(7)
      endif

      if(X.lt.1.E-37 .and. Ymax-Y.gt.1.e-4) then
        ZZ=1.-Y
        DX2=ZZ/dz
        IX1=max(DX2,0.)+4.
        IX2=min(IX1+2,NX+NZ)
        IX1=IX2-1
        if(IX1.gt.NX) then
          DX2=DX2-float(IX1-5)
        else
          IX1=1
          DX2=max(0.,(ZZ-1.+Ymax)/(dz-1.+Ymax))
        endif
        DX1=1.-DX2
      else
        if(X.lt.0.0 .or. X.gt.0.8) ierror=1
        DX2=X/dx
        IX1=max(DX2,0.)
        IX2=min(IX1+2,NX)
        IX1=IX2-1
        DX2=DX2-float(IX1-1)
        DX1=1.-DX2
      endif
c     write(6,202) 'X:',ix1,ix2,dx1,dx2
 202  format(a,2i10,1p,2E12.4)

      if(Tlog.lt.Tlogmin .or. Tlog.gt.Tlogmax) ierror=ierror+2
      DT2=(Tlog-Tlogmin)/dt
      IT1=max(DT2,0.)
      IT2=min(IT1+2,MT)
      IT1=IT2-1
      DT2=DT2-float(IT1-1)
      DT1=1.-DT2
c     write(6,202) 'T:',it1,it2,dt1,dt2

      IP(1)=(IX1-1)*MT+IT1
      IP(2)=(IX2-1)*MT+IT1
      IP(3)=(IX1-1)*MT+IT2
      IP(4)=(IX2-1)*MT+IT2
c     write(6,*) (IPTR(IP(j)),j=1,4)
      irror=0
      do j=1,4
        Pmin(j) = EOSA(IPTR(IP(j))+1,1)
        NP(j)   = IPTR(IP(j)+1)-IPTR(IP(j))-1
        DP2(j)  = (Plog-Pmin(j))/DP
        IP1(j)  = max(DP2(j),0.)
        IP2(j)  = min(IP1(j)+2,NP(j))
        IP1(j)  = IP2(j)-1
        DP2(j)  = DP2(j)-float(IP1(j)-1)
        DP1(j)  = 1.-DP2(j)
        if(DP1(j).lt.0.0 .or. DP2(j).lt.0.0) irror=4
      enddo
      ierror=ierror+irror
c     do j=1,4
c       write(6,202) 'P:',iP1(j),iP2(j),dP1(j),dP2(j)
c     enddo

      do i=1,iorder
        do j=1,4
        XX(j)=DP1(j)*EOSA(IPTR(IP(j))+IP1(j),i) + 
     &        DP2(j)*EOSA(IPTR(IP(j))+IP2(j),i)
        enddo
        if(i.le.2) then
          DS1=DT1
          DS2=DT2
        else
          DS1=min(1.,DT1)
          DS2=1.-DS1
        endif
        EOS(i)=DX1*(DS1*XX(1)+DS2*XX(3)) + DX2*(DS1*XX(2)+DS2*XX(4))

      enddo

      return
      end


      SUBROUTINE HENYEY(MODE)
      include 'parm.h'
      include 'var.h'
      logical LCONV
      DIMENSION D(MH,MH+1),CA(MH),DX(MH)
      DIMENSION corr(16)
      EQUIVALENCE (HA(1,MH+1),D(1,1)) , (D(1,MH+1),CA(1))
      COMMON/UsefulInfo/J2,Ntot,Niter
      COMMON/TimeInfo/TimeVal
      save
C
c      write(6,*) "STARTING HENYEY SUBROUTINE"
      TimeVal = TIME

      Ntot = N
      if(abs(ITMIN).gt.10 .or. mode.eq.999) then
        cf=1./max(10,ITMIN)
        cfmax=.6
        cfexp=1.2
        if(ITMIN.gt.10) cfexp=1.1
      else
c       cf=.1
        cf=1.
        cfmax=1.0
        cfexp=1.25
      endif

      ID=N
      MODE=0
      IRIT=0
      JTMAX=abs(ITMAX)
      IF(ITMIN.GT.0 .and. ITMAX.gt.0) IRIT=JTMAX/2
      ITER=0
  300 CONTINUE
c      write(6,*) "HENYEY ITERATION ",iter
      if(iter.eq.ITMAX/2) cf=.1
      cf=min(cfmax,cf*cfexp)
      ITER=ITER+1
      Niter = iter
      IF(ITER.GT.JTMAX) then
         GOTO 998
      endif
        J2=1
        J3=2
        CALL GI(J2)
c        write(6,*) 'iter',iter,'Gs: ',J2,G
c        write(6,*) 'iter',iter,'Cs: ',J2,HC        
c        write(6,*) 'iter',iter,'Ds: ',J2,HD       
c        write(6,*) 'iter',iter,'Es: ',J2,HE
        IF(ITMIN.EQ.-ITER) then
          call gid(j2)
          WRITE(6,101) 'J2=',J2,' Henyey matrices'
          WRITE(6,888) 'G:',(G(I),I=1,NG)
          do j=1,ng
            WRITE(6,888) 'C:',(HC(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'D:',(HD(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'E:',(HE(I,J),I=1,NG)
          enddo
        endif
        do I=1,MH
          do J=1,MH
            HA(I,J)=HD(I,J)
            D(I,J)=-HE(I,J)
          enddo
          D(I,MH+1)=-G(I)
          DG(I)=0.
          JG(I)=0
        enddo
        CALL GIRL (HA,MH,MH+1)
        IW=1
        do I=1,MH
          do J=1,MH+1
            HX(I,J)=D(I,J)
            HW(I,J,IW)=D(I,J)
          enddo
        enddo
  302 CONTINUE
        J1=J2
        J2=J3
        J3=J3+1
        CALL GI(J2)
c        write(6,*) 'iter',iter,'Gs: ',J2,G
c        write(6,*) 'iter',iter,'Cs: ',J2,HC        
c        write(6,*) 'iter',iter,'Ds: ',J2,HD       
c        write(6,*) 'iter',iter,'Es: ',J2,HE

        IF(ITMIN.EQ.-ITER) then
          call gid(j2)
          WRITE(6,101) 'J2=',J2,' Henyey matrices'
          WRITE(6,888) 'G:',(G(I),I=1,NG)
          do j=1,ng
            WRITE(6,888) 'C:',(HC(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'D:',(HD(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'E:',(HE(I,J),I=1,NG)
          enddo
        endif
        do I=1,MH
          do J=1,MH
            FAC=HC(I,1)*HX(1,J)
            do K=2,MH
              FAC=FAC+HC(I,K)*HX(K,J)
            enddo
            HA(I,J)=HD(I,J)+FAC
            D(I,J)=-HE(I,J)
          enddo
          FAC=HC(I,1)*HX(1,MH+1)
          do K=2,MH
            FAC=FAC+HC(I,K)*HX(K,MH+1)
          enddo
          D(I,MH+1)=-G(I)-FAC
          if(ABS(DG(I)).lt.ABS(G(I))) then
            DG(I)=G(I)
            JG(I)=J2
          endif
        enddo
        CALL GIRL (HA,MH,MH+1)
        IW=IW+1
        do J=1,MH+1
          do I=1,MH
            HX(I,J)=D(I,J)
            HW(I,J,IW)=D(I,J)
          enddo
        enddo
      IF(J2.GE.N-1) GOTO 303
      IF(IW.LT.ID) GOTO 302
        IW=0
      GOTO 302
  303 IF(ITER.GT.IRIT .and. mod(MODEL,JRIT).eq.0)
     *  WRITE(6,101) 'ITER:',ITER,' LARGEST GI:',(JG(I),DG(I),I=1,MH)
        J1=J2
        J2=J3
        J3=J2+1
        CALL GI(J2)
c        write(6,*) 'iter',iter,'Gs: ',J2,G
c        write(6,*) 'iter',iter,'Cs: ',J2,HC        
c        write(6,*) 'iter',iter,'Ds: ',J2,HD       
c        write(6,*) 'iter',iter,'Es: ',J2,HE

        IF(ITMIN.EQ.-ITER) then
          call gid(j2)
          WRITE(6,101) 'J2=',J2,' Henyey matrices'
          WRITE(6,888) 'G:',(G(I),I=1,NG)
          do j=1,ng
            WRITE(6,888) 'C:',(HC(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'D:',(HD(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'E:',(HE(I,J),I=1,NG)
          enddo
        endif
        do I=1,MH
          do J=1,MH
            FAC=HC(I,1)*HX(1,J)
            do K=2,MH
              FAC=FAC+HC(I,K)*HX(K,J)
            enddo
            D(I,J)=HD(I,J)+FAC
          enddo
          FAC=HC(I,1)*HX(1,MH+1)
          do K=2,MH
            FAC=FAC+HC(I,K)*HX(K,MH+1)
          enddo
          CA(I)=-G(I)-FAC
          G(I)=DG(I)
          DG(I)=0.
          JG(I)=0
        enddo
        CALL GIRL (D,MH,1)
  304 CONTINUE
c      write(6,*) 'cf = ',cf
c The following line IS A KLUDGE
c      cf = 1.0

      corr(1) = J2
      call invstate(21,x(j2,1),x(j2,4),hydrogen(j2),helium4(j2),Crad,   &
     &                RHO2,cP,alpha,beta,delta)
      corr(16) = RHO2

      do I=1,MH
         E(I)=CA(I)
         AFAC=E(I)*cf
         
         corr(I+1) = AFAC
         
         IF(SMAX(I).gt.0.) then
            cm39 = ABS(AFAC/X(j2,I))
            if (cm39 .gt. SMAX(I)) then
c               write(6,*) "RESCALING at iter",ITER,"corr",J2,"var",I
               AFAC = AFAC/ABS(AFAC) * SMAX(I)*X(J2,I)
            endif 
            X(J2,I)=X(J2,I)+AFAC   
            DX(I)=2.*AFAC/(X(J2,I)+VX(J2,I)+SMIN(I))
            IF(X(J2,I).lt.SMIN(I)) goto 997
         else
            AFAC=MIN(AFAC,-SMAX(I))
            AFAC=MAX(AFAC, SMAX(I))
            X(J2,I)=X(J2,I)+AFAC
            DX(I)=2.*AFAC/(ABS(X(J2,I))+ABS(VX(J2,I))+SMIN(I))
         endif
         
c         if (J2.eq.1) then
c            write(6,*)"VXS:", VX(J2,I),X(J2,I),AFAC
c         endif
         
         corr(I+5) = AFAC
         
         IF(ABS(DG(I)).lt.ABS(DX(I))) then
            DG(I)=DX(I)
            JG(I)=J2
         endif
      enddo

      corr(10) = X(J2,1)
      corr(11) = X(J2,2)
      corr(12) = X(J2,3)
      corr(13) = X(J2,4)
      corr(14) = dM(J2)
      corr(15) = zM(J2)

c      if (J2.eq.1) then
c         write(6,*) "DXS at J=",J2,":",DX
c         write(6,*) " "
c      endif
      
c        IF(ITMIN.LT.0) then
c          WRITE(6,101) 'J2=',J2,' CORRECTIONS:',(I,E(I),I=1,MH)
c        ENDIF

c.        write(6,100) "iter",ITER,"corr:",J2,(corr(ii),ii=2,16)
c. 100    format (a4,1x,i4.1,1x,a6,1x,i4.1,15G25.15)

      IF(J2.EQ.1) GOTO 305
        do I=1,MH
          CA(I)=HW(I,MH+1,IW)
          do J=1,MH
            CA(I)=CA(I)+HW(I,J,IW)*E(J)
          enddo
        enddo


        J2=J1
        J1=J1-1
        IW=IW-1

      IF(IW.GT.0 .OR. J2.EQ.1) GOTO 304
        IW=ID
      GOTO 304
  305 CONTINUE
        IF(ITER.GT.IRIT .and. mod(MODEL,JRIT).eq.0) then
          WRITE(6,101) 'ITERATION',ITER,' CORRECTIONS:',
     *    (JG(I),DG(I),I=1,MH)
          WRITE(6,101) 'ITERATION',ITER,'  NEW VALUES:',
     *    (JG(I),X(J2,i),i=1,MH)
        ENDIF
        LCONV=.FALSE.

c.        write(6,*) "On iter",iter,"cf = ",cf,"cfmax = ",cfmax

        do I=1,MH
          LCONV=LCONV .or. ABS(DG(I)) .gt. EPS(I)
c.          write(6,*) "On iter ",iter,"dXmax = ",ABS(DG(I))," EPS(I) = ",&
c.     &  EPS(I)
        enddo
        

      IF(ITER.LT.ITMIN .OR. LCONV .OR. abs(cf-cfmax).gt.0.01d0) GOTO 300
        do I=1,MH
          DG(I)=0.
          JG(I)=0
        enddo
        do J=2,N
          do I=1,MH
            E(I)=2.*(X(J,I)-VX(J,I))/(ABS(X(J,I))+ABS(VX(J,I))+SMIN(I))
            if(ABS(DG(I)).lt.ABS(E(I))) then
              DG(I)=E(I)
              JG(I)=J
            endif
          enddo
        enddo
        if(mod(MODEL,JRIT).eq.0) then
          WRITE(6,101) 'MODEL ',MODEL,' CHANGES:',(JG(I),DG(I),I=1,MH)
          WRITE(6,101) 'ITER=',ITER,' GI:',(I,G(I),I=1,MH)
        endif
      GOTO 999
 997    MODE=I
        WRITE(6,103) J2,I,X(J2,I)
      GOTO 999
 998    MODE=999
        WRITE(6,101) '--- CONVERGENCE NOT POSSIBLE: ',MODEL,' ',
     *               (JG(I),DG(I),I=1,MH),ITER
 999  CONTINUE
      CHANGE=0.
      DO I=1,MH
        CHANGE=CHANGE+ABS(DG(I))
      ENDDO
      RETURN
  101 FORMAT(1X,A,i5,A,1P,7(I5,E10.2),/,7x,8(I5,E10.2))
  103 FORMAT(/6X,'VARIABLE BELOW LIMIT:      X(',I4,',',I1,') = ',
     1       1P,E12.5)
  888 FORMAT(8x,a,1P,1x,10E10.2,10(/1X,10E10.2))
      END

      SUBROUTINE GIRL(A,N,M)
      include 'parm.h'
      DIMENSION A(1)
      SAVE
      NPM = N+M
      DO J = 1,N
        NJ = (J-1)*N
        JJ = NJ + J
        J1 = J + 1
        AMAX =ABS (A(JJ))
        JM = J
        IF( J1 .le. N ) THEN 
          DO I = J1,N
            IJ = NJ + I
            IF(ABS (A(IJ)) .gt. AMAX ) THEN
              AMAX =ABS (A(IJ))
              JM = I
            ENDIF
          ENDDO 
          IF( JM .ne. J ) THEN
            I1 = JM + NJ
            I2 = JJ
            DO I = J,NPM
              ZWI = A(I1)
              A(I1) = A(I2)
              A(I2) = ZWI
              I1 = I1 + N
              I2 = I2 + N
            ENDDO
          ENDIF
        ENDIF
        IF( A(JJ) .eq. 0.) GOTO 40
        DO I = 1,N
          IF( I .ne. J ) THEN
            IJ = NJ + I
            IK = NJ + I
            JK = JJ
            FAKTOR = - A(IJ)/A(JJ)
            DO K = J1,NPM
              JK = JK + N
              IK = IK + N
              A(IK) = A(IK) + FAKTOR * A(JK)
            ENDDO
          ENDIF
        ENDDO
        JK = JJ
        FAKTOR = 1./A(JJ)
        DO K = J1,NPM
          JK = JK + N
          A(JK) = A(JK) * FAKTOR
        ENDDO
      ENDDO
      RETURN
   40 CONTINUE
      WRITE(6,100)
  100 FORMAT(' ERROR EXIT IN SUBROUTINE GIRL')
      STOP 'girl'
       END

      subroutine GI(K2)
c
      include 'parm.h'
      include 'var.h'
      DATA GRAV/6.6704D-8/
      DATA PI  /3.14159265359d0/ 
c     DATA ELOG/ .434294d0/
      DATA XINC/.0005d0/
c      DATA XINC/.001d0/
      DATA HMASS/1.6732d-24/
      COMMON/CURRENTVALS/currentMass

      save
C
      J2=abs(K2)
      J1=J2-1
      J3=J2+1

      currentMass = zM(J2)/zM(N)
c     if(vx(j2,4).gt.5.e6) then
c       I1=-J1
c       I2=-J2
c       I3=-J3
c     else
c       I1=J1
c       I2=J2
c       I3=J3
c     endif
      if(K2.gt.0) then
        do j=1,NG
        do k=1,NG
          HC(j,k)=0.d0
          HD(j,k)=0.d0
          HE(j,k)=0.d0
        enddo
        enddo
      endif
      if(dTIME.le.0.d0) then
        brem=0.d0
      else
        brem=Cwrk/dTIME
      endif
      if (J2.eq.1) then
        zMhf2=zM(J2)-0.5*dM(J2)
        zMhf3=zM(J2)+0.5*dM(J3)
        GMdM = GRAV*zM(J2)*(dM(J2)+dM(J3))/(8.d0*PI)
        call invstate(20,x(j3,1),x(j3,4),hydrogen(j3),helium4(j3),Crad, &
     &                RHO3,cP3,alpha3,beta3,delta3)
        call invstate(21,x(j2,1),x(j2,4),hydrogen(j2),helium4(j2),Crad, &
     &                RHO2,cP,alpha,beta,delta)
        call flip(x(j2,4),hydrogen(j2),helium4(j2),helium3(j2),iflip)
        QNh1 = hydrogen(j2)*RHO2/HMASS
        QNh2 =deuterium(j2)*RHO2/2.d0/HMASS
        QNhe3 = helium3(j2)*RHO2/3.d0/HMASS
        CALL NUCRAT(RHO2,X(J2,4),QNh1,QNh2,QNhe3,DTIME,EP,              &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
        call nabla(x(j3,1),X(J3,2),X(J3,3),X(J3,4),zMhf3,               &
     &     TPNAB3,ADNAB,RADNAB,J3,IC)
        iconvect(j3) = ic 
        call nabla(X(J2,1),X(J2,2),X(J2,3),X(J2,4),zMhf2,               &
     &     TPNAB2,ADNAB,RADNAB,J2,IC)
        iconvect(j2) = ic 
        if(K2.gt.0) then
          PN=X(J2,1)*(1.d0+XINC)
          PM=X(J3,1)*(1.d0+XINC)
          call invstate(23,PN,x(j2,4),hydrogen(j2),helium4(j2),Crad,    &
     &                 RHP2,cPP,alphaP,beta,deltaP)
          call nabla(PN,X(J2,2),X(J2,3),X(J2,4),zMhf2,                  &
     &     TPNAB2P,ADNAB,RADNAB,J2,IC)
          call nabla(PM,X(J3,2),X(J3,3),X(J3,4),zMhf3,                  &
     &     TPNAB3P,ADNAB,RADNAB,J3,IC)
          RN=X(J2,2)*(1.d0+XINC)
          RM=X(J3,2)*(1.d0+XINC)
          call nabla(X(J2,1),RN,X(J2,3),X(J2,4),zMhf2,                  &
     &     TPNAB2R,ADNAB,RADNAB,J2,IC)
          call nabla(X(J3,1),RM,X(J3,3),X(J3,4),zMhf3,                  &
     &     TPNAB3R,ADNAB,RADNAB,J3,IC)
          XLN=X(J2,3)*(1.d0+XINC)
          XLM=X(J3,3)*(1.d0+XINC)
          call nabla(X(J2,1),X(J2,2),XLN,X(J2,4),zMhf2,                 &
     &     TPNAB2L,ADNAB,RADNAB,J2,IC)
          call nabla(X(J3,1),X(J3,2),XLM,X(J3,4),zMhf3,                 &
     &     TPNAB3L,ADNAB,RADNAB,J3,IC)
          TN=X(J2,4)*(1.d0-XINC)
          TM=X(J3,4)*(1.d0-XINC)
          call invstate(25,x(j2,1),TN,hydrogen(j2),helium4(j2),Crad,    &
     &                  RHT2,cPT,alphaT,beta,deltaT)
          CALL NUCRAT(RHO2,TN,QNh1,QNh2,QNhe3,DTIME,EPT,                &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
          call nabla(X(J2,1),X(J2,2),X(J2,3),TN,zMhf2,                  &
     &     TPNAB2T,ADNAB,RADNAB,J2,IC)
          call nabla(X(J3,1),X(J3,2),X(J3,3),TM,zMhf3,                  &
     &     TPNAB3T,ADNAB,RADNAB,J3,IC)
          RHD=RHO2*(1.d0+XINC)
          QNh1 = hydrogen(j2)*rhD /HMASS
          QNh2 =deuterium(j2)*rhD /2.d0/HMASS
          QNhe3 = helium3(j2)*rhD /3.d0/HMASS
          CALL NUCRAT(RHD,X(J2,4),QNh1,QNh2,QNhe3,DTIME,EPD,            &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
          dRHdP2=(RHP2-RHO2)/(XINC*X(J2,1))
          dRHdT2=(RHT2-RHO2)/(-XINC*X(J2,4))
          dEPdR=(EPD-EP)/(XINC*RHO2)
          dEPdT=(EPT-EP)/(-XINC*X(J2,4)) + dEPdR*dRHdT2
          dEPdP=dEPdR*dRHdP2
          dcPdP=(cPP-cP)/(XINC*X(J2,1))
          dcPdT=(cPT-cP)/(-XINC*X(J2,4))
          ddeldP=(deltaP-delta)/(XINC*X(J2,1))
          ddeldT=(deltaT-delta)/(-XINC*X(J2,4))
          dNAB2P=(TPNAB2P-TPNAB2)/(XINC*X(J2,1))
          dNAB2R=(TPNAB2R-TPNAB2)/(XINC*X(J2,2))
          dNAB2L=(TPNAB2L-TPNAB2)/(XINC*X(J2,3))
          dNAB2T=(TPNAB2T-TPNAB2)/(-XINC*X(J2,4))
          dNAB3P=(TPNAB3P-TPNAB3)/(XINC*X(J3,1))
          dNAB3R=(TPNAB3R-TPNAB3)/(XINC*X(J3,2))
          dNAB3L=(TPNAB3L-TPNAB3)/(XINC*X(J3,3))
          dNAB3T=(TPNAB3T-TPNAB3)/(-XINC*X(J3,4))
          HD(1,1) =-1.d0
          HD(1,2) =-4.d0*GMdM/X(J2,2)**5
          HE(1,1) = 1.d0
          HD(2,1) = dM(J2)*0.75d0/(PI*RHO2**2)*dRHdP2
          HD(2,2) = 3.d0*X(J2,2)**2
          HD(2,4) = dM(J2)*0.75d0/(PI*RHO2**2)*dRHdT2
          HD(3,1) =-dM(J2)*( dEPdP                                      &
     &             -dcPdP * (X(J2,4)-VX(J2,4))*brem                     &
     &             +delta/RHO2 *brem                                    &
     &             +ddeldP/RHO2 * (X(J2,1)-VX(J2,1))*brem               &
     &             -delta/RHO2**2 * dRHdP2 * (X(J2,1)-VX(J2,1))*brem )
          HD(3,3) = 1.d0
          HD(3,4) =-dM(J2)*( dEPdT                                      &
     &             -cP * brem                                           &
     &             -dcPdT * (X(J2,4)-VX(J2,4))*brem                     &
     &             +ddeldT/RHO2 * (X(J2,1)-VX(J2,1))*brem               &
     &             -delta/RHO2**2 * dRHdT2 * (X(J2,1)-VX(J2,1))*brem )
          HD(4,1) = 0.5d0*(TPNAB3+TPNAB2)                               &
     &                   - (X(J3,1)-X(J2,1))*0.5d0*dNAB2P
          HD(4,2) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB2R
          HD(4,3) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB2L
          HD(4,4) =-1.d0 - (X(J3,1)-X(J2,1))*0.5d0*dNAB2T
          HE(4,1) =-0.5d0*(TPNAB3+TPNAB2)                               &
     &                   - (X(J3,1)-X(J2,1))*0.5d0*dNAB3P
          HE(4,2) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB3R
          HE(4,3) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB3L
          HE(4,4) = 1.d0 - (X(J3,1)-X(J2,1))*0.5d0*dNAB3T
        endif
        G(1) = X(J3,1)-X(J2,1) + GMdM/X(J2,2)**4
        G(2) = X(J2,2)**3 - dM(J2)*0.75d0/(PI*RHO2)
        G(3) = X(J2,3) - dM(J2)*( EP - cP * (X(J2,4)-VX(J2,4))*brem     &
     &         + delta/RHO2 * (X(J2,1)-VX(J2,1))*brem )
        G(4) = (X(J3,4)-X(J2,4))                                        &
     &         -(X(J3,1)-X(J2,1))*0.5d0*(TPNAB3+TPNAB2)
c       write(6,*) j2, ' Gs ', g(1),g(2),g(3),g(4) 

        return
      endif 
      if(J2.gt.1 .and. J2.lt.N) then
        zMhf2=zM(J2)-0.5*dM(J2)
        zMhf3=zM(J2)+0.5*dM(J3)
        GMdM = GRAV*zM(J2)*(dM(J2)+dM(J3))/(8.d0*PI)
        if(K2.gt.0) then
          TPNAB2= TPNAB3
        else
          call nabla(X(J2,1),X(J2,2),X(J2,3),X(J2,4),zMhf2,             &
     &       TPNAB2,ADNAB,RADNAB,J2,IC)
        endif
        call invstate(26,x(j2,1),x(j2,4),hydrogen(j2),helium4(j2),      &
     &                Crad,RHO2,cP,alpha,beta,delta)
        call flip(x(j2,4),hydrogen(j2),helium4(j2),helium3(j2),iflip)
        call nabla(X(J3,1),X(J3,2),X(J3,3),X(J3,4),zMhf3,               &
     &     TPNAB3,ADNAB,RADNAB,J3,IC)
        iconvect(j3) = ic 
        QNh1 = hydrogen(j2)*RHO2/HMASS
        QNh2 =deuterium(j2)*RHO2/2.d0/HMASS
        QNhe3 = helium3(j2)*RHO2/3.d0/HMASS
        CALL NUCRAT(RHO2,X(J2,4),QNh1,QNh2,QNhe3,DTIME,EP,              &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
        if(K2.gt.0) then
          dNAB2P = dNAB3P
          dNAB2R = dNAB3R
          dNAB2L = dNAB3L
          dNAB2T = dNAB3T
          TPNAB2P=TPNAB3P
          TPNAB2R=TPNAB3R
          TPNAB2L=TPNAB3L
          TPNAB2T=TPNAB3T
          PN=X(J2,1)*(1.d0+XINC)
          PM=X(J3,1)*(1.d0+XINC)
          call invstate(27,PN,x(j2,4),hydrogen(j2),helium4(j2),Crad,    &
     &                  RHP2,cPP,alphaP,beta,deltaP)
          call nabla(PM,X(J3,2),X(J3,3),X(J3,4),zMhf3,                  &
     &     TPNAB3P,ADNAB,RADNAB,J3,IC)
          RM=X(J3,2)*(1.d0+XINC)
          call nabla(X(J3,1),RM,X(J3,3),X(J3,4),zMhf3,                  &
     &     TPNAB3R,ADNAB,RADNAB,J3,IC)
          XLM=X(J3,3)*(1.d0+XINC)
          call nabla(X(J3,1),X(J3,2),XLM,X(J3,4),zMhf3,                 &
     &     TPNAB3L,ADNAB,RADNAB,J3,IC)
          TN=X(J2,4)*(1.d0-XINC)
          TM=X(J3,4)*(1.d0-XINC)
          call invstate(28,x(j2,1),TN,hydrogen(j2),helium4(j2),Crad,    &
     &                  RHT2,cPT,alphaT,beta,deltaT)
          CALL NUCRAT(RHO2,TN,QNh1,QNh2,QNhe3,DTIME,EPT,                &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
          call nabla(X(J3,1),X(J3,2),X(J3,3),TM,zMhf3,                  &
     &     TPNAB3T,ADNAB,RADNAB,J3,IC)
          RHD=RHO2*(1.d0+XINC)
          QNh1 = hydrogen(j2)*rhD /HMASS
          QNh2 =deuterium(j2)*rhD /2.d0/HMASS
          QNhe3 = helium3(j2)*rhD /3.d0/HMASS
          CALL NUCRAT(RHD,X(J2,4),QNh1,QNh2,QNhe3,DTIME,EPD,            &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
          dRHdP2=(RHP2-RHO2)/(XINC*X(J2,1))
          dRHdT2=(RHT2-RHO2)/(-XINC*X(J2,4))
          dEPdR=(EPD-EP)/(XINC*RHO2)
          dEPdT=(EPT-EP)/(-XINC*X(J2,4)) + dEPdR*dRHdT2
          dEPdP=dEPdR*dRHdP2
          dcPdP=(cPP-cP)/(XINC*X(J2,1))
          dcPdT=(cPT-cP)/(-XINC*X(J2,4))
          ddeldP=(deltaP-delta)/(XINC*X(J2,1))
          ddeldT=(deltaT-delta)/(-XINC*X(J2,4))
          dNAB3P=(TPNAB3P-TPNAB3)/(XINC*X(J3,1))
          dNAB3R=(TPNAB3R-TPNAB3)/(XINC*X(J3,2))
          dNAB3L=(TPNAB3L-TPNAB3)/(XINC*X(J3,3))
          dNAB3T=(TPNAB3T-TPNAB3)/(-XINC*X(J3,4))
          HD(1,1) =-1.d0
          HD(1,2) =-4.d0*GMdM/X(J2,2)**5
          HE(1,1) = 1.d0
          HC(2,2) =-3.d0*X(J1,2)**2
          HD(2,1) = dM(J2)*0.75d0/(PI*RHO2**2)*dRHdP2
          HD(2,2) = 3.d0*X(J2,2)**2
          HD(2,4) = dM(J2)*0.75d0/(PI*RHO2**2)*dRHdT2
          HC(3,3) =-1.d0
          HD(3,1) =-dM(J2)*( dEPdP                                      &
     &             -dcPdP * (X(J2,4)-VX(J2,4))*brem                     &
     &             +delta/RHO2 *brem                                    &
     &             +ddeldP/RHO2 * (X(J2,1)-VX(J2,1))*brem               &
     &             -delta/RHO2**2 * dRHdP2 * (X(J2,1)-VX(J2,1))*brem )
          HD(3,3) = 1.d0
          HD(3,4) =-dM(J2)*( dEPdT                                      &
     &             -cP * brem                                           &
     &             -dcPdT * (X(J2,4)-VX(J2,4))*brem                     &
     &             +ddeldT/RHO2 * (X(J2,1)-VX(J2,1))*brem               &
     &             -delta/RHO2**2 * dRHdT2 * (X(J2,1)-VX(J2,1))*brem )
          HD(4,1) = 0.5d0*(TPNAB3+TPNAB2)                               &
     &                   - (X(J3,1)-X(J2,1))*0.5d0*dNAB2P
          HD(4,2) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB2R
          HD(4,3) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB2L
          HD(4,4) =-1.d0 - (X(J3,1)-X(J2,1))*0.5d0*dNAB2T
          HE(4,1) =-0.5d0*(TPNAB3+TPNAB2)                               &
     &                   - (X(J3,1)-X(J2,1))*0.5d0*dNAB3P
          HE(4,2) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB3R
          HE(4,3) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB3L
          HE(4,4) = 1.d0 - (X(J3,1)-X(J2,1))*0.5d0*dNAB3T
        endif
        G(1) = X(J3,1)-X(J2,1) + GMdM/X(J2,2)**4
        G(2) = X(J2,2)**3 - X(J1,2)**3 - dM(J2)*0.75d0/(PI*RHO2)
        G(3) = X(J2,3) - X(J1,3) - dM(J2)*( EP                          &
     &       - cP * (X(J2,4)-VX(J2,4))*brem                             &
     &       + delta/RHO2 * (X(J2,1)-VX(J2,1))*brem )
        G(4) = (X(J3,4)-X(J2,4))                                        &
     &         -(X(J3,1)-X(J2,1))*0.5d0*(TPNAB3+TPNAB2)
c       write(6,*) j2, ' Gs ', g(1),g(2),g(3),g(4) 
        return
      endif
      if (J2.eq.N) then
c$$$        if(iter.ge.abs(itmax/2).or.X(j2,3).lt.0.) then
c         call atmos(Tatm,RHOatm,Ratm,Patm,VX(J2,2),VX(J2,3),0)
c
c         if(K2.gt.0) then
c            HD(1,1) = 1.d0
c            HD(2,2) = 1.d0
c            HC(2,2) =-1.d0
c            HD(3,3) = 1.d0
c            HC(3,3) =-1.d0
c            HD(4,4) = 1.d0
c         endif

c$$$        else
          call atmos(Tatm,RHOatm,Ratm,Patm,X(J2,2),X(J2,3),0)
          if(K2.gt.0) then
            RN=X(J2,2)*(1.d0+XINC )
            call atmos(TatR,RHOatR,RatR,PatR,RN,X(J2,3),0)
            xLN=X(J2,3)*(1.d0+XINC) 
            call atmos(TatL,RHOatL,RatL,PatL,X(J2,2),xLN,0)
            dPadr = (PatR - Patm)/(XINC * X(J2,2))
            dPadL = (PatL - Patm)/(XINC * X(J2,3))
            dTadr = (TatR - Tatm)/(XINC * X(J2,2))
            dTadL = (TatL - Tatm)/(XINC * X(J2,3))
            dRadr = (RHOatR - RHOatm)/(XINC * X(J2,2))
            dRadL = (RHOatL - RHOatm)/(XINC * X(J2,3))
            datrdr = (RatR - Ratm)/(XINC * X(J2,2))
            datrdL = (RatL - Ratm)/(XINC * X(J2,3))

            HD(1,1) = 1.d0
            HD(1,2) = -dPadr
            HD(1,3) = -dPadL
            HD(2,2) = 1.d0 - dredgedr - datrdR
            HD(2,3) =      - dredgedl - datrdL
            HC(2,2) =-1.d0
            HD(3,3) = 1.d0
            HC(3,3) =-1.d0
            HD(4,4) = 1.d0
            HD(4,2) = -dTadr
            HD(4,3) = -dTadL
         endif
       endif
       G(1) = X(J2,1)-Patm
       G(2) = X(J2,2)-X(J1,2) - Ratm
       G(3) = X(J2,3)-X(J1,3)
       G(4) = X(J2,4)-Tatm
c       write(6,*) J2,' Gs ', g(1),g(2),g(3),g(4) 
        return
c      endif 
      return
      end
c
      subroutine nabla(PRES,RAD,XLUM,TEMP,XMASS,TPNAB,ADNAB,RADNAB,J,IC)
C
C ** CHECK MIXING LENGTH GRADIENT AND CONVECTIVE VELOCITIES. **
C ** KIPPENHAHN,WEIGERT,HOFMEISTER... AND TALBOT AND SMITH
C
      include 'parm.h'
      include 'var.h'
c     DATA anabcon/3.945322d+09/
      DATA GRAV/6.6704D-8/
      DATA PI  /3.14159265359d0/
      DATA SIG,RG/5.67051d-5,8.31451d7/
      data iwrite/0/
      data Trans,TRwidth/4.d5,10.d0/
      save

      anabcon=3.d0/(GRAV*PI*64.d0*SIG)
      K=abs(J)
        call invstate(31,PRES,TEMP,hydrogen(K),helium4(K),Crad,         &
     &     RHO,cP,alpha,beta,delta)
        RHOlog=log10(RHO)
        Tlog  =log10(TEMP)
        call opacity(K,Tlog,RHOlog,bkap) 
        xkap=10.d0**bkap
        radnab=anabcon/xMASS*xkap*xLUM*PRES/TEMP**4
        adnab=delta*PRES/(cP*TEMP*RHO)
       if(radnab.le.adnab) then
         IC=1
         TPNAB=TEMP/PRES*radnab
       else
         IC=0
         TPNAB=TEMP/PRES*adnab
         FACTR=max(0.d0,(vx(k,4)/Trans-1.d0)*TRwidth)
        if(FACTR.gt.1.d0) return

c$$$        WTM=Rg*RHO*TEMP/PRES
c$$$         GR = GRAV*xmass/RAD**2
c$$$         BET=1.d0
c$$$         HP=Rg*TEMP/GR/WTM/BET
c$$$         DIST=RLH*HP
c$$$C      CONVECTIVE GRADIENT
c$$$         A1=12.d0*SIG*TEMP**3/CP/RHO/RHO/DIST/XKAP
c$$$         YOU=A1*SQRT(8.d0*HP/GR/DIST/DIST/DELTA)
c$$$         E3=0.5048011d0*YOU*YOU
c$$$         A2=YOU*(4.d0/9.d0*(radnab-adnab) + 0.2373622d0*YOU*YOU)
c$$$         ddd=A2*A2+E3*E3*E3
c$$$         Wd=(A2+sqrt(ddd))**(1.d0/3.d0)
c$$$         C1=(Wd*Wd+0.7037037*Wd*YOU-E3)**2
c$$$         DCONV=max(0.d0,C1/Wd/Wd-YOU*YOU)+adnab
c$$$         if (dconv.lt. 0.d0)then
c$$$           dconv=radnab
c$$$         end if
c$$$         TPNAB=FACTR*TPNAB + (1.d0-FACTR)*TEMP/PRES*dconv
      end if
      return
      end

      subroutine GID(J2)
C+
C  The subroutine GID calculates the numerical derivatives to the
C  equations given in GI and compares them to the analytical ones.
C
C  Author:  Harold W. Yorke (JPL)
C  Date:    29-Jan-2003
C-
      include 'parm.h'
      include 'var.h'
      dimension gk(MH),gg(MH)
      data DINC/-0.0001/
c
      save
c
      j1=j2-1
      j3=j2+1
      write(6,221) '***********  GI(',j2,')=',(g(i),i=1,MH)
 221  format(a,i4,a,1p,4E12.4)
 222  format(i5,1p,4E10.2)
 223  format(5x,1p,4E10.2)
 224  format(5x,1p,4E10.2,' ********')
      do i=1,MH
        gk(i)=g(i)
      enddo
      call gi(-j2)
      write(6,221) ' numerical   GI(',j2,')=',(g(i),i=1,MH)
      if(j2.gt.1) then
      do k=1,MH
        xkeep=x(j1,k)
        x(j1,k) = x(j1,k)*(1.+DINC)
        call gi(-j2)
        x(j1,k) = xkeep
        do i=1,MH
          gg(i)=(g(i)-gk(i))/(DINC*xkeep+1.E-37)
        enddo
        write(6,222) j1,(gg(i),i=1,MH)
        write(6,222) k,(HC(i,k),i=1,MH)
        dddx=0.
        do i=1,MH
          xxxx = gg(i)
          gg(i)=.5*(HC(i,k)-gg(i))/(1.E-37+HC(i,k)+gg(i))
          HC(i,k)=xxxx
          dddx=dddx+abs(gg(i))
        enddo
        if(dddx.lt.0.2) then
          write(6,223) (gg(i),i=1,MH)
        else
          write(6,224) (gg(i),i=1,MH)
        endif
      enddo
      endif
      do k=1,MH
        xkeep = x(j2,k)
        x(j2,k) = x(j2,k)*(1.+DINC)
        call gi(-j2)
        x(j2,k) = xkeep
        do i=1,MH
          gg(i)=(g(i)-gk(i))/(DINC*xkeep+1.E-37)
        enddo
        write(6,222) j2,(gg(i),i=1,MH)
        write(6,222) k,(HD(i,k),i=1,MH)
        dddx=0.
        do i=1,MH
          xxxx = gg(i)
          gg(i)=.5*(HD(i,k)-gg(i))/(1.E-37+HD(i,k)+gg(i))
          HD(i,k)=xxxx
          dddx=dddx+abs(gg(i))
        enddo
        if(dddx.lt.0.2) then
          write(6,223) (gg(i),i=1,MH)
        else
          write(6,224) (gg(i),i=1,MH)
        endif
      enddo
      if(j2.lt.N) then
      do k=1,MH
        xkeep = x(j3,k)
        x(j3,k) = x(j3,k)*(1.+DINC)
        call gi(-j2)
        x(j3,k) = xkeep
        do i=1,MH
          gg(i)=(g(i)-gk(i))/(DINC*xkeep+1.E-37)
        enddo
        write(6,222) j3,(gg(i),i=1,MH)
        write(6,222) k,(HE(i,k),i=1,MH)
        dddx=0.
        do i=1,MH
          xxxx = gg(i)
          gg(i)=.5*(HE(i,k)-gg(i))/(1.E-37+HE(i,k)+gg(i))
          HE(i,k)=xxxx
          dddx=dddx+abs(gg(i))
        enddo
        if(dddx.lt.0.2) then
          write(6,223) (gg(i),i=1,MH)
        else
          write(6,224) (gg(i),i=1,MH)
        endif
      enddo
      endif
      call gi(j2)
      do i=1,MH
        g(i)=gk(i)
      enddo
      return
      end
      SUBROUTINE ATMOS(Tatm,RHOatm,Ratm,Patm,Rstar,Xlum,jwrite)
      include 'parm.h'
      include 'var.h'
      parameter  (MTAU=3000)
      CHARACTER*1 CS(2)
      DIMENSION TAU(MTAU),TTAU(MTAU),RHOTAU(MTAU),PTAU(MTAU),           &
     & RTAU(MTAU),ZMTAU(MTAU),ZATG(MTAU),ZRADG(MTAU),                   &
     & ZTRUG(MTAU),ICV(MTAU)
c
c     DATA Arad3/2.52197145D-15/
      DATA SIG,RG,PI/5.67051d-5,8.31451d7,3.14159265359d0/
      DATA GRAV,CC/6.6704D-8,2.99792458d10/
      DATA Z0,Z05,Z1,Z2,Z10,Z13,Z23/0.d0,0.5d0,1.d0,2.d0,10.d0,         &
     & .333333333333333333d0,.6666666666666666667d0/
      DATA CS/'*',' '/ 
      data ifirst/0/
      SAVE

c      write(6,*) "%"
c      write(6,*) "% In atmos, before calcs begin!"
c      write(6,*) "% Rstar=",Rstar,"Xlum=",Xlum
c      write(6,*) "%"

      iwrite=jwrite
      kwrite=abs(iwrite)
      ATMASS1=0.90d0*dM(N)
      ATMASS2=dM(N)
      Arad3=SIG*4.d0/(CC*3.d0)
      XX=hydrogen(N)
      YY=helium3(N)+helium4(N)
      DELTAU = .001d0
      dTAU05 = Z05*DELTAU
      ITAU23 = 0
      IDELM  = 0
      RTAU23 = Rstar
      ZTAU23 = zM(N)
      TEFF4  = Xlum/(4.d0*PI*SIG*Rstar**2)
      TEFF   = SQRT(SQRT(TEFF4))
c      write(6,*) "% Teff value is ",TEFF
      RK0    = Z0
      ZK0    = Z0
      TK0    = SQRT (SQRT(TEFF4*(0.75d0*dTAU05 + Z05)))
c      write(6,*) "% Initial TKO value is ",TK0
      Rat    = Rstar
      Zat    = zM(N)
C
C   Define atmosphere values for JK=1 at 1/2 DELTAU...
C
      Prad=Crad*Arad3*TK0**4
      Tlog=log10(TK0)
      Pout=max(Prad*Z2,Prad+10.)
c      write(6,*) "% Pout value is ",Pout
      call invstate(11,Pout,TK0,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
      RHOlog = log10(RHO)
c      write(6,*) "% Initial rho value is ",RHO
      T6 = TK0/1.0d6
      Rval = RHO/T6**3
c      write(6,*) "logRval = ",log10(Rval)
      call opacity(N,Tlog,RHOlog,bkap) 
      AKK = Z10**bkap
c      write(6,*) "% Initial kappa value is ", AKK
      GD  = DELTAU*GRAV*Zat/Rat**2
c      write(6,*) "% GD value is ", GD
      PK0 = Pout + GD/AKK
c      write(6,*) "% PK0 value is ", PK0
      call invstate(12,PK0,TK0,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
      RHOlog = log10(RHO)
c      write(6,*) "% New initial rho value is ", RHO
      call opacity(N,Tlog,RHOlog,bkap) 
      AKK = Z10**bkap
c      write(6,*) "% New initial kappa value is ", AKK
      do itr=1,30
c         write(6,*) "% "
c         write(6,*) "% ----------------------------"
c         write(6,*) "% Iteration ", itr
        G0 = (PK0-Pout)*AKK - GD
c        write(6,*) "%   PK0 value is ", PK0
c        write(6,*) "%   rhoK0 value is ", RHO
c        write(6,*) "%   kappaK0 value is ", AKK
c        write(6,*) "%   G0 value is ", G0
        PK1=PK0*1.001d0
c        write(6,*) "%   PK1 value is ", PK1
        call invstate(13,PK1,TK0,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
        RHOlog = log10(RHO)
c        write(6,*) "%   rhoK1 value is ", RHO
        T6 = TK0/1.0d6
c        write(6,*) "%   TKO = ", TK0
        Rval = RHO/(T6**3)
c        write(6,*) "%   Rval = ",Rval
c        write(6,*) "%   logRval = ",log10(Rval)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK = Z10**bkap
c        write(6,*) "%   kappaK1 value is ", AKK
        G1 = (PK1-Pout)*AKK - GD
c        write(6,*) "%   G1 value is ", G1
c        write(6,*) "%   G0 - G1 = ",G0-G1
c        write(6,*) "%   PK0 - PK1 = ",PK0 - PK1
        dGdP=(G0-G1)/(PK0-PK1)
c        write(6,*) "%   dGdP value is ", dGdP
        DELP = -G0/dGdP
        DELP = max(DELP,-Z05*PK0)
        DELP = min(DELP,.9d0*PK0)
        DELP = max(DELP,.8d0*(Prad-PK0))
c        write(6,*) "%   delP value is ", DELP
        if(itr.gt.11) DELP=Z05*DELP
        PK0 = PK0 + DELP
c.        write(6,*) "%   updated PK0 value is ", PK0
        call invstate(14,PK0,TK0,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)       
        RHOlog = log10(RHO)
c.        write(6,*) "%   updated rhoK0 value is ", RHO
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK = Z10**bkap
c.        write(6,*) "%   updated kappaK0 value is ", AKK
        if(itr.gt.2 .and. abs(DELP)/(PK0+PK1) .lt. 1.e-5) goto 3
      enddo
 3    continue
c      if(ifirst.eq.0 .or. itr.ge.30)
c     &  write(6,'(a,1p,4e12.4,a,i3,a)')                                 &
c     &     'ATMOS: First values of P,RHO,T,AKM:',PK0,RHO,TK0,AKK,       &
c     &     ' after ',itr,' iterations'

c        write(6,'(a,1p,4e12.4,a,i3,a)')                                 &
c     &     'ATMOS: First values of P,RHO,T,AKM:',PK0,RHO,TK0,AKK,       &
c     &     ' after ',itr,' iterations'
c      write(6,*) "% Done with the initial convergence loop"
c      write(6,*) "% "
c      write(6,*) "% "
      ifirst=1
      TAU(1)   = dTAU05
      PTAU(1)  = PK0
      RTAU(1)  = Z0
      ZMTAU(1) = Z0
      TTAU(1)  = TK0
      RHOTAU(1)= RHO
      call nabla(PK0,Rat,Xlum,TK0,Zat,TPNAB,ADNAB,RADNAB,N,IC)
      ZATG(1)  = ADNAB
      ZRADG(1) = RADNAB
      ZTRUG(1) = TPNAB*PK0/TK0
      ICV(1)   = IC
C
C   Done defining atmospheric values at outermost point. Now do the
C   rest of the atmosphere using Runge-Kutta integration...
C
      do JK=2,MTAU
c         write(6,*) "JK=",JK,"ZK0=",ZK0
C
C   First Runge-Kutta step: get K1 values for P (pressure), R (radius),
C   Z (mass), and T (temperature) by straightfoward extrapolation with
C   a full step of DELTAU...  R and Z are calculated from TAU = 0
C
        TAU(JK) = TAU(JK-1) + DELTAU
c        write(6,*) "JK=",JK,"deltau=",DELTAU,"ZK0=",ZK0
C
C   The following two lines determine radius and mass from center of
C   star, except when outside TAU = 2/3, where Rstar and Mstar are 
C   used.
C
        Rat = Rstar - max(Z0,RK0-RTAU23)
        Zat = zM(N) - max(Z0,ZK0-ZTAU23)
        dTAK = DELTAU/AKK
        GD   = dTAK*GRAV*Zat/Rat**2
        dPK1 = GD
        dRK1 = dTAK/RHO
        dZK1 = dTAK*4.d0*PI*Rat**2
        dTK1 = GD*TPNAB
        
        Tlog1 = log10(TK0)
        RHOlog1 = log10(RHO)

c        write(6,*) "JK=",JK,"dTAK=",dTAK
c       if(JK.eq.2)
c    &  dTK1 = SQRT (SQRT(TEFF4*(0.75d0*TAU(JK) + Z05))) - TK0
C
C   Second Runge-Kutta step: get K2 values for P, R, Z and T at 1/2 step
C
        PM2 = PK0 + Z05*dPK1
        RM2 = RK0 + Z05*dRK1
        ZM2 = ZK0 + Z05*dZK1
        TM2 = TK0 + Z05*dTK1
        Rat = Rstar - max(Z0,RM2-RTAU23)
        Zat = zM(N) - max(Z0,ZM2-ZTAU23)
        call nabla(PM2,Rat,Xlum,TM2,Zat,TPNAB,ADNAB,RADNAB,N,IC)
        call invstate(15,PM2,TM2,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
        Tlog   = log10(TM2)
        RHOlog = log10(RHO)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK = Z10**bkap
        dTAK = DELTAU/AKK
        GD   = dTAK*GRAV*Zat/Rat**2
        dPK2 = GD
        dRK2 = dTAK/RHO
        dZK2 = dTAK*4.d0*PI*Rat**2
        dTK2 = GD*TPNAB

        Tlog2 = Tlog
        RHOlog2 = RHOlog
C
C   Third Runge-Kutta step: get K3 values for P, R, Z and T at 1/2 step
C
        PM2 = PK0 + Z05*dPK2
        RM2 = RK0 + Z05*dRK2
        ZM2 = ZK0 + Z05*dZK2
        TM2 = TK0 + Z05*dTK2
        Rat = Rstar - max(Z0,RM2-RTAU23)
        Zat = zM(N) - max(Z0,ZM2-ZTAU23)
        call nabla(PM2,Rat,Xlum,TM2,Zat,TPNAB,ADNAB,RADNAB,N,IC)
        call invstate(16,PM2,TM2,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
        Tlog   = log10(TM2)
        RHOlog = log10(RHO)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK  = Z10**bkap
        dTAK = DELTAU/AKK
        GD   = dTAK*GRAV*Zat/Rat**2
        dPK3 = GD
        dRK3 = dTAK/RHO
        dZK3 = dTAK*4.d0*PI*Rat**2
        dTK3 = GD*TPNAB

        Tlog3 = Tlog
        RHOlog3 = RHOlog

C
C   Fourth Runge-Kutta step: get K4 values for P, R, Z and T at full step
C
        PM2 = PK0 + dPK3
        RM2 = RK0 + dRK3
        ZM2 = ZK0 + dZK3
        TM2 = TK0 + dTK3
        Rat = Rstar - max(Z0,RM2-RTAU23)
        Zat = zM(N) - max(Z0,ZM2-ZTAU23)
        call nabla(PM2,Rat,Xlum,TM2,Zat,TPNAB,ADNAB,RADNAB,N,IC)
        call invstate(17,PM2,TM2,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
        Tlog   = log10(TM2)
        RHOlog = log10(RHO)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK  = Z10**bkap
        dTAK = DELTAU/AKK
        GD   = dTAK*GRAV*Zat/Rat**2
c        write(6,*) "J=",JK,GD
        dPK4 = GD
        dRK4 = dTAK/RHO
        dZK4 = dTAK*4.d0*PI*Rat**2
        dTK4 = GD*TPNAB

        Tlog4 = Tlog
        RHOlog4 = RHOlog

c      write(6,*) "J=",JK
c      write(6,*) " Tlog:",Tlog1,Tlog2,Tlog3,Tlog4
c      write(6,*) " RHOlog:",RHOlog1,RHOlog2,RHOlog3,RHOlog4
c      write(6,*) "J=",JK,AK1,AK2,AK3,AK4
c      write(6,*) "J=",JK,"dz1=",dZK1,"dz2=",dZK2,"dz3=",dZK3,"dz4=",dZK4
C
C   Final Runge-Kutta step: Combine K1, K2, K3, and K4 values for P, R, Z
C   and T at the grid point JK
C
        PK0 = PK0 + (dPK1 + Z2*(dPK2 + dPK3) + dPK4)/6.d0
        RK0 = RK0 + (dRK1 + Z2*(dRK2 + dRK3) + dRK4)/6.d0
        ZK0 = ZK0 + (dZK1 + Z2*(dZK2 + dZK3) + dZK4)/6.d0
        TK0 = TK0 + (dTK1 + Z2*(dTK2 + dTK3) + dTK4)/6.d0
        Rat = Rstar - max(Z0,RK0-RTAU23)
        Zat = zM(N) - max(Z0,ZK0-ZTAU23)
        call invstate(18,PK0,TK0,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
c        write(6,*) "J=",JK,PK0,TK0,RHO
        Tlog   = log10(TK0)
        RHOlog = log10(RHO)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK = Z10**bkap
 
c        write(6,*) "J=",JK,dPk1,dPk2,dPk3,dPk4
c        write(6,*) "J=",JK,"PK0=",PK0,"TK0=",TK0

        call nabla(PK0,Rat,Xlum,TK0,Zat,TPNAB,ADNAB,RADNAB,N,IC)
        PTAU(JK)  = PK0
        RTAU(JK)  = RK0
        ZMTAU(JK) = ZK0
        TTAU(JK)  = TK0
        RHOTAU(JK)= RHO
        ZATG(JK)  = ADNAB
        ZRADG(JK) = RADNAB
        ZTRUG(JK) = TPNAB*PK0/TK0
        ICV(JK)   = IC
        if(TAU(JK).ge.Z23 .and. ITAU23.eq.0) then
          ITAU23=JK
          FAC    = (TAU(JK) - Z23)/(TAU(JK)-TAU(JK-1))
          RTAU23 = FAC*RTAU(JK-1) + (Z1-FAC)*RTAU(JK)
          ZTAU23 = FAC*ZMTAU(JK-1) + (Z1-FAC)*ZMTAU(JK)
        endif
        if(ZK0 .GE. ATMASS1 .and. IDELM.eq.0) THEN
          IDELM=JK
          FAC = (ZK0 - ATMASS1)/(ZK0-ZMTAU(JK-1))
          Tatm   = FAC*TTAU  (JK-1) + (Z1-FAC)*TTAU  (JK)
          RHOatm = FAC*RHOTAU(JK-1) + (Z1-FAC)*RHOTAU(JK)
          Patm   = FAC*PTAU  (JK-1) + (Z1-FAC)*PTAU  (JK)
        endif
        if(ZK0 .GE. ATMASS2) GO TO 75
        IF(TAU(JK) .GT. Z10) THEN
          DELTA = min(1.08d0 * DELTAU , dTAUmx)
          DELTAU = max(DELTA , 1.01d0* DELTAU )
          DELTAU = 1.01d0* DELTAU
        ELSE
          IF(TAU(JK) .GT. .1d0) DELTAU = .01d0
          IF(TAU(JK) .GT. .8d0) DELTAU = .05d0
        ENDIF
      enddo
      goto 998
75    CONTINUE
      FAC = (ZK0 - ATMASS2)/(ZK0-ZMTAU(JK-1))
      Ratm = FAC*RTAU(JK-1)+ (Z1-FAC)*RTAU(JK) - RTAU23

      if(iwrite.ne.0) then
         write(6,300) Tatm,Ratm,Patm,RHOatm,TEFF,Rstar,Xlum
 300     format('Bottom of atmos:  T=',1p,E10.2,                        &
     &        ' R=',E10.2,' P=',E10.2,' RHO=',E10.2,/,                  &
     &        'Top of atmos:  T=',E10.2,' R=',E10.2, ' L=',E10.2)                                                
      endif

c         write(6,*) 'J   C    dM       M       P      R      dR      T  &
c     &            RHO      TAU      ATG      RAD      TRU    1-beta'

      do j=JK,1,-1
c$$$         if(mod(j,kwrite).eq.0 .or. j.ge.JK-1 .or. j.le.2 .or.          &
c$$$     &        j.eq.IDELM  .or. j.eq.IDELM -1 .or.                       &
c$$$     &        j.eq.ITAU23 .or. j.eq.ITAU23-1) then
c$$$            Prad=Crad*Arad3*TTAU(j)**4
c$$$            beta=Prad/PTAU(j)
c$$$            write(6,301) j,CS(1+ICV(j)),ZMTAU(j),zM(N)-ZMTAU(j),PTAU(j),&
c$$$     &           Rstar-RTAU(j)+RTAU23,RTAU(j),TTAU(j),RHOTAU(j),        &
c$$$     &           TAU(j),ZATG(j),ZRADG(j),ZTRUG(j),beta
 301        format(1x,i4,1x,a1,1p,12E9.2)
 302        format(1x,4('*'),1x,a1,1p,12E9.2)
c$$$         endif
          if(j.eq.JK) then
            FAC = (ZMTAU(j) - ATMASS2)/(ZMTAU(j)-ZMTAU(j-1))
            PX  = FAC*PTAU(j-1) + (Z1-FAC)*PTAU(j)
            TX  = FAC*TTAU(j-1) + (Z1-FAC)*TTAU(j)
            RX  = FAC*TAU (j-1) + (Z1-FAC)*TAU (j)
            ZMX = FAC*ZMTAU(j-1) + (Z1-FAC)*ZMTAU(j)
            RHX = FAC*RHOTAU(j-1) + (Z1-FAC)*RHOTAU(j)
            ATGX= FAC*ZATG(j-1) + (Z1-FAC)*ZATG(j)
            RADX= FAC*ZRADG(j-1) + (Z1-FAC)*ZRADG(j)
            TRUX= FAC*ZTRUG(j-1) + (Z1-FAC)*ZTRUG(j)
            IC  = FAC*ICV(j-1) + (Z1-FAC)*ICV(j)
            beta= Crad*Arad3*(FAC/PTAU(j-1)*TTAU(j-1)**4                &
     &          +        (Z1-FAC)/PTAU(j  )*TTAU(j  )**4)
c            write(6,302) CS(1+IC),ZMX,zM(N)-ZMX,PX,                     &
c     &            Rstar-Ratm,Ratm+RTAU23,TX,RHX,RX,ATGX,RADX,TRUX,beta
          endif
          if(j.eq.IDELM) then
            FAC = (ZMTAU(j) - ATMASS1)/(ZMTAU(j)-ZMTAU(j-1))
            RX  = FAC*RTAU(j-1) + (Z1-FAC)*RTAU(j)
            TX  = FAC*TAU (j-1) + (Z1-FAC)*TAU (j)
            ZMX = FAC*ZMTAU(j-1) + (Z1-FAC)*ZMTAU(j)
            ATGX= FAC*ZATG(j-1) + (Z1-FAC)*ZATG(j)
            RADX= FAC*ZRADG(j-1) + (Z1-FAC)*ZRADG(j)
            TRUX= FAC*ZTRUG(j-1) + (Z1-FAC)*ZTRUG(j)
            IC  = FAC*ICV(j-1) + (Z1-FAC)*ICV(j)
            beta= Crad*Arad3*(FAC/PTAU(j-1)*TTAU(j-1)**4                &
     &          +        (Z1-FAC)/PTAU(j  )*TTAU(j  )**4)
c            write(6,302) CS(1+IC),ZMX,ZM(N)-ZMX,Patm,                   &
c     &            Rstar-RX+RTAU23,RX,Tatm,RHOatm,TX,ATGX,RADX,TRUX,beta
          endif
          if(j.eq.ITAU23) then
            FAC = (TAU(j) - Z23)/(TAU(j)-TAU(j-1))
            PX  = FAC*PTAU(j-1) + (Z1-FAC)*PTAU(j)
            TX  = FAC*TTAU(j-1) + (Z1-FAC)*TTAU(j)
            ZMX = FAC*ZMTAU(j-1) + (Z1-FAC)*ZMTAU(j)
            RHX = FAC*RHOTAU(j-1) + (Z1-FAC)*RHOTAU(j)
            ATGX= FAC*ZATG(j-1) + (Z1-FAC)*ZATG(j)
            RADX= FAC*ZRADG(j-1) + (Z1-FAC)*ZRADG(j)
            TRUX= FAC*ZTRUG(j-1) + (Z1-FAC)*ZTRUG(j)
            IC  = min(ICV(j-1),ICV(j))
            beta= Crad*Arad3*(FAC/PTAU(j-1)*TTAU(j-1)**4                &
     &          +        (Z1-FAC)/PTAU(j  )*TTAU(j  )**4)
c            write(6,302) CS(1+IC),ZMX,zM(N)-ZMX,PX,                     &
c     &                 Rstar,RTAU23,TX,RHX,Z23,ATGX,RADX,TRUX,beta
         endif
      enddo
c      endif
      taujk=TAU(JK)
      RETURN

 998  write(6,*) ' atmos: 998'
      write(6,*) "ZK0,ATMASS1,ATMASS2"
      write(6,*) ZK0,ATMASS1,ATMASS2
      write(6,*) "Tatm,RHOatm,Ratm,Patm,Rstar,Xlum,taujk,iwrite"
      write(6,*) Tatm,RHOatm,Ratm,Patm,Rstar,Xlum,taujk,iwrite
      stop 'ATMOS: 998'
      END

      SUBROUTINE ADDSUB
C+
C  THE SUBROUTINE ADDSUB DETERMINES WHETHER A GRID POINT SHOULD BE ADDED
C  OR REMOVED.
C
C  Author:  H.W. YORKE   23-SEP-02               (JPL / CALTECH)
C  Version: 12-May-2006
C  History: (see summary in CHANGES)
C  
C-
      include 'parm.h'
      include 'var.h'
      parameter(Mk=17)
      character*1 why(3)
      character*3 blank,why3
      equivalence (why,why3)
      data blank/'   '/
      data Nadd,Nsub/0,0/
      save
c
c CHECK WHETHER A POINT SHOULD BE ADDED TO THE OUTER ZONE
c
      why3=blank
      Jadd=0
      Jsub=0
      if(x(N,4).gt.Atmx) then
        do jjj=1,2
          Jadd=N
          why(1)='A'
          call add(Jadd,why)
          Nadd=Nadd+1
        enddo
      endif
c
c CHECK WHETHER A POINT SHOULD BE REMOVED FROM THE OUTER ZONE
c
      if(x(N,4).lt.Atmn .and. (dM(N-1)+dM(N))/zM(N-1) .lt. dzmax) then
        Jsub=N-1
        call sub(Jsub)
        Nsub=Nsub+1
      endif
c
c CHECK WHETHER MORE POINTS SHOULD BE ADDED CLOSE TO THE EXTERIOR
c
      if(Jadd+Jsub .eq. 0) then
        do j=N-5,N-1
          if(dM(j)/zM(j) .gt. dzmax) then
            Jadd = J
            why(1)='Z'
            call add(Jadd,why)
            Nadd=Nadd+1
            goto 10
          endif
        enddo
      endif
 10   continue
c
c CHECK WHETHER GRID POINTS SHOULD BE ADDED OR REMOVED
c
      dPL =1.E37
      Slumax=0.
      do j=1,N
        Slumax=max(Slumax,abs(x(j,3)))
      enddo
      Jadd=0
      Jsub=0
      do j=2,99999
        dP   =abs((x(j,1)-x(j-1,1))/(x(j,1)+x(j-1,1)))
        dLum =abs((x(j,3)-x(j-1,3))/Slumax)
        dX   =abs(hydrogen(j)-hydrogen(j-1))
     *       +abs(helium4 (j)-helium4 (j-1))
        if(J.le.N-2) then
          dP2  =abs((x(j+1,1)-x(j,1))/(x(j,1)+x(j+1,1)))
          dLum2=abs((x(j+1,3)-x(j,3))/Slumax)
          dX2  =abs(hydrogen(j)-hydrogen(j+1))
     *         +abs(helium4 (j)-helium4 (j+1))
          dP3  =abs((x(j+1,1)-x(j+2,1))/(x(j+2,1)+x(j+1,1)))
          dLum3=abs((x(j+1,3)-x(j+2,3))/Slumax)
          dX3  =abs(hydrogen(j+2)-hydrogen(j+1))
     *         +abs(helium4 (j+2)-helium4 (j+1))
          if((J-Jsub.gt.2 .and. J-Jadd.gt.2 .and.
     *      dLum.lt.dLmn .and. dLum2.lt.dLmn .and. dLum3.lt.dLmn .and.
     *      dP.lt.dPmn .and. dP2.lt.dPmn .and. dP3.lt.dPmn .and.
     *      dX.lt.dXmn .and. dX2.lt.dXmn .and. dX3.lt.dXmn .and.
     *      (dM(j+1)+dM(j))/zM(N) .lt. dzmax)
     *      .or. (dM(j+1)+dM(j))/zM(j) .lt. dzmin)
     *    then
            Jsub=j
            call sub(Jsub)
            Nsub=Nsub+1
          endif
        endif
        iwhy=0
        why3=blank
        if(dLum.gt.dLmx) then
          iwhy=iwhy+1
          why(iwhy)='L'
        endif
        if(dP.gt.dPmx) then
          iwhy=iwhy+1
          why(iwhy)='P'
        endif
        if(dX.gt.dXmx) then
          iwhy=iwhy+1
          why(iwhy)='X'
        endif
        if(J-Jadd.gt.2 .and. iwhy.ne.0 .and. N.lt.MJ .and.
     *     dM(J)/zM(J) .gt. 2.*dzmin .and. J-Jsub.gt.1) then
          Jadd = J
          call add(Jadd,why)
          Nadd=Nadd+1
        endif
        if(j.ge.N-1) goto 20
      enddo
 20   continue
      if(Nadd+Nsub.gt.0 .and. mod(MODEL,JRIT).eq.0 .and. JRIT.gt.1) then
        write(6,201) Nsub,Nadd,N
 201    format(' ===> A total of ',I4,' points were subtracted and ',I4,
     *         ' points were added. N=',I4)
        Nadd=0
        Nsub=0
      endif
      return
      end
      SUBROUTINE ADD(J,why)
C+
C  THE SUBROUTINE ADD ADDS A GRID POINT AT J
C
C  Author:  H.W. YORKE   30-AUG-02               (JPL / CALTECH)
C  Version: 12-May-2006
C  History: (see summary in CHANGES)
C-
      include 'parm.h'
      include 'var.h'
      character*1 why(3)
      data PI/3.14159265359d0/
      save
C
      if(N.lt.MJ .and. J.eq.N) then
        dM(J)=dM(J)*0.5d0
        dM(J+1)=dM(J)
        zM(J+1)=zM(J)
        zM(J) = zM(J)-dM(J)
        hydrogen(J+1) = hydrogen(J)
        helium3 (J+1) = helium3 (J)
        helium4 (J+1) = helium4 (J)
        deuterium(J+1)=deuterium(J)
        carbon  (J+1) = carbon  (J)
        trogen  (J+1) = trogen  (J)
        oxygen  (J+1) = oxygen  (J)
        vhydrogen(J+1) = vhydrogen(J)
        vhelium3 (J+1) = vhelium3 (J)
        vhelium4 (J+1) = vhelium4 (J)
        vdeuterium(J+1)=vdeuterium(J)
        vcarbon  (J+1) = vcarbon  (J)
        vtrogen  (J+1) = vtrogen  (J)
        voxygen  (J+1) = voxygen  (J)
C  Interpolate/extrapolate for new values of pressure at J=N,N+1
        w1 = dM(J)/(dM(J-1)+dM(J)+dM(J+1))
        dP = (x(J,1)-x(J-1,1))*w1
        x (J+1,1) = x (J,1) + dP
        x (J  ,1) = x (J,1) - dP
        dP = (vx(J,1)-vx(J-1,1))*w1
        vx(J+1,1) = vx(J,1) + dP
        vx(J  ,1) = vx(J,1) - dP
C  Use average density to estimate radius of new gridpoint at J=N
        x(J+1,2) = x(J,2)
        RHOPI  = (x(J,2)**3-x(J-1,2)**3)/(dM(J)+dM(J+1))
        x(J,2) = (x(J-1,2)**3 + dM(J)*RHOPI)**(1./3.)
        vx(J+1,2) = vx(J,2)
        RHOPI  = (vx(J,2)**3-vx(J-1,2)**3)/(dM(J)+dM(J+1))
        vx(J,2)= (vx(J-1,2)**3 + dM(J)*RHOPI)**(1./3.)
C  Interpolate for new values of luminosity at J=N
        x (J+1,3) = x(J,3)
        vx(J+1,3) = vx(J,3)
        x (J,3) = 0.5d0*(x(J-1,3)+x(J,3))
        vx(J,3) = 0.5d0*(vx(J-1,3)+vx(J,3))
C  Interpolate/extrapolate for new values of temperature at J=N,N+1
        dT = (x(J,4)-x(J-1,4))*w1
        x (J+1,4) = x (J,4) + dT
        x (J  ,4) = x (J,4) - dT
        dT = (vx(J,1)-vx(J-1,1))*w1
        vx(J+1,4) = vx(J,4) + dT
        vx(J  ,4) = vx(J,4) - dT
      endif
      if(N.lt.MJ .and. J.lt.N) then
        do i=N,J,-1
          dM      (i+1) = dM      (i)
          zM      (i+1) = zM      (i)
          hydrogen(i+1) = hydrogen(i)
          helium3 (i+1) = helium3 (i)
          helium4 (i+1) = helium4 (i)
          deuterium(i+1)=deuterium(i)
          carbon  (i+1) = carbon  (i)
          trogen  (i+1) = trogen  (i)
          oxygen  (i+1) = oxygen  (i)
          vhydrogen(i+1) = vhydrogen(i)
          vhelium3 (i+1) = vhelium3 (i)
          vhelium4 (i+1) = vhelium4 (i)
          vdeuterium(i+1)=vdeuterium(i)
          vcarbon  (i+1) = vcarbon  (i)
          vtrogen  (i+1) = vtrogen  (i)
          voxygen  (i+1) = voxygen  (i)
          do k=1,MH
            vx(i+1,k)=vx(i,k)
            x (i+1,k)=x (i,k)
          enddo
        enddo
        dM(J)=dM(J)*0.5d0
        dM(J+1)=dM(J)
        zM(J) = zM(J)-dM(J)
        if(J.gt.1) then
C  Interpolate for new values of pressure
          w1 = dM(J)/(4.d0*dM(J)+dM(J-1)+dM(J+2))
          dP = (x(J+2,1)-x(J-1,1))*w1
          x (J  ,1) = x (J  ,1) - dP
          x (J+1,1) = x (J+1,1) + dP
          dP = (vx(J+2,1)-vx(J-1,1))*w1
          vx(J  ,1) = vx(J  ,1) - dP
          vx(J+1,1) = vx(J+1,1) + dP
C  Use average density to estimate radius of new gridpoint
          RHOPI  = (x(J+2,2)**3-x(J-1,2)**3)/(dM(J)+dM(J+1))
          x(J,2) = (x(J-1,2)**3 + dM(J)*RHOPI)**(1./3.)
          RHOPI  = (vx(J+2,2)**3-vx(J-1,2)**3)/(dM(J)+dM(J+1))
          vx(J,2)= (vx(J-1,2)**3 + dM(J)*RHOPI)**(1./3.)
C  Interpolate for new values of luminosity
          x (J,3) = 0.5d0*(x(J-1,3)+x(J+1,3))
          vx(J,3) = 0.5d0*(vx(J-1,3)+vx(J+1,3))
C  Interpolate for new values of temperature
          dT = (x(J+2,4)-x(J-1,4))*w1
          x (J  ,4) = x (J  ,4) - dT
          x (J+1,4) = x (J+1,4) + dT
          dT = (vx(J+2,1)-vx(J-1,1))*w1
          vx(J  ,4) = vx(J  ,4) - dT
          vx(J+1,4) = vx(J+1,4) + dT
        else
C  Interpolate for new values of pressure at J=1
          w1 = dM(J)/(dM(J)+dM(J+1)+dM(J+2))
          dP = (x(J+2,1)-x(J+1,1))*w1
          x (J  ,1) = x (J  ,1) - dP
          x (J+1,1) = x (J+1,1) + dP
          dP = (vx(J+2,1)-vx(J+1,1))*w1
          vx(J  ,1) = vx(J  ,1) - dP
          vx(J+1,1) = vx(J+1,1) + dP
C  Use average density to estimate radius of new gridpoint at J=1
          RHOPI  = x(J+2,2)**3/(dM(J)+dM(J+1))
          x(J,2) = (dM(J)*RHOPI)**(1./3.)
          RHOPI  = vx(J+2,2)**3/(dM(J)+dM(J+1))
          vx(J,2)= (dM(J)*RHOPI)**(1./3.)
C  Interpolate for new values of luminosity at J=1
          x (J,3) = 0.5d0*x(J+1,3)
          vx(J,3) = 0.5d0*vx(J+1,3)
C  Interpolate for new values of temperature at J=1
          dT = (x(J+2,4)-x(J+1,4))*w1
          x (J  ,4) = x (J  ,4) - dT
          x (J+1,4) = x (J+1,4) + dT
          dT = (vx(J+2,1)-vx(J+1,1))*w1
          vx(J  ,4) = vx(J  ,4) - dT
          vx(J+1,4) = vx(J+1,4) + dT
        endif
      endif
      if(N.lt.MJ) then
        N=N+1
        if(mod(MODEL,JRIT).eq.0 .and. JRIT.eq.1)
     *    write(6,201) J,zM(J),dM(J+1),why
 201      format(' ===> New grid point added at J=',i4,'. M=',1p,E13.6,
     *           ' dM=',E11.4,1x,3a1)
      endif
      return
      end
      SUBROUTINE SUB(J)
C+
C  THE SUBROUTINE SUB REMOVES A GRID POINT AT J
C
C  Author:  H.W. YORKE   23-AUG-02               (JPL / CALTECH)
C  Version: 12-May-2006
C  History: (see summary in CHANGES)
C-
      include 'parm.h'
      include 'var.h'
      save
C
      if(J.lt.N .and. J.ge.1) then
        dMJ=dM(J)+dM(J+1)
        w2=dM(J)/dMJ
        w1=dM(J+1)/dMJ
        hydrogen(J) = w1*hydrogen(J) + w2*hydrogen(J+1)
        helium3 (J) = w1*helium3 (J) + w2*helium3 (J+1)
        helium4 (J) = w1*helium4 (J) + w2*helium4 (J+1)
        deuterium(J)= w1*deuterium(J)+ w2*deuterium(J+1)
        carbon  (J) = w1*carbon  (J) + w2*carbon  (J+1)
        trogen  (J) = w1*trogen  (J) + w2*trogen  (J+1)
        oxygen  (J) = w1*oxygen  (J) + w2*oxygen  (J+1)
        vhydrogen(J) = w1*vhydrogen(J) + w2*vhydrogen(J+1)
        vhelium3 (J) = w1*vhelium3 (J) + w2*vhelium3 (J+1)
        vhelium4 (J) = w1*vhelium4 (J) + w2*vhelium4 (J+1)
        vdeuterium(J)= w1*vdeuterium(J)+ w2*vdeuterium(J+1)
        vcarbon  (J) = w1*vcarbon  (J) + w2*vcarbon  (J+1)
        vtrogen  (J) = w1*vtrogen  (J) + w2*vtrogen  (J+1)
        voxygen  (J) = w1*voxygen  (J) + w2*voxygen  (J+1)
        x (J,1) = w1*x(J,1) + w2*x(J+1,1)
        x (J,2) = x(J+1,2)
        x (J,3) = x(J+1,3)
        x (J,4) = w1*x(J,4) + w2*x(J+1,4)
        vx(J,1) = w1*vx(J,1) + w2*vx(J+1,1)
        vx(J,2) = x(J+1,2)
        vx(J,3) = x(J+1,3)
        vx(J,4) = w1*vx(J,4) + w2*vx(J+1,4)
        dM(J) = dMJ
        zM(J)  = zM(J+1)
        if(J+1.lt.N-1) then
          do i=J+1,N-1
            zM      (i) = zM      (i+1)
            dM      (i) = dM      (i+1)
            hydrogen(i) = hydrogen(i+1)
            helium3 (i) = helium3 (i+1)
            helium4 (i) = helium4 (i+1)
            deuterium(i)=deuterium(i+1)
            carbon  (i) = carbon  (i+1)
            trogen  (i) = trogen  (i+1)
            oxygen  (i) = oxygen  (i+1)
            vhydrogen(i) = vhydrogen(i+1)
            vhelium3 (i) = vhelium3 (i+1)
            vhelium4 (i) = vhelium4 (i+1)
            vdeuterium(i)=vdeuterium(i+1)
            vcarbon  (i) = vcarbon  (i+1)
            vtrogen  (i) = vtrogen  (i+1)
            voxygen  (i) = voxygen  (i+1)
            do k=1,MH
              vx(i,k)=vx(i+1,k)
              x (i,k)=x (i+1,k)
            enddo
          enddo
        endif
        N=N-1
      endif
      if(JRIT.eq.1) then
        if(J.eq.1) then
          if(mod(MODEL,JRIT).eq.0) write(6,201) J,zM(J)
        else
          if(mod(MODEL,JRIT).eq.0) write(6,201) J,zM(J-1),zM(J)
        endif
      endif
 201  format(' ===> Grid point removed at J=',i4,'. M=',1p,E13.6,       &
     &     ' dM=',E11.4)
      return
      end
      SUBROUTINE GRIDMOV
C+
C  THE SUBROUTINE GRIDMOV MODIFIES THE MASS DISTRIBUTION OF GRID
C
C  Author:  H.W. YORKE   18-DEC-02               (JPL / CALTECH)
C  Version: 12-May-2006
C  History: (see summary in CHANGES)
C-
      include 'parm.h'
      include 'var.h'
      parameter(Mk=15)
      dimension WT(Mk),Zold(Mk),Fold(Mk+2),z(Mk)
      save

      Znorm = zM(N)-zM(N-Mk+1)
      Zold(1)=0.d0
      do i=2,Mk
        j=N-Mk+i
        Zold(i)=Zold(i-1)+dM(j)/Znorm
      enddo
C
C  Try to keep the bottom of the atmosphere at Tmean (Tmean is defined as
C  the average of the allowed maximum and minimum atmasphere temperatures).
      Tmean=.5*(Atmx+Atmn)
      dZ=X(N,4)/Tmean-1.d0
      if(dZ.gt.0.d0) then
        dZ=min(.5d0,dZ)*dZmax
      else
        dZ=max(-.5d0,dZ)*dZmax
      endif
      dZ=dZ*dZdt
      WX = min(.5d0,dZdt*0.05d0)
      do i=2,Mk-1
        zev  = (1.d0+dZ)*Zold(Mk-1)*float(i-1)/float(Mk-2)
        ddz  = min(0.5d0*(Zold(i+1)-Zold(i)),zev-Zold(i))
        ddz  = max(ddz,0.5d0*(Zold(i-1)-Zold(i)))
        z(i) = Zold(i) + WX*ddz
      enddo
      z(1) = 0.d0
      z(Mk)= 1.d0
      do i=2,Mk-1
        if(z(i).lt.Zold(i)) then
          WT(i)=(z(i)-Zold(i))/(Zold(i)-Zold(i-1))
        else
          WT(i)=(z(i)-Zold(i))/(Zold(i+1)-Zold(i))
        endif
      enddo
      WT(1)=0.d0
      WT(Mk)=0.d0
C  Interpolate quantities defined on mass grid (R and L)
      call interpol(WT,x (1,2),Fold,N,Mk)
      call interpol(WT,vx(1,2),Fold,N,Mk)
      call interpol(WT,x (1,3),Fold,N,Mk)
      call interpol(WT,vx(1,3),Fold,N,Mk)
C  Redefine mass grid 
      zMass = zM(N)
      do i=2,Mk
        j=N-Mk+i
        dM(j) = Znorm*(z(i)-z(i-1))
        zM(j)  = zM(j-1) + dM(j)
      enddo
      zM(N) = zMass
      if(mod(MODEL,JRIT).eq.0)
     *     write(6,201) N,(z(i)-Zold(i),i=Mk-3,Mk-1),(x(j,4),j=N-2,N)
 201  format(' GRIDMOV:  N=',i4,1p,6E11.3)
C  Redefine interpolation grid to middle of mass grid
      do i=2,Mk
        z(i)=0.5d0*(z(i-1)+z(i))
        Zold(i)=0.5d0*(Zold(i-1)+Zold(i))
      enddo
      Zold(1) = -0.5d0*dM(N-Mk+1)/Znorm
      do i=2,Mk-1
        if(z(i).lt.Zold(i)) then
          WT(i)=(z(i)-Zold(i))/(Zold(i)-Zold(i-1))
        else
          WT(i)=(z(i)-Zold(i))/(Zold(i+1)-Zold(i))
        endif
        WT(Mk)=min(0.d0,(z(Mk)-Zold(Mk))/(Zold(Mk)-Zold(Mk-1)))
      enddo
C  Interpolate quantities defined between mass gridpoints
      call interpol(WT,x (1,1),Fold,N,Mk)
      call interpol(WT,vx(1,1),Fold,N,Mk)
      call interpol(WT,x (1,4),Fold,N,Mk)
      call interpol(WT,vx(1,4),Fold,N,Mk)
      call interpol(WT,hydrogen,Fold,N,Mk)
      call interpol(WT,helium3 ,Fold,N,Mk)
      call interpol(WT,helium4 ,Fold,N,Mk)
      call interpol(WT,deuterium,Fold,N,Mk)
      call interpol(WT,carbon  ,Fold,N,Mk)
      call interpol(WT,trogen  ,Fold,N,Mk)
      call interpol(WT,oxygen  ,Fold,N,Mk)
      call interpol(WT,vhydrogen,Fold,N,Mk)
      call interpol(WT,vhelium3 ,Fold,N,Mk)
      call interpol(WT,vhelium4 ,Fold,N,Mk)
      call interpol(WT,vdeuterium,Fold,N,Mk)
      call interpol(WT,vcarbon  ,Fold,N,Mk)
      call interpol(WT,vtrogen  ,Fold,N,Mk)
      call interpol(WT,voxygen  ,Fold,N,Mk)
      return
      end
      subroutine interpol(WT,F,Fold,N,Mk)
      include 'parm.h'
      dimension WT(Mk),F(N),Fold(Mk+2)
      save

      do i=1,Mk
        j=N-Mk+i
        Fold(i)=F(j)
      enddo
      do i=2,Mk-1
        if(WT(i).gt.0.d0) then
          W3=WT(i)
          W2=1.d0-W3
          W1=0.d0
        else
          W1=-WT(i)
          W2=1.d0-W1
          W3=0.d0
        endif
        j=N-Mk+i
        F(j) = W1*Fold(i-1)+W2*Fold(i)+W3*Fold(i+1)
      enddo
      W1=-WT(Mk)
      W2=1.d0-W1
      F(N) = W1*Fold(Mk-1) + W2*Fold(Mk)
      return
      end
      SUBROUTINE MASSFLUX(time,Zflux)
C+
C  THE SUBROUTINE MASSFLUX SPECIFIES THE CURRENT MASS FLUX ONTO
C  THE STAR (OR MASS LOSS FROM THE STAR IF ZFLUX IS NEGATIVE)
C
C  Author:  H.W. YORKE   18-DEC-02               (JPL / CALTECH)
C  Version: 12-May-2006
C  History: (see summary in CHANGES)
C-
      include 'parm.h'
      parameter (MF=20)
      dimension tim(MF),flx(MF)
      data ifirst/0/
      save

      if(ifirst.eq.0) then
        ifirst=1
        do i=1,MF
          tim(i)=-1.
          call comrd
          READ (5,*,err=10,end=10) TIM(i),FLX(i)
          if(tim(i).lt.0.d0) goto 10
        enddo
        i=MF
 10     NF=i-1
      endif
      do i=1,NF
        if(time.lt.tim(i)) goto 20
      enddo
      i=NF+1
 20   continue
      if(i.eq.1) then
        Zflux=flx(1)
      else
        if(i.gt.NF) then
          Zflux=flx(NF)
        else
          dt=(time-tim(i-1))/(tim(i)-tim(i-1))
          Zflux=flx(i-1) + dt * (flx(i)-flx(i-1))
        endif
      endif
      return
      end

      subroutine addmass
      include 'parm.h'
      include 'var.h'
      save
c
      call massflux(time,Zflux)
      dM(N) = dM(N) + Zflux*dtime
      if(dM(N).le.0.d0) stop 'addmass'
c     do j=1,N
c       dM(j) = dM(j)*(1.+1.d-6)
c       zM(j) = zM(j)*(1.+1.d-6)
c     enddo
      zM(N)  = zM(N-1) + dM(N)
      Zmass  = zM(N)
      return
      end
