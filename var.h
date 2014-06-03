      COMMON/ GIC / X(MJ,MH), VX(MJ,MH),                                &
     &     zM(MJ),dM(MJ),Zmass,Zflux,                                   &
     &     TIME,DTIME,FACTIM,DTMIN,DTMAX,CHGMIN,CHGMAX,CHANGE,          &
     &     Atmx,Atmn,dLmx,dLmn,dPmx,dPmn,dXmn,dXmx,dZmax,dZmin,         &
     &     dZdt,dTAUmx,RLH,taujk,taufactor,Crad,Cwrk,xH2,               &
     &     JRIT,NRIT,MODEL,N,NG,METHD,NPROB,NMOD,iUNIT,jUNIT,NREC,NATM, &
     &     TMAX, TWRT
      COMMON/  HY  /HA(MH,2*MH+1),HW(MH,MH+1,MJ),HX(MH,MH+1),           &
     &              E(MH),G(MH),HC(MH,MH),HD(MH,MH),HE(MH,MH)
      COMMON/ CORR /DG(MH),SMIN(MH),SMAX(MH),EPS(MH),                   &
     &              JG(MH),ITMIN,ITMAX,ITER
      COMMON/helium/hydrogen(MJ),helium3(MJ),helium4(MJ),               &
     &              deuterium(MJ),carbon(MJ),trogen(MJ),oxygen(MJ),     &
     &              vhydrogen(MJ),vhelium3(MJ),vhelium4(MJ),            &
     &              vdeuterium(MJ),vcarbon(MJ),vtrogen(MJ),voxygen(MJ)
      COMMON/XOPAC/ZKAP(14,60)
      common/ahycopacity/rhyr(23),thyr(28),auxchyop(23,28)
      common/ahecopacity/rhel(25),thel(32),auxcheop(25,32)
      COMMON/MIX/ICONVECT(MJ) 
	COMMON/ORBVALS/ecc,aval,orbn,qval,zstar,efac,frequency
