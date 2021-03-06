	  SUBROUTINE CREEP(DECRA,DESWA,STATEV,SERD,EC,ESW,P,QTILD,
     1 TEMP,DTEMP,PREDEF,DPRED,TIME,DTIME,CMNAME,LEXIMP,LEND,
     2 COORDS,NSTATV,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
	  REAL STRAIN
	  REAL NOMR
	  REAL DENOMR
	  REAL CRINC
	  REAL NOMA
	  REAL NOMB
C
      DIMENSION DECRA(5),DESWA(5),STATEV(*),PREDEF(*),DPRED(*),
     1 TIME(3),COORDS(*),EC(2),ESW(2)
C
C     DEFINE CONSTANTS
C
	  CA = 21.818d0
	  CB = 432.647d0
	  CC = 10.093d0
	  CD = 24.249d0
	  CE = 493.213d0
	  CF = 17.3194d0
	  CG = 2.4305d0
	  CH = 60.5657d0
	  CI = 1.8360d0
C
C     Update current STRAIN level
      STRAIN = EC(1) + 0.001
C     Calculation of creep increment
      NOMR = (QTILD-CB*EXP(-CC*STRAIN) + CE*EXP(-CF*STRAIN) - CH*EXP(CI*STRAIN))
	  DENOMR = (CA*EXP(-CC*STRAIN) - CD*EXP(-CF*STRAIN) + CG*EXP(CI*STRAIN))
	  CRINC = 10**(NOMR/DENOMR)
	  DECRA(1) = DTIME*CRINC
	  IF(LEXIMP.EQ.1) THEN
	  NOMA = (CB*CC*EXP(-CC*STRAIN) - CE*CF*EXP(-CF*STRAIN) - CH*CI*EXP(CI*STRAIN))*DENOMR
	  NOMB = (-CA*CC*EXP(-CC*STRAIN) + CD*CF*EXP(-CF*STRAIN)+CG*CI*EXP(CI*STRAIN))*NOMR
	  DECRA(2) = DTIME*LOG(10.)*CRINC*((NOMA-NOMB)/(DENOMR)**2)
	  DECRA(5) = DTIME*LOG(10.)*CRINC*(1/(DENOMR))
	  END IF     
C
      RETURN
      END SUBROUTINE CREEP