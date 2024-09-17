* Encoding: UTF-8.
/* MEMORE for SPSS Version 3.Beta*/.
/* Copyright 2021 */.
/* by Amanda Kay Montoya */.
/* akmontoya.com*/.
/* Documentation available online at akmontoya.com */.
/* or by email to akmontoya@ucla.edu */.


preserve. 
set printback=off.

/* Permission is hereby granted, free of charge, to any person obtaining a copy of this */.
/* software and associated documentation files (the "Software"), to use the software */.
/* in this form.  Distribution after modification is prohibited, as is its use for any */.  
/* commercial purpose without authorization. This software should not be posted or */.
/* stored on any webpage, server, or directory accessible to the public whether free */.
/* or for a charge unless written permission has been granted by the copyright holder.*/.
/* The copyright holder requests that this software be distributed by directing users */.
/* to akmontoya.com where the latest release of the software and documentation is */.
/* archived and  can be downloaded.*/.

/* THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, */.
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF */.
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT */.
/* IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, */.
/*  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT */.
/* OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE */.
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE */.

/* The above text should be included in any distribution of the software */.

preserve. 
set printback=off.

define CDFINVT (p = !charend('/') /df = !charend('/')). 
compute p0=-.322232431088.
compute p1 = -1.
compute p2 = -.342242088547.
compute p3 = -.0204231210245.
compute p4 = -.0000453642210148.
compute q0 = .0993484626060.
compute q1 = .588581570495.
compute q2 = .531103462366.
compute q3 = .103537752850.
compute q4 = .0038560700634.
compute ppv = !p.
  do if (!p > .5).
    compute ppv = 1-!p.
  end if.
  compute y5=sqrt(-2*ln(ppv)).
  compute xp=y5+((((y5*p4+p3)*y5+p2)*y5+p1)*y5+p0)/((((y5*q4+q3)*y5+q2)*y5+q1)*y5+q0).
  do if (!p <= .5).
    compute xp = -xp.
  end if.
compute toutput = sqrt(!df*(exp((!df-(5/6))*(xp**2)/(!df-(2/3)+.1/!df)**2)-1)). 
!enddefine.  

define CHOOSE(r = !charend('/') /k = !charend('/')). 
COMPUTE r = trunc(!r). 
COMPUTE k = trunc(!k). 
COMPUTE rfact = 1. 
COMPUTE rmkfact = 1. 
COMPUTE kfact = 1. 
LOOP z = 1 TO r. 
COMPUTE rfact = rfact*z. 
END LOOP. 
LOOP z = 1 TO (r-k). 
COMPUTE rmkfact = rmkfact*z. 
END LOOP. 
LOOP z = 1 TO k. 
COMPUTE kfact = kfact*z. 
END LOOP. 
COMPUTE rchoosek = rfact/(rmkfact*kfact). 
!enddefine. 

define SOBEL(a = !charend('/') /sea = !charend('/') /b = !charend('/') /seb = !charend('/')). 
    COMPUTE product = !a*!b. 
    COMPUTE sobse = sqrt((!a**2)*(!seb**2)+(!b**2)*(!sea**2)). 
    COMPUTE sobelZ = product/sobse. 
    COMPUTE sobelp = 2*cdfnorm((-1)*abs(sobelZ)). 
    COMPUTE sobelres = {product, sobse, sobelZ, sobelp}. 
!enddefine. 

define DICHOT(modcount = !charend('/') /dat = !charend('/')). 
    COMPUTE dich = MAKE(!modcount, 3, -999). 
    LOOP q = 1 to !modcount. 
        COMPUTE uniqdes = DESIGN(!dat(:,q)). 
        COMPUTE dich(q,1) = (ncol(uniqdes) = 2). 
        DO IF (dich(q,1) = 1). 
            COMPUTE dichsort = MAKE(nrow(!dat), 1, -999). 
            COMPUTE dichgrad = GRADE(!dat(:,q)). 
            COMPUTE dichsort(dichgrad,1) = !dat(:,q). 
            COMPUTE dich(q,2) = dichsort(1,1). 
            COMPUTE dich(q,3) = dichsort(N,1). 
       END IF. 
   END LOOP.
!enddefine.

define PROBRES (coef = !charend('/') /values = !charend('/') /semat = !charend('/') /df = !charend('/')). 
    COMPUTE pbresmat = {!values, MAKE(nrow(!values), 6, -999)}. 
    COMPUTE values = {MAKE(nrow(!values),1, 1), !values}.
    COMPUTE pbresmat(:,ncol(!values)+1) = values*!coef. 
    COMPUTE pbresmat(:,ncol(!values)+2) = sqrt(diag(values*!semat*t(values))). 
    COMPUTE pbresmat(:,ncol(!values)+3) = pbresmat(:,ncol(!values)+1)/pbresmat(:,ncol(!values)+2).  
    COMPUTE pbresmat(:,ncol(!values)+4) = 2*(1-tcdf(abs(pbresmat(:,ncol(!values)+3)), !df)).
    CDFINVT p = temp/ df = !df. 
    COMPUTE tcritb = toutput.
    COMPUTE pbresmat(:,(ncol(!values)+5):(ncol(!values)+6)) = {pbresmat(:,ncol(!values)+1)-tcritb*pbresmat(:,ncol(!values)+2), pbresmat(:,ncol(!values)+1)+tcritb*pbresmat(:,ncol(!values)+2)}.
!enddefine. 

define CENTERD (centdat = !charend('/')). 
          COMPUTE centdat = !centdat. 
          DICHOT modcount = ncol(centdat) /dat = centdat. 
          COMPUTE avgs = csum(centdat)/nrow(centdat). 
          DO IF (ncol(centdat) > 1). 
              COMPUTE centmean = mdiag(avgs(1,1:(ncol(centdat)))&*t(1-dich(:,1)*(center=2))). 
          ELSE. 
              COMPUTE centmean = avgs(1,1:(ncol(centdat)))&*t(1-dich(:,1)*(center=2)). 
          END IF. 
          COMPUTE outdat = centdat - MAKE(nrow(centdat), ncol(centdat),1)*centmean. 
!enddefine. 

define JNprobe (coefone = !charend('/') /coeftwo = !charend('/') /seone = !charend('/') /setwo = !charend('/') /cov = !charend('/') /critt = !charend('/') /dfJN = !charend('/')). 
    COMPUTE coefone = !coefone. 
    COMPUTE coeftwo = !coeftwo. 
    COMPUTE seOne = !seOne. 
    COMPUTE seTwo = !seTwo. 
    COMPUTE cov = !cov. 
    COMPUTE critt = !critt. 
    COMPUTE dfJN = !dfJN. 
    COMPUTE cquad = (coefOne**2) - (critt**2)*seOne. 
    COMPUTE bquad = 2*coefOne*coefTwo - 2*(critt**2)*cov. 
    COMPUTE aquad = (coefTwo**2) - (critt**2)*seTwo. 
    DO IF ((bquad**2 - 4*cquad*aquad) >= 0). 
        COMPUTE JNsoln = {(-1*bquad + sqrt(bquad**2 - 4*cquad*aquad))/(2*aquad); (-1*bquad - sqrt(bquad**2 - 4*cquad*aquad))/(2*aquad)}. 
        COMPUTE Solngrad = grade(JNsoln).  
        COMPUTE JNsoln(Solngrad,1) = JNsoln(:,1).
        COMPUTE Pcntabv = csum(({moddat(:,1), moddat(:,1)} - make(N, 2, 1)*MDIAG(JNsoln)) > 0)/N*100. 
        COMPUTE NumJN = 2-rsum(Pcntabv = 100)-rsum(Pcntabv = 0). 
        COMPUTE Toohigh = rsum(Pcntabv = 0). 
        COMPUTE Toolow = rsum(Pcntabv = 100). 
        DO IF (Toolow = 1). 
            COMPUTE JNsoln = JNsoln(2,1). 
            COMPUTE Pcntabv = Pcntabv(2). 
        ELSE IF (Toohigh = 1). 
            COMPUTE JNsoln = JNsoln(1,1). 
            COMPUTE Pcntabv = Pcntabv(1). 
        END IF.  
    ELSE IF ((bquad**2 - 4*cquad*aquad) < 0). 
        COMPUTE NumJN = 0. 
    END IF.  
    DO IF (NumJN > 0). 
        COMPUTE JNWcomb = MAKE(20+numJN,2,1).
        COMPUTE MinW = MMIN(moddat(:,1)). 
        COMPUTE MaxW = MMAX(moddat(:,1)). 
        COMPUTE Range = MaxW - minW. 
        LOOP i = 1 TO 20. 
            COMPUTE JNWcomb(i,2) = MinW+(Range)/19*(i-1). 
        END LOOP.  
        COMPUTE JNWcomb(21:(20+numJN),2) = JNsoln.
        COMPUTE JNgrad = grade(JNWcomb(:,2)).
        COMPUTE JNWcomb(JNgrad,2) = JNWcomb(:,2). 
        COMPUTE JNres = {JNWcomb(:,2), MAKE(nrow(JNWcomb), 6, -999)}.
        COMPUTE JNres(:,2) = JNWcomb*{coefOne;coefTwo }. 
        COMPUTE JNres(:,3) = sqrt(diag(JNWcomb*{seOne , cov ; cov , seTwo }*t(JNWcomb))). 
        COMPUTE JNres(:,4) = JNres(:,2)/JNres(:,3).  
        COMPUTE JNres(:,5) = 2*(1-tcdf(abs(JNres(:,4)), dfJN )).
        COMPUTE JNres(:,6:7) = {JNres(:,2)-critt*JNres(:,3), JNres(:,2)+critt*JNres(:,3)}.
    END IF. 
*Outputs JNres: (20 row matrix with all probed points). 
*            JNsoln: list of all the JN solutions.
*            Pcntabv: percent of sample above each of the JN solutions.  
!enddefine. 

DEFINE MEMORE (Y = !charend('/') /M = !charend('/') !default(xxxxxxx)/W = !charend('/') !default(xxxxxxx)/Conf = !charend('/') !default(95) /mc = !charend('/') !default(0) 
   /samples = !charend('/') !default(5000) /normal = !charend('/') !default(0) /decimals=!charend('/') !default(F10.4) /save = !charend('/') !default(0) /bc = !charend('/') !default(0)
   /seed = !charend('/') !default(random) /contrast = !charend('/') !default(0) /xmint = !charend('/') !default(1) /serial = !charend('/') !default(0) /model = !charend('/') !default(1) /jn = !charend('/') !default(0) 
   /plot = !charend('/') !default(0) /quantile = !charend('/') !default(0) /center = !charend('/') !default(0) /wmodval1 = !charend('/') !default(999.99) /wmodval2 = !charend('/') !default(999.99)
   /wmodval3 = !charend('/') !default(999.99)). 
set mxloop = 100000000.
set seed = !seed. 
   
matrix. 
COMPUTE runnotes = MAKE(31,1,0). 
COMPUTE criterr = 0.  
COMPUTE model = !model. 
COMPUTE modelmt2 = {1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12; 13; 14; 15; 16; 17; 18}. 
COMPUTE validtst = (model = modelmt2). 
COMPUTE validmod = ANY(validtst). 

DO IF (validmod <> 1). 
   compute criterr = 1. 
   compute runnotes(23,1) = 23. 
END IF. 

COMPUTE w = !quote(!w). 
COMPUTE m = !quote(!m). 

DO IF ((w = "xxxxxxx") and (model <> 1)). 
   compute criterr = 1. 
   COMPUTE wnames = {"xxxxxxx"}.
   COMPUTE mnames = {"xxxxxxx"}.
   COMPUTE runnotes(21,1) = 21. 
END IF. 
DO IF ((m = "xxxxxxx") and ((model = 1) OR (model >= 4))). 
   compute criterr = 1. 
   COMPUTE mnames = {"xxxxxxx"}.
   COMPUTE wnames = {"xxxxxxx"}.
   compute runnotes(22,1) = 22. 
END IF. 

DO IF (criterr <> 1). 
   GET ydat / variables = !Y / names = ynames /missing = OMIT. 
   do if (model = 1).  
      GET mdat / variables = !M /names = mnames /missing = OMIT. 
      GET data / Variables = !M !Y  / Names = namevec /missing = OMIT. 
      GET fulldat / Variables = !M !Y /missing = 999. 
      COMPUTE wnames = {"xxxxxxx"}. 
   else if ((model = 2) OR (model = 3)). 
      GET wdat / variables = !W /names = wnames /missing = OMIT. 
      GET data / Variables = !W !Y  / Names = namevec /missing = OMIT. 
     GET fulldat / Variables = !W !Y /missing = 999. 
      COMPUTE mnames = {"xxxxxxx"}. 
   else if (model >= 4). 
      GET mdat / variables = !M /names = mnames /missing = OMIT. 
      GET wdat / variables = !W /names = wnames /missing = OMIT. 
      GET data / Variables = !W !M !Y  / Names = namevec /missing = OMIT. 
      GET fulldat / Variables = !W !M !Y /missing = 999. 
   end if. 
   COMPUTE missing = nrow(fulldat) - nrow(data). 
end if. 

COMPUTE mc = (!mc=1). 
COMPUTE serial = (!serial = 1). 
COMPUTE jn = (!jn = 1). 
COMPUTE plot = (!plot = 1). 
COMPUTE quantile = (!quantile = 1). 
COMPUTE contrast = (!contrast = 1). 
COMPUTE normal = (!normal = 1). 
COMPUTE xmint = (!xmint = 1). 
COMPUTE bc = (!bc = 1).

DO IF ((model = 1) AND (plot = 1)). 
    COMPUTE plot = 0. 
    COMPUTE runnotes(30,1)=30. 
END IF. 

COMPUTE center = !center. 
COMPUTE centervals = {0, 1, 2}. 
DO IF (1-any(center = centervals)). 
    COMPUTE runnotes(28,1) = 28. 
    COMPUTE center = 0. 
END IF. 

DO IF ((model >= 4) AND (serial = 1)). 
    COMPUTE serial = 0. 
    COMPUTE runnotes(25,1) = 25. 
END IF. 

DO IF (((model = 2) OR (model = 3)) AND (ncol(wnames) > 1) AND (jn = 1)). 
      COMPUTE jn = 0. 
      COMPUTE runnotes(16,1) = 16. 
END IF.  

DO IF ((model = 17) AND (xmint = 0)).
    COMPUTE criterr = 1. 
    COMPUTE runnotes(27,1) = 27.
END IF. 

COMPUTE apathmod = 0. 
COMPUTE bpathmod = 0. 
COMPUTE cppthmd = 0. 
COMPUTE dpathmod = 0. 

COMPUTE atvec = {4, 5, 6, 7, 9, 10, 12, 15}. 
COMPUTE btvec = {4, 5, 6, 8, 9, 11, 13, 16}. 
COMPUTE dtvec = {4, 5, 7, 8, 10, 11, 14, 17}.
COMPUTE cptvec = {4, 6, 7, 8, 12, 13, 14, 18}.  

DO IF (any(model = atvec)). 
    COMPUTE apathmod = 1. 
END IF. 

DO IF (any(model = btvec)).
    COMPUTE bpathmod = 1. 
END IF. 

DO IF (any(model = cptvec)).
    COMPUTE cppthmd = 1. 
END IF. 

DO IF (any(model = dtvec)). 
    COMPUTE dpathmod = 1. 
END IF. 

DO IF (xmint = 0).
    COMPUTE dpathmod = 0. 
END IF. 

COMPUTE anymod = apathmod+bpathmod+cppthmd+dpathmod. 
COMPUTE anymod = (anymod > 0). 

!let !toomany=0.
!do !i !in (!M).
  !do !j = 1 !to !length(!i).
    !if ((!j > 8) !and (!toomany = 0)) !then.
      compute criterr = 1.
      compute runnotes(11,1) = 11.
      !let !toomany = 1.
    !ifend.
  !doend.
!doend.
!do !i !in (!Y).
  !do !j = 1 !to !length(!i).
    !if ((!j > 8) !and (!toomany = 0)) !then.
      compute criterr = 1.
      compute runnotes(11,1) = 11.
      !let !toomany = 1.
    !ifend.
  !doend.
!doend.
!do !i !in (!W).
  !do !j = 1 !to !length(!i).
    !if ((!j > 8) !and (!toomany = 0)) !then.
      compute criterr = 1.
      compute runnotes(11,1) = 11.
      !let !toomany = 1.
    !ifend.
  !doend.
!doend.

DO IF (criterr <> 1). 
   DO IF (missing > 0). 
      COMPUTE runnotes(1,1) =1. 
   END IF. 
   DO IF (ncol(ydat) <> 2). 
      COMPUTE runnotes(2, 1) = 2.  
      COMPUTE criterr = 1. 
   END IF. 
END IF.
COMPUTE Mcount = ncol(mnames). 
COMPUTE Wcount = ncol(wnames). 
DO IF (any(wnames = "xxxxxxx")). 
    COMPUTE Wcount = 0. 
END IF. 

DO IF (((model = 2) OR (model = 3)) OR (anymod > 0)). 
   COMPUTE wmodval1 = !concat("{", !wmodval1, "}"). 
   COMPUTE wmodval2 = !concat("{", !wmodval2, "}").
   COMPUTE wmodval3 = !concat("{", !wmodval3, "}").
   COMPUTE yesval = MAKE(3,1, -999). 
   COMPUTE length = MAKE(3,1, -999).
   COMPUTE modvmat = MAKE(3, Wcount, -999). 
   LOOP i = 1 to 3.
      DO IF (i = 1). 
         COMPUTE wmodval = wmodval1. 
      ELSE IF (i = 2). 
         COMPUTE wmodval = wmodval2. 
      ELSE IF (i = 3). 
         COMPUTE wmodval = wmodval3. 
      END IF. 
      COMPUTE yesval(i,1) = ANY(wmodval - 999.99). 

      COMPUTE length(i,1) = (ncol(wmodval) = Wcount). 

      DO IF ((yesval(i,1) = 1)AND(length(i,1) = 1)). 
         COMPUTE modvmat(i,:) = wmodval. 
      ELSE IF ((yesval(i,1) = 1)AND(length(i,1) = 0)). 
         COMPUTE runnotes(18,1) = 18. 
         COMPUTE setswv = i-1. 
         BREAK. 
      ELSE. 
         COMPUTE setswv = i-1. 
         BREAK. 
      END IF.
      COMPUTE setswv = 3.  
   END LOOP. 

    DO IF (setswv > 0). 
       COMPUTE modvmat = modvmat(1:setswv, :). 
   ELSE. 
       COMPUTE modvmat = MAKE(1,1, -999). 
   END IF. 
END IF. 

DO IF (((Mcount = 0) OR (Mcount = 1) OR mod(Mcount,2) > 0) OR (Mcount >20)) AND (criterr <> 1) AND ((model =1) OR (model >= 4)). 
   COMPUTE runnotes(6, 1) = 6. 
   COMPUTE criterr = 1. 
END IF. 

DO IF (((model = 2) OR (model = 3)) AND ((Wcount = 0) OR (Wcount > 5))). 
   COMPUTE runnotes(15,1) = 15. 
   COMPUTE criterr = 1. 
END IF. 

DO IF ((model >= 4) AND (Wcount <> 1)). 
   COMPUTE runnotes(9,1) = 9. 
   COMPUTE criterr = 1. 
END IF. 
DO IF (criterr = 0). 

DO IF ((model = 1) OR (model >= 4)). 
   LOOP i = 1 TO Mcount. 
      DO IF ((mnames(1,i) = ynames(1,1)) OR (mnames(1,i) = ynames(1,2))). 
      COMPUTE runnotes(8,1) = 8. 
      COMPUTE criterr = 1. 
      END IF. 
   END LOOP. 
END IF. 

DO IF ((model = 2) OR (model = 3) OR (model >= 4)). 
   LOOP i = 1 TO Wcount. 
      DO IF ((wnames(1,i) = ynames(1,1)) OR (wnames(1,i) = ynames(1,2))). 
      COMPUTE runnotes(8,1) = 8. 
      COMPUTE criterr = 1. 
      END IF. 
   END LOOP. 
   END IF. 
END IF. 

DO IF ((model >= 4) AND (Wcount > 0)). 
    LOOP i = 1 TO Mcount. 
      DO IF (mnames(1,i) = wnames(1,1)). 
      COMPUTE runnotes(8,1) = 8. 
      COMPUTE criterr = 1. 
      END IF. 
   END LOOP. 
END IF. 

DO IF ((serial = 1) AND (mc = 1)). 
   COMPUTE runnotes(12,1) = 12. 
   COMPUTE mc = 0. 
END IF. 

DO IF (criterr <> 1). 
   COMPUTE zero = MAKE (nrow(data),1,0). 
   LOOP i = 1 TO (ncol(data)-1). 
      LOOP j = i+1 TO (ncol(data)). 
         COMPUTE diff = data(:,i) - data(:,j). 
         COMPUTE copy = csum(diff = zero). 
         DO IF (copy = nrow(data)). 
            COMPUTE copyname = {namevec(1,i), namevec(1,j)}. 
         BREAK. 
         END IF.  
      END LOOP. 
      DO IF (copy = nrow(data)). 
         BREAK. 
      END IF. 
   END LOOP. 

   DO IF (copy = nrow(data)). 
      COMPUTE runnotes(7,1) = 7. 
      COMPUTE criterr = 1. 
   END IF. 

END IF. 



DO IF (!samples = 0). 
   COMPUTE samples = 5000. 
   COMPUTE mc = 1. 
ELSE. 
    COMPUTE samples = abs(trunc(!samples))*(abs(trunc(!samples)) >= 1000) + 5000*(abs(trunc(!samples)) < 1000). 
END IF. 
DO IF (samples <> !samples). 
    COMPUTE runnotes(3, 1) = 3. 
END IF.

COMPUTE Conf = !Conf. 
DO IF (!Conf < 50 OR !Conf > 99.99). 
   COMPUTE Conf = 95. 
   COMPUTE runnotes (5, 1) = 5. 
END IF. 

DO IF (mc = 1 AND bc = 1).  
   COMPUTE runnotes(31,1) = 31. 
   COMPUTE bc = 0. 
END IF.  

DO IF (contrast = 1 AND (Mcount/2) = 1) AND (model = 1). 
   COMPUTE runnotes(10,1) = 10. 
END IF. 

DO IF ((contrast = 1) AND (serial = 1) AND (Mcount > 6) AND (model = 1)). 
   COMPUTE contrast = 0. 
   COMPUTE runnotes(20,1) = 20. 
END IF. 

DO IF ((contrast = 1) AND ((apathmod = 1) OR (bpathmod = 1))). 
    COMPUTE contrast = 0. 
    COMPUTE runnotes(26,1) = 26. 
END IF.  

DO IF (((model = 1) OR (model >= 4)) AND serial = 1) AND ((Mcount < 4) OR (Mcount > 10)). 
   COMPUTE runnotes(13,1) = 13.
   COMPUTE criterr = 1. 
END IF. 

DO IF ((model = 1) AND (center > 0)). 
   COMPUTE runnotes(24,1) = 24.
END IF. 


DO IF (criterr = 0). 

    COMPUTE N = nrow(data).

   DO IF (model = 3). 
      COMPUTE intcount = make(Wcount,1,-999). 
      LOOP i = 1 TO Wcount. 
         CHOOSE r = Wcount /k = i.   
         COMPUTE intcount(i,1) = rchoosek. 
      END LOOP. 
      COMPUTE Nfail = Wcount + csum(intcount)+2. 
   ELSE IF (Model = 2).
       COMPUTE Nfail = Wcount+2. 
   ELSE IF (Model = 1). 
      COMPUTE Nfail = mcount+2.
   ELSE IF (Model >= 4). 
      COMPUTE Nfail = 2*mcount+3.
   END IF. 
   DO IF (N < Nfail). 
      COMPUTE runnotes(19,1) = 19. 
      COMPUTE criterr = 1. 
   END IF. 
END IF. 

DO IF (criterr = 0). 
    DO IF ((model = 1) OR (model >= 4)). 
       COMPUTE Mpairs = Mcount/2.  
       COMPUTE mnamemat = reshape(mnames, Mpairs, 2).
       COMPUTE transmat = {1,1/2;-1,1/2}. 
       COMPUTE Tmat = make(ncol(data), ncol(data), 0).   
       DO IF (model >= 4). 
            COMPUTE Tmat(1:Wcount, 1:Wcount) = IDENT(Wcount). 
       END IF. 
       LOOP i = (1+Wcount) TO (2*Mpairs+1+Wcount) BY 2. 
          COMPUTE Tmat(i:(i+1), i:(i+1)) = transmat. 
       END LOOP. 
       COMPUTE dataT = data*Tmat. 
       COMPUTE select = {1,3,5,7,9,11, 13, 15, 17, 19, 21}.
       DO IF (model >= 4). 
           COMPUTE select = {1:Wcount, select + Wcount}. 
       END IF.  
       DO IF (xmint = 0). 
          COMPUTE dataT = dataT(: , select(1:(Mpairs+1+Wcount))).  
       ELSE IF (xmint = 1). 
          COMPUTE dataT = dataT(:,1:(ncol(dataT)-1)).  
          LOOP j = Wcount+2 TO (ncol(dataT)-1) by 2. 
               COMPUTE summean = csum(dataT(:,j))/N.  
               COMPUTE dataT(:,j) = (dataT(:,j) - summean).
          END LOOP. 
       END IF.

       DO IF (anymod = 1). 
           COMPUTE moddat = dataT(:,1:Wcount).
           DO IF (center > 0). 
               CENTERD centdat = moddat. 
               COMPUTE dataT(:,1:Wcount) = outdat. 
               COMPUTE moddat = outdat.
            END IF.
            LOOP i = (Wcount+1) to (ncol(dataT)-1) by (1+xmint). 
                DO IF (bpathmod = 1). 
                    COMPUTE moddat = {moddat, moddat(:,1)&*dataT(:,i)}.
                END IF. 
                DO IF (dpathmod = 1). 
                    COMPUTE moddat = {moddat, moddat(:,1)&*dataT(:,i+1)}.
                END IF. 
            END LOOP. 
       END IF. 

      
    ELSE IF (model = 2). 
       COMPUTE Mpairs = 0. 
       COMPUTE moddat = data(:,1:Wcount).
       COMPUTE tempvec = MAKE(Wcount,1, 0). 
       COMPUTE tempvec = {tempvec; 1; -1}. 
       COMPUTE transmat = {IDENT(Wcount+2, Wcount), tempvec}. 
       COMPUTE dataT = data*transmat. 
       DO IF (center > 0). 
          CENTERD centdat = moddat. 
          COMPUTE dataT(:,1:Wcount) = outdat. 
          COMPUTE moddat = outdat.
       END IF. 
    ELSE IF (model = 3). 
       COMPUTE Mpairs = 0. 
       COMPUTE moddat = data(:,1:Wcount).
       DO IF (center > 0). 
           CENTERD centdat = moddat. 
           COMPUTE moddat = outdat. 
       END IF. 
       DO IF (Wcount > 1). 
            LOOP h = 1 to Wcount - 1. 
             LOOP i = 1 TO Wcount-h. 
                CHOOSE r = Wcount-i /k = h. 
                LOOP j = csum(intcount(1:h,1))-rchoosek+1 TO csum(intcount(1:h,1)). 
                  COMPUTE moddat = {moddat, moddat(:,i)&*moddat(:,j)}.
                END LOOP. 
             END LOOP. 
            END LOOP.
       END IF. 
        COMPUTE tempvec = MAKE(ncol(moddat),1, 0). 
       COMPUTE tempvec = {tempvec; 1; -1}. 
       COMPUTE transmat = {IDENT(ncol(moddat)+2, ncol(moddat)), tempvec}. 
       COMPUTE data = {moddat, data(:,(ncol(data)-1):ncol(data))}.
       COMPUTE dataT = data*transmat. 
    END IF. 

    COMPUTE N = nrow(data).
    COMPUTE alpha = (1-.01*Conf). 
    COMPUTE temp = alpha/2. 

    DO IF ((model = 1) OR (model >= 4)).
          COMPUTE aresmat = MAKE (Mpairs*(1+apathmod), 6, -999). 
          COMPUTE ades = make(N,1, 1).
          DO IF (apathmod = 1). 
              COMPUTE ades = {ades, dataT(:,1:Wcount)}. 
              COMPUTE amodsum = make(mpairs, 7, -999).
          END IF. 
          CDFINVT p = temp/ df = (N-ncol(ades)). 
          COMPUTE tcrita = toutput. 
          COMPUTE tcritc = toutput.
          COMPUTE counterj = 1. 
          COMPUTE sem3aall = MAKE(ncol(ades)*mpairs, ncol(ades), -999). 
          COMPUTE aWcount = (apathmod = 1)*Wcount. 
          LOOP j = 1 TO Mpairs. 
             COMPUTE colj = (1+xmint)*j + Wcount - xmint. 
             COMPUTE rowj = j*(1+aWcount)-aWcount.
             DO IF (det(t(ades)*ades)=0). 
                 COMPUTE criterr = 1. 
                 COMPUTE runnotes(29,1) = 29. 
                 BREAK. 
             END IF. 
             COMPUTE avec = inv(t(ades)*ades)*t(ades)*dataT(:,colj). 
             COMPUTE M3pred = ades*avec. 
             COMPUTE M3ssr = csum((dataT(:,colj)-M3pred)&**2). 
             COMPUTE M3sst = csum((dataT(:,colj) - csum(dataT(:,colj))/N)&**2).  
             COMPUTE M3Rsq = 1-M3ssr/m3sst.  
             COMPUTE Rsqmat = {0, M3Rsq}. 
             COMPUTE M3r = sqrt(mmax(Rsqmat)).    
             COMPUTE M3df1 = ncol(ades) - 1.  
             COMPUTE M3df2 = (N - ncol(ades)). 
             COMPUTE M3msr = M3ssr/M3df2. 
            DO IF (apathmod = 1). 
                 COMPUTE M3F = m3df2*m3rsq/(m3df1*(1-m3rsq)).  
                 COMPUTE M3p = 1-FCDF(M3F, M3df1, m3df2). 
                 COMPUTE amodsum(j,:) =  {M3r, m3rsq, m3msr, m3F, m3df1, m3df2, m3p}. 
            END IF. 
             COMPUTE sem3amat = (m3msr*inv(t(ades)*ades)). 
             COMPUTE sem3aall(rowj:(rowj+aWcount), :) =  sem3amat.
             COMPUTE sem3a = (diag(sem3amat))&**(1/2).
             COMPUTE aresmat(rowj:(rowj+aWcount),1) = avec. 
             COMPUTE aresmat(rowj:(rowj+aWcount),2) = sem3a.  
             COMPUTE aresmat(rowj:(rowj+aWcount), 3) = aresmat(rowj:(rowj+aWcount),1) &/ aresmat(rowj:(rowj+aWcount), 2). 
             COMPUTE aresmat(rowj:(rowj+aWcount), 4) = 2*(1-tcdf(abs(aresmat(rowj:(rowj+aWcount),3)), m3df2)). 
             COMPUTE aresmat(rowj:(rowj+aWcount), 5) = aresmat(rowj:(rowj+aWcount), 1) - tcrita*aresmat(rowj:(rowj+aWcount),2). 
             COMPUTE aresmat(rowj:(rowj+aWcount), 6) = aresmat(rowj:(rowj+aWcount),1) + tcrita*aresmat(rowj:(rowj+aWcount),2). 
          END LOOP.  
          COMPUTE ccols = ncol(dataT).
          COMPUTE cols = 0. 
          COMPUTE cdes = make(N,1, 1).  
          DO IF (anymod > 0). 
              COMPUTE cdes = {cdes, dataT(:,1:Wcount)}. 
          END IF. 
          DO IF (det(t(cdes)*cdes)=0). 
                 COMPUTE criterr = 1. 
                 COMPUTE runnotes(29,1) = 29. 
                 COMPUTE cvec = MAKE(ncol(cdes), 1, -999).
          ELSE.  
              COMPUTE cvec = inv(t(cdes)*cdes)*t(cdes)*dataT(:,ccols). 
          END IF. 
          COMPUTE M4pred = cdes*cvec. 
          COMPUTE M4ssr = csum((dataT(:,ccols)-M4pred)&**2). 
          COMPUTE M4sst = csum((dataT(:,ccols) - csum(dataT(:,ccols))/N)&**2). 
          COMPUTE M4Rsq = 1-M4ssr/m4sst. 
          COMPUTE Rsqmat = {0, M4Rsq}.  
          COMPUTE M4r = sqrt(mmax(Rsqmat)).     
          COMPUTE M4df1 = ncol(cdes) - 1. 
          COMPUTE M4df2 = (N - ncol(cdes)). 
          COMPUTE M4msr = M4ssr/M4df2.   
          DO IF (anymod > 0). 
                 COMPUTE M4F = m4df2*m4rsq/(m4df1*(1-m4rsq)). 
                 COMPUTE M4p = 1-FCDF(M4F, M4df1, m4df2).
                 COMPUTE cmodsum =  {M4r, m4rsq, m4msr, m4F, m4df1, m4df2, m4p}. 
          END IF. 
          COMPUTE sem4cmat = (m4msr*inv(t(cdes)*cdes)). 
          COMPUTE sem4c = (diag(sem4cmat))&**(1/2).
          COMPUTE cresmat = MAKE(1+Wcount, 6, -999). 
          COMPUTE cresmat(1:(1+Wcount),1) = cvec.  
          COMPUTE cresmat(1:(1+Wcount),2) = sem4c. 
          COMPUTE cresmat(1:(1+Wcount), 3) = cresmat(1:(1+Wcount),1) &/ cresmat(1:(1+Wcount), 2). 
          COMPUTE cresmat(1:(1+Wcount), 4) = 2*(1-tcdf(abs(cresmat(1:(1+Wcount),3)), m4df2)). 
          COMPUTE cresmat(1:(1+Wcount), 5) = cresmat(1:(1+Wcount), 1) - tcritc*cresmat(1:(1+Wcount),2). 
          COMPUTE cresmat(1:(1+Wcount), 6) = cresmat(1:(1+Wcount),1) + tcritc*cresmat(1:(1+Wcount),2).    
          DO IF (serial = 1). 
                DO IF (xmint = 1). 
                   COMPUTE cols = mpairs**2 -1. 
                   COMPUTE serres = make((mpairs**2-1), 6, 0). 
                ELSE IF (xmint = 0). 
                   COMPUTE cols = ((mpairs - 1)**2 + 3*(mpairs-1))/2. 
                   COMPUTE serres = make(cols, 6, 0). 
                END IF. 
                COMPUTE smodsum = make(mpairs-1, 7, -999). 
                COMPUTE start = 1. 
                COMPUTE counterj = 1. 
              LOOP j = (2+xmint) TO (ncol(dataT)-1) by (1+xmint). 
                COMPUTE serdes = {make(N,1,1), dataT(:, 1:(j-1))}. 
                DO IF (det(t(serdes)*serdes)=0). 
                     COMPUTE criterr = 1. 
                     COMPUTE runnotes(29,1) = 29. 
                     COMPUTE M2modbs = MAKE(ncol(serdes), 1, -999).
                ELSE.  
                      COMPUTE M2modbs = inv(t(serdes)*serdes)*t(serdes)*dataT(:,j). 
                END IF.     
                COMPUTE M2pred = serdes*m2modbs. 
                COMPUTE M2ssr = csum((dataT(:,j)-M2pred)&**2). 
                COMPUTE M2sst = csum((dataT(:,j) - csum(dataT(:,j))/N)&**2). 
                COMPUTE M2Rsq = 1-M2ssr/m2sst. 
                COMPUTE Rsqmat = {0, M2Rsq}. 
                COMPUTE M2r = sqrt(mmax(Rsqmat)).      
                COMPUTE M2msr = M2ssr/(N - ncol(serdes)). 
                COMPUTE M2df1 = ncol(serdes) - 1. 
                COMPUTE M2df2 = (N - ncol(serdes)). 
                DO IF (m2rsq = 1). 
                    COMPUTE M2F = -999. 
                ELSE. 
                    COMPUTE M2F = m2df2*m2rsq/(m2df1*(1-m2rsq)). 
                END IF. 
                DO IF (det(t(serdes)*serdes)=0). 
                     COMPUTE sem2bmat = MAKE(ncol(serdes), ncol(serdes), 1). 
                ELSE.  
                     COMPUTE sem2bmat = (m2msr*inv(t(serdes)*serdes)). 
                END IF.  
                COMPUTE M2p = 1-FCDF(M2F, M2df1, m2df2). 
                COMPUTE sem2b = (diag(sem2bmat))&**(1/2). 
                COMPUTE smodsum(counterj,:) = {M2r, m2rsq, m2msr, m2F, m2df1, m2df2, m2p}. 
                DO IF (xmint = 1). 
                   COMPUTE end = (counterj+1)**2 - 1. 
                ELSE IF (xmint = 0).
                   COMPUTE end = (counterj**2 + 3*counterj)/2. 
                END IF. 
                CDFINVT p = temp /df = M2df2. 
                COMPUTE sercritt = toutput. 
                COMPUTE serres(start:end,1) = M2modbs. 
                COMPUTE serres(start:end, 2) = sem2b. 
                COMPUTE serres(start:end, 3) = serres(start:end, 1) &/ serres(start:end, 2). 
                COMPUTE serres(start:end, 4) = 2*(1-tcdf(abs(serres(start:end,3)), m2df2)). 
                COMPUTE serres(start:end, 5) = serres(start:end, 1) - sercritt*serres(start:end,2). 
                COMPUTE serres(start:end, 6) = serres(start:end, 1) + sercritt*serres(start:end,2). 
                COMPUTE aresmat(counterj+1,:) = serres(start,:). 
                COMPUTE start = end +1. 
                COMPUTE counterj = counterj + 1. 
             END LOOP.  
          END IF. 

    END IF.  

    /*for all models*/. 
    COMPUTE bcpdes = {make(N,1,1), dataT(:,(1+anymod):(ncol(dataT)-1))}. 
    DO IF (cppthmd = 1). 
        COMPUTE bcpdes = {make(N,1,1), dataT(:,1:(ncol(dataT)-1))}.
    END IF. 
    DO IF ((bpathmod = 1) OR (dpathmod = 1)). 
        COMPUTE bcpdes = {bcpdes, moddat(:,2:ncol(moddat))}. 
    END IF. 
    DO IF (det(t(bcpdes)*bcpdes)=0). 
        COMPUTE criterr = 1. 
        COMPUTE runnotes(29,1) = 29. 
        COMPUTE  bcpvec = MAKE(ncol(bcpdes), 1, -999).
    ELSE.  
        COMPUTE bcpvec = inv(t(bcpdes)*bcpdes)*t(bcpdes)*dataT(:,ncol(dataT)).
    END IF.  
    COMPUTE ypred = bcpdes*bcpvec. 
    COMPUTE ssr = csum((dataT(:,ncol(dataT)) - ypred)&**2). 
    COMPUTE sst = csum((dataT(:,ncol(dataT)) - csum(dataT(:,ncol(dataT)))/N)&**2).
    COMPUTE msr = ssr/(N-ncol(bcpdes)). 
    COMPUTE rsqfull = 1- ssr/sst. 
    COMPUTE Rsqmat = {0, rsqfull}. 
    COMPUTE rfull = sqrt(mmax(Rsqmat)).    
    COMPUTE df1 = (ncol(bcpdes)-1). 
    COMPUTE df2 = (N - ncol(bcpdes)). 
    COMPUTE Ffull = df2*Rsqfull/((df1)*(1-rsqfull)). 
    COMPUTE pfull =1- FCDF(Ffull, df1, df2). 
    DO IF (det(t(bcpdes)*bcpdes)=0). 
        COMPUTE sebcpmat = MAKE(ncol(bcpdes), ncol(bcpdes), 1). 
    ELSE.  
        COMPUTE sebcpmat = (msr*inv(t(bcpdes)*bcpdes)). 
    END IF.
    COMPUTE sebcp = (diag(sebcpmat))&**(1/2).  
    COMPUTE modsumr = {Rfull, Rsqfull, MSR, Ffull, df1, df2, pfull}. 
    
    CDFINVT p = temp/ df = df2. 
    COMPUTE tcritb = toutput. 
    COMPUTE tcritcp = toutput. 
    COMPUTE tcritd = toutput. 
    
    COMPUTE LCII = rnd((1-.01*Conf)/2*samples). 
    COMPUTE UCII = trunc((1-((1-.01*Conf)/2))*samples)+1. 
    DO IF ((LCII  < 1) OR (UCII > samples)). 
       COMPUTE runnotes(4, 1) = 4.  
       COMPUTE criterr = 1. 
       COMPUTE LCII = 1. 
       COMPUTE UCII = samples. 
    END IF. 

    DO IF ((Model <> 2) AND (Model <> 3) AND (CRITERR = 0)). 
          COMPUTE bresmat = MAKE(Mpairs*(bpathmod+1), 6, -999). 
          COMPUTE dresmat = MAKE(Mpairs*(dpathmod+1), 6, -999). 
          COMPUTE cpresmat = MAKE(cppthmd +1, 6, -999). 
          COMPUTE serind = 0. 
          LOOP i = 2 to Mpairs. 
                CHOOSE r = Mpairs /k = i. 
                COMPUTE serind = serind + rchoosek. 
          END LOOP. 
          COMPUTE indres = MAKE((Mpairs*(1+(apathmod = 1))*(1+(bpathmod = 1))+1+serind*(serial=1)),1,0). 
          DO IF (normal = 1). 
             COMPUTE normres = MAKE((Mpairs*(1+(apathmod = 1))*(1+(bpathmod = 1))+serind*(serial=1)), 4, 0). 
          END IF. 
          COMPUTE cpindx = 1+cppthmd.  
          COMPUTE cpresmat(1:cpindx,1) = bcpvec(1:cpindx,1). 
          COMPUTE cpresmat(1:cpindx,2) = sebcp(1:cpindx,1). 
          COMPUTE cpresmat(1:cpindx,3) = bcpvec(1:cpindx,1)/sebcp(1:cpindx,1).
          COMPUTE cpresmat(1:cpindx,4) = 2*(1-tcdf(abs(cpresmat(1:cpindx,3)), df2)).
          COMPUTE cpresmat(1:cpindx,5) = bcpvec(1:cpindx,1)-tcritcp*sebcp(1:cpindx,1). 
          COMPUTE cpresmat(1:cpindx,6) = bcpvec(1:cpindx,1)+tcritcp*sebcp(1:cpindx,1). 
          /*MC Setup*/. 
          DO IF (mc = 1). 
              COMPUTE mcsamps = samples.  
              COMPUTE randsamp = sqrt(-2*ln(uniform(mcsamps,Mpairs*(bpathmod+1))))&*cos((2*3.14159265358979)*uniform(mcsamps,Mpairs*(bpathmod+1))).
              COMPUTE MCres = MAKE((Mpairs*(1+(apathmod = 1))*(1+(bpathmod = 1))+1+serind*(serial=1)), 4, 0). 
              COMPUTE bindices = MAKE(Mpairs*(1+bpathmod), 1, 0). 
              LOOP i = 1 to Mpairs. 
                  COMPUTE bindices(i,1) = 2+cppthmd+(1+xmint)*(i-1). 
              END LOOP.  
              DO IF (bpathmod=1). 
                  COMPUTE seq = {1;2;3;4;5;6;7;8;9;10}. 
                  COMPUTE bindices((mpairs+1):(2*mpairs), 1) = seq(1:Mpairs,1)+(1 +cppthmd+Mpairs*(1+xmint)). 
              END IF. 
              COMPUTE MCcorr = sebcpmat(bindices,bindices). 
              COMPUTE rndnb = randsamp*chol(MCcorr).  
              COMPUTE rndna = sqrt(-2*ln(uniform(mcsamps,Mpairs*(1+apathmod))))&*cos((2*3.14159265358979)*uniform(mcsamps,Mpairs*(1+apathmod))).
              COMPUTE mcsave = make(samples, Mpairs*(2+apathmod+bpathmod), 0). 
              COMPUTE totsav = Mpairs*(2+apathmod+bpathmod). 
              COMPUTE mcind = make(samples, nrow(indres), 0). 
          END IF. 

          **Loop for simple effects. 
          COMPUTE counteri = 1. 
          LOOP i = 2+cppthmd TO (1+cppthmd+mpairs*(1+xmint)) by (1+xmint). 
             COMPUTE bpath = bcpvec(i,1). 
             COMPUTE sebpath = sebcp(i,1). 
             COMPUTE tbpath = bpath/sebpath. 
             COMPUTE pbpath = 2*(1-tcdf(abs(tbpath), df2)).
             COMPUTE lcib = bpath-tcritb*sebpath. 
             COMPUTE ucib = bpath+tcritb*sebpath. 
             COMPUTE bresmat(counteri, :) = {bpath, sebpath, tbpath, pbpath, lcib, ucib}. 
    
             DO IF ((xmint = 1) AND (counteri <= nrow(dresmat))). 
                COMPUTE dpath = bcpvec((i+1), 1). 
                COMPUTE sedpath = sebcp((i+1),1). 
                COMPUTE tdpath = dpath/sedpath. 
                COMPUTE pdpath = 2*(1-tcdf(abs(tdpath),df2)). 
                COMPUTE lcid = dpath-tcritd*sedpath. 
                COMPUTE ucid = dpath+tcritd*sedpath. 
                COMPUTE dresmat(counteri, :) = {dpath, sedpath, tdpath, pdpath, lcid, ucid}. 
             END IF. 

            /*Calculating Indirects in Sample*/.
             COMPUTE aindx = counteri*(1+apathmod)-apathmod -((1+apathmod)*mpairs)*(counteri > mpairs).
             COMPUTE amat = aresmat(aindx:(aindx+(apathmod = 1)),1). 
             COMPUTE indirect = amat*bresmat(counteri,1). 
             COMPUTE length = (1+apathmod)*(1+bpathmod). 
             COMPUTE indres((1+length/2*(counteri>mpairs) +(counteri-1-mpairs*(counteri > mpairs))*length):(1+length/2*(counteri>mpairs) +(counteri-1-mpairs*(counteri > mpairs))*length+apathmod),1) = indirect. 

             /* Normal theory tests*/. 
             DO IF (normal = 1). 
                SOBEL a = aresmat(aindx,1) /sea = aresmat(aindx,2) /b = bresmat(counteri,1) /seb = bresmat(counteri,2).
                COMPUTE normres(((apathmod+1)*counteri-apathmod),:) = sobelres. 
                DO IF (apathmod = 1). 
                    SOBEL a = aresmat(aindx+1,1) /sea = aresmat(aindx+1,2) /b = bresmat(counteri,1) /seb = bresmat(counteri,2).
                    COMPUTE normres(((apathmod+1)*counteri),:) = sobelres. 
                END IF. 
             END IF. 


            /*Monte Carlo Confidence Interval*/.
             DO IF (mc = 1). 
                 COMPUTE ones = MAKE(mcsamps, 1, 1). 
                 COMPUTE asamp = rndna(:,counteri:(counteri + apathmod))&*(ones*t(aresmat(aindx:(aindx+apathmod),2)))+(ones*t(aresmat(aindx:(aindx + apathmod),1)))). 
                 COMPUTE bsamp = rndnb(:,counteri)+bresmat(counteri,1).  
                 COMPUTE absamp = asamp&*(bsamp*MAKE(1,ncol(asamp),1)). 
                 COMPUTE mcsave(:,(1+(counteri-1)*(1+apathmod)):(counteri*(1+apathmod))) = asamp.  
                 COMPUTE mcsave(:,(1+Mpairs*(1+apathmod)+(counteri-1)*(bpathmod+1))) = bsamp.  
                 COMPUTE mcind(:,(1+length*(counteri-1)):(length*(counteri-1)+1+apathmod)) = absamp. 
                 LOOP j = 1 to ncol(absamp). 
                     COMPUTE mcgrad = grade(absamp(:,j)). 
                     COMPUTE mcsort = absamp(:,j). 
                     COMPUTE mcsort(mcgrad) = absamp(:,j).
                     COMPUTE MCLLCI = mcsort(LCII). 
                     COMPUTE MCULCI = mcsort(UCII). 
                     COMPUTE seMC = sqrt(csum((mcsort(:,1)-(csum(mcsort(:,1))/mcsamps))&**2)/(mcsamps-1)). 
                     COMPUTE MCres(j+(1+apathmod)*(counteri-1),:) = {indirect(j,1), seMC, MCLLCI, MCULCI}.  
                  END LOOP. 
             END IF. 
             COMPUTE counteri = counteri + 1. 
          END LOOP. 
          
          **Loop for interactions. 
         DO IF ((bpathmod = 1) OR (dpathmod = 1)). 
             LOOP i = (2+cppthmd+mpairs*(1+xmint)) TO nrow(bcpvec) by (bpathmod+dpathmod). 
                 DO IF (bpathmod = 1). 
                     COMPUTE bpath = bcpvec(i,1). 
                     COMPUTE sebpath = sebcp(i,1). 
                     COMPUTE tbpath = bpath/sebpath. 
                     COMPUTE pbpath = 2*(1-tcdf(abs(tbpath), df2)).
                     COMPUTE lcib = bpath-tcritb*sebpath. 
                     COMPUTE ucib = bpath+tcritb*sebpath. 
                     COMPUTE bresmat(counteri, :) = {bpath, sebpath, tbpath, pbpath, lcib, ucib}. 
                 END IF. 
                 
                 DO IF (dpathmod = 1). 
                    COMPUTE dpath = bcpvec((i+bpathmod), 1). 
                    COMPUTE sedpath = sebcp((i+bpathmod),1). 
                    COMPUTE tdpath = dpath/sedpath. 
                    COMPUTE pdpath = 2*(1-tcdf(abs(tdpath),df2)). 
                    COMPUTE lcid = dpath-tcritd*sedpath. 
                    COMPUTE ucid = dpath+tcritd*sedpath. 
                    COMPUTE dresmat(counteri, :) = {dpath, sedpath, tdpath, pdpath, lcid, ucid}. 
                 END IF. 
                DO IF (bpathmod = 1). 
                    /*Calculating Indirects in Sample*/.
                     COMPUTE aindx = counteri*(1+apathmod)-apathmod -((1+apathmod)*mpairs)*(counteri > mpairs).
                     COMPUTE amat = aresmat(aindx:(aindx+(apathmod = 1)),1). 
                     COMPUTE indirect = amat*bresmat(counteri,1). 
                     COMPUTE length = (1+apathmod)*(1+bpathmod). 
                     COMPUTE indres((1+length/2*(counteri>mpairs) +(counteri-1-mpairs*(counteri > mpairs))*length):(1+length/2*(counteri>mpairs) +(counteri-1-mpairs*(counteri > mpairs))*length+apathmod),1) = indirect. 
                     
                     /* Normal theory tests*/. 
                     DO IF (normal = 1). 
                        SOBEL a = aresmat(aindx,1) /sea = aresmat(aindx,2) /b = bresmat(counteri,1) /seb = bresmat(counteri,2).
                        COMPUTE normres(((apathmod+1)*counteri-apathmod),:) = sobelres. 
                        DO IF (apathmod = 1). 
                            SOBEL a = aresmat(aindx+1,1) /sea = aresmat(aindx+1,2) /b = bresmat(counteri,1) /seb = bresmat(counteri,2).
                            COMPUTE normres(((apathmod+1)*counteri),:) = sobelres. 
                        END IF. 
                     END IF. 
          
                    /*Monte Carlo Confidence Interval*/.
                    DO IF (mc = 1). 
                         COMPUTE ones = MAKE(mcsamps, 1, 1).  
                         COMPUTE asamp = rndna(:,(counteri-mpairs):(counteri-mpairs + apathmod))&*(ones*t(aresmat(aindx:(aindx+apathmod),2)))+(ones*t(aresmat(aindx:(aindx + apathmod),1)))).
                         COMPUTE bsamp = rndnb(:,counteri)+bresmat(counteri,1). 
                         COMPUTE absamp = asamp&*(bsamp*MAKE(1,ncol(asamp),1)). 
                         COMPUTE mcsave(:,(2+Mpairs*(1+apathmod)+(counteri-mpairs-1)*(bpathmod+1))) = bsamp. 
                         COMPUTE mcind(:,(1+length*(counteri-mpairs-1)+length/2):(length*(counteri-mpairs-1)+length)) = absamp. 
                         
                         LOOP j = 1 to ncol(absamp). 
                             COMPUTE mcgrad = grade(absamp(:,j)). 
                             COMPUTE mcsort = absamp(:,j). 
                             COMPUTE mcsort(mcgrad) = absamp(:,j).
                             COMPUTE MCLLCI = mcsort(LCII). 
                             COMPUTE MCULCI = mcsort(UCII). 
                             COMPUTE seMC = sqrt(csum((mcsort(:,1)-(csum(mcsort(:,1))/mcsamps))&**2)/(mcsamps-1)).
                             COMPUTE MCres(j+(1+apathmod)*(counteri-1),:) = {indirect(j,1), seMC, MCLLCI, MCULCI}.  
                          END LOOP. 
                     END IF.
                 END IF. 
                 COMPUTE counteri = counteri + 1. 
              END LOOP.   
          END IF.      
          
      /*serial mediation indirect paths*/.

        DO IF (serial = 1). 
          COMPUTE counter = 1. 
          COMPUTE indse = MAKE(serind*(serial=1), Mpairs+1, 1). 
          LOOP j = 1 TO Mpairs-1. 
             LOOP m = 1 TO Mpairs-j.
                COMPUTE step1 = aresmat(m,1).
                CHOOSE r = (Mpairs-m) /k = j. 
                COMPUTE indse(counter:(counter+rchoosek-1),1) = MAKE(rchoosek, 1, aresmat(m,2)**2).
                COMPUTE  indse(counter:(counter+rchoosek-1),2:(Mpairs+1)) = MAKE(rchoosek, Mpairs, aresmat(m,1)**2).
                LOOP l = m to (Mpairs - j).  
                   CHOOSE r = (Mpairs-l-1) /k = (j-1). 
                   DO IF (xmint = 1). 
                      COMPUTE srindx2 = (l**2-1)+2*m. 
                   ELSE IF (xmint = 0). 
                      COMPUTE srindx2 = l*(l+1)/2 + m. 
                   END IF.  
                   COMPUTE step2 = step1*serres(srindx2,1).
                   COMPUTE indse(counter:(counter+rchoosek-1),2) =   indse(counter:(counter+rchoosek-1),2)*(serres(srindx2,2)**2).
                   COMPUTE indse(counter:(counter+rchoosek-1),{1,3:(Mpairs+1)}) = indse(counter:(counter+rchoosek-1),{1,3:(Mpairs+1)})*(serres(srindx2,1)**2).
                   DO IF (j > 1).                  
                      LOOP h = (l+1) to Mpairs-j+1. 
                         CHOOSE r = (Mpairs - h -1) /k = (j-2).
                         DO IF (xmint = 1). 
                            COMPUTE srindx3 = (h**2-1)+2*(l+1). 
                         ELSE IF (xmint = 0). 
                            COMPUTE srindx3 = h*(h+1)/2 + l+1. 
                         END IF. 
                         COMPUTE step3 = step2*serres(srindx3,1).
                         COMPUTE indse(counter:(counter+rchoosek-1),3) =   indse(counter:(counter+rchoosek-1),3)*(serres(srindx3,2)**2).
                         COMPUTE indse(counter:(counter+rchoosek-1),{1:2,4:(Mpairs+1)}) = indse(counter:(counter+rchoosek-1),{1:2,4:(Mpairs+1)})*(serres(srindx3,1)**2).
                         DO IF (j > 2). 
                            LOOP i = h+1 to Mpairs-j+2. 
                               DO IF (xmint = 1). 
                                  COMPUTE srindx4 = ((i**2-1)+2*(h+1)). 
                               ELSE IF (xmint = 0). 
                                  COMPUTE srindx4 = i*(i+1)/2 + h+1. 
                               END IF. 
                            COMPUTE step4 = step3*serres(srindx4,1). 
                            COMPUTE indse(counter,4) =   indse(counter,4)*(serres(srindx4,2)**2).
                            COMPUTE indse(counter,{1:3,5:(Mpairs+1)}) = indse(counter,{1:3,5:(Mpairs+1)})*(serres(srindx4,1)**2).
                            DO IF (j > 3). 
                               DO IF (xmint = 1). 
                                  COMPUTE srindx5 = 23. 
                               ELSE IF (xmint = 0). 
                                  COMPUTE srindx5 = 14. 
                               END IF.  
                             COMPUTE step5 = step4*serres(srindx5,1).
                             COMPUTE indse(counter,5) =   indse(counter,5)*(serres(srindx5,2)**2).
                               COMPUTE indse(counter,{1:4,6:(Mpairs+1)}) = indse(counter,{1:4,6:(Mpairs+1)})*(serres(srindx5,1)**2).
                            COMPUTE indres(mpairs+counter,1) = step5*bresmat(5,1).
                            COMPUTE indse(counter, 6) = indse(counter, 6)*(bresmat(5,2)**2).
                            COMPUTE indse(counter,1:5) = indse(counter,1:5)*(bresmat(5,1)**2). 
                            COMPUTE counter = counter+1.
                            ELSE. 
                                 DO IF (((srindx4 < 4) AND (xmint = 1))OR((srindx4 < 3) AND (xmint = 0))).
                                     COMPUTE bindx = 2. 
                                  ELSE IF (((srindx4 < 9) AND (xmint = 1))OR((srindx4 < 6) AND (xmint = 0))).
                                     COMPUTE bindx = 3. 
                                  ELSE IF (((srindx4 < 16) AND (xmint = 1))OR((srindx4 < 10) AND (xmint = 0))).
                                     COMPUTE bindx = 4. 
                                  ELSE. 
                                     COMPUTE bindx = 5. 
                                  END IF. 
                                  COMPUTE indres(mpairs+counter,1) = step4*bresmat(bindx,1).
                                  COMPUTE indse(counter, 5) = indse(counter, 5)*(bresmat(bindx,2)**2).
                                  COMPUTE indse(counter,1:4) = indse(counter,1:4)*(bresmat(bindx,1)**2). 
                                  DO IF (Mpairs > 4). 
                                     COMPUTE indse(counter, 6) = 0.
                                  END IF.
                                  COMPUTE counter = counter+1.
                             END IF. 
                            END LOOP. 
                            /*end iloop */.
                        ELSE. 
                            DO IF (((srindx3 < 4) AND (xmint = 1))OR((srindx3 < 3) AND (xmint = 0))).
                                COMPUTE bindx = 2. 
                            ELSE IF (((srindx3 < 9) AND (xmint = 1))OR((srindx3 < 6) AND (xmint = 0))).
                                COMPUTE bindx = 3. 
                            ELSE IF (((srindx3 < 16) AND (xmint = 1))OR((srindx3 < 10) AND (xmint = 0))).
                                COMPUTE bindx = 4. 
                            ELSE. 
                                COMPUTE bindx = 5. 
                            END IF.   
                            COMPUTE indres(mpairs+counter,1) = step3*bresmat(bindx,1).
                            COMPUTE indse(counter, 4) = indse(counter, 4)*(bresmat(bindx,2)**2).
                            COMPUTE indse(counter,1:3) = indse(counter,1:3)*(bresmat(bindx,1)**2). 
                            DO IF (Mpairs > 3). 
                               COMPUTE indse(counter, 5:(Mpairs+1)) = MAKE(1, ncol({5:(Mpairs+1)}), 0).
                            END IF.
                            COMPUTE counter = counter+1.
                         END IF. 
                      END LOOP. 
                      /*end hloop */. 
                   ELSE. 
                   DO IF (((srindx2 < 4) AND (xmint = 1))OR((srindx2 < 3) AND (xmint = 0))).
                          COMPUTE bindx = 2. 
                   ELSE IF (((srindx2 < 9) AND (xmint = 1))OR((srindx2 < 6) AND (xmint = 0))).
                          COMPUTE bindx = 3. 
                   ELSE IF (((srindx2 < 16) AND (xmint = 1))OR((srindx2 < 10) AND (xmint = 0))).
                          COMPUTE bindx = 4. 
                    ELSE. 
                          COMPUTE bindx = 5. 
                    END IF. 
                      COMPUTE indres(mpairs+counter,1) = step2*bresmat(bindx,1).
                      COMPUTE indse(counter, 3) = indse(counter, 3)*(bresmat(bindx,2)**2). 
                      COMPUTE indse(counter,1:2) = indse(counter,1:2)*(bresmat(bindx,1)**2).
                      DO IF (Mpairs > 2). 
                         COMPUTE indse(counter, 4:(Mpairs+1)) = MAKE(1, ncol({4:(Mpairs+1)}), 0).
                      END IF. 
                      COMPUTE counter = counter+1.
                   END IF. 
                END LOOP. 
                /*end lloop*/. 
             END LOOP. 
             /*end mloop*/. 
          END LOOP. 
          /*end jloop*/. 


          DO IF (normal = 1). 
             COMPUTE serlind = indres((mpairs+1):(nrow(indres)-1)). 
             COMPUTE serialse = sqrt(rsum(indse)). 
             COMPUTE sobelZ = serlind&/serialse. 
             COMPUTE sobelp = 2*cdfnorm((-1)*abs(sobelZ)). 
             COMPUTE serindn = {serlind, serialse, sobelZ, sobelp}. 
             COMPUTE normres((mpairs+1):(mpairs+serind),:) = serindn.
           END IF. 
       
       END IF. 
     
       /*Total monte carlo*/.
        DO IF (mc = 1). 
          COMPUTE mcind(:,ncol(mcind)) = rsum(mcind). 
          COMPUTE mcsort = mcind(:,ncol(mcind)). 
          COMPUTE mcgrad = grade(mcsort). 
          COMPUTE mcsort(mcgrad) = mcsort. 
          COMPUTE MCLLCI = mcsort(LCII). 
          COMPUTE MCULCI = mcsort(UCII). 
          COMPUTE seMC = sqrt(csum((mcsort(:,1)-(csum(mcsort(:,1))/mcsamps))&**2)/(mcsamps-1)).
          COMPUTE MCres(nrow(MCres), :) = {csum(indres), seMC, MCLLCI, MCULCI}.
    
          DO IF ((contrast = 1) AND (Mpairs >1)). 
                COMPUTE npairs = Mpairs*(Mpairs-1)/2. 
                COMPUTE contres = MAKE(npairs, 4,0). 
                COMPUTE contsamp = MAKE(samples, npairs, 0). 
                COMPUTE contsort = contsamp. 
                COMPUTE counter = 1. 
                LOOP i = 1 TO Mpairs-1. 
                   LOOP j = i+1 TO Mpairs. 
                   COMPUTE contsamp(:,counter) = mcind(:,i) - mcind(:,j). 
                   COMPUTE contres(counter, 1) = indres(i,1) - indres(j,1). 
                   COMPUTE contres(counter, 2) = sqrt(csum((contsamp(:,counter)-(csum(contsamp(:,counter))/mcsamps))&**2)/(mcsamps-1)).
                   COMPUTE contgrad = grade(contsamp(:,counter)). 
                   COMPUTE contsort(contgrad,counter) = contsamp(:,counter). 
                   COMPUTE contres(counter, 3) = contsort(LCII, counter). 
                   COMPUTE contres(counter,4) = contsort(UCII, counter). 
                   COMPUTE counter = counter+1. 
                   END LOOP. 
                END LOOP. 
         END IF.
    END IF.
    
      COMPUTE badboot = 0. 
         /*Bootstrapping*/.
      DO IF (mc <> 1). 
         *Check the determinant for however many models there are. 
         COMPUTE detcheck = make(Mpairs+2, 1, -999). 
         *Bootsamp contains indirect effects (1, 2, or 4 depending on moderation) X Mpairs, 1 total indirect effect, serial indirect effects (if applicable).
         COMPUTE Bootsamp = make(samples, Mpairs*(1+(apathmod = 1))*(1+(bpathmod = 1))+1+serind*(serial=1), 0).  
         DO IF (serial = 1). 
               DO IF (xmint = 1). 
                  COMPUTE bootcol = mpairs**2-1. 
               ELSE IF (xmint = 0). 
                  COMPUTE bootcol = ((mpairs - 1)**2 + 3*(mpairs-1))/2. 
               END IF. 
         ELSE IF (serial = 0). 
            COMPUTE bootcol = -999. 
         END IF.  
         *Bootsave just contains path coefficients. 
         DO IF (serial = 0). 
             COMPUTE totsav = nrow(aresmat) + nrow(bcpvec) + nrow(cresmat). 
         ELSE IF (serial = 1). 
             COMPUTE totsav =  nrow(serres) +1 + nrow(bcpvec) + nrow(cresmat). 
         END IF. 
         COMPUTE Bootsave = make(samples, totsav,0). 
         LOOP i = 1 TO samples.
            LOOP k = 1 TO 10000. 
                COMPUTE counterj = 1. 
                COMPUTE sortvar = trunc(uniform(N,1)*N)+1.
                COMPUTE bootdat = dataT(sortvar(:,1),:).  
                DO IF (xmint = 1).
                  LOOP j = 1 TO Mpairs. 
                      COMPUTE summean = csum(bootdat(:,2*j+(1-(anymod=0))))/N.
                      COMPUTE bootdat(:,2*j+(1-(anymod=0))) = bootdat(:,2*j+(1-(anymod=0))) - summean.
                  END LOOP.
                END IF. 
                DO IF ((bpathmod = 1) OR (dpathmod = 1)). 
                       COMPUTE btmoddat = bootdat(:,1:Wcount).
                       DO IF (center > 0). 
                              CENTERD centdat = bootdat(:,1:Wcount). 
                              COMPUTE btmoddat = outdat. 
                        END IF.
                        LOOP l = (Wcount+1) to (ncol(bootdat)-1) by (1 + xmint). 
                            DO IF (bpathmod = 1). 
                                COMPUTE btmoddat = {btmoddat, btmoddat(:,1)&*bootdat(:,l)}.
                            END IF. 
                            DO IF (dpathmod = 1). 
                                COMPUTE btmoddat = {btmoddat, btmoddat(:,1)&*bootdat(:,l+1)}.
                            END IF. 
                        END LOOP. 
               END IF.
               
               *bootdes for b-model. 
               COMPUTE bootdes = {make(N,1,1), bootdat(:,(1+anymod):(ncol(bootdat)-1))}.
               DO IF (cppthmd = 1). 
                    COMPUTE bootdes = {make(N,1,1), bootdat(:,1:(ncol(bootdat)-1))}. 
               END IF.
               DO IF (bpathmod = 1) OR (dpathmod = 1). 
                    COMPUTE bootdes = {bootdes, btmoddat(:,2:ncol(btmoddat))}. 
               END IF. 
               COMPUTE detcheck(1,1) = (det(t(bootdes)*bootdes)=0).
               COMPUTE counterj = counterj+1.  
               COMPUTE bootcdes = make(N,1, 1).
               DO IF (anymod > 0). 
                    COMPUTE bootcdes = {bootcdes, bootdat(:,1:Wcount)}. 
               END IF.
               COMPUTE detcheck(counterj,1) = (det(t(bootcdes)*bootcdes)=0). 
                      COMPUTE counterj = counterj + 1. 
               LOOP j = 1 TO Mpairs. 
                      COMPUTE bootades = make(N,1, 1).
                      DO IF ((serial = 1) and (j > 1)).     
                          COMPUTE bootades = {make(N,1,1), bootdat(:, 1:((j-1)*(1+xmint)))}.
                      END IF. 
                      DO IF (apathmod = 1). 
                                COMPUTE bootades = {bootades, bootdat(:,1:Wcount)}. 
                      END IF.
                      COMPUTE detcheck(counterj,1) = (det(t(bootades)*bootades)=0). 
                      COMPUTE counterj = counterj + 1. 
               END LOOP.                                                                    
               COMPUTE badboot = badboot+(k = 2). 
           END LOOP IF (csum(detcheck(:,1)) = 0). 
           COMPUTE bootbeta = inv(t(bootdes)*bootdes)*t(bootdes)*bootdat(:,ncol(bootdat)). 
           COMPUTE bootavec = MAKE(mpairs*(1+apathmod), 1, -999).  
           COMPUTE bootbvec = MAKE(mpairs*(1+bpathmod),1,-999). 
           COMPUTE bootdvec = MAKE(mpairs*(1+dpathmod), 1, -999). 
           COMPUTE bootcvec = MAKE(2 - (anymod = 0), 1, -999).  
           COMPUTE bootcpvc = bootbeta(1:(1+cppthmd), 1). 
           DO IF (serial = 1). 
               COMPUTE bootser = MAKE(bootcol, 1, -999).
           END IF. 
           COMPUTE start = 1.  
           COMPUTE counterj = 1. 
           
           *Total effect bootstrap model. 
           COMPUTE bootcdes = make(N,1,1). 
           DO IF (anymod > 0). 
                  COMPUTE bootcdes = {bootcdes, bootdat(:,1:Wcount)}. 
           END IF.              
           COMPUTE bootcvec= inv(t(bootcdes)*bootcdes)*t(bootcdes)*bootdat(:,ncol(bootdat)).
        
           *a-path bootstrap model. 
           LOOP j = 1 TO Mpairs.  
              COMPUTE bootades = make(N,1, 1).
                      DO IF (apathmod = 1). 
                                COMPUTE bootades = {bootades, bootdat(:,1:Wcount)}. 
                      END IF. 
               COMPUTE colj = (1+xmint)*j + Wcount - xmint. 
               COMPUTE rowj = j*(1+aWcount)-aWcount.
               COMPUTE bootavec(rowj:(rowj+aWcount),1) = inv(t(bootades)*bootades)*t(bootades)*bootdat(:,colj).
               COMPUTE bootbvec(j,1) = bootbeta(2+cppthmd+(j-1)*(1+xmint), 1). 
               DO IF (xmint = 1). 
                  COMPUTE bootdvec(j,1) = bootbeta(3+cppthmd+(j-1)*(1+xmint), 1). 
 
                  DO IF (dpathmod = 1). 
                      COMPUTE bootdvec(j+mpairs,1) = bootbeta(2+bpathmod+mpairs*(1+xmint)+cppthmd+(j-1)*(bpathmod+xmint), 1). 
                  END IF. 
               END IF. 
               DO IF (bpathmod = 1).  
                  COMPUTE bootbvec(j+mpairs,1) = bootbeta(2+mpairs*(1+xmint)+cppthmd+(j-1)*(1+dpathmod), 1). 
               END IF.   
               
               COMPUTE bootb = bootbvec(j,1).
               DO IF (bpathmod = 1). 
                      COMPUTE bootb = {bootb, bootbvec(j+mpairs,1)}. 
               END IF. 
               DO IF ((serial = 1) AND (j>1)).
                  COMPUTE bootades = {make(N,1,1), bootdat(:, 1:((j-1)*(1+xmint)))}.
                  DO IF (xmint = 1).
                        COMPUTE end = (j)**2 - 1. 
                  ELSE IF (xmint = 0).
                        COMPUTE end = ((j-1)**2 + 3*(j-1))/2. 
                  END IF.
                  COMPUTE bootser(start:end,1) = inv(t(bootades)*bootades)*t(bootades)*bootdat(:,((j-1)*(1+xmint))+1)). 
                  COMPUTE bootavec(j,1) = bootser(start,1). 
                  COMPUTE start = end +1. 
               END IF. 
               COMPUTE aindx =  (j*(2-(apathmod = 0))-(apathmod=1)). 
               COMPUTE amat = bootavec(aindx:(aindx+(apathmod = 1))). 
               COMPUTE indmat = t(bootb)*t(amat). 
               COMPUTE indvec = RESHAPE(indmat, 1, (1+(apathmod = 1))*(1+(bpathmod = 1))).  
               COMPUTE bootsamp(i,counterj:(counterj-1+(1+(apathmod = 1))*(1+(bpathmod = 1)))) = indvec.
               COMPUTE counterj = counterj + (1+(apathmod = 1))*(1+(bpathmod = 1)). 
            
           END LOOP. 
           
            /*Serial loop goes here*/. 
              DO IF (serial = 1). 
               COMPUTE counter = 1. 
               LOOP j = 1 TO Mpairs-1. 
                  LOOP m = 1 TO Mpairs-j.
                     COMPUTE step1 = bootavec(m,1).
                     LOOP l = m to (Mpairs - j). 
                         DO IF (xmint = 1). 
                              COMPUTE srindx2 = (l**2-1)+2*m. 
                         ELSE IF (xmint = 0). 
                              COMPUTE srindx2 = l*(l+1)/2 + m. 
                         END IF.
                        COMPUTE step2 = step1*bootser(srindx2,1).
                        DO IF (j > 1).                  
                           LOOP h = (l+1) to Mpairs-j+1. 
                              DO IF (xmint = 1). 
                                 COMPUTE srindx3 = (h**2-1)+2*(l+1). 
                              ELSE IF (xmint = 0). 
                                 COMPUTE srindx3 = h*(h+1)/2 + l+1. 
                              END IF. 
                              COMPUTE step3 = step2*bootser(srindx3,1).
                              DO IF (j > 2). 
                                 LOOP o = h+1 to Mpairs-j+2. 
                                    DO IF (xmint = 1). 
                                       COMPUTE srindx4 = ((o**2-1)+2*(h+1)). 
                                    ELSE IF (xmint = 0). 
                                       COMPUTE srindx4 = o*(o+1)/2 + h+1. 
                                    END IF. 
                                 COMPUTE step4 = step3*bootser(srindx4,1). 
                                 DO IF (j > 3). 
                                    DO IF (xmint = 1). 
                                       COMPUTE srindx5 = 23. 
                                    ELSE IF (xmint = 0). 
                                       COMPUTE srindx5 = 14. 
                                    END IF.  
                                 COMPUTE step5 = step4*bootser(srindx5,1).
                                 COMPUTE bootsamp(i, mpairs + counter) = step5*bootbvec(5,1).
                                 COMPUTE counter = counter+1.
                                 ELSE. 
                                      DO IF (((srindx4 < 4) AND (xmint = 1))OR((srindx4 < 3) AND (xmint = 0))).
                                          COMPUTE bindx = 2. 
                                       ELSE IF (((srindx4 < 9) AND (xmint = 1))OR((srindx4 < 6) AND (xmint = 0))).
                                          COMPUTE bindx = 3. 
                                       ELSE IF (((srindx4 < 16) AND (xmint = 1))OR((srindx4 < 10) AND (xmint = 0))).
                                          COMPUTE bindx = 4. 
                                       ELSE. 
                                          COMPUTE bindx = 5. 
                                       END IF. 
                                      COMPUTE bootsamp(i, mpairs + counter) = step4*bootbvec(bindx,1).
                                     COMPUTE counter = counter+1.
                                  END IF. 
                                 END LOOP. 
                                /*end oloop*/. 
                             ELSE. 
                                 DO IF (((srindx3 < 4) AND (xmint = 1))OR((srindx3 < 3) AND (xmint = 0))).
                                     COMPUTE bindx = 2. 
                                 ELSE IF (((srindx3 < 9) AND (xmint = 1))OR((srindx3 < 6) AND (xmint = 0))).
                                     COMPUTE bindx = 3. 
                                 ELSE IF (((srindx3 < 16) AND (xmint = 1))OR((srindx3 < 10) AND (xmint = 0))).
                                     COMPUTE bindx = 4. 
                                 ELSE. 
                                     COMPUTE bindx = 5. 
                                 END IF.  
                                COMPUTE bootsamp(i, mpairs + counter) = step3*bootbvec(bindx,1).
                                COMPUTE counter = counter+1.
                              END IF. 
                           END LOOP. 
                           /*end hloop*/. 
                        ELSE. 
                           DO IF (((srindx2 < 4) AND (xmint = 1))OR((srindx2 < 3) AND (xmint = 0))).
                                  COMPUTE bindx = 2. 
                           ELSE IF (((srindx2 < 9) AND (xmint = 1))OR((srindx2 < 6) AND (xmint = 0))).
                                  COMPUTE bindx = 3. 
                           ELSE IF (((srindx2 < 16) AND (xmint = 1))OR((srindx2 < 10) AND (xmint = 0))).
                                  COMPUTE bindx = 4. 
                            ELSE. 
                                  COMPUTE bindx = 5. 
                            END IF. 
                           COMPUTE bootsamp(i, mpairs + counter) = step2*bootbvec(bindx,1).
                           COMPUTE counter = counter+1.  
                        END IF. 
                     END LOOP. 
                     /*end lloop*/. 
                  END LOOP. 
                  /*end mloop*/. 
               END LOOP. 
               /*end jloop*/. 
            END IF. 
            /*serial = 1*/.        
            
            /*Save bootstrap estimates*/. 
            DO IF ((serial = 0) AND (xmint = 1)). 
                COMPUTE bootsave(i, :) = {t(bootcvec), t(bootavec),t(bootcpvc), t(bootbvec), t(bootdvec)}.
            ELSE IF ((serial = 0) AND (xmint = 0)). 
                COMPUTE bootsave(i, :) = {t(bootcvec), t(bootavec),t(bootcpvc), t(bootbvec)}.
            ELSE IF ((serial = 1) AND (xmint = 1)). 
                COMPUTE bootsave(i, :) = {t(bootcvec), bootavec(1,1), t(bootser),t(bootcpvc), t(bootbvec), t(bootdvec)}.
            ELSE IF ((serial = 1) AND (xmint = 0)). 
                COMPUTE bootsave(i, :) = {t(bootcvec), bootavec(1,1), t(bootser),t(bootcpvc), t(bootbvec)}.
            END IF.             
            
     END LOOP. 
         DO IF (badboot >0). 
            COMPUTE runnotes(14,1) = 14. 
         END IF.  
         DO IF ((contrast = 1) AND (Mpairs >1)). 
            COMPUTE npairs = (ncol(bootsamp)-1)*(ncol(bootsamp)-2)/2. 
            COMPUTE contres = MAKE(npairs, 4,0). 
            COMPUTE contsamp = MAKE(samples, npairs, 0). 
            COMPUTE counter = 1. 
            LOOP i = 1 TO ncol(bootsamp)-2. 
               LOOP j = i+1 TO ncol(bootsamp)-1. 
               COMPUTE contsamp(:,counter) = bootsamp(:,i) - bootsamp(:,j). 
               COMPUTE contres(counter, 1) = indres(i,1) - indres(j,1). 
               COMPUTE counter = counter+1. 
               END LOOP. 
            END LOOP. 
         END IF. 

         COMPUTE bootsamp(:,ncol(bootsamp)) = rsum(bootsamp(:,1:(ncol(bootsamp)-1))).    
         COMPUTE indres(nrow(indres),1) = csum(indres).
         COMPUTE bootsort = bootsamp. 
         COMPUTE seboots = MAKE(nrow(indres), 1, 0).  
         COMPUTE bccires = MAKE(nrow(indres, 4, -999). 
         COMPUTE BootLLCI = MAKE(1,ncol(bootsamp),0). 
         COMPUTE BootULCI = MAKE(1,ncol(bootsamp),0). 
         COMPUTE zalpha2 = sqrt(-2*ln(alpha/2)).
         COMPUTE zalpha2 = (zalpha2+((((zalpha2*p4+p3)*zalpha2+p2)*zalpha2+p1)*zalpha2+p0)/((((zalpha2*q4+q3)*zalpha2+q2)*zalpha2+q1)*zalpha2+q0)).
         LOOP i = 1 TO ncol(bootsamp). 
            COMPUTE bootgrad = grade(bootsamp(:,i)). 
            COMPUTE bootsort(bootgrad,i) = bootsamp(:,i). 
            COMPUTE seboots(i,1) = sqrt(csum((bootsort(:,i)-(csum(bootsort(:,i))/samples))&**2)/(samples-1)).
            DO IF (bc = 1). 
                COMPUTE bccires(1,i) = csum(bootsamp(:,i)<indres(i,1))/samples.
                COMPUTE bccires(2,i) = bccires(1,i). 
                DO IF (bccires(1,i) > .5). 
                   COMPUTE bccires(2,i) = 1-bccires(1,i). 
                END IF. 
                COMPUTE bccires(3,i) = sqrt(-2*ln(bccires(2,i))). 
                COMPUTE bccires(4,i) = bccires(3,i)+((((bccires(3,i)*p4+p3)*bccires(3,i)+p2)*bccires(3,i)+p1)*bccires(3,i)+p0)/((((bccires(3,i)*q4+q3)*bccires(3,i)+q2)*bccires(3,i)+q1)*bccires(3,i)+q0).
                DO IF (bccires(1,i) <= .5). 
                   COMPUTE bccires(4,i) = -bccires(4,i). 
                END IF. 
                COMPUTE BCLLII = (cdfnorm(2*bccires(4,i)-zalpha2))*samples.
                COMPUTE BCUCII = (cdfnorm(2*bccires(4,i)+zalpha2))*samples.
                COMPUTE LCII = rnd(BCLLII). 
                COMPUTE UCII = trunc(BCUCII)+1. 
                DO IF (LCII < 1 OR UCII > samples). 
                   COMPUTE runnotes(4, 1) = 4.  
                   COMPUTE criterr = 1. 
                   COMPUTE LCII = 1. 
                   COMPUTE UCII = samples.
                END IF. 
                COMPUTE BootLLCI(1,i) = bootsort(LCII,i). 
                COMPUTE BootULCI(1,i) = bootsort(UCII,i). 
            END IF.  
         END LOOP. 
         DO IF (bc <>1). 
             COMPUTE BootLLCI = bootsort(LCII, :). 
             COMPUTE BootULCI = bootsort(UCII, :).
         END IF. 
         COMPUTE BootCI = {t(bootllci),t(bootulci)}. 
         COMPUTE bootres = {indres, seboots, bootci}.

         DO IF  (contrast = 1) AND (Mpairs >1). 
                 COMPUTE bccicont = MAKE(4,ncol(contsamp), 0).
                 COMPUTE contsort = contsamp. 
                 COMPUTE ContLLCI = MAKE(1,ncol(contsamp),0). 
                 COMPUTE ContULCI = MAKE(1,ncol(contsamp),0). 
                 LOOP i = 1 TO ncol(contsamp). 
                    COMPUTE contgrad = grade(contsamp(:,i)). 
                    COMPUTE contsort(contgrad,i) = contsamp(:,i). 
                    COMPUTE contres(i,2) = sqrt(csum((contsort(:,i)-(csum(contsort(:,i))/samples))&**2)/(samples-1)).  
                    DO IF (bc = 1). 
                       COMPUTE bccicont(1,i) = csum(contsamp(:,i)<contres(i,1))/samples. 
                       COMPUTE bccicont(2,i) = bccicont(1,i). 
                       DO IF (bccicont(1,i) > .5). 
                          COMPUTE bccicont(2,i) = 1-bccicont(1,i). 
                       END IF. 
                       COMPUTE bccicont(3,i) = sqrt(-2*ln(bccicont(2,i))). 
                       COMPUTE bccicont(4,i) = bccicont(3,i)+((((bccicont(3,i)*p4+p3)*bccicont(3,i)+p2)*bccicont(3,i)+p1)*bccicont(3,i)+p0)/((((bccicont(3,i)*q4+q3)*bccicont(3,i)+q2)*bccicont(3,i)+q1)*bccicont(3,i)+q0).
                       DO IF (bccicont(1,i) <= .5). 
                          COMPUTE bccicont(4,i) = -bccicont(4,i). 
                       END IF. 
                       COMPUTE CBCLLII = (cdfnorm(2*bccicont(4,i)-zalpha2))*samples.
                       COMPUTE CBCUCII = (cdfnorm(2*bccicont(4,i)+zalpha2))*samples.
                       COMPUTE LCII = rnd(CBCLLII). 
                       COMPUTE UCII = trunc(CBCUCII)+1. 
                       DO IF (LCII < 1 OR UCII > samples). 
                          COMPUTE runnotes(4, 1) = 4.  
                          COMPUTE criterr = 1. 
                          COMPUTE LCII = 1. 
                          COMPUTE UCII = samples.
                       END IF. 
                       COMPUTE ContLLCI(1,i) = contsort(LCII,i). 
                       COMPUTE ContULCI(1,i) = contsort(UCII,i).
                 END IF. 
                 END LOOP. 
                 DO IF (bc <>1). 
                     COMPUTE ContLLCI = contsort(LCII, :). 
                     COMPUTE ContULCI = contsort(UCII, :). 
                 END IF. 
                 COMPUTE ContCI = {t(contllci),t(contulci)}. 
                 COMPUTE contres(:,3:4) = contCI.  
         END IF. *DO IF  (contrast = 1) AND (Mpairs >1).
    END IF. 
    *close MC <> 1. 
END IF. 
*close Not 2 and Not 3.  

    DO IF ((((model = 2) OR (model = 3)) OR (anymod > 0)) AND (criterr = 0)). 
       DICHOT modcount = Wcount /dat = moddat. 
       DO IF ((quantile = 0) OR (plot = 1)). 
          COMPUTE modmeans = csum(moddat(:,1:Wcount)/N).
          DO IF (Wcount = 1). 
             COMPUTE dmeans = modmeans. 
          ELSE. 
             COMPUTE dmeans = MDIAG(modmeans). 
          END IF. 
          COMPUTE meansmat = make(N, Wcount, 1)*dmeans.
          COMPUTE modsds = sqrt(csum((moddat(:,1:Wcount)-meansmat)&**2)/(N-1)). 
          COMPUTE dimmc = 3**(Wcount - csum(dich(:,1)))*2**(csum(dich(:,1))). 
          COMPUTE modcomb = MAKE(dimmc, Wcount, -999). 
          COMPUTE counter = 1. 
          COMPUTE last = dimmc.
          LOOP i = 1 to Wcount. 
             DO IF (dich(i,1) = 0). 
                LOOP j = 1 to counter.  
                COMPUTE modcomb((last*(j-1)+1):(last*j),i) = {MAKE(last/3, 1, modmeans(1,i)-modsds(1,i));  MAKE(last/3, 1, modmeans(1,i)); MAKE(last/3, 1, modmeans(1,i)+modsds(1,i))}.
                END LOOP.
                COMPUTE last = last/3. 
                COMPUTE counter = counter*3. 
             ELSE. 
                LOOP j = 1 to counter.  
                 COMPUTE modcomb((last*(j-1)+1):(last*j),i) = {MAKE(last/2, 1, dich(i,2));  MAKE(last/2, 1, dich(i,3))}.
                END LOOP.
                COMPUTE last = last/2. 
                COMPUTE counter = counter*2. 
             END IF. 
          END LOOP. 
          COMPUTE plotdat = modcomb. 
       END IF. 
       
       DO IF (quantile = 1). 
          COMPUTE perctls = {.10, .25, .50, .75, .90}. 
          COMPUTE pctindx = rnd(perctls*N).  
          DO IF (pctindx(1,1)  < 1). 
             COMPUTE pctindx(1,1) = 1. 
          END IF. 
          DO IF (pctindx(1,5) > N). 
             COMPUTE pctindx(1,5) = N. 
          END IF. 
          COMPUTE modtemp = data(:,1:Wcount). 
          COMPUTE modsort = MAKE(N, Wcount, -999). 
          LOOP i = 1 TO Wcount. 
             COMPUTE modgrad = grade(modtemp(:,i)). 
             COMPUTE modsort(modgrad,i) = moddat(:,i). 
          END LOOP. 
          COMPUTE pctnum = modsort(pctindx, :).  
          COMPUTE dimmc = 5**(Wcount - csum(dich(:,1)))*2**(csum(dich(:,1))). 
          COMPUTE modcomb = MAKE(dimmc, Wcount, -999). 
          COMPUTE counter = 1. 
          COMPUTE last = dimmc. 
          LOOP i = 1 to Wcount. 
             DO IF (dich(i,1) = 0). 
                LOOP j = 1 to counter. 
                COMPUTE modcomb((last*(j-1)+1):(last*j),i)= {MAKE(last/5, 1, pctnum(1,i)); MAKE(last/5, 1, pctnum(2,i)); 
                                 MAKE(last/5, 1, pctnum(3,i)); MAKE(last/5, 1, pctnum(4,i)); MAKE(last/5, 1, pctnum(5,i))}.
                END LOOP.
                COMPUTE last = last/5. 
                COMPUTE counter = counter*5. 
             ELSE. 
                LOOP j = 1 to counter.  
                 COMPUTE modcomb((last*(j-1)+1):(last*j),i) = {MAKE(last/2, 1, dich(i,2));  MAKE(last/2, 1, dich(i,3))}.
                END LOOP.
                COMPUTE last = last/2. 
                COMPUTE counter = counter*2. 
             END IF. 
          END LOOP.
          COMPUTE plotdat = modcomb. 
       END IF. 
    END IF.      
        DO IF (anymod > 0).
            /*probe conditional effect of X on Y*/. 
            /*V3: with multiple moderators, need to account for interactions*/. 
            PROBRES coef = cvec /values = modcomb /semat = sem4cmat /df = m4df2. 
            COMPUTE XYgWres = pbresmat. 
            DO IF (setswv > 0). 
                PROBRES coef = cvec /values = modvmat /semat = sem4cmat /df = m4df2. 
                COMPUTE XYgWvres = pbresmat. 
            END IF. 
            
            DO IF (quantile = 0). 
                 COMPUTE plotdat = {XYgWres(:,1:(Wcount+1))}.
                 DO IF (apathmod = 1).
                     COMPUTE plotdat = {plotdat, MAKE(nrow(XYgWres), Mpairs, -999)}.
                 END IF. 
                 COMPUTE plotdes = modcomb.  
                 COMPUTE plotdes = {MAKE(nrow(plotdes),1, 1), plotdes}.
             ELSE. 
                 COMPUTE plotdes = plotdat. 
                 COMPUTE plotdes = {MAKE(nrow(plotdes),1, 1), plotdes}.
                 COMPUTE plotdat = {plotdat, MAKE(nrow(plotdat), 1+apathmod*Mpairs, -999)}.
                 COMPUTE plotdat(:,Wcount+1) = plotdes*cvec.
            END IF. 
            
            DO IF (jn = 1).  
                 JNprobe coefOne = cvec(1,1) /coefTwo = cvec(2,1) /seOne = sem4cmat(1,1) /seTwo = sem4cmat(2,2) /cov = sem4cmat(2,1) /critt = tcritc /dfJN = M4df2. 
                 compute cJNres = JNres. 
                 compute cJNsoln = JNsoln. 
                 compute cPcntabv = Pcntabv. 
             END IF. 
           
       END IF. 
        
        DO IF (apathmod = 1). 
            COMPUTE XMgWres = MAKE(dimmc*mpairs, Wcount+6, -999). 
            DO IF (setswv > 0). 
                COMPUTE XMgWvres = MAKE(setswv*mpairs, Wcount+6, -999). 
            END IF. 
            DO IF (jn = 1). 
                compute aJNres = MAKE(22*mpairs, 7, -999). 
                compute aJNsoln = MAKE(mpairs, 2, -999). 
                compute aPcntabv = MAKE(mpairs, 2, -999). 
                compute aNumJN = MAKE(mpairs, 1, -999). 
            END IF. 
            LOOP j = 1 to Mpairs. 
                 COMPUTE rowj = j*(1+Wcount)-Wcount.
                 COMPUTE avec = aresmat(rowj:(rowj+Wcount), 1). 
                 COMPUTE ase = sem3aall(rowj:(rowj+wcount), :). 
                 /*V3: for serial mediation need a vector of degrees of freedom*/.
                 PROBRES coef = avec /values = modcomb /semat = ase /df = m3df2. 
                 COMPUTE XMgWres((j*dimmc-(dimmc-1)):(j*dimmc), :) = pbresmat. 
                 DO IF (setswv > 0). 
                     PROBRES coef = avec /values = modvmat /semat = ase /df = m3df2. 
                     COMPUTE XMgWvres((j*setswv-(setswv-1)):(j*setswv), :) = pbresmat. 
                 END IF. 
                 DO IF (quantile = 0). 
                     COMPUTE plotdat(:,Wcount+1+j) = XMgWres((j*dimmc-(dimmc-1)):(j*dimmc), 2).
                 ELSE.  
                     COMPUTE plotdat(:,Wcount+1+j) = plotdes*avec.
                END IF. 
                DO IF (jn = 1).  
                     JNprobe coefOne = avec(1,1) /coefTwo = avec(2,1) /seOne = ase(1,1) /seTwo = ase(2,2) /cov = ase(2,1) /critt = tcrita /dfJN = M3df2.
                     compute aNumJN(j,:) = NumJN.
                     DO IF (NumJN > 0). 
                         compute aJNres((1+22*(j-1)):(20+NumJN+22*(j-1)),:) = JNres. 
                         compute aJNsoln(j,1:NumJN) = t(JNsoln). 
                         compute aPcntabv(j,1:NumJN) = Pcntabv. 
                     END IF.  
                 END IF.
            END LOOP. 
        END IF. 

       DO IF (bpathmod = 1). 
            *probe conditional effect of M on Y. 
            COMPUTE MYgWres = MAKE(dimmc*mpairs, Wcount+6, -999). 
            DO IF (setswv > 0). 
                COMPUTE MYgWvres = MAKE(setswv*mpairs, Wcount+6, -999). 
            END IF. 
            COMPUTE bplotdat = MAKE(dimmc*3*Mpairs, 3, -999).   
            DO IF (jn = 1). 
                compute bJNres = MAKE(22*mpairs, 7, -999). 
                compute bJNsoln = MAKE(mpairs, 2, -999). 
                compute bPcntabv = MAKE(mpairs, 2, -999). 
                compute bNumJN = MAKE(mpairs, 1, -999). 
            END IF. 
            LOOP j = 1 to Mpairs. 
                COMPUTE bvec = {bresmat(j, 1); bresmat( j + mpairs, 1)}. 
                COMPUTE firstloc = 2+cppthmd+(j-1)*(1+xmint). 
                COMPUTE secloc =2+cppthmd+mpairs*(xmint+1)+(j-1)*(1+dpathmod).
                COMPUTE bse = {sebcpmat(firstloc, firstloc), sebcpmat(firstloc, secloc); sebcpmat(secloc, firstloc), sebcpmat(secloc, secloc)}. 
                PROBRES coef = bvec /values = modcomb /semat = bse /df = df2.
                COMPUTE MYgWres((j*dimmc-(dimmc-1)):(j*dimmc), :) = pbresmat.  
                DO IF (setswv > 0). 
                     PROBRES coef = bvec /values = modvmat /semat = bse /df = df2.  
                     COMPUTE MYgWvres((j*setswv-(setswv-1)):(j*setswv), :) = pbresmat. 
                 END IF. 
                 COMPUTE bpvec = {bcpvec(1:(1+cppthmd),1); bvec}. 
                 COMPUTE Mdmean = csum(dataT(:,Wcount+j+(j-1)*xmint))/N. 
                 COMPUTE Mdsd = sqrt(csum((dataT(:,Wcount+j+(j-1)*xmint)-Mdmean)&**2)/(N-1)). 
                 COMPUTE Mdvec = {Mdmean-Mdsd; Mdmean; Mdmean+Mdsd}. 
                 COMPUTE longW = MAKE(3,1,plotdes(1,2)). 
                 LOOP k = 2 to dimmc.
                     COMPUTE longW = { longW ; MAKE(3,1,plotdes(k,2))}.
                 END LOOP. 
                 COMPUTE Mdvecrep = Mdvec.
                 LOOP k = 2 to dimmc. 
                     COMPUTE Mdvecrep = {Mdvecrep; Mdvec}. 
                 END LOOP. 
                 COMPUTE bplotdes = MAKE(dimmc*3,1,1). 
                 DO IF (cppthmd = 1).
                     COMPUTE bplotdes = {bplotdes, longW}.
                 END IF. 
                 COMPUTE bplotdes = {bplotdes, Mdvecrep}. 
                 COMPUTE bplotdes = {bplotdes, Mdvecrep&*longW}. 
                 COMPUTE bplotdat((1+(j-1)*dimmc*3):(j*dimmc*3), : ) = {longW, mdvecrep, bplotdes*bpvec}.  
                 DO IF (jn = 1).  
                     JNprobe coefOne = bvec(1,1) /coefTwo = bvec(2,1) /seOne = bse(1,1) /seTwo = bse(2,2) /cov = bse(2,1) /critt = tcritb /dfJN = df2.
                     compute bNumJN(j,:) = NumJN.
                     DO IF (NumJN > 0).  
                         compute bJNres((1+22*(j-1)):(20+NumJN+22*(j-1)),:) = JNres. 
                         compute bJNsoln(j,1:NumJN) = t(JNsoln). 
                         compute bPcntabv(j,1:NumJN) = Pcntabv. 
                     END IF.  
                 END IF. 
            END LOOP. 
        END IF.
        DO IF (cppthmd = 1). 
            DO IF (xmint = 1). 
                COMPUTE modcomb2 = modcomb. 
                DO IF (setswv > 0). 
                    COMPUTE modcomb2 = {modcomb;modvmat}. 
                END IF. 
                COMPUTE modcmbcp = {modcomb2, MAKE(dimmc+setswv, mpairs*(1+dpathmod), -999)}. 
                COMPUTE cpvec = {cpresmat(:,1); MAKE(mpairs*(1+dpathmod), 1, -999)}.  
                COMPUTE locindx = MAKE(mpairs*(1+dpathmod), 1, -999). 
                LOOP j = 1 to Mpairs. 
                        COMPUTE mmoddes = {MAKE(N, 1, 1), moddat(:, 1:Wcount)}. 
                        COMPUTE locindx(j) = 1+cppthmd+j*2. 
                        COMPUTE mmodvec = inv(t(mmoddes)*mmoddes)*t(mmoddes)*bcpdes(:,locindx(j)). 
                        COMPUTE mapred = mmoddes*mmodvec. 
                        COMPUTE ssr = csum((bcpdes(:,locindx(j)) - mapred)&**2). 
                        COMPUTE sst = csum((bcpdes(:,locindx(j)) - csum(bcpdes(:,locindx(j)))/N)&**2).
                        COMPUTE msr = ssr/(N-ncol(mmoddes)). 
                        COMPUTE mmoddf2 = (N - ncol(mmoddes)). 
                        COMPUTE semmodmt = (msr*inv(t(mmoddes)*mmoddes)). 
                        PROBRES coef = mmodvec /values = modcomb2 /semat = semmodmt /df = mmoddf2. 
                        COMPUTE modcmbcp(:,Wcount+j) = pbresmat(:,1+Wcount).
                        COMPUTE cpvec(2+j) = bcpvec(locindx(j)). 
                        DO IF (dpathmod = 1). 
                            COMPUTE locindx(j+mpairs) =2+cppthmd+mpairs*(xmint+1)+(j-1)*(1+bpathmod)+bpathmod. 
                            COMPUTE cpvec(2+j+mpairs) = bcpvec(locindx(j+mpairs)).
                            COMPUTE modcmbcp(:,Wcount+mpairs+j) = modcmbcp(:,1)&*modcmbcp(:,Wcount+j). 
                        END IF. 
                 END LOOP. 
                 COMPUTE locindx = {t(1:(1+Wcount)); locindx}.  
                 COMPUTE cpse = sebcpmat(locindx, locindx). 
            ELSE. 
                COMPUTE cpvec = cpresmat(:,1). 
                COMPUTE cpse = sebcpmat(1:(1+Wcount), 1:(1+Wcount)). 
                COMPUTE modcomb2 = modcomb. 
                DO IF (setswv > 0). 
                    COMPUTE modcomb2 = {modcomb;modvmat}. 
                END IF. 
                COMPUTE modcmbcp  = modcomb2. 
            END IF. 
            PROBRES coef = cpvec /values = modcmbcp /semat = cpse /df = df2. 
            DO IF (xmint = 1). 
                COMPUTE XYgWcMrs = {pbresmat(:,1:(Wcount+Mpairs)),pbresmat(:,(ncol(modcmbcp)+1):ncol(pbresmat))}. 
            ELSE IF (xmint =0). 
                COMPUTE XYgWcMrs = pbresmat.  
           END IF. 
            DO IF (setswv>0). 
                COMPUTE XYgWcMv =  XYgWcMrs((dimmc+1):(dimmc+setswv), :). 
                COMPUTE XYgWcMrs = XYgWcMrs(1:dimmc, :). 
            END IF. 
            COMPUTE cppltdat = {XYgWcMrs(:,1), XYgWcMrs(:,2+MPairs*xmint)}. 
            DO IF (jn = 1).  
                 JNprobe coefOne = cpvec(1,1) /coefTwo = cpvec(2,1) /seOne = cpse(1,1) /seTwo = cpse(2,2) /cov = cpse(2,1) /critt = tcritcp /dfJN = df2. 
                 compute cpJNres = JNres. 
                 compute cpJNsoln = JNsoln. 
                 compute cpPcntabv = Pcntabv. 
             END IF.
            
        END IF. 
       
        DO IF (dpathmod = 1). 
            COMPUTE MAYgWres = MAKE(dimmc*mpairs, Wcount+6, -999). 
            DO IF (setswv > 0). 
                COMPUTE MAYgWvrs = MAKE(setswv*mpairs, Wcount+6, -999). 
            END IF.
            COMPUTE dplotdat = MAKE(dimmc*3*Mpairs, 3, -999). 
             DO IF (jn = 1). 
                compute dJNres = MAKE(22*mpairs, 7, -999). 
                compute dJNsoln = MAKE(mpairs, 2, -999). 
                compute dPcntabv = MAKE(mpairs, 2, -999). 
                compute dNumJN = MAKE(mpairs, 1, -999). 
            END IF. 
            LOOP j = 1 to Mpairs. 
                COMPUTE dvec = {dresmat(j, 1); dresmat( j + mpairs, 1)}.
                COMPUTE firstloc = 3+cppthmd+(j-1)*(2). 
                COMPUTE secloc = 3+cppthmd+(mpairs-1)*(xmint+1)+(j)*(1+bpathmod).
                COMPUTE dse = {sebcpmat(firstloc, firstloc), sebcpmat(firstloc, secloc); sebcpmat(secloc, firstloc), sebcpmat(secloc, secloc)}. 
                PROBRES coef = dvec /values = modcomb /semat = dse /df = df2.
                COMPUTE MAYgWres((j*dimmc-(dimmc-1)):(j*dimmc), :) = pbresmat. 
                DO IF (setswv > 0). 
                     PROBRES coef = dvec /values = modvmat /semat = dse /df = df2. 
                     COMPUTE MAYgWvrs((j*setswv-(setswv-1)):(j*setswv), :) = pbresmat. 
                 END IF.        
                 COMPUTE dpvec = {bcpvec(1:(1+cppthmd),1); dvec}. 
                 COMPUTE MAmean = csum(dataT(:,Wcount+j+1+(j-1)*xmint))/N. 
                 COMPUTE MAsd = sqrt(csum((dataT(:,Wcount+j+1+(j-1)*xmint)-MAmean)&**2)/(N-1)).  
                 COMPUTE MAvec = {MAmean-MAsd; MAmean; MAmean+MAsd}. 
                 
                 COMPUTE longW = MAKE(3,1,plotdes(1,2)); 
                 LOOP k = 2 to dimmc. 
                     COMPUTE longW = { longW ; MAKE(3,1,plotdes(k,2))}.
                 END LOOP. 
                 COMPUTE MAvecrep = MAvec.
                 LOOP k = 2 to dimmc. 
                     COMPUTE MAvecrep = {MAvecrep; MAvec}. 
                 END LOOP. 
                 COMPUTE dplotdes = MAKE(dimmc*3,1,1). 
                 DO IF (cppthmd = 1).
                     COMPUTE dplotdes = {dplotdes, longW}.
                 END IF. 
                 COMPUTE dplotdes = {dplotdes, MAvecrep}. 
                 COMPUTE dplotdes = {dplotdes, MAvecrep&*longW}. 
                 COMPUTE dplotdat((1+(j-1)*dimmc*3):(j*dimmc*3), : ) = {longW, MAvecrep, dplotdes*dpvec}.  
                 DO IF (jn = 1).  
                     JNprobe coefOne = dvec(1,1) /coefTwo = dvec(2,1) /seOne = dse(1,1) /seTwo = dse(2,2) /cov = dse(2,1) /critt = tcritb /dfJN = df2.
                     compute dNumJN(j,:) = NumJN.
                     DO IF (NumJN > 0). 
                         compute dJNres((1+22*(j-1)):(20+NumJN+22*(j-1)),:) = JNres. 
                         compute dJNsoln(j,1:NumJN) = t(JNsoln). 
                         compute dPcntabv(j,1:NumJN) = Pcntabv. 
                     END IF.  
                 END IF. 
            END LOOP. 
        END IF. 
        
*Conditional Indirect Effects. 
DO IF (((apathmod = 1) OR (bpathmod = 1))) AND (criterr = 0). 
    COMPUTE immres = MAKE(mpairs, 4, -999). 
    COMPUTE cindres = MAKE((dimmc+setswv)*mpairs, 5, -999). 
    COMPUTE modmat = {MAKE(1, dimmc, 1); t(modcomb)}. 
    DO IF (setswv > 0). 
          COMPUTE modmat = {modmat, ({MAKE(1, setswv, 1); t(modvmat)})}. 
    END IF.
    COMPUTE cindres(:,1) = RESHAPE(MAKE(mpairs, 1, 1)*modmat(2,:), (dimmc+setswv)*mpairs, 1).  
        LOOP i = 1 to Mpairs. 
            DO IF (apathmod = 1). 
                    DO IF (mc = 0). 
                        COMPUTE asamps = bootsave(:, (1+2*i):(2+2*i)). 
                    ELSE IF (mc = 1). 
                        COMPUTE asamps = mcsave(:,(-1+2*i):(2*i)). 
                    END IF. 
                    COMPUTE casamps = asamps*modmat.                     
                    COMPUTE condas = {XMgWres((1+(i-1)*dimmc):(i*dimmc),2)}.
                    DO IF (setswv > 0). 
                        COMPUTE condas = {condas; XMgWvres((1+(i-1)*setswv):(i*setswv),2)}. 
                    END IF. 
            ELSE. 
                DO IF (mc=0). 
                    COMPUTE casamps = bootsave(:, 2+i).  
                ELSE IF (mc=1). 
                    COMPUTE casamps = mcsave(:,i). 
                END IF. 
                COMPUTE casamps = casamps*MAKE(1, dimmc+setswv, 1). 
                COMPUTE condas = MAKE((dimmc+setswv), 1, aresmat(i,1)). 
            END IF. 
            
            DO IF (bpathmod = 1). 
                DO IF (mc=0). 
                    COMPUTE bsamps = {bootsave(:,3+(1+apathmod)*mpairs+cppthmd+i), bootsave(:,3+(2+apathmod)*mpairs+cppthmd+i)}. 
                ELSE IF (mc=1). 
                    COMPUTE bsamps = mcsave(:,(-1+(1+apathmod)*mpairs+2*i):((1+apathmod)*mpairs+2*i)). 
                END IF. 
                COMPUTE cbsamps = bsamps*modmat. 
                COMPUTE condbs = {MYgWres((1+(i-1)*dimmc):(i*dimmc),2)}.
                DO IF (setswv > 0). 
                        COMPUTE condbs = {condbs; MYgWvres((1+(i-1)*setswv):(i*setswv),2)}. 
                END IF. 
             ELSE. 
                 DO IF (mc=0). 
                     COMPUTE cbsamps = bootsave(:, 3+(1+apathmod)*mpairs+cppthmd+i)*MAKE(1, dimmc+setswv, 1).
                 ELSE IF (mc=1). 
                     COMPUTE cbsamps = mcsave(:, (1+apathmod)*mpairs+i). 
                      COMPUTE cbsamps = cbsamps*MAKE(1, dimmc+setswv, 1).  
                 END IF.   
                 COMPUTE condbs = MAKE((dimmc+setswv), 1, bresmat(i,1)). 
             END IF.  
            COMPUTE cindres((1+(i-1)*(dimmc+setswv)):(i*(dimmc+setswv)),2) = condas&*condbs. 
            COMPUTE cabsamps = casamps&*cbsamps. 
            COMPUTE scabsamps = cabsamps. 
            DO IF (dimmc = 2).  
                COMPUTE immsamp = {cabsamps(:,2) - cabsamps(:,1)}. 
                COMPUTE immres(i,1) = condas(2,1)*condbs(2,1) - condas(1,1)*condbs(1,1). 
                COMPUTE sampgrad = grade(immsamp).  
                COMPUTE immsamp(sampgrad) = immsamp. 
                COMPUTE immres(i,2) = sqrt(csum((immsamp-(csum(immsamp)/samples))&**2)/(samples-1)).
                COMPUTE immres(i,3) = immsamp(LCII, :).  
                COMPUTE immres(i,4) = immsamp(UCII, :).
            END IF. 
            LOOP j = 1 TO ncol(cabsamps). 
                COMPUTE sampgrad = grade(cabsamps(:,j)). 
                COMPUTE scabsamps(sampgrad,j) = cabsamps(:,j). 
                COMPUTE cindres((i-1)*(dimmc+setswv)+j,3) = sqrt(csum((scabsamps(:,j)-(csum(scabsamps(:,j))/samples))&**2)/(samples-1)).
             END LOOP. 
             COMPUTE cindres((1+(i-1)*(dimmc+setswv)):(i*(dimmc+setswv)),4) = t(scabsamps(LCII, :)). 
             COMPUTE cindres((1+(i-1)*(dimmc+setswv)):(i*(dimmc+setswv)),5) = t(scabsamps(UCII, :)).
        END LOOP. 
END IF. 

DO IF ((model = 2) OR (model = 3)) AND (criterr = 0). 
   COMPUTE wnamemat = {t(wnames), MAKE(Wcount, 1, " ")}. 
   COMPUTE modres = MAKE(ncol(bcpdes), 6, -999). 
   COMPUTE modres(:,1) = bcpvec. 
   COMPUTE modres(:,2) = sebcp. 
   COMPUTE modres(:,3) = modres(:,1)/modres(:,2).  
   COMPUTE modres(:,4) = 2*(1-tcdf(abs(modres(:,3)), df2)).
   COMPUTE modres(:,5:6) = {modres(:,1)-tcritb*modres(:,2), modres(:,1)+tcritb*modres(:,2)}. 



   /*Conditional effect of X on Y at values of W*/. 

   COMPUTE XYgWres = {modcomb, MAKE(dimmc, 6, -999)}. 
      DO IF (model = 3). 
         LOOP h = 1 to Wcount - 1. 
            LOOP i = 1 TO Wcount-h. 
               CHOOSE r = Wcount-i /k = h. 
               LOOP j = csum(intcount(1:h,1))-rchoosek+1 TO csum(intcount(1:h,1)). 
                 COMPUTE modcomb = {modcomb, modcomb(:,i)&*modcomb(:,j)}.
                 DO IF (setswv > 0). 
                 COMPUTE modvmat = {modvmat, modvmat(:,i)&*modvmat(:,j)}. 
                 END IF. 
               END LOOP. 
            END LOOP. 
          END LOOP.
      END IF. 
      PROBRES coef = bcpvec /values = modcomb /semat = sebcpmat /df = df2. 
      COMPUTE XYgWres = pbresmat. 

      DO IF (setswv > 0). 
      COMPUTE wvres = MAKE(setswv, Wcount+6, -999).
      COMPUTE wvres(:,1:Wcount) = modvmat(:,1:Wcount). 
      COMPUTE modvmat= {MAKE(setswv,1,1), modvmat}. 
      COMPUTE wvres(:,Wcount+1) = modvmat*bcpvec.
      COMPUTE wvres(:,Wcount+2) = sqrt(diag(modvmat*sebcpmat*t(modvmat))). 
      COMPUTE wvres(:,Wcount+3) = wvres(:,Wcount+1)/wvres(:,Wcount+2).  
      COMPUTE wvres(:,Wcount+4) = 2*(1-tcdf(abs(wvres(:,Wcount+3)), df2)).
      COMPUTE wvres(:,(Wcount+5):(Wcount+6)) = {wvres(:,Wcount+1)-tcritb*wvres(:,Wcount+2), wvres(:,Wcount+1)+tcritb*wvres(:,Wcount+2)}.
      END IF. 
 

      DO IF (csum(dich(:,1)) > 0) AND (jn = 1). 
         COMPUTE jn = 0. 
         COMPUTE runnotes(17,1) = 17. 
      END IF. 

      DO IF (jn = 1).  
         COMPUTE cquad = (bcpvec(1,1)**2) - (tcritb**2)*sebcpmat(1,1). 
         COMPUTE bquad = 2*bcpvec(1,1)*bcpvec(2,1) - 2*(tcritb**2)*sebcpmat(1,2). 
         COMPUTE aquad = (bcpvec(2,1)**2) - (tcritb**2)*sebcpmat(2,2). 
         DO IF ((bquad**2 - 4*cquad*aquad) >= 0). 
            COMPUTE JNsoln = {(-1*bquad + sqrt(bquad**2 - 4*cquad*aquad))/(2*aquad); (-1*bquad - sqrt(bquad**2 - 4*cquad*aquad))/(2*aquad)}. 
            COMPUTE Solngrad = grade(JNsoln). 
            COMPUTE JNsoln(Solngrad,1) = JNsoln(:,1).
            COMPUTE pcntabv = csum(({moddat, moddat} - make(N, 2, 1)*MDIAG(JNsoln)) > 0)/N*100. 
            COMPUTE numJN = 2-rsum(pcntabv = 100)-rsum(pcntabv = 0). 
            COMPUTE toohigh = rsum(pcntabv = 0). 
            COMPUTE toolow = rsum(pcntabv = 100). 
            DO IF (toolow = 1). 
               COMPUTE JNsoln = JNsoln(2,1). 
               COMPUTE pcntabv = pcntabv(2). 
            ELSE IF (toohigh = 1). 
               COMPUTE JNsoln = JNsoln(1,1). 
               COMPUTE pcntabv = pcntabv(1). 
            END IF.  
         ELSE IF ((bquad**2 - 4*cquad*aquad) < 0). 
            COMPUTE numJN = 0. 
         END IF.  
         DO IF (numJN > 0). 
            COMPUTE JNWcomb = MAKE(20+numJN,2,1).
            COMPUTE MinW = MMIN(moddat). 
            COMPUTE MaxW = MMAX(moddat). 
            COMPUTE Range = MaxW - minW. 
            LOOP i = 1 TO 20. 
               COMPUTE JNWcomb(i,2) = MinW+(Range)/19*(i-1). 
            END LOOP. 
            COMPUTE JNWcomb(21:(20+numJN),2) = JNsoln. 
            COMPUTE JNgrad = grade(JNWcomb(:,2)). 
            COMPUTE JNWcomb(JNgrad,2) = JNWcomb(:,2). 
            COMPUTE JNres = {JNWcomb(:,2), MAKE(nrow(JNWcomb), 6, -999)}.
            COMPUTE JNres(:,2) = JNWcomb*bcpvec. 
            COMPUTE JNres(:,3) = sqrt(diag(JNWcomb*sebcpmat*t(JNWcomb))). 
            COMPUTE JNres(:,4) = JNres(:,2)/JNres(:,3).  
            COMPUTE JNres(:,5) = 2*(1-tcdf(abs(JNres(:,4)), df2)).
            COMPUTE JNres(:,6:7) = {JNres(:,2)-tcritb*JNres(:,3), JNres(:,2)+tcritb*JNres(:,3)}.
         END IF. 
      END IF. 

 
   
   /*Probing Effect of Ws on Y*/. 
   COMPUTE prbmsum = MAKE(2,7, -999). 
   COMPUTE prbmdres = MAKE(2*ncol(bcpdes), 6, -999). 
   LOOP i = 1 TO 2. 
         COMPUTE probevec = inv(t(bcpdes)*bcpdes)*t(bcpdes)*data(:,ncol(data)-2+i).
         COMPUTE prbypred = bcpdes*probevec. 
         COMPUTE probessr = csum((data(:,ncol(data)-2+i) - prbypred)&**2). 
         COMPUTE probesst = csum((data(:,ncol(data)-2+i) - csum(data(:,ncol(data)-2+i))/N)&**2).
         COMPUTE probemsr = probessr/(N-ncol(bcpdes)). 
         COMPUTE prsqfull = 1- probessr/probesst. 
         COMPUTE Rsqmat = {0, prsqfull}. 
         COMPUTE prbrfull = sqrt(mmax(Rsqmat)).    
         COMPUTE probedf1 = (ncol(bcpdes)-1). 
         COMPUTE probedf2 = (N - ncol(bcpdes)). 
         COMPUTE prbFfull = probedf2*pRsqfull/((probedf1)*(1-prsqfull)). 
         COMPUTE prbpfull =1- FCDF(prbFfull, probedf1, probedf2). 
         COMPUTE seprbmat = (probemsr*inv(t(bcpdes)*bcpdes)).
         COMPUTE seprb = (diag(seprbmat))&**(1/2). 
         COMPUTE prbmsum(i,:) = {prbrfull, prsqfull, probemsr, prbFfull, probedf1, probedf2, prbpfull}. 
         
         COMPUTE prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),1) = probevec. 
         COMPUTE prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),2) = seprb. 
         COMPUTE prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),3) = prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),1)/prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),2).  
         COMPUTE prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),4) = 2*(1-tcdf(abs(prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),3)), df2)).
         COMPUTE prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),5) = prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),1)-tcritb*prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),2).
         COMPUTE prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),6) = prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),1)+tcritb*prbmdres((ncol(bcpdes)*(i-1)+1):(ncol(bcpdes)*i),2).
         
   END LOOP. 

    DO IF (plot = 1). 
          DO IF (quantile = 0). 
             COMPUTE plotdat = {XYgWres(:,1:(Wcount+1)), MAKE(nrow(XYgWres), 2, -999)}.
             COMPUTE plotdes = modcomb.  
             COMPUTE plotdes = {MAKE(nrow(plotdes),1, 1), plotdes}.
          ELSE. 
             COMPUTE plotdes = plotdat. 
             COMPUTE plotdat = {plotdat, MAKE(nrow(plotdat), 3, -999)}.
             DO IF (model = 3). 
                LOOP h = 1 to Wcount - 1. 
                   LOOP i = 1 TO Wcount-h. 
                      CHOOSE r = Wcount-i /k = h. 
                      LOOP j = csum(intcount(1:h,1))-rchoosek+1 TO csum(intcount(1:h,1)). 
                        COMPUTE plotdes = {plotdes, plotdes(:,i)&*plotdes(:,j)}.
                      END LOOP. 
                   END LOOP. 
                 END LOOP.
             END IF. 
             COMPUTE plotdes = {MAKE(nrow(plotdes),1, 1), plotdes}. 
             COMPUTE plotdat(:,Wcount+1) = plotdes*bcpvec.
         END IF. 
    
       /*Plot of Y for each Condition by Ws*/.  
       COMPUTE plotdat(:,Wcount+2) = plotdes*prbmdres(1:ncol(bcpdes),1).
       COMPUTE plotdat(:,Wcount+3) = plotdes*prbmdres((ncol(bcpdes)+1):(2*ncol(bcpdes)),1).
    END IF. 

END IF. 



END IF. 

print /title = "*********************** MEMORE Procedure for SPSS Version 3.0 ************************".
print /title = "                           Written by Amanda Montoya       ".
print /title = "                    Documentation available at akmontoya.com ".
print /title = "**************************************************************************************".

DO IF (criterr = 0).

    print model /title = "Model:". 
    COMPUTE varrlabs = {'Y = '}. 
    DO IF ((Model <> 1) AND (Wcount = 1)). 
        COMPUTE varrlabs = {varrlabs, 'W = '}. 
    ELSE IF ((Model <> 1) AND (Wcount > 1)). 
        COMPUTE varrlabs = {varrlabs, 'W1 = ', 'W2 = ', 'W3 = ', 'W4 = ', 'W5 = '}.
    END IF. 
    DO IF (((Model = 1) OR (Model >= 4)) AND (Mpairs = 1)). 
        COMPUTE varrlabs = {varrlabs, 'M = '}. 
    ELSE IF (((Model = 1) OR (Model >= 4)) AND (Mpairs > 1)). 
        COMPUTE varrlabs = {varrlabs, 'M1 = ', 'M2 = ', 'M3 = ', 'M4 = ', 'M5 = ', 'M6 = ', 'M7 = ', 'M8 = ', 'M9 = ', 'M10 = '}.
        COMPUTE temp9 = {"(M1)", "(M2)", "(M3)", "(M4)", "(M5)", "(M6)", "(M7)", "(M8)", "(M9)", "(M10)"}. 
    END IF. 
    
    DO IF (model = 1). 
        print {ynames; mnamemat} /title = "Variables: " /rnames = varrlabs  /format = a8.  
    ELSE IF ((Model = 2) OR (Model = 3)). 
        print {ynames; wnamemat} /title = "Variables: " /rnames = varrlabs  /format = a8. 
    END IF. 
    DO IF (Model >= 4). 
        COMPUTE wnamemat = {t(wnames), MAKE(Wcount, 1, " ")}. 
        print {ynames; wnamemat; mnamemat} /title = "Variables: " /rnames = varrlabs  /format = a8. 
    END IF. 

    DO IF ((model = 1) OR (model >= 4)). 
          COMPUTE compname = MAKE(((1+xmint+dpathmod+bpathmod)*Mpairs+1),9, 0).   
          COMPUTE compname(1,:) = {' ', ynames(1,1), ' - ', ynames(1,2), ' ', ' ', ' ', ' ', ' '}.
          LOOP j = 1 TO Mpairs.
             COMPUTE compname((1+j),:) ={' ', mnamemat(j,1), ' - ', mnamemat(j,2), ' ', ' ', ' ', ' ', ' '}. 
               DO IF (xmint = 1). 
                   COMPUTE compname((Mpairs+1+j),:) = { '(', mnamemat(j,1), ' + ', mnamemat(j,2), ')', '/2', ' ', ' ', 'Centered'}. 
               END IF. 
          END LOOP. 
          COMPUTE temp1 = {'M1diff = ','M2diff = ','M3diff = ', 'M4diff = ', 'M5diff = ','M6diff = ','M7diff = ','M8diff = ', 'M9diff = ', 'M10diff = '}. 
          COMPUTE temp2 = {'M1avg  = ','M2avg  = ','M3avg  = ', 'M4avg  = ', 'M5avg  = ','M6avg  = ','M7avg  = ','M8avg  = ', 'M9avg  = ', 'M10avg  = '}.
          COMPUTE temp6 = {'M1diff','M2diff','M3diff', 'M4diff', 'M5diff','M6diff','M7diff','M8diff', 'M9diff', 'M10diff'}. 
          COMPUTE temp7 = {'M1avg','M2avg','M3avg', 'M4avg', 'M5avg','M6avg','M7avg','M8avg', 'M9avg', 'M10avg'}.
          DO IF (xmint = 1). 
                DO IF (Mpairs = 1). 
                   COMPUTE temprnam = {'Ydiff = ', 'Mdiff = ', 'Mavg = '}. 
                ELSE. 
                   COMPUTE temprnam = {'Ydiff = ', temp1(1,1:Mpairs), temp2(1,1:Mpairs)}. 
                END IF. 
          ELSE IF (xmint = 0). 
                DO IF (Mpairs = 1). 
                   COMPUTE temprnam = {'Ydiff = ', 'Mdiff = '}. 
                ELSE. 
                   COMPUTE temprnam = {'Ydiff = ', temp1(1,1:Mpairs)}. 
                END IF. 
         END IF. 
         COMPUTE temp3 = {'Int1  = ', 'Int2  = ', 'Int3  = ', 'Int4  = ', 'Int5  = ', 'Int6  = ', 'Int7  = ', 'Int8  = ', 'Int9  = ', 'Int10 = '}.
         COMPUTE temp4 = {'Int11 = ', 'Int12 = ', 'Int13 = ', 'Int14 = ', 'Int15 = ', 'Int16 = ', 'Int17 = ', 'Int18 = ', 'Int19 = ', 'Int20 = '}.
         COMPUTE temp5 = {'Int21 = ', 'Int22 = ', 'Int23 = ', 'Int24 = ', 'Int25 = ', 'Int26 = ', 'Int27 = ', 'Int28 = ', 'Int29 = ', 'Int30 = '}.
         COMPUTE temprnam = {temprnam, temp3, temp4, temp5}.  
         COMPUTE counter = (1+xmint)*Mpairs + 2. 
         DO IF (bpathmod = 1). 
             LOOP j = 1 TO mpairs.  
                 COMPUTE compname(counter,:) = {wnamemat(1,1), '*', '(', mnamemat(j,1), ' - ', mnamemat(j,2), ')', ' ', ' '}. 
                 COMPUTE counter = counter + 1. 
             END LOOP. 
         END IF.
         DO IF (dpathmod = 1). 
             LOOP j = 1 TO mpairs. 
                 COMPUTE compname(counter,:) = {wnamemat(1,1), '*', '(', mnamemat(j,1), ' + ', mnamemat(j,2), ')', '/2', 'Centered'}.
                 COMPUTE counter = counter + 1. 
             END LOOP. 
         END IF.  
    END IF. 

    DO IF ((model = 2) OR ((model = 3) and (Wcount = 1))). 
          COMPUTE compname = MAKE(1,7, " "). 
          COMPUTE compname(1,1:4) = {' ', ynames(1,1), ' - ', ynames(1,2)}.
          COMPUTE temprnam = {'Ydiff = '}. 
    ELSE IF ((model = 3) AND (Wcount > 1)). 
          COMPUTE nint = ncol(moddat) - Wcount.
          COMPUTE compname = MAKE(1+nint, 10, " ").
          COMPUTE compname(1,1:4) = {' ', ynames(1,1), ' - ', ynames(1,2)}. 
          COMPUTE temp1 = {'Int1  = ', 'Int2  = ', 'Int3  = ', 'Int4  = ', 'Int5  = ', 'Int6  = ', 'Int7  = ', 'Int8  = ', 'Int9  = ', 'Int10 = '}.
          COMPUTE temp2 = {'Int11 = ', 'Int12 = ', 'Int13 = ', 'Int14 = ', 'Int15 = ', 'Int16 = ', 'Int17 = ', 'Int18 = ', 'Int19 = ', 'Int20 = '}.
          COMPUTE temp3 = {'Int21 = ', 'Int22 = ', 'Int23 = ', 'Int24 = ', 'Int25 = ', 'Int26 = ', 'Int27 = ', 'Int28 = ', 'Int29 = ', 'Int30 = '}.
          COMPUTE temprnam = {'Ydiff = ', temp1, temp2, temp3}.
          COMPUTE intnames = MAKE(csum(intcount), 9, " "). 
          COMPUTE intnames(1:Wcount, 1) = t(wnames). 
          COMPUTE counter = 1. 
             LOOP h = 1 to Wcount - 1. 
                   LOOP i = 1 TO Wcount-h. 
                      CHOOSE r = Wcount-i /k = h. 
                      LOOP j = csum(intcount(1:h,1))-rchoosek+1 TO csum(intcount(1:h,1)). 
                        COMPUTE intnames(Wcount+counter,1:((2*h)+1)) = {intnames(i,1), " x ", intnames(j,1:(2*h-1))}. 
                         COMPUTE counter = counter +1. 
                      END LOOP. 
                   END LOOP. 
                 END LOOP. 
          COMPUTE compname(2:(nint+1), 2:10) = intnames((Wcount+1):(Wcount+nint), :). 
    END IF. 
    print compname /title = "Computed Variables:" /rnames = temprnam /format = a8. 
    
    
    print N /title = "Sample Size:". 
    do if (!quote(!seed) <> "random"). 
    print !seed /title = "Seed:". 
    end if. 
    DO IF ((model = 1) OR (model >= 4)). 
       print {"Ydiff =" , ynames(1,1), ' - ', ynames(1,2)} /title = "**************************************************************************************" /rlabels = "Outcome:" /format = A8.
       COMPUTE collab = {"Effect", "SE", "t", "p", "LLCI", "ULCI"}. 
        DO IF (anymod=1). 
            print cmodsum /title = "Model Summary" /clabels = "R", "R-sq", "MSE", "F" , "df1" , "df2", "p" /format = !decimals. 
        END IF. 
       print cresmat /title = "Model" /rnames = {"constant", "W"} /cnames = collab /format = !decimals. 
       print M4df2 /title = "Degrees of freedom for all regression coefficient estimates:".
       DO IF (anymod > 0). 
                   COMPUTE coeflabs = {"Effect" , "SE", "t", "p", "LLCI", "ULCI"}.
                   COMPUTE XYgWlabs = {wnames, coeflabs}.  
                   COMPUTE condnam = MAKE(3, 3, "                  "). 
                   COMPUTE condnam(:,1) = {"Focal:"; "Outcome:"; "Mod:"}. 
                   COMPUTE condnam(:, 2) = {"'X'"; "Ydiff"; wnamemat(1,1)}. 
                   COMPUTE condnam(:, 3) = {"(X)"; "(Y)"; "(W)"}. 
                   print condnam /title = "Conditional Effect of Focal Predictor on Outcome at values of Moderator(s)" /format = A8. 
                   print XYgWres /cnames = XYgWlabs /format = !decimals /title = " ". 
                DO IF (csum(dich(:,1)) <> Wcount). 
                    DO IF (quantile = 1). 
                       print /title = "     Values for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentile.".
                    ELSE IF (quantile = 0). 
                       print /title = "     Values for quantitative moderators are the mean and plus/minus one SD from the mean.".
                    END IF. 
                END IF. 
                DO IF (csum(dich(:,1)) > 0). 
                   print /title = "     Values for dichotomous moderators are the two values of the moderator.".
                END IF. 
                DO IF (setswv >0).  
                print condnam /title = "Conditional Effect of Focal Predictor on Outcome at requested values of Moderator(s)" /format = A8. 
                print XYgWvres /title = " " /cnames = XYgWlabs /format = !decimals.  
                END IF.    
                print M4df2 /title = "Degrees of freedom for all conditional effects:".
       END IF. 
       DO IF (jn = 1). 
            print /title = "--------------------------------------------------------------------------------------".
            print /title = "                    JOHNSON-NEYMAN PROCEDURE: Total Effect Model". 
            DO IF (ncol(cJNsoln) <> 0). 
            print /title = "Moderator value(s) defining Johnson-Neyman significance region(s) and percent of ".
            print {cJNsoln, t(cpcntabv)} /title "observed data above value:" /clabels = "Value", "% Abv" /format = !decimals /space = 0. 
            print cJNRes /title = "Conditional Total Effect of 'X' on Y at values of moderator" /cnames = XYgWlabs /format = !decimals.
            print M4df2 /title = "Degrees of freedom for all conditional effects:".
            ELSE IF (ncol(cJNsoln) = 0). 
            print /title = "There are no statistically significant transition points within the observed range of data.". 
            END IF. 
        END IF. 

   COMPUTE alabs = {"a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10"}.
   COMPUTE start = 1. 
   LOOP j = 1 to Mpairs. 
      DO IF (Mpairs = 1). 
         print {"Mdiff = " , mnamemat(j,1), ' - ', mnamemat(j,2)} /title = "**************************************************************************************" /rlabels = "Outcome:" /format = A8.
      ELSE. 
         print {temp1(1,j) , mnamemat(j,1), ' - ', mnamemat(j,2)} /title = "**************************************************************************************" /rlabels = "Outcome:" /format = A8.
      END IF. 
      DO IF ((serial = 1) AND (j > 1)). 
         print smodsum(j-1,:) /title = "Model Summary" /clabels = "R", "R-sq", "MSE", "F" , "df1" , "df2", "p" /format = !decimals. 
         DO IF (xmint = 1). 
            COMPUTE m2labs = {"constant", "M1diff", "M1avg", "M2diff", "M2avg", "M3diff", "M3avg", "M4diff", "M4avg"}. 
            COMPUTE end = j**2 - 1. 
        ELSE IF (xmint = 0). 
            COMPUTE m2labs = {"constant", "M1diff", "M2diff", "M3diff", "M4diff"}. 
            COMPUTE end = ((j-1)**2 + 3*(j-1))/2.  
         END IF. 
         print serres(start:end,:) /title "Model" /rnames = m2labs /clabels = "coeff" , "SE", "t", "p", "LLCI", "ULCI" /format = !decimals.
         COMPUTE start = end+1. 
         print smodsum(j-1,6) /title = "Degrees of freedom for all regression coefficient estimates:".
      ELSE. 
          DO IF (apathmod = 1). 
              print amodsum(j,:) /title = "Model Summary" /clabels = "R", "R-sq", "MSE", "F" , "df1" , "df2", "p" /format = !decimals. 
              COMPUTE rowj = j*(1+Wcount)-Wcount.
              print aresmat(rowj:(rowj+Wcount),:) /title = "Model" /rnames = {"constant", "W"} /cnames = collab /format = !decimals. 
          ELSE. 
              print aresmat(j,:) /title = "Model" /rnames = {"constant"} /cnames = collab /format = !decimals. 
          END IF. 

          print M3df2 /title = "Degrees of freedom for all regression coefficient estimates:".
          DO IF (apathmod = 1). 
               COMPUTE coeflabs = {"Effect" , "SE", "t", "p", "LLCI", "ULCI"}.
               COMPUTE XMgWlabs = {wnames, coeflabs}.  
               DO IF (mpairs = 1). 
                   COMPUTE condnam(:, 2) = {"'X'"; "Mdiff"; wnamemat(1,1)}. 
               ELSE. 
                   COMPUTE condnam(:, 2) = {"'X'"; temp6(1,j); wnamemat(1,1)}. 
               END IF. 
               COMPUTE condnam(:, 3) = {"(X)"; "(M)"; "(W)"}. 
               print condnam /title = "Conditional Effect of Focal Predictor on Outcome at values of Moderator(s)" /format = A8. 
               print XMgWres((j*dimmc-(dimmc-1)):(j*dimmc), :) /title = " " /cnames = XYgWlabs /format = !decimals. 
                DO IF (csum(dich(:,1)) <> Wcount). 
                    DO IF (quantile = 1). 
                       print /title = "     Values for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentile.".
                    ELSE IF (quantile = 0). 
                       print /title = "     Values for quantitative moderators are the mean and plus/minus one SD from the mean.".
                    END IF. 
                END IF. 
                DO IF (csum(dich(:,1)) > 0). 
                   print /title = "     Values for dichotomous moderators are the two values of the moderator.".
                END IF. 
                DO IF (setswv >0).  
                    print condnam /title = "Conditional Effect of Focal Predictor on Outcome at requested values of Moderator(s)" /format = A8. 
                    print XMgWvres((j*setswv-(setswv-1)):(j*setswv), :) /title = " " /cnames = XYgWlabs /format = !decimals.  
                END IF.    
              print M3df2 /title = "Degrees of freedom for all conditional effects:".
              
              DO IF (jn = 1). 
                    print /title = "--------------------------------------------------------------------------------------".
                    print /title = "                    JOHNSON-NEYMAN PROCEDURE: Mediator Model". 
                    DO IF (Mpairs = 1). 
                         print {"Mdiff = " , mnamemat(j,1), ' - ', mnamemat(j,2)} /title = " " /rlabels = "Outcome:" /format = A8.
                      ELSE. 
                         print {temp1(1,j) , mnamemat(j,1), ' - ', mnamemat(j,2)} /title = " " /rlabels = "Outcome:" /format = A8.
                      END IF. 
                    DO IF (aNumJN(j,1) <> 0). 
                        print /title = "Moderator value(s) defining Johnson-Neyman significance region(s) and percent of ".
                        print {t(aJNsoln(j,1:aNumJN(j,1))), t(apcntabv(j,1:aNumJN(j,1)))} /title "observed data above value:" /clabels = "Value", "% Abv" /format = !decimals /space = 0. 
                        print aJNres((1+22*(j-1)):(20+aNumJN(j,1)+22*(j-1)),:) /title = "Conditional Effect of 'X' on M at values of moderator" /cnames = XMgWlabs /format = !decimals.
                        print M3df2 /title = "Degrees of freedom for all conditional effects:".
                    ELSE IF (aNumJN(j,1) = 0). 
                        print /title = "There are no statistically significant transition points within the observed range of data.". 
                    END IF. 
                END IF. 
              
              
          END IF. 
      END IF. 
   END LOOP. 

END IF. 

print {"Ydiff =" , ynames(1,1), ' - ', ynames(1,2)} /title = "**************************************************************************************" /rlabels = "Outcome:" /format = A8.
print modsumr /title = "Model Summary" /clabels = "R", "R-sq", "MSE", "F" , "df1" , "df2", "p" /format = !decimals. 

DO IF ((model = 1) OR (model >= 4)). 
      DO IF (xmint = 1). 
      COMPUTE modres = {cpresmat;bresmat;dresmat}. 
      ELSE IF (xmint = 0). 
       COMPUTE modres = {cpresmat;bresmat}. 
      END IF. 
   COMPUTE bdlabs = {"M1diff", "M2diff", "M3diff", "M4diff", "M5diff", "M6diff", "M7diff", "M8diff", "M9diff", "M10diff"}.
   COMPUTE bslabs = {"M1avg", "M2avg", "M3avg", "M4avg", "M5avg", "M6avg", "M7avg", "M8avg", "M9avg", "M10avg"}.
   COMPUTE intlabs1 = {"Int1", "Int2", "Int3", "Int4", "Int5", "Int6", "Int7", "Int8", "Int9", "Int10"}. 
   COMPUTE intlabs2 = {"Int11", "Int12", "Int13", "Int14", "Int15", "Int16", "Int17", "Int18", "Int19", "Int20"}. 
   COMPUTE intlabs3 = {"Int21", "Int22", "Int23", "Int24", "Int25", "Int26", "Int27", "Int28", "Int29", "Int30"}. 
   COMPUTE intlabs = {intlabs1, intlabs2, intlabs3}. 
   DO IF (Mpairs = 1). 
      COMPUTE modlabs = {"constant"}.
      DO IF (cppthmd = 1). 
          COMPUTE modlabs = {modlabs, "W"}.
      END IF. 
      COMPUTE modlabs = {modlabs, "Mdiff"}. 
      COMPUTE intcount = 1. 
      DO IF (bpathmod = 1). 
          COMPUTE modlabs = {modlabs, intlabs(1, intcount)}.
          COMPUTE intcount = intcount + 1. 
      END IF. 
      COMPUTE modlabs = {modlabs, "Mavg", intlabs(1, intcount)}. 
      COMPUTE bdlabs = {"Mdiff"}.  
   ELSE. 
      COMPUTE modlabs = {"constant"}.
      DO IF (cppthmd = 1). 
          COMPUTE modlabs = {modlabs, "W"}. 
      END IF. 
      COMPUTE modlabs = {modlabs, bdlabs(1, 1:Mpairs)}.
      COMPUTE intcount = 1.
      DO IF (bpathmod = 1). 
          COMPUTE modlabs = {modlabs, intlabs(1, intcount:(intcount+Mpairs-1))}. 
          COMPUTE intcount = intcount + Mpairs. 
      END IF. 
      COMPUTE modlabs = {modlabs, bslabs(1, 1:Mpairs), intlabs(1, intcount:(intcount + Mpairs-1))}. 
   END IF. 
ELSE IF ((Model = 2) OR (model = 3)). 
   COMPUTE wlabs = wnames.
   COMPUTE intlabs1 = {"Int1", "Int2", "Int3", "Int4", "Int5", "Int6", "Int7", "Int8", "Int9", "Int10"}. 
   COMPUTE intlabs2 = {"Int11", "Int12", "Int13", "Int14", "Int15", "Int16", "Int17", "Int18", "Int19", "Int20"}. 
   COMPUTE intlabs3 = {"Int21", "Int22", "Int23", "Int24", "Int25", "Int26", "Int27", "Int28", "Int29", "Int30"}. 
   COMPUTE intlabs = {intlabs1, intlabs2, intlabs3}. 
   COMPUTE modlabs = {"constant"; t(wlabs)}.
   DO IF ((model = 3) AND (Wcount >1)). 
      COMPUTE modlabs = {modlabs; t(intlabs(1,1:nint))}. 
   END IF. 
END IF. 
print modres /title "Model" /rnames = modlabs /clabels = "coeff" , "SE", "t", "p", "LLCI", "ULCI" /format = !decimals.
print df2 /title = "Degrees of freedom for all regression coefficient estimates:".
DO IF (cppthmd = 1).
     print /title = "--------------------------------------------------------------------------------------". 
    COMPUTE XYgWcMlb = {wnames}. 
    DO IF (xmint = 1). 
        DO IF (Mpairs = 1). 
            COMPUTE bsmpairs = "Mavg". 
        ELSE. 
            COMPUTE bsmpairs = bslabs(1:Mpairs). 
        END IF. 
        COMPUTE XYgWcMlb = {XYgWcMlb, bsmpairs}. 
    END IF. 
    COMPUTE XYgWcMlb = {XYgWcMlb, coeflabs}. 
    COMPUTE condnam(:, 2) = {"'X'"; "Ydiff"; wnamemat(1,1)}.  
    print condnam /title = "Conditional Effect of Focal Predictor on Outcome at values of Moderator(s)" /format = A8. 
    print XYgWcMrs /title = " " /cnames = XYgWcMlb /format = !decimals.
    DO IF (csum(dich(:,1)) <> Wcount). 
        DO IF (quantile = 1). 
           print /title = "     Values for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentile.".
        ELSE IF (quantile = 0). 
           print /title = "     Values for quantitative moderators are the mean and plus/minus one SD from the mean.".
        END IF. 
    END IF. 
    DO IF (csum(dich(:,1)) > 0). 
        print /title = "     Values for dichotomous moderators are the two values of the moderator.".
    END IF. 
    DO IF (xmint = 1). 
          print /title = "     Values for mediator averages are the conditional values based on the values of the moderator.".
    END IF. 
    DO IF (setswv >0).  
         print condnam /title = "Conditional Effect of Focal Predictor on Outcome at requested values of Moderator(s)" /format = A8. 
         print XYgWcMv /title = " " /cnames = XYgWcMlb /format = !decimals.  
         DO IF (xmint = 1). 
             print /title = "     Values for mediator averages (Mavg) are the conditional values based on the values of the moderator.".
        END IF.
    END IF.
    
    DO IF (jn = 1). 
            print /title = "                    JOHNSON-NEYMAN PROCEDURE: Direct Effect". 
            DO IF (ncol(cpJNsoln) <> 0). 
            print /title = "Moderator value(s) defining Johnson-Neyman significance region(s) and percent of ".
            print {cpJNsoln, t(cppcntabv)} /title "observed data above value:" /clabels = "Value", "% Abv" /format = !decimals /space = 0. 
            print cpJNRes /title = "Conditional Direct Effect of 'X' on Y at values of moderator" /cnames = XYgWlabs /format = !decimals.
            ELSE IF (ncol(cJNsoln) = 0). 
            print /title = "There are no statistically significant transition points within the observed range of data.". 
            END IF. 
    END IF. 
    print /title = "--------------------------------------------------------------------------------------".
END IF. 

DO IF (bpathmod = 1). 
    LOOP j = 1 to Mpairs. 
         DO IF (mpairs = 1). 
               COMPUTE condnam(:, 2) = {"Mdiff"; "Ydiff"; wnamemat(1,1)}. 
         ELSE. 
               COMPUTE condnam(:, 2) = {temp6(1,j); "Ydiff"; wnamemat(1,1)}. 
         END IF. 
         COMPUTE condnam(:, 3) = {"(M)"; "(Y)"; "(W)"}. 
         print condnam /title = "Conditional Effect of Focal Predictor on Outcome at values of Moderator(s)" /format = A8. 
          print MYgWres((j*dimmc-(dimmc-1)):(j*dimmc), :) /title = " " /cnames = XYgWlabs /format = !decimals. 
          DO IF (csum(dich(:,1)) <> Wcount). 
              DO IF (quantile = 1). 
                 print /title = "     Values for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentile.".
              ELSE IF (quantile = 0). 
                 print /title = "     Values for quantitative moderators are the mean and plus/minus one SD from the mean.".
              END IF. 
          END IF. 
          DO IF (csum(dich(:,1)) > 0). 
              print /title = "     Values for dichotomous moderators are the two values of the moderator.".
          END IF. 
          DO IF (setswv >0).  
               print condnam /title = "Conditional Effect of Focal Predictor on Outcome at requested values of Moderator(s)" /format = A8. 
               print MYgWvres((j*setswv-(setswv-1)):(j*setswv), :) /title = " " /cnames = XYgWlabs /format = !decimals.  
          END IF.       
          
         DO IF (jn = 1). 
                    print condnam /title = "JOHNSON-NEYMAN PROCEDURE: b-paths (Mdiff --> Y)" /format = A8. 
                    DO IF (bNumJN(j,1) <> 0). 
                    print /title = "Moderator value(s) defining Johnson-Neyman significance region(s) and percent of ".
                    print {t(bJNsoln(j,1:bNumJN(j,1))), t(bpcntabv(j,1:bNumJN(j,1)))} /title "observed data above value:" /clabels = "Value", "% Abv" /format = !decimals /space = 0. 
                    print bJNres((1+22*(j-1)):(20+bNumJN(j,1)+22*(j-1)),:) /title = "Conditional Effect of M on Y at values of moderator" /cnames = XYgWlabs /format = !decimals.
                    ELSE IF (bNumJN(j,1) = 0). 
                    print /title = "There are no statistically significant transition points within the observed range of data.". 
                    END IF. 
        END IF.     
        print /title = "--------------------------------------------------------------------------------------". 
    END LOOP. 
END IF. 

DO IF (dpathmod = 1). 
     print /title = "--------------------------------------------------------------------------------------".      
    LOOP j = 1 to Mpairs. 
         DO IF (mpairs = 1). 
               COMPUTE condnam(:, 2) = {"Mavg"; "Ydiff"; wnamemat(1,1)}. 
         ELSE. 
               COMPUTE condnam(:, 2) = {temp7(1,j); "Ydiff"; wnamemat(1,1)}. 
         END IF. 
         COMPUTE condnam(:, 3) = {"(XM)"; "(Y)"; "(W)"}. 
         print condnam /title = "Conditional Effect of Focal Predictor on Outcome at values of Moderator(s)" /format = A8. 
          print MAYgWres((j*dimmc-(dimmc-1)):(j*dimmc), :) /title = " " /cnames = XYgWlabs /format = !decimals. 
          DO IF (csum(dich(:,1)) <> Wcount). 
              DO IF (quantile = 1). 
                 print /title = "     Values for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentile.".
              ELSE IF (quantile = 0). 
                 print /title = "     Values for quantitative moderators are the mean and plus/minus one SD from the mean.".
              END IF. 
          END IF. 
          DO IF (csum(dich(:,1)) > 0). 
              print /title = "     Values for dichotomous moderators are the two values of the moderator.".
          END IF. 
          DO IF (setswv >0).  
               print condnam /title = "Conditional Effect of Focal Predictor on Outcome at requested values of Moderator(s)" /format = A8. 
               print MAYgWvrs((j*setswv-(setswv-1)):(j*setswv), :) /title = " " /cnames = XYgWlabs /format = !decimals.  
          END IF.    
         
        DO IF (jn = 1). 
                    print condnam /title = "JOHNSON-NEYMAN PROCEDURE: d-paths (Mavg --> Y)" /format = A8. 
                    DO IF (dNumJN(j,1) <> 0). 
                    print /title = "Moderator value(s) defining Johnson-Neyman significance region(s) and percent of ".
                    print {t(dJNsoln(j,1:dNumJN(j,1))), t(dpcntabv(j,1:dNumJN(j,1)))} /title "observed data above value:" /clabels = "Value", "% Abv" /format = !decimals /space = 0. 
                    print dJNres((1+22*(j-1)):(20+dNumJN(j,1)+22*(j-1)),:) /title = "Conditional Effect of Mavg on Y at values of moderator" /cnames = XYgWlabs /format = !decimals.
                    ELSE IF (dNumJN(j,1) = 0). 
                    print /title = "There are no statistically significant transition points within the observed range of data.". 
                    END IF. 
        END IF.      
                
    print /title = "--------------------------------------------------------------------------------------".                        
    END LOOP. 
END IF. 

DO IF ((cppthmd =1) OR (bpathmod = 1) OR (dpathmod = 1)).
    print df2 /title = "Degrees of freedom for all conditional effects:".
END IF. 

DO IF ((model = 1) OR (model > 3)). 

   COMPUTE collab = {"Effect", "SE", "t", "df", "p", "LLCI", "ULCI"}. 
   DO IF (anymod = 0). 
       print /title = "************************* TOTAL, DIRECT, AND INDIRECT EFFECTS *************************" .
       COMPUTE cresmat = {cresmat(:, 1:3), MAKE(nrow(cresmat), 1, M4df2), cresmat(:, 4:6)}. 
       print cresmat /title = "Total effect of X on Y" /cnames = collab /format = !decimals. 
   ELSE IF (anymod > 0). 
          print /title = "******************* CONDITIONAL TOTAL, DIRECT, AND INDIRECT EFFECTS *******************" .
       DO IF (setswv >0). 
           COMPUTE XYgWres = {XygWres; XYgWvres}. 
       END IF. 
       COMPUTE XYgWlabs = {XYgWlabs(1:(ncol(XYgWlabs)-3)), "df", XYgWlabs((ncol(XYgWlabs)-2):ncol(XYgWlabs))}. 
       COMPUTE XYgWres = {XYgWres(:,1:(ncol(XYgWres) - 3)), MAKE(nrow(XYgWres), 1, M4df2), XYgWres(:, (ncol(XYgWres)-2):ncol(XYgWres))}. 
       print XYgWres /cnames = XYgWlabs /format = !decimals /title = "Conditional Total Effect of X on Y at values of the Moderator(s)".    
       DO IF (csum(dich(:,1)) <> Wcount). 
           DO IF (quantile = 1). 
                 print /title = "     Values for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentile.".
           ELSE IF (quantile = 0). 
                 print /title = "     Values for quantitative moderators are the mean and plus/minus one SD from the mean.".
           END IF. 
       END IF. 
       DO IF (csum(dich(:,1)) > 0). 
              print /title = "     Values for dichotomous moderators are the two values of the moderator.".
       END IF. 
       DO IF (setswv >0). 
           print /title = "     Requested values for moderators included in table above.".
       END IF. 
   END IF.    
   
   DO IF (cppthmd = 0).      
       COMPUTE cpresmat = {cpresmat(1:3), df2, cpresmat(4:6)}.
       print cpresmat /title = "Direct effect of X on Y" /cnames = collab /format = !decimals. 
   ELSE IF (cppthmd = 1). 
      DO IF (setswv >0). 
           COMPUTE XYgWcMrs = {XYgWcMrs; XYgWcMv}. 
       END IF. 
       COMPUTE XYgWcMrs = {XYgWcMrs(:,1:(ncol(XYgWcMrs) - 3)), MAKE(nrow(XYgWcMrs), 1, df2), XYgWcMrs(:, (ncol(XYgWcMrs)-2):ncol(XYgWcMrs))}. 
       COMPUTE XYgWcMlb = {XYgWcMlb(1:(ncol(XYgWcMlb) - 3)), "df", XYgWcMlb((ncol(XYgWcMlb) - 2):ncol(XYgWcMlb))}. 
       print XYgWcMrs /cnames = XYgWcMlb /format = !decimals /title = "Conditional Direct Effect of X on Y at values of the Moderator(s)".    
       DO IF (csum(dich(:,1)) <> Wcount). 
           DO IF (quantile = 1). 
                 print /title = "     Values for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentile.".
           ELSE IF (quantile = 0). 
                 print /title = "     Values for quantitative moderators are the mean and plus/minus one SD from the mean.".
           END IF. 
       END IF. 
       DO IF (csum(dich(:,1)) > 0). 
              print /title = "     Values for dichotomous moderators are the two values of the moderator.".
       END IF. 
       DO IF (setswv >0). 
           print /title = "     Requested values for moderators included in table above.".
       END IF. 
       DO IF (xmint = 1). 
             print /title = "     Values for mediator averages (Mavg) are the conditional values based on the values of the moderator.".
        END IF.
   END IF. 
   
   DO IF (mc = 1). 
      COMPUTE indlabs = {"Effect", "MCSE", "MCLLCI", "MCULCI"}. 
      COMPUTE indres = MCres. 
   ELSE. 
      COMPUTE indlabs = {"Effect", "BootSE", "BootLLCI", "BootULCI"}. 
      COMPUTE indres = bootres. 
   END IF.  

   COMPUTE mlab = {"Ind1", "Ind2", "Ind3", "Ind4", "Ind5", "Ind6", "Ind7", "Ind8", "Ind9", "Ind10"}. 
   COMPUTE mlab2 = {"Ind11", "Ind12", "Ind13", "Ind14", "Ind15", "Ind16", "Ind17", "Ind18", "Ind19", "Ind20"}. 
   COMPUTE mlab3 = {"Ind21", "Ind22", "Ind23", "Ind24", "Ind25", "Ind26", "Ind27", "Ind28", "Ind29", "Ind30"}. 
   COMPUTE mlab4 = {"Ind31", "Ind32", "Ind33", "Ind34", "Ind35", "Ind36", "Ind37", "Ind38", "Ind39", "Ind40"}. 
   COMPUTE mlab = {mlab, mlab2, mlab3, mlab4}. 


   DO IF ((apathmod = 0) AND (bpathmod = 0)). 
       COMPUTE m2lab = {mlab(1,1:(nrow(indres)-1)), 'Total'}. 
       DO IF (Mpairs = 1). 
           print indres(1,:) /title = "Indirect Effect of X on Y through M" /rnames = m2lab / cnames = indlabs /format = !decimals. 
       ELSE. 
       print indres /title = "Indirect Effect of X on Y through M" /rnames = m2lab / cnames = indlabs /format = !decimals. 
           END IF. 
       DO IF (normal = 1). 
           print normres /title = "Normal Theory Tests for Indirect Effect" /rnames = mlab /clabels = "Effect", "SE", "Z", "p" /format = !decimals. 
       END IF.
   ELSE IF ((apathmod =1) OR (bpathmod = 1)). 
       COMPUTE m2lab = {mlab(1,1:mpairs), 'Total'}. 
       COMPUTE indlabsW = {wnamemat(1,1), indlabs}. 
       COMPUTE condnam1 = {"Ind:", " ", " "}. 
       COMPUTE condnam2 = {"Med:", "  ", "  "}.
       COMPUTE condnam = {condnam1; condnam2}. 
       LOOP i = 1 to Mpairs. 
              DO IF (Mpairs = 1). 
                  COMPUTE condnam(2,2) = {"Mdiff"}. 
                  COMPUTE condnam(2,3) = {"(M)"}. 
                  COMPUTE condnam(1,2) = {"Ind1"}. 
              ELSE.
                  COMPUTE condnam(2,2) = temp6(1, i). 
                  COMPUTE condnam(2,3) = temp9(1, i). 
                  COMPUTE condnam(1,2) = mlab(1,i). 
              END IF. 
              print condnam  /format = A9 /title = "Conditional Indirect Effect of X on Y through Mediator at values of the Moderator".    
              print cindres((1+(i-1)*(dimmc+setswv)):(i*(dimmc+setswv)),:) /cnames = indlabsW /format = !decimals /title = " ". 
              DO IF (csum(dich(:,1)) <> Wcount). 
                  DO IF (quantile = 1). 
                      print /title = "     Values for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentile.".
                  ELSE IF (quantile = 0). 
                      print /title = "     Values for quantitative moderators are the mean and plus/minus one SD from the mean.".
                  END IF. 
              END IF. 
              DO IF (csum(dich(:,1)) > 0). 
                  print /title = "     Values for dichotomous moderators are the two values of the moderator.".
              END IF. 
              DO IF (setswv >0). 
                  print /title = "     Requested values for moderators included in table above.".
              END IF. 
       END LOOP. 
  END IF. 
     COMPUTE indkey = make((ncol(m2lab)-1), (3+2*Mpairs), " "). .
   LOOP i = 1 to mpairs.
      COMPUTE indkey(i,1:5) = {"'X'" , "->", bdlabs(1,i), "->", "Ydiff"}. 
   END LOOP.  
   DO IF (serial = 1). 
       COMPUTE counter = 1. 
           LOOP j = 1 TO Mpairs-1. 
               LOOP m = 1 TO Mpairs-j.
                    COMPUTE step1 = {"'X'", "->", bdlabs(1,m)}. 
                              LOOP l = m+1 to (Mpairs). 
                                 COMPUTE step2 = {step1, "->", bdlabs(1,l)}.
                                 DO IF (j > 1).                  
                                    LOOP h = (l+1) to Mpairs. 
                                       COMPUTE step3 = {step2, "->", bdlabs(1,h)}.
                                       DO IF (j > 2). 
                                          LOOP o = h+1 to Mpairs.  
                                             COMPUTE step4 = {step3, "->", bdlabs(1,o)}. 
                                             DO IF (j > 3). 
                                             COMPUTE step5 = {step4, "->", bdlabs(1,5), "->", "YDiff" }.
                                             COMPUTE indkey(mpairs+counter,1:13) = step5. 
                                             COMPUTE counter = counter+1.
                                            ELSE. 
                                                COMPUTE indkey(mpairs+counter,1:11) = {step4, "->", "YDiff" }. 
                                                 COMPUTE counter = counter+1.
                                              END IF. 
                                          END LOOP IF (j > 3). 
                                          *end oloop. 
                                      ELSE. 
                                             COMPUTE indkey(mpairs+counter,1:9) = {step3, "->", "YDiff" }. 
                                         COMPUTE counter = counter+1.
                                       END IF. 
                                    END LOOP IF (j > 3). 
                                    /*end hloop*/. 
                                 ELSE. 
                                    COMPUTE indkey(mpairs+counter,1:7) = {step2, "->", "YDiff" }. 
                                    COMPUTE counter = counter+1.
                                 END IF. 
                              END LOOP IF (j > 3). 
                              /*end lloop*/. 
                           END LOOP IF (j > 3). 
                           /*end mloop*/. 
                        END LOOP. 
                        /*end jloop*/. 
   ELSE. 
       COMPUTE indkey = indkey(:, 1:5). 
   END IF.  
   print indkey /title = "Indirect Key" /rnames =m2lab /format = A8. 
   
   DO IF  (contrast = 1) AND (Mpairs >1).
   COMPUTE contlab1 = {'(C1)', '(C2)', '(C3)', '(C4)', '(C5)', '(C6)', '(C7)', '(C8)', '(C9)', '(C10)'}.
   COMPUTE contlab2 = {'(C11)', '(C12)', '(C13)', '(C14)', '(C15)', '(C16)', '(C17)', '(C18)', '(C19)', '(C20)'}.
   COMPUTE contlab3 =  {'(C21)', '(C22)', '(C23)', '(C24)', '(C25)', '(C26)', '(C27)', '(C28)', '(C29)', '(C30)'}.
   COMPUTE contlab4 = {'(C31)', '(C32)', '(C33)', '(C34)', '(C35)', '(C36)', '(C37)', '(C38)', '(C39)', '(C40)'}.
   COMPUTE contlab5 =  {'(C41)', '(C42)', '(C43)', '(C44)', '(C45)'}.
   COMPUTE contlab = {contlab1, contlab2, contlab3, contlab4}. 
   print contres /title = "Pairwise Contrasts Between Specific Indirect Effects" /rnames = contlab /cnames = indlabs /format = !decimals. 
   COMPUTE contkey = MAKE(npairs, 3, 0). 
   COMPUTE counter = 1. 
   LOOP i = 1 to nrow(indres)-2. 
   LOOP j = i+1 to nrow(indres)-1. 
   COMPUTE contkey(counter,:) = {mlab(1,i), ' - ', mlab(1,j)}. 
   COMPUTE counter = counter+1. 
   END LOOP. 
   END LOOP. 
   print contkey /title = "Contrast Key:" /rnames = contlab /format = A8. 
   END IF. 
   
    DO IF (anymod > 0).  
          print /title = "******************************** INDICES OF MODERATION ********************************" .
          
          COMPUTE cresmat= {cresmat(2:(1+Wcount), 1:3), M4df2, cresmat(2:(1+Wcount), 4:6)}.
          print cresmat /title = "Test of Moderation of the Total Effect" /rnames = {"W"} /cnames = collab /format = !decimals. 

          DO IF (cppthmd = 1). 
              COMPUTE cpresmat = {cpresmat(2:(1+Wcount),1:3), df2, cpresmat(2:(1+Wcount),4:6)}.
              print cpresmat /title = "Test of Moderation of the Direct Effect" /rnames = {"W"} /cnames = collab /format = !decimals. 
          END IF. 

          DO IF ((apathmod = 1) OR (bpathmod = 1)).
              DO IF (csum(dich(:,1)) > 0). 
                   print immres /rnames = mlab/cnames = indlabs /format = !decimals /title = "Index of Moderated Mediation for each Indirect Effect.". 
              END IF. 
              DO IF (csum(dich(:,1)) = 0). 
                  DO IF ((apathmod = 1) AND (bpathmod = 1)). 
                       print /title = "The INDEX OF MODERATED MEDIATION is not generated for this model because".
                       print /title = "the indirect effect is a non-linear function of the moderator.".
                  ELSE. 
                       COMPUTE indices = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20}. 
                       DO IF (mc=0). 
                           COMPUTE immres = bootres(indices(1:Mpairs),:).
                       ELSE IF (mc=1). 
                           COMPUTE immres = MCres(indices(1:Mpairs),:).
                       END IF. 
                           print immres /rnames = mlab/cnames = indlabs /format = !decimals /title = "Index of Moderated Mediation for each Indirect Effect.". 
                  END IF. 
              END IF. 
         END IF. 
        
        DO IF (plot = 1). 
            print /title = "************************************ PLOTS *******************************************" .
            
            print /title = "Data for visualizing conditional effect of X on Y at values of W" /space = 0.
            DO IF (apathmod = 1). 
                print /title = "Data for visualizing conditional effect of X on M at values of W".
            END IF. 
            print/title = "Paste text below into a SPSS syntax window and execute to produce plot."/space=0.
            DO IF (center = 1). 
            print /title = "Note: All moderator values have been centered." /space = 0. 
            ELSE IF (center = 2). 
            print /title = "Note: All continuous moderator values have been centered." /space = 0. 
            END IF. 
            !let !line0 = !concat("DATA LIST FREE/", !w, " ", "YdiffHAT").
            DO IF (apathmod = 1). 
                DO IF (Mpairs = 1).  
                    !let !line0 = !concat(!line0, " " ,"MdiffHAT").
                    print /title = !quote(!concat(!line0, "."))/space=1.
                END IF.
                DO IF (Mpairs > 1). 
                    !let !line0 = !concat("DATA LIST FREE/", !w, " ", "YdiffHAT", " ").
                    !let !line0 = !concat(!line0, "M1diffHAT").
                END IF. 
                DO IF (Mpairs >= 2). 
                    !let !line0 = !concat(!line0, " ","M2diffHAT").
                    DO IF (Mpairs = 2). 
                        print /title = !quote(!concat(!line0, "."))/space=1.
                    END IF. 
                END IF. 
                DO IF (Mpairs >= 3). 
                    !let !line0 = !concat(!line0, " ","M3diffHAT").
                    DO IF (Mpairs = 3). 
                        print /title = !quote(!concat(!line0, "."))/space=1.
                    END IF. 
                END IF. 
                DO IF (Mpairs >= 4). 
                    !let !line0 = !concat(!line0, " ","M4diffHAT").
                    DO IF (Mpairs = 4). 
                        print /title = !quote(!concat(!line0, "."))/space=1.
                    END IF. 
                END IF. 
                DO IF (Mpairs >= 5). 
                    !let !line0 = !concat(!line0, " ","M5diffHAT").
                    DO IF (Mpairs = 5). 
                        print /title = !quote(!concat(!line0, "."))/space=1.
                    END IF. 
                END IF.             
                DO IF (Mpairs >= 6). 
                    !let !line0 = !concat(!line0, " ", "M6diffHAT").
                    DO IF (Mpairs = 6). 
                        print /title = !quote(!concat(!line0, "."))/space=1.
                    END IF. 
                END IF.             
                DO IF (Mpairs >= 7). 
                    !let !line0 = !concat(!line0," ", "M7diffHAT").
                    DO IF (Mpairs = 7). 
                        print /title = !quote(!concat(!line0, "."))/space=1.
                    END IF. 
                END IF. 
                DO IF (Mpairs >= 8). 
                    !let !line0 = !concat(!line0, " ","M8diffHAT").
                    DO IF (Mpairs = 8). 
                        print /title = !quote(!concat(!line0, "."))/space=1.
                    END IF. 
                END IF.             
                DO IF (Mpairs >= 9). 
                    !let !line0 = !concat(!line0, " ","M9diffHAT").
                    DO IF (Mpairs = 9). 
                        print /title = !quote(!concat(!line0, "."))/space=1.
                    END IF. 
                END IF.             
                DO IF (Mpairs >= 10). 
                    !let !line0 = !concat(!line0, " ","M10diffHAT").
                    DO IF (Mpairs = 10). 
                        print /title = !quote(!concat(!line0, "."))/space=1.
                    END IF. 
                END IF. 
            END IF. 
            !let !line1 = "BEGIN DATA.".
            !let !line6 = "END DATA.".
            
            print /title = !quote(!line1)/space=1.
            print plotdat /title = " "/format = !decimals/space=0.
            print /title = !quote(!line6)/space=1.
            
            !let !line7 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "YdiffHAT", ".").
            !let !line8 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "MdiffHAT", ".").
            !let !line9 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "M1diffHAT", ".").
            !let !line10 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "M2diffHAT", ".").
            !let !line11 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "M3diffHAT", ".").
            !let !line12 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "M4diffHAT", ".").
            !let !line13 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "M5diffHAT", ".").
            !let !line14 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "M6diffHAT", ".").
            !let !line15 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "M7diffHAT", ".").
            !let !line16 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "M8diffHAT", ".").
            !let !line17 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "M9diffHAT", ".").
            !let !line18 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "M10diffHAT", ".").

            print /title = !quote(!line7)/space=1.
            DO IF (apathmod = 1). 
                DO IF (Mpairs=1). 
                    print /title = !quote(!line8)/space=1.
                END IF. 
                DO IF (Mpairs > 1). 
                    print /title = !quote(!line9)/space=1.
                    DO IF (Mpairs >= 2). 
                        print /title = !quote(!line10)/space=1.
                        DO IF (Mpairs >= 3). 
                            print /title = !quote(!line11)/space=1. 
                            DO IF (Mpairs >= 4). 
                                print /title = !quote(!line12)/space=1.
                                DO IF (Mpairs >= 5). 
                                    print /title = !quote(!line13)/space=1.
                                    DO IF (Mpairs >= 6). 
                                        print /title = !quote(!line14)/space=1.
                                        DO IF (Mpairs >= 7). 
                                            print /title = !quote(!line15)/space=1.
                                            DO IF (Mpairs >= 8). 
                                                print /title = !quote(!line16)/space=1.
                                                DO IF (Mpairs >= 9). 
                                                    print /title = !quote(!line17)/space=1.
                                                    DO IF (Mpairs = 10). 
                                                        print /title = !quote(!line18)/space=1.
                                                    END IF. 
                                                END IF. 
                                            END IF. 
                                        END IF. 
                                    END IF. 
                                END IF. 
                            END IF. 
                        END IF.
                    END IF. 
                END IF. 
            END IF. 
            
            DO IF (bpathmod = 1). 
                print /title = "--------------------------------------------------------------------------------------".
                print /title = "Data for visualizing conditional effect of M (Mdiff) on Y at values of W.".
                print/title = "Paste text below into a SPSS syntax window and execute to produce plot."/space=0.
                DO IF (center = 1). 
                print /title = "Note: All moderator values have been centered." /space = 0. 
                ELSE IF (center = 2). 
                print /title = "Note: All continuous moderator values have been centered." /space = 0. 
                END IF. 
                LOOP i = 1 to Mpairs. 
                    print mnames(1, (2*i-1):(2*i)) /title = "The following section applies to the mediators defined by:" /format = A8. 
                    !let !line0 = !concat("DATA LIST FREE/", !w, " ", "Mdiff", " ", "YdiffHAT", ".").                   
                    print /title = !quote(!line0)/space=1. 
                    print /title = !quote(!line1)/space=1.
                    print bplotdat((1+(i-1)*(3*(3-csum(dich(:,1))))):(i*(3*(3-csum(dich(:,1))))),:) /title = " "/format = !decimals/space=0.
                    print /title = !quote(!line6)/space=1.       
                    !let !line7 = !concat("GRAPH/SCATTERPLOT = ", "Mdiff", " WITH ", "YdiffHAT ", "BY ", !HEAD(!w),".").
                    print /title = !quote(!line7)/space=1. 
                END LOOP. 
            END IF.     
            
            DO IF (cppthmd = 1). 
                print /title = "--------------------------------------------------------------------------------------".
                print /title = "Data for visualizing conditional direct effect of X on Y at values of W.".
                print/title = "Paste text below into a SPSS syntax window and execute to produce plot."/space=0.
                DO IF (center = 1). 
                print /title = "Note: All moderator values have been centered." /space = 0. 
                ELSE IF (center = 2). 
                print /title = "Note: All continuous moderator values have been centered." /space = 0. 
                END IF. 
                DO IF (xmint = 1). 
                    print /title = "Note: All mediator averages are conditioned on their predicted value based on W." /space = 0. 
                END IF. 
                !let !line0 = !concat("DATA LIST FREE/", !w, " ", "YdiffHAT", ".").
                print /title = !quote(!line0)/space=1. 
                print /title = !quote(!line1)/space=1.
                print cppltdat /title = " "/format = !decimals/space=0.
                print /title = !quote(!line6)/space=1.
                !let !line7 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "YdiffHAT.").
                print /title = !quote(!line7)/space=1. 
            END IF. 
            
            
            DO IF (dpathmod = 1). 
                print /title = "--------------------------------------------------------------------------------------".
                print /title = "Data for visualizing differential conditional effect of M (Mavg) on Y at values of W.".
                print/title = "Paste text below into a SPSS syntax window and execute to produce plot."/space=0.
                DO IF (center = 1). 
                print /title = "Note: All moderator values have been centered." /space = 0. 
                ELSE IF (center = 2). 
                print /title = "Note: All continuous moderator values have been centered." /space = 0. 
                END IF. 
                LOOP i = 1 to Mpairs. 
                    print mnames(1, (2*i-1):(2*i)) /title = "The following section applies to the mediators defined by:" /format = A8. 
                    !let !line0 = !concat("DATA LIST FREE/", !w, " ", "Mavg", " ", "YdiffHAT", ".").
                    
                    !let !line1 = "BEGIN DATA.".
                    !let !line6 = "END DATA.".
                    
                    print /title = !quote(!line0)/space=1. 
                    print /title = !quote(!line1)/space=1.
                    print dplotdat((1+(i-1)*9):(i*9),:) /title = " "/format = !decimals/space=0.
                    print /title = !quote(!line6)/space=1.
                    !let !line7 = !concat("GRAPH/SCATTERPLOT = ", "Mavg", " WITH ", "YdiffHAT ", "BY ", !HEAD(!w),".").
                    print /title = !quote(!line7)/space=1. 
                END LOOP. 
            END IF.     
            
            
            /*GRAPH/SCATTERPLOT=Y2 WITH Y1 BY J1 BY K2 (NAME) /PANEL ROWVAR=J2 COLVAR = K1*/.
        END IF. 
         
         
         
   END IF. 
   


   
ELSE IF ((Model = 2) OR (Model = 3)). 
       print /title = "**************************************************************************************" .
       COMPUTE coeflabs = {"Effect" , "SE", "t", "p", "LLCI", "ULCI"}.
       COMPUTE XYgWlabs = {wnames, coeflabs}.  
       print XYgWres /title = "Conditional Effect of 'X' on Y at values of moderator(s)" /cnames = XYgWlabs /format = !decimals. 
       print df2 /title = "Degrees of freedom for all conditional effects:".
    DO IF (csum(dich(:,1)) <> Wcount). 
        DO IF (quantile = 1). 
           print /title = "Values for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentile.".
        ELSE IF (quantile = 0). 
           print /title = "Values for quantitative moderators are the mean and plus/minus one SD from the mean.".
        END IF. 
    END IF. 
    DO IF (csum(dich(:,1)) > 0). 
       print /title = "Values for dichotomous moderators are the two values of the moderator.".
    END IF. 
    DO IF (setswv >0). 
    print /title = "--------------------------------------------------------------------------------------".
    print wvres /title = "Conditional Effect of 'X' on Y at requested values of modederator(s)" /cnames = XYgWlabs /format = !decimals. 
    print df2 /title = "Degrees of freedom for all conditional effects:".
    
    END IF. 
    
    DO IF (jn = 1). 
        print /title = "****************************** JOHNSON-NEYMAN PROCEDURE *******************************". 
        DO IF (numJN <> 0). 
        print /title = "Moderator value(s) defining Johnson-Neyman significance region(s) and percent of ".
        print {JNsoln, t(pcntabv)} /title "observed data above value:" /clabels = "Value", "% Abv" /format = !decimals /space = 0. 
        print JNRes /title = "Conditional Effect of 'X' on Y at values of moderator" /cnames = XYgWlabs /format = !decimals.
        print df2 /title = "Degrees of freedom for all conditional effects:".
        ELSE IF (numJN = 0). 
        print /title = "There are no statistically significant transition points within the observed range of data.". 
        END IF. 
    
    END IF. 
    
    print /title = "**************************************************************************************" .
       print /title = "Conditional Effect of Moderator(s) on Y in each Condition". 
       print ynames(1,1) /title = "Condition 1 Outcome:" /format = A8.
       print prbmsum(1,:) /title = "Model Summary" /clabels = "R", "R-sq", "MSE", "F" , "df1" , "df2", "p" /format = !decimals. 
       print prbmdres(1:(nrow(prbmdres)/2), :) /title "Model" /rnames = modlabs /clabels = "coeff" , "SE", "t", "p", "LLCI", "ULCI" /format = !decimals.
       print probedf2 /title = "Degrees of freedom for all conditional effects:".
    print /title = "--------------------------------------------------------------------------------------".
       print ynames(1,2) /title = "Condition 2 Outcome:" /format = A8.
       print prbmsum(2,:) /title = "Model Summary" /clabels = "R", "R-sq", "MSE", "F" , "df1" , "df2", "p" /format = !decimals. 
       print prbmdres((nrow(prbmdres)/2+1):nrow(prbmdres), :) /title "Model" /rnames = modlabs /clabels = "coeff" , "SE", "t", "p", "LLCI", "ULCI" /format = !decimals.
       print probedf2 /title = "Degrees of freedom for all conditional effects:".

DO IF (plot = 1).         
        print /title = "**************************************************************************************" .
        
        print /title = "Data for visualizing conditional effect of X on Y.".
        print/title = "Paste text below into a SPSS syntax window and execute to produce plot."/space=0.
        DO IF (center = 1). 
        print /title = "Note: All moderator values have been centered." /space = 0. 
        ELSE IF (center = 2). 
        print /title = "Note: All continuous moderator values have been centered." /space = 0. 
        END IF. 
        !let !line0 = !concat("DATA LIST FREE/", !w, " ", "YdiffHAT", " ", !HEAD(!Y), "HAT", !TAIL(!Y), "HAT", ".").
        !let !line1 = "BEGIN DATA.".
        !let !line6 = "END DATA.".
        
        print /title = !quote(!line0)/space=1.
        print /title = !quote(!line1)/space=1.
        print plotdat /title = " "/format = !decimals/space=0.
        print /title = !quote(!line6)/space=1.
        
        !let !line7 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "YdiffHAT").
        !let !line8 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", !HEAD(!Y), "HAT").
        !let !line9 = !concat("GRAPH/SCATTERPLOT = ",  !HEAD(!w), " WITH", !TAIL(!Y), "HAT").
        !let !line7 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", "YdiffHAT").
        !let !line8 = !concat("GRAPH/SCATTERPLOT = ", !HEAD(!w), " WITH ", !HEAD(!Y), "HAT").
        !let !line9 = !concat("GRAPH/SCATTERPLOT = ",  !HEAD(!w), " WITH", !TAIL(!Y), "HAT").
        
        DO IF (Wcount = 1). 
        print /title = !quote(!concat(!line7, "."))/space=1.
        print /title = !quote(!concat(!line8, "."))/space=0.
        print /title = !quote(!concat(!line9, "."))/space=0.
        END IF. 
        
        
        !let !line7 = !concat(!line7, " BY ", !HEAD(!TAIL(!w))).
        !let !line8 = !concat(!line8, " BY ", !HEAD(!TAIL(!w))).
        !let !line9 = !concat(!line9, " BY ", !HEAD(!TAIL(!w))).
        
        DO IF (Wcount = 2). 
        print /title = !quote(!concat(!line7, "."))/space=1.
        print /title = !quote(!concat(!line8, "."))/space=0.
        print /title = !quote(!concat(!line9, "."))/space=0.
        END IF. 
        
        !let !line10 = !concat(!line7, " /PANEL ROWVAR = ", !HEAD(!TAIL(!TAIL(!w)))).
        !let !line11 = !concat(!line8, " /PANEL ROWVAR = ", !HEAD(!TAIL(!TAIL(!w)))).
        !let !line12 = !concat(!line9, " /PANEL ROWVAR = ", !HEAD(!TAIL(!TAIL(!w)))).
        
        DO IF (Wcount = 3). 
        print /title = !quote(!concat(!line10, "."))/space=1.
        print /title = !quote(!concat(!line11, "."))/space=0.
        print /title = !quote(!concat(!line12, "."))/space=0.
        END IF.
        
        !let !line10 = !concat(!line10, " COLVAR = ", !HEAD(!TAIL(!TAIL(!TAIL(!w))))).
        !let !line11 = !concat(!line11, " COLVAR = ", !HEAD(!TAIL(!TAIL(!TAIL(!w))))).
        !let !line12 = !concat(!line12, " COLVAR = ", !HEAD(!TAIL(!TAIL(!TAIL(!w))))).
        
        DO IF (Wcount = 4). 
        print /title = !quote(!concat(!line10, "."))/space=1.
        print /title = !quote(!concat(!line11, "."))/space=0.
        print /title = !quote(!concat(!line12, "."))/space=0.
        END IF.
        
        !let !line7 = !concat(!line7, " ", "BY", !TAIL(!TAIL(!TAIL(!TAIL(!w)))), " ", "(NAME)", " /PANEL ROWVAR = ", !HEAD(!TAIL(!TAIL(!w))), " COLVAR = ", !HEAD(!TAIL(!TAIL(!TAIL(!w))))).
        !let !line8 = !concat(!line8, " ", "BY", !TAIL(!TAIL(!TAIL(!TAIL(!w)))), " ", "(NAME)", " /PANEL ROWVAR = ", !HEAD(!TAIL(!TAIL(!w))), " COLVAR = ", !HEAD(!TAIL(!TAIL(!TAIL(!w))))).
        !let !line9 = !concat(!line9, " ", "BY", !TAIL(!TAIL(!TAIL(!TAIL(!w)))), " ", "(NAME)", " /PANEL ROWVAR = ", !HEAD(!TAIL(!TAIL(!w))), " COLVAR = ", !HEAD(!TAIL(!TAIL(!TAIL(!w))))).
        
        DO IF (Mcount = 5). 
        print /title = !quote(!concat(!line7, "."))/space=1.
        print /title = !quote(!concat(!line8, "."))/space=0.
        print /title = !quote(!concat(!line9, "."))/space=0.
        END IF. 
        
        /*GRAPH/SCATTERPLOT=Y2 WITH Y1 BY J1 BY K2 (NAME) /PANEL ROWVAR=J2 COLVAR = K1*/.
    
END IF. 



END IF. 
    

    



END IF. 
print /title = "**************************** ANALYSIS NOTES AND WARNINGS *****************************". 
LOOP i = 1 to nrow(runnotes).
   DO IF (runnotes(i,1) = 1). 
      PRINT missing /title = "NOTE: Some cases were deleted due to missing data. The number of cases was:". 
  ELSE IF (runnotes(i,1) = 2). 
      PRINT /title = "ERROR: Two Y variables are needed.".
    ELSE IF (runnotes(i,1) = 3). 
        PRINT /title = "NOTE: An invalid number of samples was provided.". 
   ELSE IF (runnotes(i,1) = 4). 
   PRINT /title = "ERROR: The number of samples specified is insufficient for desired confidence.".
       PRINT /title = "       Please increase number of samples or decrease confidence." /space = 0.
      PRINT {Conf, Samples} /title = "Current Confidence & Samples:" /space = 0.
   ELSE IF (runnotes(i,1) = 5).
      PRINT /title = "NOTE: The confidence specified was not between 50 and 99.99. Level of confidence". 
      PRINT {95} /title =  "       was adjusted to:" /space = 0.
   ELSE IF (runnotes(i,1) = 6). 
      PRINT /title = "ERROR: An even number of variables in M list is required. ". 
      PRINT /title = "       No more than 20 M variables can be specified." /space = 0. 
   ELSE IF (runnotes(i,1) = 7). 
      PRINT copyname /title = "ERROR: Two of the specified variables are copies. The variable names are:" /format = A8.
   ELSE IF(runnotes(i,1) = 8). 
      PRINT /title = "ERROR: All specified variables must be unique. No variables may be the same in W, M, and Y". 
   ELSE IF (runnotes(i,1) = 9). 
      PRINT /title = "ERROR: At least one and no more than one W variable can be specified in the W list for Model 4 - 18.".
   ELSE IF (runnotes(i,1) = 10). 
      PRINT /title = "NOTE: Contrast command was specified with only 1 pair of mediators.".
      PRINT /title = "       No contrasts calculated." /space = 0.
   ELSE IF (runnotes(i,1) = 11). 
      PRINT /title = "ERROR: All variable names must have 8 characters or fewer.".
   ELSE IF (runnotes(i,1) = 12). 
      PRINT /title = "NOTE: Monte Carlo confidence intervals are not available for serial mediation". 
   ELSE IF (runnotes(i,1) = 13). 
      PRINT /title = "ERROR: The serial mediation model must have between 2 and 5 pairs of mediators.".
   ELSE IF (runnotes(i,1) = 14). 
      PRINT badboot /title = "NOTE: Some bootstrap samples had to be replaced.  The number of such replacements was:".
   ELSE IF (runnotes(i,1) = 15). 
      PRINT /title = "ERROR: At least one W variable, and no more than five W variables can be specified in the W list for Models 2 and 3.".
   ELSE IF (runnotes(i,1) = 16). 
      PRINT /title = "ERROR: Johnson-Neyman procedure not available for models with more than one moderator. ".
   ELSE IF (runnotes(i,1) = 17). 
      PRINT /title = "ERROR: Johnson-Neyman procedure not available for models with a dichotomous moderator. ".
   ELSE IF (runnotes(i,1) = 18). 
      PRINT /title = "ERROR: You must specify a value to probe at for each moderator in the w list.".
      PRINT /title = "       Wmodval lists should be the same length as w list." /space = 0.
   ELSE IF (runnotes(i,1) = 19). 
      PRINT /title = "ERROR: Sample size is not large enough. Model is unidentified.".
   ELSE IF (runnotes(i,1) = 20). 
      PRINT /title = "ERROR: Contrasts not available for serial mediation with more than 3 pairs of mediators.".
   ELSE IF (runnotes(i,1) = 21). 
      PRINT /title = "ERROR: Models 2, 3, and 4 require a variable name in the W argument.".
   ELSE IF (runnotes(i,1) = 22). 
      PRINT /title = "ERROR: Models 1 and 4 require at least two variable names in the M argument.".
   ELSE IF (runnotes(i,1) = 23). 
      PRINT /title = "ERROR: Invalid model number. Valid Models are numbered 1 through 18. Please see documentation for model descriptions".
   ELSE IF (runnotes(i,1) = 24). 
      PRINT /title = "NOTE: Centering command has no effect for Model 1.".
   ELSE IF (runnotes(i,1) = 25). 
      PRINT /title = "NOTE: Moderated mediation and serial option not available together. Parallel model estimated.".
   ELSE IF (runnotes(i,1) = 26). 
      PRINT /title = "ERROR: Contrasts not available for models with moderated indirect effects.".
   ELSE IF (runnotes(i, 1) = 27). 
       PRINT /title = "ERROR: Model 17 only involves moderation of the XM interaction, xmint cannot be set to zero for this model.".
       PRINT /title = "       Try Model = 1, xmint = 0 for a model with no moderation and no XM interaction.".
   ELSE IF (runnotes(i, 1) = 28). 
       PRINT /title = "ERROR: Invalid entry for CENTER option. Please select 0 (no centering), 1 (all moderators centered), ". 
       PRINT /title = "       or 2 (only continuous moderators centered). No centering conducted.".
   ELSE IF (runnotes(i, 1) = 29). 
       PRINT /title = "ERROR: Non-invertible design matrix. Results not produced. ". 
       PRINT /title = "       Check that variables have non-zero variance and are not perfectly correlated with each other.".    
   ELSE IF (runnotes(i,1) = 30). 
      PRINT /title = "NOTE: Plotting option only available for models with moderators. No plots generated.".
   END IF. 
END LOOP. 

DO IF (criterr = 0). 
    DO IF ((model = 1) OR (model > 3)). 
       DO IF (mc <>1 AND bc = 1). 
          print /title = "Bootstrap confidence interval method used: Bias corrected.".
       ELSE IF (mc <>1 AND bc <>1). 
          print /title = "Bootstrap confidence interval method used: Percentile bootstrap.".
       END IF.  
       DO IF (mc = 1). 
       print samples /title = "Number of samples for Monte Carlo condifidence intervals:".
       ELSE. 
       print samples /title = "Number of bootstrap samples for bootstrap confidence intervals:".
       END IF. 
       DO IF (xmint = 1). 
       COMPUTE centvars = MAKE(Mpairs, 6,0). 
          LOOP j = 1 TO Mpairs. 
             COMPUTE centvars(j,:) = { '(', mnamemat(j,1), ' + ', mnamemat(j,2), ')', '/2' }. 
          END LOOP. 
          print centvars /title = "The following variables were mean centered prior to analysis:" /format = A8. 
       END IF. 
       print conf /title = "Level of confidence for all confidence intervals in output:" /format = F10.2. 
       
       COMPUTE savelab = make(1, totsav,0). 
       COMPUTE blabs = {"b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10"}. 
       COMPUTE bwlabs = {"b1.1", "b2.1", "b3.1", "b4.1", "b5.1", "b6.1", "b7.1", "b8.1", "b9.1", "b10.1"}.
       COMPUTE bwintlabs = {"b1.2", "b2.2", "b3.2", "b4.2", "b5.2", "b6.2", "b7.2", "b8.2", "b9.2", "b10.2"}.
       COMPUTE awlabs = {"a1.0", "a2.0", "a3.0", "a4.0", "a5.0", "a6.0", "a7.0", "a8.0", "a9.0", "a10.0"}.
       COMPUTE awintlabs = {"a1.1", "a2.1", "a3.1", "a4.1", "a5.1", "a6.1", "a7.1", "a8.1", "a9.1", "a10.1"}.
       COMPUTE dlabs = {"d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9", "d10"}. 
       COMPUTE dwlabs = {"d1.1", "d2.1", "d3.1", "d4.1", "d5.1", "d6.1", "d7.1", "d8.1", "d9.1", "d10.1"}.
       COMPUTE dwintlabs = {"d1.2", "d2.2", "d3.2", "d4.2", "d5.2", "d6.2", "d7.2", "d8.2", "d9.2", "d10.2"}. 
       DO IF (mc = 0). 
           COMPUTE savelab(1,1) = {"c0"}. 
           COMPUTE counter = 2.
       ELSE IF (mc = 1). 
           COMPUTE counter = 1. 
       END IF. 
           DO IF ((anymod = 1) AND (mc=0)). 
               COMPUTE savelab(1,counter) = {"c1"}. 
               COMPUTE counter = counter+1. 
           END IF. 
           DO IF (serial = 0). 
              DO IF (apathmod = 0). 
                  COMPUTE savelab(1,counter:(counter+mpairs-1)) =  alabs(1, 1:Mpairs). 
                  COMPUTE counter = counter+Mpairs. 
              ELSE. 
                  COMPUTE alabtemp = {t(awlabs(1,1:Mpairs)), t(awintlabs(1,1:Mpairs))}. 
                  COMPUTE alabtemp = RESHAPE(alabtemp, 1, Mpairs*2). 
                  COMPUTE savelab(1, counter:(counter+2*mpairs-1)) = alabtemp. 
                  COMPUTE counter = counter+2*Mpairs. 
              END IF. 
              DO IF (mc = 0). 
                  DO IF (cppthmd = 0). 
                      COMPUTE savelab(1, counter) = {"cp"}. 
                      COMPUTE counter = counter + 1. 
                  ELSE. 
                      COMPUTE savelab(1, counter:(counter+1)) = {"cp0", "cp1"}. 
                      COMPUTE counter = counter + 2. 
                  END IF. 
              END IF. 
              DO IF (bpathmod = 0). 
                  COMPUTE savelab(1, counter:(counter+mpairs-1)) = blabs(1, 1:Mpairs). 
                  COMPUTE counter = counter+Mpairs. 
              ELSE. 
                  COMPUTE savelab(1, counter:(counter+mpairs-1)) = bwlabs(1, 1:Mpairs). 
                  COMPUTE counter = counter+Mpairs. 
                  COMPUTE savelab(1, counter:(counter+mpairs-1)) = bwintlabs(1, 1:Mpairs). 
                  COMPUTE counter = counter+Mpairs.
              END IF. 
              DO IF ((xmint = 1) AND (mc=0)). 
                  DO IF (dpathmod = 0). 
                      COMPUTE savelab(1, counter:(counter+mpairs-1)) = dlabs(1, 1:Mpairs). 
                      COMPUTE counter = counter+Mpairs. 
                  ELSE. 
                      COMPUTE savelab(1, counter:(counter+mpairs-1)) = dwlabs(1, 1:Mpairs). 
                      COMPUTE counter = counter+Mpairs. 
                      COMPUTE savelab(1, counter:(counter+mpairs-1)) = dwintlabs(1, 1:Mpairs). 
                      COMPUTE counter = counter+Mpairs.
                  END IF.  
             END IF. 
            ELSE. 
              COMPUTE savelab(1, 1:2) = {'c', 'a1'}. 
              DO IF (xmint = 1). 
                 COMPUTE serlab = {'a2', 'f1', 'g1', 'a3', 'f2', 'g2', 'f3', 'g3', 'a4', 'f4', 'g4', 'f5', 'g5', 'f6', 'g6', 'a5', 'f7', 'g7', 'f8', 'g8', 'f9', 'g9', 'f10', 'g10'}. 
                 COMPUTE savelab(1, 3:(1+mpairs**2)) = serlab(1, 1:(mpairs**2-1)). 
                 COMPUTE savelab(1, 2+mpairs**2) = {'cp'}. 
                 COMPUTE savelab(1, (3+mpairs**2):(2+mpairs**2+mpairs)) = blabs(1, 1:Mpairs). 
                 COMPUTE savelab(1, (3+mpairs**2+mpairs):(2+mpairs**2+2*mpairs)) = dlabs(1, 1:Mpairs). 
              ELSE IF (xmint = 0). 
                 COMPUTE serlab = {'a2', 'f1', 'a3', 'f2', 'f3', 'a4', 'f4',  'f5', 'f6','a5', 'f7', 'f8',  'f9', 'f10'}.  
                 COMPUTE savelab(1, 3:(2+cols)) = serlab(1, 1:cols). 
                 COMPUTE savelab(1, 3+cols) = {'cp'}. 
                 COMPUTE savelab(1, (4+cols):(3+cols+Mpairs)) = blabs(1, 1:Mpairs). 
              END IF. 
    
           END IF. 
       DO IF (mc = 1 AND !save = 1). 
          COMPUTE savelab = {savelab, "TotalInd"}. 
          SAVE mcsave /outfile =* /names = savelab.
       ELSE IF (mc <> 1 AND !save = 1).  
          SAVE bootsave /outfile =* /names = savelab.  
       END IF. 
    END IF. 
    
    
    DO IF ((model <> 1) AND (center > 0)). 
        DICHOT modcount = Wcount /dat = moddat. 
        DO IF (Wcount-(center=2)*csum(dich(:,1)) > 0). 
            COMPUTE centvars = MAKE(1, Wcount-(center=2)*csum(dich(:,1)), -999). 
            COMPUTE indx = 1. 
            LOOP i = 1 to Wcount. 
                DO IF ((center = 1) OR (dich(i,1) = 0)). 
                COMPUTE centvars(1,indx) = wnames(1,i). 
                COMPUTE indx = indx + 1. 
                END IF. 
            END LOOP. 
            print centvars /title = "The following variables were mean centered prior to analysis:" /format = A8. 
        ELSE. 
            print /title = "No variables were mean centered prior to analysis". 
        END IF. 
    END IF. 

END IF. 

end matrix. 
!ENDDEFINE. 
restore. 