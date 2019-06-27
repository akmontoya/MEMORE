* Encoding: UTF-8.
/* MEMORE for SPSS Version 2.1*/.
/* Copyright 2019 */.
/* by Amanda Kay Montoya */.
/* akmontoya.com*/.
/* Documentation available online at akmontoya.com */.
/* or by email to akmontoya@ucla.edu */.


preserve. 
set printback=off.

/* Permission is hereby granted, free of charge, to any person obtaining a copy of this software */.
/* and associated documentation files (the "Software"), to use the software in this form.  Distribution */.
/* after modification is prohibited, as is its use for any commercial purpose without authorization */.  
/* This software should not be posted or stored on any webpage, server, or directory accessible to */.
/* the public whether free or for a charge unless written permission has been granted by the copyright */.
/* holder.  The copyright holder requests that this software be distributed by directing users to */.
/* afhayes.com where the latest release of the software and documentation is archived and */.
/* can be downloaded.

/* THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF */.
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT */.
/* IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, */.
/*  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT */.
/* OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE */.
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE */.

/* The above text should be included in any distribution of the software */.


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

DEFINE MEMORE (Y = !charend('/') /M = !charend('/') !default(xxxxxxx)/W = !charend('/') !default(xxxxxxx)/Conf = !charend('/') !default(95) /mc = !charend('/') !default(0) 
   /samples = !charend('/') !default(5000) /normal = !charend('/') !default(0) /bc = !charend('/') !default(0) /decimals=!charend('/') !default(F10.4) /save = !charend('/') !default(0)
   /seed = !charend('/') !default(random) /contrast = !charend('/') !default(0) /xmint = !charend('/') !default(1) /serial = !charend('/') !default(0) /model = !charend('/') !default(1) /jn = !charend('/') !default(0) 
   /plot = !charend('/') !default(0) /quantile = !charend('/') !default(0) /center = !charend('/') !default(0) /wmodval1 = !charend('/') !default(999.99) /wmodval2 = !charend('/') !default(999.99)
   /wmodval3 = !charend('/') !default(999.99)). 
set mxloop = 100000000.
set seed = !seed. 
   
matrix. 
COMPUTE runnotes = MAKE(24,1,0). 
COMPUTE criterr = 0.  
COMPUTE model = !model. 
COMPUTE modelmat = MAKE(3,1,1)*model. 
COMPUTE modelmt2 = {1; 2; 3}. 
COMPUTE validtst = (Modelmat = modelmt2). 
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
DO IF ((m = "xxxxxxx") and (model = 1)). 
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
   end if. 
   COMPUTE missing = nrow(fulldat) - nrow(data). 
end if. 

COMPUTE mc = (!mc=1). 
COMPUTE serial = (!serial = 1). 
COMPUTE jn = (!jn = 1). 
COMPUTE plot = (!plot = 1). 
COMPUTE quantile = (!quantile = 1). 
COMPUTE center = (!center = 1). 
COMPUTE contrast = (!contrast = 1). 
COMPUTE normal = (!normal = 1). 
COMPUTE xmint = (!xmint = 1). 

DO IF (((model = 2) OR (model = 3)) AND (ncol(wnames) > 1) AND (jn = 1)). 
      COMPUTE jn = 0. 
      COMPUTE runnotes(16,1) = 16. 
END IF.  



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

DO IF ((model =2) OR (model = 3)). 
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

DO IF (((Mcount = 0) OR (Mcount = 1) OR mod(Mcount,2) > 0) OR (Mcount >20)) AND (criterr <> 1) AND (model =1). 
   COMPUTE runnotes(6, 1) = 6. 
   COMPUTE criterr = 1. 
END IF. 

DO IF (((model = 2) OR (model = 3)) AND ((Wcount = 0) OR (Wcount > 5))). 
   COMPUTE runnotes(15,1) = 15. 
   COMPUTE criterr = 1. 
END IF. 

DO IF (criterr = 0). 

DO IF (model = 1). 
   LOOP i = 1 TO Mcount. 
      DO IF ((mnames(1,i) = ynames(1,1)) OR (mnames(1,i) = ynames(1,2))). 
      COMPUTE runnotes(8,1) = 8. 
      COMPUTE criterr = 1. 
      END IF. 
   END LOOP. 
END IF. 

DO IF ((model = 2) OR (model = 3)). 
   LOOP i = 1 TO Wcount. 
      DO IF ((wnames(1,i) = ynames(1,1)) OR (wnames(1,i) = ynames(1,2))). 
      COMPUTE runnotes(8,1) = 8. 
      COMPUTE criterr = 1. 
      END IF. 
   END LOOP. 
   END IF. 
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

COMPUTE bc = trunc(!bc). 
DO IF (mc = 1 AND bc = 1).  
   COMPUTE runnotes(9,1) = 9. 
   COMPUTE bc = 0. 
END IF. 

DO IF (contrast = 1 AND (Mcount/2) = 1) AND (model = 1). 
   COMPUTE runnotes(10,1) = 10. 
END IF. 

DO IF ((contrast = 1) AND (serial = 1) AND (Mcount > 6) AND (model = 1)). 
   COMPUTE contrast = 0. 
   COMPUTE runnotes(20,1) = 20. 
END IF. 

DO IF (model = 1 AND serial = 1) AND ((Mcount < 4) OR (Mcount > 10)). 
   COMPUTE runnotes(13,1) = 13.
   COMPUTE criterr = 1. 
END IF. 

DO IF ((model = 1) AND (center = 1)). 
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
   END IF. 
   DO IF (N < Nfail). 
      COMPUTE runnotes(19,1) = 19. 
      COMPUTE criterr = 1. 
   END IF. 
END IF. 

DO IF (criterr = 0). 
DO IF (model = 1). 
   COMPUTE Mpairs = Mcount/2.  
   COMPUTE mnamemat = reshape(mnames, Mpairs, 2).
   COMPUTE transmat = {1,1/2;-1,1/2}. 
   COMPUTE Tmat = make(ncol(data), ncol(data), 0).  
   LOOP i = 1 TO (2*Mpairs+1) BY 2. 
      COMPUTE Tmat(i:(i+1), i:(i+1)) = transmat. 
   END LOOP. 
   COMPUTE dataT = data*Tmat.
   COMPUTE select = {1,3,5,7,9,11, 13, 15, 17, 19, 21}.
   DO IF (xmint = 0). 
      COMPUTE dataT = dataT(: , select(1:(Mpairs+1))). 
   ELSE IF (xmint = 1). 
      COMPUTE dataT = dataT(:,1:(ncol(dataT)-1)). 
   END IF. 
ELSE IF (model = 2). 
   COMPUTE Mpairs = 0. 
   COMPUTE moddat = data(:,1:Wcount).
   COMPUTE tempvec = MAKE(Wcount,1, 0). 
   COMPUTE tempvec = {tempvec; 1; -1}. 
   COMPUTE transmat = {IDENT(Wcount+2, Wcount), tempvec}. 
   COMPUTE dataT = data*transmat. 
   DO IF (center = 1). 
      COMPUTE centmean = mdiag(csum(dataT)/N). 
      COMPUTE dataT(:,1:Wcount) = dataT(:,1:Wcount) - MAKE(N,Wcount,1)*centmean(1:Wcount, 1:Wcount). 
      COMPUTE moddat = dataT(:,1:Wcount). 
   END IF. 
ELSE IF (model = 3). 
   COMPUTE Mpairs = 0. 
   COMPUTE moddat = data(:,1:Wcount).
   DO IF (center = 1). 
      DO IF (Wcount = 1). 
         COMPUTE centmean = csum(moddat)/N. 
      ELSE.
         COMPUTE centmean = mdiag(csum(moddat)/N). 
      END IF. 
      COMPUTE moddat = moddat - MAKE(N,Wcount,1)*centmean.
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
DO IF (model = 1). 
      COMPUTE aresmat = MAKE (Mpairs, 6, 0). 
      COMPUTE ghostdes = make(N,1, 1).
      CDFINVT p = temp/ df = (N-1). 
      COMPUTE tcrita = toutput. 
      COMPUTE tcritc = toutput.
      COMPUTE counterj = 1. 
      LOOP j = 1 TO (ncol(dataT)-1) by (1+xmint). 
         DO IF (xmint = 1). 
            COMPUTE summean = csum(dataT(:,j+1))/N. 
            COMPUTE dataT(:,j+1) = (dataT(:,j+1) - summean).
         END IF.    
         COMPUTE apath = inv(t(ghostdes)*ghostdes)*t(ghostdes)*dataT(:,j). 
         COMPUTE seMdiff = sqrt((csum((dataT(:,j)-(csum(dataT(:,j))/N))&**2)/(N-1))/N).
         COMPUTE tapath = apath/seMdiff.  
         COMPUTE dfapath = N-1. 
         COMPUTE papath = 2*(1-tcdf(abs(tapath), dfapath)). 
         COMPUTE lcia = apath-tcrita*seMdiff. 
         COMPUTE ucia = apath+tcrita*seMdiff. 
         COMPUTE aresmat(counterj,:) = {apath, seMdiff, tapath, papath, lcia, ucia}. 
         COMPUTE counterj = counterj + 1. 
      END LOOP. 

      COMPUTE ccols = ncol(dataT).
      COMPUTE cols = 0. 
      COMPUTE cpath = inv(t(ghostdes)*ghostdes)*t(ghostdes)*dataT(:,ccols). 
      COMPUTE seYdiff = sqrt(csum((dataT(:,ccols)-(csum(dataT(:,ccols))/N))&**2)/(N-1)/N).
      COMPUTE tcpath = cpath/seYdiff. 
      COMPUTE dfcpath = N-1. 
      COMPUTE pcpath = 2*(1-tcdf(abs(tcpath), N-1)). 
      COMPUTE lcic = cpath-tcritc*seYdiff. 
      COMPUTE ucic = cpath+tcritc*seYdiff. 
      COMPUTE cresmat = {cpath, seYdiff, tcpath, pcpath, lcic, ucic}. 

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
            COMPUTE M2modbs = inv(t(serdes)*serdes)*t(serdes)*dataT(:,j). 
            COMPUTE M2pred = serdes*m2modbs. 
            COMPUTE M2ssr = csum((dataT(:,j)-M2pred)&**2). 
            COMPUTE M2sst = csum((dataT(:,j) - csum(dataT(:,j))/N)&**2). 
            COMPUTE M2Rsq = 1-M2ssr/m2sst. 
            COMPUTE M2r = sqrt(M2Rsq).    
            COMPUTE M2msr = M2ssr/(N - ncol(serdes)). 
            COMPUTE M2df1 = ncol(serdes) - 1. 
            COMPUTE M2df2 = (N - ncol(serdes)). 
            COMPUTE M2F = m2df2*m2rsq/(m2df1*(1-m2rsq)). 
            COMPUTE M2p = 1-FCDF(M2F, M2df1, m2df2). 
            COMPUTE sem2bmat = (m2msr*inv(t(serdes)*serdes)). 
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

/*for all models. 
COMPUTE bcpdes = {make(N,1,1), dataT(:,1:(ncol(dataT)-1))}.
COMPUTE bcpvec = inv(t(bcpdes)*bcpdes)*t(bcpdes)*dataT(:,ncol(dataT)).
COMPUTE ypred = bcpdes*bcpvec. 
COMPUTE ssr = csum((dataT(:,ncol(dataT)) - ypred)&**2). 
COMPUTE sst = csum((dataT(:,ncol(dataT)) - csum(dataT(:,ncol(dataT)))/N)&**2).
COMPUTE msr = ssr/(N-ncol(bcpdes)). 
COMPUTE rsqfull = 1- ssr/sst. 
COMPUTE rfull = sqrt(rsqfull). 
COMPUTE df1 = (ncol(bcpdes)-1). 
COMPUTE df2 = (N - ncol(bcpdes)). 
COMPUTE Ffull = df2*Rsqfull/((df1)*(1-rsqfull)). 
COMPUTE pfull =1- FCDF(Ffull, df1, df2). 
COMPUTE sebcpmat = (msr*inv(t(bcpdes)*bcpdes)). 
COMPUTE sebcp = (diag(sebcpmat))&**(1/2). 
COMPUTE modsumr = {Rfull, Rsqfull, MSR, Ffull, df1, df2, pfull}. 

CDFINVT p = temp/ df = df2. 
COMPUTE tcritb = toutput. 
COMPUTE tcritcp = toutput. 
COMPUTE tcritd = toutput. 

COMPUTE LCII = rnd((1-.01*Conf)/2*samples). 
COMPUTE UCII = trunc((1-((1-.01*Conf)/2))*samples)+1. 
DO IF (LCII  < 1 OR UCII > samples). 
   COMPUTE runnotes(4, 1) = 4.  
   COMPUTE criterr = 1. 
   COMPUTE LCII = 1. 
   COMPUTE UCII = samples. 
END IF. 

DO IF (model = 1). 
      COMPUTE bresmat = MAKE(Mpairs, 6, 0). 
      COMPUTE dresmat = MAKE(Mpairs, 6, -999). 
      COMPUTE serind = 0. 
      LOOP i = 2 to Mpairs. 
            CHOOSE r = Mpairs /k = i. 
            COMPUTE serind = serind + rchoosek. 
      END LOOP. 
      COMPUTE indres = MAKE(Mpairs+1+serind*(serial=1),1,0). 
      DO IF (normal = 1). 
         COMPUTE normres = MAKE(Mpairs, 4, 0). 
      END IF. 
      COMPUTE cppath = bcpvec(1,1). 
      COMPUTE secppath = sebcp(1,1). 
      COMPUTE tcppath = cppath/secppath.
      COMPUTE pcppath = 2*(1-tcdf(abs(tcppath), df2)).
      COMPUTE lcicp = cppath-tcritcp*secppath. 
      COMPUTE ucicp = cppath+tcritcp*secppath. 
      COMPUTE cpresmat = {cppath, secppath, tcppath, pcppath, lcicp, ucicp}. 
      /*MC Setup. 
      DO IF (mc = 1). 
      COMPUTE mcsamps = samples. 
      COMPUTE randsamp = sqrt(-2*ln(uniform(mcsamps,Mpairs)))&*cos((2*3.14159265358979)*uniform(mcsamps,Mpairs)).
      COMPUTE MCres = MAKE (Mpairs+1, 4, 0). 
      COMPUTE MCcorr = MAKE(Mpairs, Mpairs, 1). 
      COMPUTE counteri = 1. 
      LOOP i = 2 TO (nrow(sebcpmat)-xmint) by (1+xmint). 
         COMPUTE MCcorr(counteri,counteri) = sebcpmat(i,i).
         DO IF ((Mpairs > 1) AND (counteri <> Mpairs)). 
         COMPUTE counterj = counteri + 1. 
         LOOP j = (i+1+xmint) TO (nrow(sebcpmat)-xmint) by (1+xmint).
            COMPUTE MCcorr(counteri,counterj) = sebcpmat(i,j). 
            COMPUTE MCcorr(counterj,counteri) = sebcpmat(i,j). 
            COMPUTE counterj = counterj + 1. 
         END LOOP. 
         END IF. 
         COMPUTE counteri = counteri + 1. 
      END LOOP.
      COMPUTE rndnb = randsamp*chol(MCcorr).  
      COMPUTE rndna = sqrt(-2*ln(uniform(mcsamps,Mpairs)))&*cos((2*3.14159265358979)*uniform(mcsamps,Mpairs)).
      COMPUTE mcsave = make(samples, 3*Mpairs+1, 0). 
      COMPUTE mcsave2 = make(samples, Mpairs, 0). 
      END IF. 

      COMPUTE counteri = 1. 
      LOOP i = 2 TO (nrow(bcpvec)-xmint) by (1+xmint). 
         COMPUTE bpath = bcpvec(i,1). 
         COMPUTE sebpath = sebcp(i,1). 
         COMPUTE tbpath = bpath/sebpath. 
         COMPUTE pbpath = 2*(1-tcdf(abs(tbpath), df2)).
         COMPUTE lcib = bpath-tcritb*sebpath. 
         COMPUTE ucib = bpath+tcritb*sebpath. 
         COMPUTE bresmat(counteri, :) = {bpath, sebpath, tbpath, pbpath, lcib, ucib}. 

         DO IF (xmint = 1). 
            COMPUTE dpath = bcpvec((i+1), 1). 
            COMPUTE sedpath = sebcp((i+1),1). 
            COMPUTE tdpath = dpath/sedpath. 
            COMPUTE pdpath = 2*(1-tcdf(abs(tdpath),df2)). 
            COMPUTE lcid = dpath-tcritd*sedpath. 
            COMPUTE ucid = dpath+tcritd*sedpath. 
            COMPUTE dresmat(counteri, :) = {dpath, sedpath, tdpath, pdpath, lcid, ucid}. 
         END IF. 

         /*Normal Tests. 
         COMPUTE indirect = aresmat(counteri,1)*bresmat(counteri,1). 
         COMPUTE indres(counteri,1) = indirect.  
         DO IF (normal = 1). 
            COMPUTE sobseab = sqrt((aresmat(counteri,1)**2)*(bresmat(counteri,2)**2)+(bresmat(counteri,1)**2)*(aresmat(counteri,2)**2)). 
            COMPUTE sobelZ = indirect/sobseab.
            COMPUTE sobelp = 2*cdfnorm((-1)*abs(sobelZ)). 
            COMPUTE normres(counteri,:) = {indirect, sobseab, sobelZ, sobelp}. 
         END IF. 

         /*Monte Carlo Confidence Interval*/.

         DO IF (mc = 1). 
         COMPUTE asamp = rndna(:,counteri)*aresmat(counteri,2)+aresmat(counteri,1). 
         COMPUTE bsamp = rndnb(:,counteri)+bresmat(counteri,1). 
         COMPUTE absamp = asamp&*bsamp. 
         COMPUTE mcgrad = grade(absamp). 
         COMPUTE mcsort = absamp. 
         COMPUTE mcsort(mcgrad) = absamp.
         COMPUTE mcsave(:,(3*counteri-2):(3*counteri)) = {asamp, bsamp,absamp}. 
         COMPUTE mcsave2(:,counteri) = absamp. 
         COMPUTE MCLLCI = mcsort(LCII). 
         COMPUTE MCULCI = mcsort(UCII). 
         COMPUTE seMC = sqrt(csum((mcsort(:,1)-(csum(mcsort(:,1))/mcsamps))&**2)/(mcsamps-1)).
         COMPUTE MCres(counteri,:) = {indirect, seMC, MCLLCI, MCULCI}.  
         END IF. 
         COMPUTE counteri = counteri + 1. 
      END LOOP. 
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
         COMPUTE serialse = rsum(indse). 
         COMPUTE sobelZ = serlind&/serialse. 
         COMPUTE sobelp = 2*cdfnorm((-1)*abs(sobelZ)). 
         COMPUTE serindn = {serlind, serialse, sobelZ, sobelp}. 
         COMPUTE normres = {normres; serindn}.
       END IF. 

   END IF. 
      /*Total monte carlo*/.

    DO IF (mc = 1). 
      COMPUTE mcsave(:,3*Mpairs+1) = rsum(mcsave2(:,:)). 
      COMPUTE mcsort = mcsave(:,3*Mpairs+1). 
      COMPUTE mcgrad = grade(mcsort). 
      COMPUTE mcsort(mcgrad) = mcsort. 
      COMPUTE MCLLCI = mcsort(LCII). 
      COMPUTE MCULCI = mcsort(UCII). 
      COMPUTE seMC = sqrt(csum((mcsort(:,1)-(csum(mcsort(:,1))/mcsamps))&**2)/(mcsamps-1)).
      COMPUTE MCres(Mpairs+1, :) = {csum(indres), seMC, MCLLCI, MCULCI}. 

      DO IF ((contrast = 1) AND (Mpairs >1)). 
            COMPUTE npairs = Mpairs*(Mpairs-1)/2. 
            COMPUTE contres = MAKE(npairs, 4,0). 
            COMPUTE contsamp = MAKE(samples, npairs, 0). 
            COMPUTE contsort = contsamp. 
            COMPUTE counter = 1. 
            LOOP i = 1 TO Mpairs-1. 
               LOOP j = i+1 TO Mpairs. 
               COMPUTE contsamp(:,counter) = mcsave2(:,i) - mcsave2(:,j). 
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
         COMPUTE detcheck = make(((Mpairs-1)*(serial=1) +1), 1, -999). 
         COMPUTE Bootsamp = make(samples, Mpairs+1+serind*(serial=1), 0). 
         DO IF (serial = 1). 
               DO IF (xmint = 1). 
                  COMPUTE bootcol = mpairs**2-1. 
               ELSE IF (xmint = 0). 
                  COMPUTE bootcol = ((mpairs - 1)**2 + 3*(mpairs-1))/2. 
               END IF. 
         ELSE IF (serial = 0). 
            COMPUTE bootcol = -999. 
         END IF. 
         COMPUTE Bootsave = make(samples, 3*Mpairs+3+(bootcol+serind-mpairs+1)*(serial=1),0). 
         COMPUTE indtemp = make(samples,Mpairs+serind*(serial=1),0).
         LOOP i = 1 TO samples.
            LOOP k = 1 TO 10000. 
               COMPUTE sortvar = trunc(uniform(N,1)*N)+1.
               COMPUTE bootdat = dataT(sortvar(:,1),:).  
               DO IF (xmint = 1).
                  LOOP j = 1 TO Mpairs. 
                  COMPUTE summean = csum(bootdat(:,2*j))/N.
                  COMPUTE bootdat(:,2*j) = (bootdat(:,2*j) - summean).
                  END LOOP.
               END IF. 
               COMPUTE bootdes = {make(N,1,1), bootdat(:,1:(ncol(bootdat)-1))}. 
               COMPUTE detcheck(1,1) = (det(t(bootdes)*bootdes)=0).
               DO IF (serial = 1). 
                  COMPUTE counterj = 2. 
                  LOOP j = (2+xmint) TO (ncol(dataT)-1) by (1+xmint). . 
                          COMPUTE bootadat = {make(N,1,1), bootdat(:, 1:(j-1))}. 
                          COMPUTE detcheck(counterj,1) = (det(t(bootadat)*bootadat)=0). 
                  COMPUTE counterj = counterj + 1. 
                  END LOOP. 
               END IF.                  
               COMPUTE badboot = badboot+(k = 2). 
            END LOOP IF (csum(detcheck(:,1)) = 0). 
            COMPUTE bootbeta = inv(t(bootdes)*bootdes)*t(bootdes)*bootdat(:,ncol(bootdat)). 
            COMPUTE bootavec = MAKE(mpairs, 1, -999).  
            COMPUTE bootbvec = MAKE(mpairs,1,-999). 
            DO IF (serial = 1). 
               COMPUTE bootser = make(bootcol, 1, -999).
            END IF. 
            COMPUTE start = 1.  
            COMPUTE counterj = 1. 
            LOOP j = 1 TO (ncol(dataT)-1) by (1+xmint).  
               COMPUTE bootavec(counterj,1) = inv(t(ghostdes)*ghostdes)*t(ghostdes)*bootdat(:,j).
               COMPUTE bootb = bootbeta(j+1,1).
               COMPUTE bootbvec(counterj,1) = bootb. 
               DO IF ((serial = 1)AND(j>1)).
                  COMPUTE serdes = {make(N,1,1), bootdat(:, 1:(j-1))}.  
                  DO IF (xmint = 1).
                        COMPUTE end = (counterj)**2 - 1. 
                  ELSE IF (xmint = 0).
                        COMPUTE end = ((counterj-1)**2 + 3*(counterj-1))/2. 
                  END IF.
                  COMPUTE bootser(start:end,1) = inv(t(serdes)*serdes)*t(serdes)*bootdat(:,j). 
                  COMPUTE bootavec(counterj,1) = bootser(start,1). 
                 COMPUTE  bootsave(i,(2*counterj):(2*counterj+1)) = {bootb, bootavec(counterj,1)*bootb}. 
                  COMPUTE start = end +1. 
               END IF. 
                  COMPUTE bootsamp(i,counterj) = bootb*bootavec(counterj,1).
                 DO IF ((serial = 0) OR (j = 1)). 
                  COMPUTE bootsave(i,(3*counterj-2):(3*counterj)) = {bootavec(counterj,1), bootb, bootavec(counterj,1)*bootb}.
                END IF. 
                  COMPUTE indtemp(i,counterj) = (bootavec(counterj,1)*bootb).  
            COMPUTE counterj = counterj + 1. 
           END LOOP. 

            /*Serial loop goes here*/. 
              DO IF (serial = 1). 
               COMPUTE bootsave(i, (4+2*(mpairs-1)):(3+2*(mpairs-1) + bootcol)) = t(bootser). 
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
                                 COMPUTE indtemp(i,mpairs+counter) = step5*bootbvec(5,1).
                                 COMPUTE bootsamp(i, mpairs + counter) = step5*bootbvec(5,1).
                                 COMPUTE bootsave(i, (1+2*mpairs + bootcol +counter)) = step5*bootbvec(5,1).
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
                                      COMPUTE indtemp(i,mpairs+counter) = step4*bootbvec(bindx,1).
                                      COMPUTE bootsamp(i, mpairs + counter) = step4*bootbvec(bindx,1).
                                      COMPUTE bootsave(i, (1+2*mpairs + bootcol +counter)) = step4*bootbvec(bindx,1).
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
                                COMPUTE indtemp(i,mpairs+counter) = step3*bootbvec(bindx,1).
                                COMPUTE bootsamp(i, mpairs + counter) = step3*bootbvec(bindx,1).
                                COMPUTE bootsave(i, (1+2*mpairs + bootcol +counter)) = step3*bootbvec(bindx,1).
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
                           COMPUTE indtemp(i,mpairs+counter) = step2*bootbvec(bindx,1).
                           COMPUTE bootsamp(i, mpairs + counter) = step2*bootbvec(bindx,1).
                           COMPUTE bootsave(i, (1+2*mpairs + bootcol +counter)) = step2*bootbvec(bindx,1).
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

            /*direct effect*/. 
            COMPUTE bootsave(i,ncol(bootsave)-1) = bootbeta(1,1). 
            /*total effect*/.
            COMPUTE bootsave(i,ncol(bootsave)) = rsum(indtemp(i,:))+bootbeta(1,1). 
         END LOOP. 

         DO IF (badboot >0). 
            COMPUTE runnotes(14,1) = 14. 
         END IF. 

         DO IF ((contrast = 1) AND (Mpairs >1)). 
            COMPUTE npairs = ncol(indtemp)*(ncol(indtemp)-1)/2. 
            COMPUTE contres = MAKE(npairs, 4,0). 
            COMPUTE contsamp = MAKE(samples, npairs, 0). 
            COMPUTE counter = 1. 
            LOOP i = 1 TO ncol(indtemp)-1. 
               LOOP j = i+1 TO ncol(indtemp). 
               COMPUTE contsamp(:,counter) = indtemp(:,i) - indtemp(:,j). 
               COMPUTE contres(counter, 1) = indres(i,1) - indres(j,1). 
               COMPUTE counter = counter+1. 
               END LOOP. 
            END LOOP. 
         END IF. 

         COMPUTE bootsamp(:,ncol(bootsamp)) = rsum(bootsamp(:,1:(ncol(bootsamp)-1))).
         COMPUTE bootsave(:,(ncol(bootsave)-2)) =  rsum(bootsamp(:,1:(ncol(bootsamp)-1))).

         COMPUTE indres(nrow(indres),1) = csum(indres). 
         COMPUTE bootsort = bootsamp. 
         COMPUTE seboots = MAKE(Mpairs+1+serind*(serial=1), 1, 0). 
         COMPUTE bccires = MAKE(4,Mpairs+1+serind*(serial=1), 0). 
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
      END IF. 

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
         END IF. 


ELSE IF ((model = 2) OR (model = 3)). 
   COMPUTE wnamemat = {t(wnames), MAKE(Wcount, 1, " ")}. 
   COMPUTE modres = MAKE(ncol(bcpdes), 6, -999). 
   COMPUTE modres(:,1) = bcpvec. 
   COMPUTE modres(:,2) = sebcp. 
   COMPUTE modres(:,3) = modres(:,1)/modres(:,2).  
   COMPUTE modres(:,4) = 2*(1-tcdf(abs(modres(:,3)), df2)).
   COMPUTE modres(:,5:6) = {modres(:,1)-tcritb*modres(:,2), modres(:,1)+tcritb*modres(:,2)}. 

   /*Conditional effect of X on Y at values of W*/. 
   COMPUTE dich = MAKE(Wcount, 3, -999). 
   LOOP i = 1 to Wcount. 
      COMPUTE uniqdes = DESIGN(moddat(:,i)). 
      COMPUTE dich(i,1) = (ncol(uniqdes) = 2). 
      DO IF (dich(i,1) = 1). 
            COMPUTE dichsort = MAKE(N, 1, -999).
            COMPUTE dichgrad = grade(moddat(:,i)). 
            COMPUTE dichsort(dichgrad,1) = moddat(:,i).
            COMPUTE dich(i,2) = dichsort(1,1). 
            COMPUTE dich(i,3) = dichsort(N,1).
      END IF. 
   END LOOP. 

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
   END IF. 
   COMPUTE XYgWres = {modcomb, MAKE(nrow(modcomb), 6, -999)}. 
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
      COMPUTE modcomb = {MAKE(nrow(modcomb),1, 1), modcomb}.
      COMPUTE XYgWres(:,Wcount+1) = modcomb*bcpvec. 
      COMPUTE XYgWres(:,Wcount+2) = sqrt(diag(modcomb*sebcpmat*t(modcomb))). 
      COMPUTE XYgWres(:,Wcount+3) = XYgWres(:,Wcount+1)/XYgWres(:,Wcount+2).  
      COMPUTE XYgWres(:,Wcount+4) = 2*(1-tcdf(abs(XYgWres(:,Wcount+3)), df2)).
      COMPUTE XYgWres(:,(Wcount+5):(Wcount+6)) = {XYgWres(:,Wcount+1)-tcritb*XYgWres(:,Wcount+2), XYgWres(:,Wcount+1)+tcritb*XYgWres(:,Wcount+2)}.

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
         COMPUTE prbrfull = sqrt(prsqfull). 
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
 

print /title = "********************* MEMORE Procedure for SPSS Version 2.1 ***********************".
print /title = "                           Written by Amanda Montoya       ".
print /title = "                    Documentation available at akmontoya.com ".
print /title = "**************************************************************************************".
DO IF (criterr = 0).

print model /title = "Model:". 
DO IF (Model = 1 AND Mpairs = 1) .
COMPUTE varrlabs = {'Y = ', 'M = '}. 
ELSE IF (Model <>1 AND Wcount = 1). 
COMPUTE varrlabs = {'Y = ', 'W = '}. 
ELSE IF (Model = 1 AND Mpairs > 1). 
COMPUTE varrlabs = {'Y = ', 'M1 = ', 'M2 = ', 'M3 = ', 'M4 = ', 'M5 = ', 'M6 = ', 'M7 = ', 'M8 = ', 'M9 = ', 'M10 = '}.
ELSE IF (Model <>1 AND Wcount > 1). 
COMPUTE varrlabs = {'Y = ', 'W1 = ', 'W2 = ', 'W3 = ', 'W4 = ', 'W5 = '}.
END IF. 

DO IF (model = 1). 
print {ynames; mnamemat} /title = "Variables: " /rnames = varrlabs  /format = a8. 
ELSE IF (Model <> 1). 
print {ynames; wnamemat} /title = "Variables: " /rnames = varrlabs  /format = a8. 
END IF. 

DO IF (model = 1). 
      COMPUTE compname = MAKE(((1+xmint)*Mpairs+1),7, 0). 
      COMPUTE compname(1,:) = {' ', ynames(1,1), ' - ', ynames(1,2), ' ', ' ', ' '}.
      LOOP j = 1 TO Mpairs.
         COMPUTE compname((1+j),:) ={' ', mnamemat(j,1), ' - ', mnamemat(j,2), ' ', ' ', ' '}. 
           DO IF (xmint = 1). 
               COMPUTE compname((Mpairs+1+j),:) = { '(', mnamemat(j,1), ' + ', mnamemat(j,2), ')', '/2', 'Centered'}. 
         END IF. 
      END LOOP. 
      COMPUTE temp1 = {'M1diff = ','M2diff = ','M3diff = ', 'M4diff = ', 'M5diff = ','M6diff = ','M7diff = ','M8diff = ', 'M9diff = ', 'M10diff = '}. 
      COMPUTE temp2 = {'M1avg  = ','M2avg  = ','M3avg  = ', 'M4avg  = ', 'M5avg  = ','M6avg  = ','M7avg  = ','M8avg  = ', 'M9avg  = ', 'M10avg  = '}.
      DO IF (xmint = 1). 
            DO IF (Mpairs = 1). 
               COMPUTE temprnam = {'Ydiff = ', 'Mdiff = ', 'Mavg = '}. 
            ELSE. 
               compute temprnam = {'Ydiff = ', temp1(1,1:Mpairs), temp2(1,1:Mpairs)}. 
            END IF. 
      ELSE IF (xmint = 0). 
            DO IF (Mpairs = 1). 
               COMPUTE temprnam = {'Ydiff = ', 'Mdiff = '}. 
            ELSE. 
               compute temprnam = {'Ydiff = ', temp1(1,1:Mpairs)}. 
            END IF. 
     END IF. 
ELSE IF ((model = 2) OR (Wcount = 1)). 
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
DO IF (model = 1). 
   print {"Ydiff =" , ynames(1,1), ' - ', ynames(1,2)} /title = "**************************************************************************************" /rlabels = "Outcome:" /format = A8.
   COMPUTE collab = {"Effect", "SE", "t", "p", "LLCI", "ULCI"}. 
   print cresmat /title = "Model" /rnames = {"'X'"} /cnames = collab /format = !decimals. 
   print dfcpath /title = "Degrees of freedom for all regression coefficient estimates:".
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
            COMPUTE m2labs = {"'X'", "M1diff", "M1avg", "M2diff", "M2avg", "M3diff", "M3avg", "M4diff", "M4avg"}. 
            COMPUTE end = j**2 - 1. 
        ELSE IF (xmint = 0). 
            COMPUTE m2labs = {"'X'", "M1diff", "M2diff", "M3diff", "M4diff"}. 
            COMPUTE end = ((j-1)**2 + 3*(j-1))/2.  
         END IF. 
         print serres(start:end,:) /title "Model" /rnames = m2labs /clabels = "coeff" , "SE", "t", "p", "LLCI", "ULCI" /format = !decimals.
         COMPUTE start = end+1. 
         print smodsum(j-1,6) /title = "Degrees of freedom for all regression coefficient estimates:".
      ELSE. 
      print aresmat(j,:) /title = "Model" /rnames = {"'X'"} /cnames = collab /format = !decimals. 
      print dfapath /title = "Degrees of freedom for all regression coefficient estimates:".
      END IF. 
   END LOOP. 

END IF. 
print {"Ydiff =" , ynames(1,1), ' - ', ynames(1,2)} /title = "**************************************************************************************" /rlabels = "Outcome:" /format = A8.
print modsumr /title = "Model Summary" /clabels = "R", "R-sq", "MSE", "F" , "df1" , "df2", "p" /format = !decimals. 

DO IF (model = 1). 
      DO IF (xmint = 1). 
      COMPUTE modres = {cpresmat;bresmat;dresmat}. 
      ELSE IF (xmint = 0). 
       COMPUTE modres = {cpresmat;bresmat}. 
      END IF. 
   COMPUTE bdlabs = {"M1diff", "M2diff", "M3diff", "M4diff", "M5diff", "M6diff", "M7diff", "M8diff", "M9diff", "M10diff"}.
   COMPUTE bslabs = {"M1avg", "M2avg", "M3avg", "M4avg", "M5avg", "M6avg", "M7avg", "M8avg", "M9avg", "M10avg"}.
   DO IF (Mpairs = 1). 
      COMPUTE modlabs = {"'X'", "Mdiff", "Mavg"}. 
      COMPUTE bdlabs = {"Mdiff"}. 
   ELSE. 
      COMPUTE modlabs = {"'X'", bdlabs(1, 1:Mpairs), bslabs(1, 1:Mpairs)}. 
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
DO IF (model = 1). 
   print /title = "************************* TOTAL, DIRECT, AND INDIRECT EFFECTS *************************" .
   COMPUTE collab = {"Effect", "SE", "t", "df", "p", "LLCI", "ULCI"}. 
   COMPUTE cresmat = {cresmat(1:3), dfcpath, cresmat(4:6)}. 
   COMPUTE cpresmat = {cpresmat(1:3), df2, cpresmat(4:6)}.
   print cresmat /title = "Total effect of X on Y" /cnames = collab /format = !decimals. 
   print cpresmat /title = "Direct effect of X on Y" /cnames = collab /format = !decimals. 
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
   COMPUTE m2lab = {mlab(1,1:(nrow(indres)-1)), 'Total'}. 

   DO IF (Mpairs = 1). 
   print indres(1,:) /title = "Indirect Effect of X on Y through M" /rnames = m2lab / cnames = indlabs /format = !decimals. 
   ELSE. 
   print indres /title = "Indirect Effect of X on Y through M" /rnames = m2lab / cnames = indlabs /format = !decimals. 
   END IF. 

   DO IF (normal = 1). 
   print normres /title = "Normal Theory Tests for Indirect Effect" /rnames = mlab /clabels = "Effect", "SE", "Z", "p" /format = !decimals. 
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
ELSE IF ((Model = 2) OR (Model = 3)). 
   print /title = "**************************************************************************************" .
   COMPUTE coeflabs = {"Effect" , "SE", "t", "p", "LLCI", "ULCI"}.
   COMPUTE XYgWlabs = {wnames, coeflabs}.  
   print XYgWres /title = "Conditional Effect of 'X' on Y at values of moderator(s)" /cnames = XYgWlabs /format = !decimals. 
   print df2 /title = "Degrees of freedom for all conditional effects:".
DO IF (quantile = 1). 
   print /title = "Values for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentile.".
ELSE IF (quantile = 0). 
   print /title = "Values for quantitative moderators are the mean and plus/minus one SD from the mean.".
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
print /title = "****************************** ANALYSIS NOTES AND WARNINGS *******************************". 
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
      PRINT /title = "ERROR: An even number of variables in M list is required.". 
   ELSE IF (runnotes(i,1) = 7). 
      PRINT copyname /title = "ERROR: Two of the specified variables are copies. The variable names are:" /format = A8.
   ELSE IF(runnotes(i,1) = 8). 
      PRINT /title = "ERROR: All specified variables must be unique. No variables may be the same in W, M, and Y". 
   ELSE IF (runnotes(i,1) = 9). 
      PRINT /title = "NOTE: Both Monte Carlo Confidence Interval and Bias-Correction Bootstrap ". 
      PRINT /title = "       Confidence Interval were selected. Monte Carlo CI was calculated." /space = 0.
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
      PRINT /title = "ERROR: Johnson-Neyman procedure not available for models with a dichtomous moderator. ".
   ELSE IF (runnotes(i,1) = 18). 
      PRINT /title = "ERROR: You must specify a value to probe at for each moderator in the m list.".
      PRINT /title = "       Wmodval lists should be the same length as w list." /space = 0.
   ELSE IF (runnotes(i,1) = 19). 
      PRINT /title = "ERROR: Sample size is not large enough. Model is unidentified.".
   ELSE IF (runnotes(i,1) = 20). 
      PRINT /title = "ERROR: Contrasts not available for serial mediation with more than 3 pairs of mediators.".
   ELSE IF (runnotes(i,1) = 21). 
      PRINT /title = "ERROR: Models 2 and 3 require a variable name in the W argument.".
   ELSE IF (runnotes(i,1) = 22). 
      PRINT /title = "ERROR: Model 1 requires at least two variable names in the M argument.".
   ELSE IF (runnotes(i,1) = 23). 
      PRINT /title = "ERROR: Invalid model number. Please select Model 1, 2, or 3.".
   ELSE IF (runnotes(i,1) = 24). 
      PRINT /title = "NOTE: Centering command has no effect for Model 1.".
   END IF. 
END LOOP. 

DO IF (criterr = 0). 
DO IF (model = 1). 
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
   COMPUTE blabs = {"b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10"}. 

   COMPUTE savelab = make(1, 3*Mpairs+3+(cols+serind-mpairs+1)*(serial=1),0).
   DO IF (serial = 0). 
      LOOP i = 1 to Mpairs. 
         COMPUTE savelab(1,3*i) = mlab(1,i). 
         DO IF Mpairs = 1. 
         COMPUTE savelab(1,3*i-1) = {'b'}. 
         COMPUTE savelab(1,3*1-2) = {'a'}. 
         ELSE. 
         COMPUTE savelab(1,3*i-1) =blabs(1,i). 
         COMPUTE savelab(1,3*i-2) = alabs(1,i). 
         END IF. 
      END LOOP. 
      
      COMPUTE savelab(1, (ncol(savelab) - 2):(ncol(savelab))) = {"TotalInd", "Direct", "Total"}. 
    ELSE. 
      COMPUTE savelab(1, 1:3) = {'a1', 'b1', 'Ind1'}. 
      LOOP i = 2 to Mpairs. 
         COMPUTE savelab(1,(2*i):(2*i+1))= {blabs(1,i), mlab(1,i)}. 
      END LOOP. 
      DO IF (xmint = 1). 
         COMPUTE serlab = {'a2', 'd1', 'f1', 'a3', 'd2', 'f2', 'd3', 'f3', 'a4', 'd4', 'f4', 'd5', 'f5', 'd6', 'f6', 'a5', 'd7', 'f7', 'd8', 'f8', 'd9', 'f9', 'd10', 'f10'}. 
         COMPUTE savelab(1, (2+2*mpairs):(2*mpairs+mpairs**2)) = serlab(1, 1:(mpairs**2-1)). 
         COMPUTE savelab(1, (2*mpairs+mpairs**2+1):(2*mpairs+mpairs**2+serind)) = mlab(1,(nrow(indres) - serind):(nrow(indres)-1)). 
      ELSE IF (xmint = 0). 
         COMPUTE serlab = {'a2', 'd1', 'a3', 'd2', 'd3', 'a4', 'd4',  'd5', 'd6','a5', 'd7', 'd8',  'd9', 'd10'}.  
         COMPUTE savelab(1,(2+2*mpairs):(1+2*mpairs+cols)) = serlab(1, 1:cols). 
         COMPUTE savelab(1, (2+2*mpairs+cols):(1+2*mpairs+cols+serind)) = mlab(1,(nrow(indres) - serind):(nrow(indres)-1)). 
       END IF. 
   COMPUTE savelab(1, (ncol(savelab) - 2):(ncol(savelab))) = {"TotalInd", "Direct", "Total"}. 
   END IF. 
   DO IF (mc = 1 AND !save = 1). 
      COMPUTE savelab = {savelab, "TotalInd"}. 
      SAVE mcsave /outfile =* /names = savelab.
   ELSE IF (mc <> 1 AND !save = 1).  
      SAVE bootsave /outfile =* /names = savelab. 
   END IF. 
END IF. 

print conf /title = "Level of confidence for all confidence intervals in output:" /format = F10.2. 
DO IF (((model = 2) OR (model = 3)) AND (center = 1)). 
    COMPUTE centvars = wnames. 
    print centvars /title = "The following variables were mean centered prior to analysis:" /format = A8. 
END IF. 

END IF. 

end matrix. 
!ENDDEFINE. 
restore. 



COMMENT BOOKMARK;LINE_NUM=613;ID=3.
COMMENT BOOKMARK;LINE_NUM=776;ID=2.
COMMENT BOOKMARK;LINE_NUM=848;ID=6.
COMMENT BOOKMARK;LINE_NUM=1373;ID=4.
COMMENT BOOKMARK;LINE_NUM=1736;ID=5.
