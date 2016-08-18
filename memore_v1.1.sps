/* MEMORE for SPSS Version 1.1*/.
/* Copyright 2016 */.
/* by Amanda Kay Montoya */.
/* www.afhayes.com*/.
/* Documentation available online at afhayes.com */.
/* or by email to montoya.29@osu.edu */.


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

DEFINE MEMORE (Y = !charend('/') /M = !charend('/') /Conf = !charend('/') !default(95) /mc = !charend('/') !default(0) 
   /samples = !charend('/') !default(5000) /normal = !charend('/') !default(0) /bc = !charend('/') !default(0) /decimals=!charend('/') !default(F10.4) /save = !charend('/') !default(0)
   /seed = !charend('/') !default(random) /contrast = !charend('/') !default(0) /serial = !charend('/') !default(0)). 
set mxloop = 100000000.
set seed = !seed. 
   
matrix. 
COMPUTE runnotes = MAKE(14,1,0). 
COMPUTE criterr = 0.  
GET data/ Variables = !M !Y / Names = namevec /missing = OMIT. 
GET mdat / variables = !M / names = mnames /missing = OMIT. 
GET ydat / variables = !Y / names = ynames /missing = OMIT. 
GET fulldat / Variables = !M !Y /missing = 999. 
COMPUTE missing = nrow(fulldat) - nrow(data). 

COMPUTE mc = (!mc=1). 
COMPUTE serial = (!serial = 1). 

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

DO IF (missing > 0). 
   COMPUTE runnotes(1,1) =1. 
END IF. 
DO IF (ncol(ydat) <> 2). 
   COMPUTE runnotes(2, 1) = 2.  
   COMPUTE criterr = 1. 
END IF. 
COMPUTE Mcount = ncol(mnames). 
DO IF (((Mcount = 0) OR (Mcount = 1) OR mod(mcount,2) > 0) OR (Mcount >20) AND (criterr <> 1)). 
   COMPUTE runnotes(6, 1) = 6. 
   COMPUTE criterr = 1. 
END IF. 

DO IF (criterr = 0). 
LOOP i = 1 TO Mcount. 
   DO IF ((mnames(1,i) = ynames(1,1)) OR (mnames(1,i) = ynames(1,2))). 
   COMPUTE runnotes(8,1) = 8. 
   COMPUTE criterr = 1. 
   END IF. 
END LOOP. 
END IF. 

DO IF ((serial = 1) AND (mc = 1)). 
   COMPUTE runnotes(12,1) = 12. 
   COMPUTE mc = 0. 
END IF. 

COMPUTE zero = MAKE (nrow(data),1,0). 
LOOP i = 1 TO (ncol(data)-1) BY 2. 
   COMPUTE diff = data(:,i) - data(:,i+1). 
   COMPUTE copy = csum(diff = zero). 
   DO IF (copy = nrow(data)). 
      COMPUTE copyname = {namevec(1,i), namevec(1,i+1)}. 
      BREAK. 
   END IF.  
END LOOP. 
DO IF (copy = nrow(data)). 
   COMPUTE runnotes(7,1) = 7. 
   COMPUTE criterr = 1. 
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
DO IF (mc = 1 AND trunc(!bc) = 1).  
   COMPUTE runnotes(9,1) = 9. 
   COMPUTE bc = 0. 
END IF. 

DO IF (!contrast = 1 AND (Mcount/2) =1). 
   COMPUTE runnotes(10,1) = 10. 
END IF. 

DO IF (serial = 1 AND ((Mcount >4)OR (Mcount <4))). 
   COMPUTE runnotes(13,1) = 13.
   COMPUTE criterr = 1. 
END IF. 

DO IF (criterr = 0). 

COMPUTE Mpairs = Mcount/2. 
COMPUTE mnamemat = reshape(mnames, Mpairs, 2).
COMPUTE N = nrow(data). 
COMPUTE transmat = {1,1/2;-1,1/2}. 
COMPUTE Tmat = make(ncol(data), ncol(data), 0). 
LOOP i = 1 TO (2*Mpairs+1) BY 2. 
   COMPUTE Tmat(i:(i+1), i:(i+1)) = transmat. 
END LOOP. 

COMPUTE dataT = data*Tmat.
COMPUTE aresmat = MAKE (Mpairs, 7, 0). 
COMPUTE ghostdes = make(N,1, 1).
COMPUTE alpha = (1-.01*Conf). 
COMPUTE temp = alpha/2. 
CDFINVT p = temp/ df = (N-1). 
COMPUTE tcrita = toutput. 
COMPUTE tcritc = toutput.

LOOP j = 1 TO Mpairs. 
   COMPUTE summean = csum(dataT(:,2*j))/N. 
   COMPUTE dataT(:,2*j) = (dataT(:,2*j) - summean).    
   COMPUTE apath = inv(t(ghostdes)*ghostdes)*t(ghostdes)*dataT(:,(2*j-1)).
   COMPUTE seMdiff = sqrt((csum((dataT(:,2*j-1)-(csum(dataT(:,2*j-1))/N))&**2)/(N-1))/N).
   COMPUTE tapath = apath/seMdiff.  
   COMPUTE dfapath = N-1. 
   COMPUTE papath = 2*(1-tcdf(abs(tapath), dfapath)). 
   COMPUTE lcia = apath-tcrita*seMdiff. 
   COMPUTE ucia = apath+tcrita*seMdiff. 
   COMPUTE aresmat(j,:) = {apath, seMdiff, tapath, dfapath, papath, lcia, ucia}. 
END LOOP. 

COMPUTE cpath = inv(t(ghostdes)*ghostdes)*t(ghostdes)*dataT(:,(2*Mpairs+1)). 
COMPUTE seYdiff = sqrt((csum((dataT(:,2*Mpairs+1)-(csum(dataT(:,2*Mpairs+1))/N))&**2)/(N-1))/N).
COMPUTE tcpath = cpath/seYdiff. 
COMPUTE dfcpath = N-1. 
COMPUTE pcpath = 2*(1-tcdf(abs(tcpath), N-1)). 
COMPUTE lcic = cpath-tcritc*seYdiff. 
COMPUTE ucic = cpath+tcritc*seYdiff. 
COMPUTE cresmat = {cpath, seYdiff, tcpath, dfcpath, pcpath, lcic, ucic}. 

DO IF (serial = 1). 
   COMPUTE serres = make(3, 7, 0). 
   COMPUTE serdes = {make(N,1,1), dataT(:, 1:(ncol(dataT)-4))}. 
   COMPUTE M2modbs = inv(t(serdes)*serdes)*t(serdes)*dataT(:,(ncol(dataT)-3)). 
   COMPUTE M2pred = serdes*m2modbs. 
   COMPUTE M2ssr = csum((dataT(:,ncol(dataT)-3)-M2pred)&**2). 
   COMPUTE M2sst = csum((dataT(:,(ncol(dataT)-3)) - csum(dataT(:,(ncol(dataT)-3)))/N)&**2). 
   COMPUTE M2Rsq = 1-M2ssr/m2sst. 
   COMPUTE M2r = sqrt(M2Rsq).    
   COMPUTE M2msr = M2ssr/(N - ncol(serdes)). 
   COMPUTE M2df1 = ncol(serdes) - 1. 
   COMPUTE M2df2 = (N - ncol(serdes)). 
   COMPUTE M2F = m2df2*m2rsq/(m2df1*(1-m2rsq)). 
   COMPUTE M2p = 1-FCDF(M2F, M2df1, m2df2). 
   COMPUTE sem2bmat = (m2msr*inv(t(serdes)*serdes)). 
   COMPUTE sem2b = (diag(sem2bmat))&**(1/2). 
   COMPUTE m2modsum = {M2r, m2rsq, m2msr, m2F, m2df1, m2df2, m2p}. 
   
   CDFINVT p = temp /df = M2df2. 
   COMPUTE sercritt = toutput. 
   COMPUTE serres(1:3,1) = M2modbs. 
   COMPUTE serres(1:3, 2) = sem2b. 
   COMPUTE serres(1:3, 3) = serres(1:3, 1) &/ serres(1:3, 2). 
   COMPUTE serres(1:3, 4) = make(3,1, m2df2). 
   COMPUTE serres(1:3, 5) = 2*(1-tcdf(abs(serres(1:3,3)), m2df2)). 
   COMPUTE serres(1:3,6) = serres(1:3, 1) - sercritt*serres(1:3,2). 
   COMPUTE serres(1:3, 7) = serres(1:3, 1) + sercritt*serres(1:3,2). 
   COMPUTE aresmat(2,:) = serres(1,:). 
END IF. 

COMPUTE bcpdes = {make(N,1,1), dataT(:,1:(ncol(dataT)-2))}.
COMPUTE bcpvec = inv(t(bcpdes)*bcpdes)*t(bcpdes)*dataT(:,(ncol(dataT)-1)).
COMPUTE ypred = bcpdes*bcpvec. 
COMPUTE ssr = csum((dataT(:,(ncol(dataT)-1)) - ypred)&**2). 
COMPUTE sst = csum((dataT(:,(ncol(dataT)-1)) - csum(dataT(:,(ncol(dataT)-1)))/N)&**2).
COMPUTE msr = ssr/(N-ncol(bcpdes)). 
COMPUTE rsqfull = 1- ssr/sst. 
COMPUTE rfull = sqrt(rsqfull). 
COMPUTE df1 = (ncol(bcpdes)-1). 
COMPUTE df2 = (N - ncol(bcpdes)). 
COMPUTE Ffull = df2*Rsqfull/((df1)*(1-rsqfull)). 
COMPUTE pfull =1- FCDF(Ffull, df1, df2). 
COMPUTE sebcpmat = (msr*inv(t(bcpdes)*bcpdes)). 
COMPUTE sebcp = (diag(sebcpmat))&**(1/2). 

COMPUTE bresmat = MAKE(Mpairs, 7, 0). 
COMPUTE dresmat = MAKE(Mpairs, 7, 0). 
COMPUTE indres = MAKE(Mpairs+1+(serial=1),1,0). 
DO IF (!normal = 1). 
   COMPUTE normres = MAKE(Mpairs, 4, 0). 
END IF. 

CDFINVT p = temp/ df = df2. 
COMPUTE tcritb = toutput. 
COMPUTE tcritcp = toutput. 
COMPUTE tcritd = toutput. 

COMPUTE cppath = bcpvec(1,1). 
COMPUTE secppath = sebcp(1,1). 
COMPUTE tcppath = cppath/secppath.
COMPUTE pcppath = 2*(1-tcdf(abs(tcppath), df2)).
COMPUTE lcicp = cppath-tcritcp*secppath. 
COMPUTE ucicp = cppath+tcritcp*secppath. 
COMPUTE cpresmat = {cppath, secppath, tcppath, df2, pcppath, lcicp, ucicp}. 

COMPUTE LCII = rnd((1-.01*Conf)/2*samples). 
COMPUTE UCII = trunc((1-((1-.01*Conf)/2))*samples)+1. 
DO IF (LCII  < 1 OR UCII > samples). 
   COMPUTE runnotes(4, 1) = 4.  
   COMPUTE criterr = 1. 
   COMPUTE LCII = 1. 
   COMPUTE UCII = samples. 
END IF. 

/*MC Setup. 
DO IF (mc = 1). 
COMPUTE mcsamps = samples. 
COMPUTE randsamp = sqrt(-2*ln(uniform(mcsamps,Mpairs)))&*cos((2*3.14159265358979)*uniform(mcsamps,Mpairs)).
COMPUTE MCres = MAKE (Mpairs+1, 4, 0). 
COMPUTE MCcorr = MAKE(Mpairs, Mpairs, 1). 
LOOP i = 1 TO Mpairs. 
   COMPUTE MCcorr(i,i) = sebcpmat(2*i,2*i).
   DO IF ((Mpairs > 1) AND (i <> Mpairs)). 
   LOOP j = (i+1) TO Mpairs. 
   COMPUTE MCcorr(i,j) = sebcpmat(2*i,(2*j)). 
   COMPUTE MCcorr(j,i) = sebcpmat(2*i,2*j). 
   END LOOP. 
   END IF. 
END LOOP.
COMPUTE rndnb = randsamp*chol(MCcorr).  
COMPUTE rndna = sqrt(-2*ln(uniform(mcsamps,Mpairs)))&*cos((2*3.14159265358979)*uniform(mcsamps,Mpairs)).
COMPUTE mcsave = make(samples, 3*Mpairs+1, 0). 
COMPUTE mcsave2 = make(samples, Mpairs, 0). 
END IF. 

LOOP i = 1 TO Mpairs. 
   COMPUTE bpath = bcpvec(2*i,1). 
   COMPUTE sebpath = sebcp(2*i,1). 
   COMPUTE tbpath = bpath/sebpath. 
   COMPUTE pbpath = 2*(1-tcdf(abs(tbpath), df2)).
   COMPUTE lcib = bpath-tcritb*sebpath. 
   COMPUTE ucib = bpath+tcritb*sebpath. 
   COMPUTE bresmat(i, :) = {bpath, sebpath, tbpath, df2, pbpath, lcib, ucib}. 

   COMPUTE dpath = bcpvec((2*i+1), 1). 
   COMPUTE sedpath = sebcp((2*i+1),1). 
   COMPUTE tdpath = dpath/sedpath. 
   COMPUTE pdpath = 2*(1-tcdf(abs(tdpath),df2)). 
   COMPUTE lcid = dpath-tcritd*sedpath. 
   COMPUTE ucid = dpath+tcritd*sedpath. 
   COMPUTE dresmat(i, :) = {dpath, sedpath, tdpath, df2, pdpath, lcid, ucid}. 

   /*Normal Tests. 
   COMPUTE indirect = aresmat(i,1)*bresmat(i,1). 
   COMPUTE indres(i,1) = indirect.  
   DO IF (!normal = 1). 
      COMPUTE sobseab = sqrt((aresmat(i,1)**2)*(bresmat(i,2)**2)+(bresmat(i,1)**2)*(aresmat(i,2)**2)). 
      COMPUTE sobelZ = indirect/sobseab.
      COMPUTE sobelp = 2*cdfnorm((-1)*abs(sobelZ)). 
      COMPUTE normres(i,:) = {indirect, sobseab, sobelZ, sobelp}. 
   END IF. 

   /*Monte Carlo Confidence Interval*/
   DO IF (mc = 1). 
   COMPUTE asamp = rndna(:,i)*aresmat(i,2)+aresmat(i,1). 
   COMPUTE bsamp = rndnb(:,i)+bresmat(i,1). 
   COMPUTE absamp = asamp&*bsamp. 
   COMPUTE mcgrad = grade(absamp). 
   COMPUTE mcsort = absamp. 
   COMPUTE mcsort(mcgrad) = absamp.
   COMPUTE mcsave(:,(3*i-2):(3*i)) = {asamp, bsamp,absamp}. 
   COMPUTE mcsave2(:,i) = absamp. 
   COMPUTE MCLLCI = mcsort(LCII). 
   COMPUTE MCULCI = mcsort(UCII). 
   COMPUTE seMC = sqrt(csum((mcsort(:,1)-(csum(mcsort(:,1))/mcsamps))&**2)/(mcsamps-1)).
   COMPUTE MCres(i,:) = {indirect, seMC, MCLLCI, MCULCI}.  
   END IF. 

END LOOP. 
 
/*serial mediation indirect paths*/
DO IF (serial = 1).
   COMPUTE indres(mpairs+1, 1) = aresmat(1,1)*serres(2,1)*bresmat(2,1).
   DO IF (!normal = 1)). 
   COMPUTE indirect = aresmat(1,1)*serres(2,1)*bresmat(2,1). 
   COMPUTE sobseab = sqrt((aresmat(1,1)**2)*(serres(2,1)**2)*(bresmat(2,2)**2)+(aresmat(1,1)**2)*(bresmat(2,1)**2)*(serres(2,2)**2)+(serres(2,1)**2)*(bresmat(2,1)**2)*(aresmat(2,1)**2)). 
   COMPUTE sobelZ = indirect/sobseab. 
   COMPUTE sobelp = 2*cdfnorm((-1)*abs(sobelZ)). 
   COMPUTE serind = {indirect, sobseab, sobelZ, sobelp}. 
   COMPUTE normres = {normres; serind}.

   END IF. 
END IF. 

/*Total monte carlo*/
DO IF (mc = 1). 
COMPUTE mcsave(:,3*Mpairs+1) = rsum(mcsave2(:,:)). 
COMPUTE mcsort = mcsave(:,3*Mpairs+1). 
COMPUTE mcgrad = grade(mcsort). 
COMPUTE mcsort(mcgrad) = mcsort. 
COMPUTE MCLLCI = mcsort(LCII). 
COMPUTE MCULCI = mcsort(UCII). 
COMPUTE seMC = sqrt(csum((mcsort(:,1)-(csum(mcsort(:,1))/mcsamps))&**2)/(mcsamps-1)).
COMPUTE MCres(Mpairs+1, :) = {csum(indres), seMC, MCLLCI, MCULCI}. 

DO IF ((!contrast = 1) AND (Mpairs >1)). 
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
   /*Bootstrapping*/
DO IF (mc <> 1). 
   COMPUTE detcheck = make(((Mpairs-1)*(serial=1) +1), 1, -999). 
   COMPUTE Bootsamp = make(samples, Mpairs+1+(serial=1), 0). 
   COMPUTE Bootsave = make(samples, 3*Mpairs+3+2*(serial=1),0). 
   COMPUTE indtemp = make(samples,Mpairs+(serial=1),0).
   LOOP i = 1 TO samples.
      LOOP k = 1 TO 10000.  
      COMPUTE sortvar = trunc(uniform(N,1)*N)+1.
      COMPUTE bootdat = dataT(sortvar(:,1),:). 
      LOOP j = 1 TO Mpairs. 
      COMPUTE summean = csum(bootdat(:,2*j))/N.
      COMPUTE bootdat(:,2*j) = (bootdat(:,2*j) - summean).
      END LOOP.
      COMPUTE bootdes = {make(N,1,1), bootdat(:,1:(ncol(bootdat)-2))}. 
      COMPUTE detcheck(1,1) = (det(t(bootdes)*bootdes)=0).
      DO IF (serial = 1).
      LOOP j = 2 TO Mpairs. 
         COMPUTE bootadat = bootdes(:,1:(2*j-1)). 
         COMPUTE detcheck(j,1) = (det(t(bootadat)*bootadat)=0). 
      END LOOP. 
      END IF. 
      COMPUTE badboot = badboot+(k = 2). 
      END LOOP IF (csum(detcheck(:,1)) = 0). 
      
      COMPUTE bootbeta = inv(t(bootdes)*bootdes)*t(bootdes)*bootdat(:,(ncol(bootdat)-1)). 
      LOOP j = 1 TO Mpairs. 
      COMPUTE boota = inv(t(ghostdes)*ghostdes)*t(ghostdes)*bootdat(:,(2*j-1)).
      COMPUTE bootb = bootbeta(2*j,1).

      DO IF ((serial = 1)AND(j>1)).
         COMPUTE bootadat = bootdes(:,1:(2*j-1)). 
         COMPUTE bootserb = inv(t(bootadat)*bootadat)*t(bootadat)*bootdes(:,2*j). 
         COMPUTE boota = bootserb(1,1). 
         COMPUTE bootd = bootserb(2,1). 
         COMPUTE bootsamp(i,j+1) = bootsave(i,1)*bootd*bootb. 
         COMPUTE bootsave(i,3*j+2) = bootsave(i,1)*bootd*bootb.
         COMPUTE bootsave(i,3*j+1) = bootd. 
         COMPUTE indtemp (i,j+1) = bootsave(i,1)*bootd*bootb.
      END IF. 

      COMPUTE bootsamp(i,j) = bootb*boota. 
      COMPUTE bootsave(i,(3*j-2):(3*j)) = {boota, bootb, boota*bootb}. 
      COMPUTE indtemp(i,j) = (boota*bootb). 
      END LOOP. 
      COMPUTE bootsave(i,ncol(bootsave)-1) = bootbeta(1,1). 
      COMPUTE bootsave(i,ncol(bootsave)) = rsum(indtemp(i,:))+bootbeta(1,1). 
   END LOOP. 

   DO IF (badboot >0). 
      COMPUTE runnotes(14,1) = 14. 
   END IF. 
      


   DO IF ((!contrast = 1) AND (Mpairs >1)). 
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
   COMPUTE seboots = MAKE(Mpairs+1+(serial=1), 1, 0). 
   COMPUTE bccires = MAKE(4,Mpairs+1+(serial=1), 0). 
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

   DO IF  (!contrast = 1) AND (Mpairs >1). 
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

END IF. 

print/title = "******************* MEMORE Procedure for SPSS Version 1.1 *********************".
print/title = "                           Written by Amanda Montoya       ".
print/title = "                    Documentation available at afhayes.com ".
print /title = "********************************************************************************".
DO IF (criterr = 0). 
DO IF Mpairs = 1. 
COMPUTE varrlabs = {'Y = ', 'M = '}. 
ELSE. 
COMPUTE varrlabs = {'Y = ', 'M1 = ', 'M2 = ', 'M3 = ', 'M4 = ', 'M5 = ', 'M6 = ', 'M7 = ', 'M8 = ', 'M9 = ', 'M10 = '}.
END IF. 
print {ynames; mnamemat} /title = "Variables: " /rnames = varrlabs  /format = a8. 
COMPUTE compname = MAKE((2*Mpairs+1),7, 0). 
COMPUTE compname(1,:) = {' ', ynames(1,1), ' - ', ynames(1,2), ' ', ' ', ' '}. 
LOOP j = 1 TO Mpairs. 
   COMPUTE compname((1+j),:) ={' ', mnamemat(j,1), ' - ', mnamemat(j,2), ' ', ' ', ' '}. 
   COMPUTE compname((Mpairs+1+j),:) = { '(', mnamemat(j,1), ' + ', mnamemat(j,2), ')', '/2', 'Centered'}. 
END LOOP. 
Compute temp1 = {'M1diff = ','M2diff = ','M3diff = ', 'M4diff = ', 'M5diff = ','M6diff = ','M7diff = ','M8diff = ', 'M9diff = ', 'M10diff = '}. 
COMPUTE temp2 = {'M1avg  = ','M2avg  = ','M3avg  = ', 'M4avg  = ', 'M5avg  = ','M6avg  = ','M7avg  = ','M8avg  = ', 'M9avg  = ', 'M10avg  = '}.
DO IF (Mpairs = 1). 
COMPUTE temprnam = {'Ydiff = ', 'Mdiff = ', 'Mavg = '}. 
ELSE. 
compute temprnam = {'Ydiff = ', temp1(1,1:Mpairs), temp2(1,1:Mpairs)}. 
END IF. 
print compname /title = "Computed Variables:" /rnames = temprnam /format = a8. 
print N /title = "Sample Size:". 
do if (!quote(!seed) <> "random"). 
print !seed /title = "Seed:". 
end if. 
print {"Ydiff =" , ynames(1,1), ' - ', ynames(1,2)} /title = "********************************************************************************" /rlabels = "Outcome:" /format = A8.
COMPUTE collab = {"Effect", "SE", "t", "df", "p", "LLCI", "ULCI"}. 
print cresmat /title = "Model" /rnames = {"'X'"} /cnames = collab /format = !decimals. 
COMPUTE alabs = {"a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10"}.
LOOP j = 1 to Mpairs. 
DO IF (serial = 1 AND j >1). 
   BREAK. 
END IF. 
DO IF (Mpairs = 1). 
print {"Mdiff = " , mnamemat(j,1), ' - ', mnamemat(j,2)} /title = "********************************************************************************" /rlabels = "Outcome:" /format = A8.
ELSE. 
print {temp1(1,j) , mnamemat(j,1), ' - ', mnamemat(j,2)} /title = "********************************************************************************" /rlabels = "Outcome:" /format = A8.
END IF. 
print aresmat(j,:) /title = "Model" /rnames = {"'X'"} /cnames = collab /format = !decimals. 
END LOOP. 
DO IF (serial = 1). 
print {"M2diff =" , mnamemat(2,1), ' - ', mnamemat(2,2)} /title = "********************************************************************************" /rlabels = "Outcome:" /format = A8.
print m2modsum /title = "Model Summary" /clabels = "R", "R-sq", "MSE", "F" , "df1" , "df2", "p" /format = !decimals. 
COMPUTE m2labs = {"'X'", "M1diff", "M1avg"}. 
print serres /title "Model" /rnames = m2labs /clabels = "coeff" , "SE", "t", "df", "p", "LLCI", "ULCI" /format = !decimals.
END IF. 
print {"Ydiff =" , ynames(1,1), ' - ', ynames(1,2)} /title = "********************************************************************************" /rlabels = "Outcome:" /format = A8.
COMPUTE modsumr = {Rfull, Rsqfull, MSR, Ffull, df1, df2, pfull}. 
print modsumr /title = "Model Summary" /clabels = "R", "R-sq", "MSE", "F" , "df1" , "df2", "p" /format = !decimals. 

COMPUTE modres = {cpresmat;bresmat;dresmat}. 
COMPUTE bdlabs = {"M1diff", "M2diff", "M3diff", "M4diff", "M5diff", "M6diff", "M7diff", "M8diff", "M9diff", "M10diff"}.
COMPUTE bslabs = {"M1avg", "M2avg", "M3avg", "M4avg", "M5avg", "M6avg", "M7avg", "M8avg", "M9avg", "M10avg"}.
DO IF (Mpairs = 1). 
COMPUTE modlabs = {"'X'", "Mdiff", "Mavg"}. 
ELSE. 
COMPUTE modlabs = {"'X'", bdlabs(1, 1:Mpairs), bslabs(1, 1:Mpairs)}. 
END IF. 
print modres /title "Model" /rnames = modlabs /clabels = "coeff" , "SE", "t", "df", "p", "LLCI", "ULCI" /format = !decimals.
print /title = "********************** TOTAL, DIRECT, AND INDIRECT EFFECTS **********************" .
COMPUTE collab = {"Effect", "SE", "t", "df", "p", "LLCI", "ULCI"}. 
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
COMPUTE m2lab = {mlab(1,1:(nrow(indres)-1)), 'Total'}. 

DO IF (Mpairs = 1). 
print indres(1,:) /title = "Indirect Effect of X on Y through M" /rnames = m2lab / cnames = indlabs /format = !decimals. 
ELSE. 
print indres /title = "Indirect Effect of X on Y through M" /rnames = m2lab / cnames = indlabs /format = !decimals. 
END IF. 

DO IF (!normal = 1). 
print normres /title = "Normal Theory Tests for Indirect Effect" /rnames = mlab /clabels = "Effect", "SE", "Z", "p" /format = !decimals. 
END IF.

COMPUTE indkey = make((ncol(m2lab)-1), 7, 0). 
LOOP i = 1 to (ncol(m2lab)-1).
   COMPUTE indkey(i,:) = {"X" , "->", bdlabs(1,i), "->", "Ydiff", " ", " "}. 
END LOOP. 
DO IF (serial = 1). 
   COMPUTE indkey(3,:) = {"X" , "->", bdlabs(1,1), "->", bdlabs(1,2), "->", "Ydiff"}. 
END IF. 
print indkey /title = "Indirect Key" /rnames =mlab /format = A8. 

DO IF  (!contrast = 1) AND (Mpairs >1).
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
END IF. 

print /title = "************************* ANALYSIS NOTES AND WARNINGS **************************". 
LOOP i = 1 to nrow(runnotes).
   DO IF (runnotes(i,1) = 1). 
      PRINT missing /title = "NOTE: Some cases were deleted due to missing data. The number of cases was:". 
  ELSE IF (runnotes(i,1) = 2 ). 
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
      PRINT /title = "ERROR: All specified variables must be unique. No variables may be the same in M and Y". 
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
      PRINT /title = "ERROR: The serial mediation model must have 2 pairs of mediators.".
   ELSE IF (runnotes(i,1) = 14). 
      PRINT badboot /title = "NOTE: Some bootstrap samples had to be replaced.  The number of such replacements was:".
   END IF. 
END LOOP. 

DO IF (criterr = 0). 
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
print conf /title = "Level of confidence for all confidence intervals in output:" /format = F10.2. 
COMPUTE centvars = MAKE(Mpairs, 6,0). 
LOOP j = 1 TO Mpairs. 
   COMPUTE centvars(j,:) = { '(', mnamemat(j,1), ' + ', mnamemat(j,2), ')', '/2' }. 
END LOOP. 
print centvars /title = "The following variables were mean centered prior to analysis:" /format = A8. 

COMPUTE blabs = {"b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10"}. 

COMPUTE savelab = make(1, 3*Mpairs,0).
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

DO IF (mc = 1 AND !save = 1). 
   COMPUTE savelab = {savelab, "TotalInd"}. 
   SAVE mcsave /outfile =* /names = savelab.
ELSE IF (mc <> 1 AND !save = 1). 
   DO IF (serial = 1 and mpairs > 1). 
   COMPUTE savelab = {savelab, "a3", "Ind3"}.
   END IF. 
   COMPUTE savelab = {savelab, "TotalInd", "cprime", "c"}.
   SAVE bootsave /outfile =* /names = savelab. 
END IF. 

END IF. 

end matrix. 
!ENDDEFINE. 
restore. 
COMMENT BOOKMARK;LINE_NUM=346;ID=1.
COMMENT BOOKMARK;LINE_NUM=423;ID=2.
