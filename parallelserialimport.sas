PROC IMPORT OUT= WORK.parallelserial 
            DATAFILE= "C:\QRCLabData\Amanda\MEMORE\MEMOREgit\MEMORE\ParallelSerial.sav" 
            DBMS=SPSS REPLACE;

RUN;


data WORK.parallelserial;set WORK.parallelserial;dich=(M1S >=1);run;



PROC IMPORT OUT= WORK.Parents 
            DATAFILE= "C:\QRCLabData\Amanda\memoreandmplusresults\ParentsMLM2 Wide.sav" 
            DBMS=SPSS REPLACE;

RUN;
C:\QRCLabData\Amanda\memoreandmplusresults

