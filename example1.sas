
%include 'gformula3.sas ';
%include 'gformula_rct_270918.sas';

%macro create_sample(eof=0 ) ;
%let condition = ;
%let condition = diab or cens  or dead  ;
%if &eof = 1   %then %let condition = cens ;



**SAMPLE Data;
    data sample(drop = i j ahbp1-ahbp8 aact1-aact8 acont1 - acont8);

        call streaminit(5027);

        do i=1 to 10000;
            baseage = int( 35 + 25*rand('uniform'));

            array ahbp(8);
            array aact(8);
            array acont(8);

            do j=1 to 8;
                ahbp(j) = (0.2>rand('uniform'));
                acont(j) = rand('normal');
                if j > 1 & ahbp(j-1) = 1 then ahbp(j) = 1 ;

        aact(j)=(0.7>rand('uniform'));
        if aact(j)=1 then do;
                   aact(j) = int(exp(3.5+0.4*(rand('normal'))));
        end;
                end;

            do j=3 to 8  until ( &condition   ) ;
                id=i;
                time=j-3;
                

                hbp     = ahbp(j);
                hbp_l1  = ahbp(j-1);
                hbp_l2  = ahbp(j-2);
                

                act     = aact(j);
                act_l1  = aact(j-1);
                act_l2  = aact(j-2);


                cont = acont(j);
                cont_l1 = acont(j-1);
                cont_l2 = acont(j-2);
                

              diab = ( (j/500) >rand('uniform'));
                cens  = (0.05>rand('uniform'));
                dead      = (0.05>rand('uniform'));
               

                %if &eof = 1 %then if j < 8 then diab = . ;;
                if j = 8 then conteof = rand('normal') ;
                output;

                end;
           end;




    run;


data sample ;
set sample; 
call streaminit(1234);
if cens=1 then do;
       diab= .;
       dead= .;
end;
else do;
   if dead=1 then diab= .;
end;
treat = rand('bernouli',0.6);
run;

data sample;
set sample;
rename time=visit;
run;

proc means data=sample;
title 'Means of SAMPLE data';
run;
%mend ;


%let interva1 = intno=1 ,  nintvar=1,
    intlabel='All subjects exercise at least 30 minutes per day in all intervals',
    intvar1 = act, inttype1 = 2, intmin1=30, intpr1=1, inttimes1 = 0 1 2 3 4 5 ;


%let intervb1 = intno=1 ,  nintvar=1,
    intlabel='All subjects exercise at least 45 minutes per day in all intervals',
    intvar1 = act, inttype1 = 2, intmin1=45, intpr1=1, inttimes1 = 0 1 2 3 4 5 ;

%create_sample(eof = 0 ); 

proc datasets library = work nolist ;
save sample ;
quit;
ods graphics off ;
 
 
%gformula_rct(
data= sample,
id=id,
time=visit,
/* time=time, */
timepoints = 6,
outc= diab,
outctype=binsurv,
comprisk =  dead  ,
arm = treat ,

timeptype= concat, 
timeknots = 1 2 3 4 5,


seed= 9458, 

hazardratio = 1,
bootstrap_hazard = 1,
forrct = 1 ,
/* intcomp =  , */
nsimul = , /* use same number in each simulated data set for each arm */
nsamples = 100,
survdata = mytest,
intervname = myintervt,
resultsdata = myresultst,
 rungraphs = 0,

/* for first arm */
ncova=2,
numinta=1 ,
runnca = 1 ,
refinta = 1,
fixedcova = baseage,
cova1  = hbp,    cova1otype  = 2, cova1ptype = lag1bin,
cova2  = act,    cova2otype  = 4, cova2ptype = cumavg,

/* for second arm */
ncovb = 2 ,
numintb = 1,
runncb = 1,
refintb = 1,
fixedcovb= baseage ,
covb1  = hbp,    covb1otype  = 2, covb1ptype = lag1bin,
covb2  = act,    covb2otype  = 4, covb2ptype = cumavg,
 
usespline = 0 ,
graphfilea=gfilea.pdf,
graphfileb=gfileb.pdf

);


ods graphics on ; 
