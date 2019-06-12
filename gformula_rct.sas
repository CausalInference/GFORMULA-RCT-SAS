




%macro gformula_rct(
    data=,           
    id=,      
    time=,       
    timeptype = ,
    timeknots = ,
    timeinc = 1,
    timefuncgen = ,   
    timepoints=,    
    interval=1,      

    outc=,          
    outctype=binsurv,  
    outcinteract=,   
    outcwherem = (1=1) ,  
    outcwherenosim=(1=0),  
    outcnosimelsemacro =,  
    comprisk=,        
    compriskinteract=,   
    compriskwherem = (1=1) , 
    compriskwherenosim=(1=0),   
    comprisknosimelsemacro =, 
  
    arm = , /* variable for separating the data into two arms, need to have values {0,1} for  {a,b},   */
    numinta= 0 ,
    refinta = 0 ,
    runnca = 1,
    fixedcova=,      
    ncova=,          
    usebetadataa = 0 ,
    betadataa = ,
    cova1=,cova1otype=1,cova1ptype=,cova1mtype= all ,cova1cumint= ,cova1skip=-1,cova1inc=0,cova1knots=,cova1interact=,cova1wherem=(1=1),cova1wherenosim=(1=0),
    cova1nosimelsemacro=, cova1class=, cova1classelse=, cova1addvars=,cova1genmacro=,cova1modusermacro =,cova1moddatausermacro =,cova1setblvar =,
    cova1simusermacro = ,cova1barray = ,cova1sarray =,cova1randomvisitp=,cova1visitpmaxgap=9e10,cova1visitpwherem=(1=1),cova1visitpcount = ,

/** for second call **/
 
   numintb= 0 ,
   refintb = 0 ,
   runncb = 1,
   fixedcovb=, 
   ncovb= ,
   usebetadatab = 0 ,
   betadatab = ,
   covb1=,covb1otype=1,covb1ptype=,covb1mtype= all ,covb1cumint= ,covb1skip=-1,covb1inc=0,covb1knots=,covb1interact=,covb1wherem=(1=1),covb1wherenosim=(1=0),
   covb1nosimelsemacro=, covb1class=, covb1classelse=, covb1addvars=,covb1genmacro=,covb1modusermacro =,covb1moddatausermacro =,covb1setblvar =,
   covb1simusermacro = ,covb1barray = ,covb1sarray =,covb1randomvisitp=,covb1visitpmaxgap=9e10,covb1visitpwherem=(1=1),covb1visitpcount = ,

    /*other override options used for each gformula call */

    wherevars=,        /* list of variables referenced in any of the cov#wherem conditions, change JGY*/ 
    keepsimuldata=,    /*list of variables not created by a given ptype or otype that will be needed in simulated data set new change JGY*/
    equalitiessimuldata=,/*user defined macro that equates pre-baseline simulated vars to observed new change JGY*/
    eventaddvars=, /*list of variables to be added to event predictor list new change JGY*/
     
    compriskaddvars=,
    
    censladdvars=,
    

 
    
    simuldata=,     /* data set to store simulated data */
    resultsdata=,   /* data set to store results */
    survdata=,       /*data set to store cumulative survival probabilities at each time point under all interventions*/
    outputs = yes,  /* whether to print regression results */
    print_stats = 1 ,
    check_cov_models = 0, /* create data set for difference of mean of observed covs and mean of simulated covs under natural 
                              course */
    print_cov_means = 0,  /* print out tables of comparison of observed and simulated variables */
    covmeandata =  ,
    save_raw_covmean = 0,
    observed_surv = ,
    intervname = ,
    
   
    seed= 7834,          /* random numbers seed */
    nsamples=50,    /* number of bootstrap samples (default is 50) */
    nsimul=,        /* size (# subjects) of simulated sample, default is sample size */
    nparam=,        /* size (# subjects) of parameter sample, default is sample size */
    hazardratio=0 , /* calculate the hazard ratio for two interventions. This will increase the run time for the macro*/
    forrct = 0 ,
    intcomp =    ,  /* needs to be a list with two numbers for comparing the hazard ratio when hazardratio = 1 */
    bootstrap_hazard = 0, /* when running bootstrap samples also include calculation of hazard ratios */
    hazardname =  , /* name of data set to hold hazard ratio when runnig bootstraps in parts, data will be saved in savelib. */
   
    sample_start = 0 ,  /* first sample to use in bootstraps (can be 0 for original data ) */
    sample_end = -1,     /* last sample to use in bootstraps (should be at most equal to nsamples) */
    savelib = work,  /* location for saving intermediate results for chunks */


    rungraphs = 0 ,
    title1a=,
    title2a=,
    title3a=,
    titledataa= ,
    graphfilea=gfilea.pdf ,
    tsizea=1 ,
    title1b=,
    title2b=,
    title3b=,
    titledatab= ,
    graphfileb=gfileb.pdf ,
    tsizeb=1 ,
    runnc = 1 ,
    weight = ,
    printlogstats = 1,
    usespline = 1 ,
    checkaddvars = 1,
    minimalistic = no, /* only keep results for outcome variables */
    testing = no /* keep each simulated data set for each intervention. will be named savelib.simulated&intno */
    )  /parmbuff;

 

    %let suffixlist=  otype  ptype  mtype   cumint  skip inc knots interact wherem wherenosim
     nosimelsemacro class classelse addvars genmacro modusermacro  moddatausermacro  setblvar simusermacro barray  
     sarray  randomvisitp visitpmaxgap visitpwherem  visitpcount ;


    %do i = 1 %to &ncova ;
        %let covlist = cova&i ;        
            %do j = 1 %to 25 ; 
                %let word = %scan(&suffixlist,&j) ;                
                %let covlist = &covlist cova&i.&word  ;
               
            %end;             
            %local &covlist ;
            %let cova&i.otype=1  ;
            %let cova&i.ptype=  ;
            %let cova&i.cumint=   ;    
            %let cova&i.mtype= all  ;
            %let cova&i.knots= ;     
            %let cova&i.skip=-1 ;     
            %let cova&i.inc=0 ;    
            %let cova&i.interact= ;  
            %let cova&i.wherem= (1=1) ;                   
            %let cova&i.wherenosim=(1=0) ; 
            %let cova&i.nosimelsemacro = ; 
            %let cova&i.class= ;    
            %let cova&i.classelse= ;     
            %let cova&i.addvars = ;
            %let cova&i.genmacro = ; 
            %let cova&i.modusermacro = ; 
            %let cova&i.moddatausermacro = ;  
            %let cova&i.setblvar = ;  
            %let cova&i.simusermacro = ;  
            %let cova&i.barray = ;   
            %let cova&i.sarray = ; 
            %let cova&i.randomvisitp = ;  
            %let cova&i.visitpmaxgap=9e10 ;     
            %let cova&i.visitpwherem=(1=1) ; 
            %let cova&i.visitpcount=  ;  
      %end;
 
 
     
    %do i = 1 %to &ncovb ;
        %let covlist = covb&i ;        
            %do j = 1 %to 25 ; 
                %let word = %scan(&suffixlist,&j) ;                
                %let covlist = &covlist covb&i.&word  ;               
            %end;             
            %local &covlist ;
 
            %let covb&i.otype=1  ;
            %let covb&i.ptype=  ;
            %let covb&i.cumint=   ;    
            %let covb&i.mtype= all  ;
            %let covb&i.knots= ;     
            %let covb&i.skip=-1 ;     
            %let covb&i.inc=0 ;    
            %let covb&i.interact= ;  
            %let covb&i.wherem= (1=1) ;                   
            %let covb&i.wherenosim=(1=0) ; 
            %let covb&i.nosimelsemacro = ; 
            %let covb&i.class= ;    
            %let covb&i.classelse= ;     
            %let covb&i.addvars = ;
            %let covb&i.genmacro = ; 
            %let covb&i.modusermacro = ; 
            %let covb&i.moddatausermacro = ;  
            %let covb&i.setblvar = ;  
            %let covb&i.simusermacro = ;  
            %let covb&i.barray = ;   
            %let covb&i.sarray = ; 
            %let covb&i.randomvisitp = ;  
            %let covb&i.visitpmaxgap=9e10 ;     
            %let covb&i.visitpwherem=(1=1) ; 
            %let covb&i.visitpcount=  ;  
    %end;
    
 

 %let mycount = %sysfunc(countw(&syspbuff,','));
  
 %do i = 1 %to &mycount ;  
    %let mydef = %scan(&syspbuff,&i,%str((),));   
    %let myvar = %scan(&mydef,1,%str(=));
    %let myval = %scan(&mydef,2,%str(=));
    %let &myvar = &myval ;
 %end;


  %if &hazardratio = 0 %then %do;
         %let forrct = 0 ;
         %let bootstrap_hazard = 0 ;
  %end; 

  %put &hazardratio &forrct &bootstrap_hazard &nsamples ;
  %put ;
  %if &nsamples = 0 %then %let bootstrap_hazard = 0 ;


 %if &ncova > 0 and &ncovb > 0 %then %do;
    %if %bquote(&survdata) = %then %let survdata = survholder ;
 %end;
    
 %if &ncova > 0 %then %modelA ;
 %if &ncovb > 0 %then %modelB;

 %if &hazardratio=1 and &forrct = 1 %then %hazardrct ;

 

 %let genresults = 0;
 %if &numinta = 0 and &numintb = 0 %then %let genresults = 1;
 %if &numinta = 0 and &numintb > 0 and &runncb = 1 %then %let genresults = 1 ;
 %if &numinta > 0 and &runnca = 1 and &numintb = 0  %then %let genresults = 1;
 %if &numinta > 0 and &numintb > 0 %then %let genresults = 1;

 %if &genresults = 1 %then %results2;
 %if &genresults = 0 %then %do ;
    %put NO COMPARISON DONE DUE TO NO MATCHING INTERVENTIONS ;
    %put( numinta = &numinta , runnca = &runnca , startint = %eval(1-&runnca) ) ( numintb = &numintb , runncb = &runncb , startint = %eval(1 - &runncb)) ;
%end; 

 /*****/
     
%mend ;
%macro modelA ;

 data dataa ;
 set &data(where = (&arm = 0));
 run;


     
    %let simuldataa = ;
    %let resultsdataa= _resultsA_ ;
    %let survdataa = ;
    %let covmeandataa = ;
    %let observed_surva = ;
    %let intervnamea = ;
    

    %if %bquote(&simuldata)^= %then %let simuldataa = &simuldata.a ;          
     
    %if %bquote(&survdata)^=  %then %let survdataa = &survdata.a  ;
    %if %bquote(&covmeandata)^=  %then %let covmeandataa = &covmeandata.a ;      
    %if %bquote(&observed_surv)^=  %then %let observed_surva = &observed_surv.a; 
    %if %bquote(&intervname)^=  %then %let intervnamea = &intervname.a ;
    
    
    %if &numinta > 0 %then %do;
        %do i = 1 %to &numinta ;
            %let interv&i = &&interva&i ;
            %put interv&i for modelA ::: &&interv&i ;
        %end;
    %end;
 
 %gformula(
    data= dataa,           
    id= &id,      
    time= &time,       
    timeptype = &timeptype ,
    timeknots = &timeknots,
    timeinc = &timeinc,
    timefuncgen = &timefuncgen,   
    timepoints= &timepoints,    
    interval=&interval,      

    outc= &outc,          
    outctype=&outctype,  
    outcinteract=&outcinteract,   
    outcwherem = &outcwherem ,  
    outcwherenosim= &outcwherenosim,  
    outcnosimelsemacro = &outcnosimelsemacro,  
    comprisk= &comprisk,        
    compriskinteract= &compriskinteract,   
    compriskwherem = &compriskwherem , 
    compriskwherenosim= &compriskwherenosim,   
    comprisknosimelsemacro = &comprisknosimelsemacro, 
  
     
    numint= &numinta ,
    fixedcov= &fixedcova,      
    ncov= &ncova,          

    
    %do i = 1 %to &ncova ;
        cov&i = &&cova&i ,
        %do j = 1 %to 25 ;
            %let word = %scan(&suffixlist,&j);
            cov&i.&word = &&cova&i.&word ,
        %end;
    %end;
    /*other override options used for each gformula call */

    wherevars= &wherevars,        
    keepsimuldata= &keepsimuldata,     
    equalitiessimuldata= &equalitiessimuldata, 
    eventaddvars= &eventaddvars,  
     
    compriskaddvars= &compriskaddvars,
    
     
 
    usebetadata =&usebetadataa,
    betadata= &betadataa  ,      
    simuldata=  &simuldataa,     
    resultsdata= &resultsdataa,    
    survdata= &survdataa,       
    outputs = &outputs,   
    print_stats = &print_stats ,
    check_cov_models = &check_cov_models,  
    print_cov_means = &print_cov_means,  
    covmeandata = &covmeandataa ,
    save_raw_covmean = &save_raw_covmean,
    observed_surv = &observed_surva,
    intervname = &intervnamea,
    runnc = &runnca ,
    refint = &refinta,
    seed= &seed,          
    nsamples=&nsamples,     
    nsimul=&nsimul,         
    nparam= &nparam,        
    hazardratio=&hazardratio ,  
     
    intcomp =   &intcomp ,  /* needs to be a list with two numbers for comparing the hazard ratio when hazardratio = 1 */
    bootstrap_hazard = &bootstrap_hazard,  
    hazardname = &hazardname.a ,  
   
    sample_start = &sample_start ,  
    sample_end = &sample_end,      
    savelib = &savelib,   


    rungraphs = &rungraphs ,
    title1= &title1a,
    title2= &title2a,
    title3= &title3a,
    titledata= &titledataa ,
    graphfile=&graphfilea ,
    tsize=&tsizea ,     
    weight = &weight,
    printlogstats = &printlogstats,
    usespline = &usespline ,
    checkaddvars = &checkaddvars,
    minimalistic = &minimalistic, 
    testing = &testing  
    ) 

    %if &forrct = 1 %then %do;
        
            proc datasets library = work nolist ;

            change %do intno = %eval(1-&runnca) %to &numinta ; hazard_all_&intno = hazard_all_&intno._a %end ; ;
            quit;
   %end;



%mend;

%macro modelB ;

 
 data datab ;
 set &data(where = (&arm = 1));
 run;

    %let simuldatab = ;
    %let resultsdatab= _resultsB_  ;
    %let survdatab = ;
    %let covmeandatab = ;
    %let observed_survb = ;
    %let intervnameb = ;
    

    %if %bquote(&simuldata)^= %then %let simuldatab = &simuldata.b ;               
    %if %bquote(&survdata)^=  %then %let survdatab = &survdata.b  ;
    %if %bquote(&covmeandata)^=  %then %let covmeandatab = &covmeandata.b ;      
    %if %bquote(&observed_surv)^=  %then %let observed_survb = &observed_surv.b; 
    %if %bquote(&intervname)^=  %then %let intervnameb = &intervname.b ;
    
    
    %if &numintb > 0 %then %do;
        %do i = 1 %to &numintb ;
            %let interv&i = &&intervb&i ;
            %put interv&i for modelB ::: &&interv&i ;
        %end;
    %end;

 %gformula(
    data= datab,           
    id= &id,      
    time= &time,       
    timeptype = &timeptype ,
    timeknots = &timeknots,
    timeinc = &timeinc,
    timefuncgen = &timefuncgen,   
    timepoints= &timepoints,    
    interval=&interval,      

    outc= &outc,          
    outctype=&outctype,  
    outcinteract=&outcinteract,   
    outcwherem = &outcwherem ,  
    outcwherenosim= &outcwherenosim,  
    outcnosimelsemacro = &outcnosimelsemacro,  
    comprisk= &comprisk,        
    compriskinteract= &compriskinteract,   
    compriskwherem = &compriskwherem , 
    compriskwherenosim= &compriskwherenosim,   
    comprisknosimelsemacro = &comprisknosimelsemacro, 
  
     
    numint= &numintb ,
    fixedcov= &fixedcovb,      
    ncov= &ncovb,          

    
    %do i = 1 %to &ncovb ;
        cov&i = &&covb&i ,
        %do j = 1 %to 25 ;
            %let word = %scan(&suffixlist,&j);
            cov&i.&word = &&covb&i.&word ,
        %end;
    %end;
    /*other override options used for each gformula call */

    wherevars= &wherevars,        
    keepsimuldata= &keepsimuldata,     
    equalitiessimuldata= &equalitiessimuldata, 
    eventaddvars= &eventaddvars,  
     
    compriskaddvars= &compriskaddvars,
    
     
 
    usebetadata =&usebetadatab,
    betadata= &betadatab ,      
    simuldata=  &simuldatab,     
    resultsdata= &resultsdatab,    
    survdata= &survdatab,       
    outputs = &outputs,   
    print_stats = &print_stats ,
    check_cov_models = &check_cov_models,  
    print_cov_means = &print_cov_means,  
    covmeandata = &covmeandatab ,
    save_raw_covmean = &save_raw_covmean,
    observed_surv = &observed_survb,
    intervname = &intervnameb,
    runnc = &runncb ,
    refint = &refintb,
    seed= &seed,          
    nsamples=&nsamples,     
    nsimul=&nsimul,         
    nparam= &nparam,        
    hazardratio=&hazardratio ,  
     
    intcomp =   &intcomp ,  /* needs to be a list with two numbers for comparing the hazard ratio when hazardratio = 1 */
    bootstrap_hazard = &bootstrap_hazard,  
    hazardname = &hazardname.a ,  
   
    sample_start = &sample_start ,  
    sample_end = &sample_end,      
    savelib = &savelib,   


    rungraphs = &rungraphs ,
    title1= &title1b,
    title2= &title2b,
    title3= &title3b,
    titledata= &titledatab ,
    graphfile=&graphfileb ,
    tsize=&tsizeb ,  
    weight = &weight,
    printlogstats = &printlogstats,
    usespline = &usespline ,
    checkaddvars = &checkaddvars,
    minimalistic = &minimalistic, 
    testing = &testing  
    ) ;

    
    
    %if &forrct = 1 %then %do;
        
            proc datasets library = work nolist  ;
            change %do intno = %eval(1-&runncb) %to &numintb ; hazard_all_&intno = hazard_all_&intno._b %end ; ;
            quit;
   %end;



%mend ;

%macro hazardrct ;
   %local numint start intno ; 

   %if &runnca = 1 and &runncb = 1 %then %let start = 0;
   %else %let start = 1 ;
  
   %let numint = %sysfunc(min(&numinta, &numintb)); 
 %put _local_ ;
   %do intno = &start %to &numint ;
        data both ;
        set hazard_all_&intno._a (rename = (int = int_a)) hazard_all_&intno._b (rename = (int = int_b)) ;
        run;

        proc sort data = both ;
        by sample ;
        run;

      
        data both ;
        set both ;
 
        if &outc = 1 then event = 1 ;
        else event = 0 ;
        if int_a = &intno then int = 0;
        else int = 1 ;
        run;

 
        ods select none ;
        proc phreg data = both   ;
        ods output ParameterEstimates=_inthr_rct&intno._  ; 
        model newtime*event(0) =  int / rl  ;
        by sample ;
        run;
        ods select all ;

        proc datasets library = work nolist ;
        delete both ;
        quit;
    %end;
%mend;


%macro createhazard ;
    /* this version will replace the createhazard sub-macro supplied by the gformula2.sas code */
    /* modified version of createhazard from gformula macro code.  */
    %local firstint secondint intno runphreg nintcomp;

      %let nintcomp = %numargs(&intcomp);

      %let runphreg = 0 ;

      %if &nintcomp = 2 %then %do;     
        %let firstint = %sysfunc(compress(%scan(&intcomp,1)));
        %let secondint = %sysfunc(compress(%scan(&intcomp,2)));    
        %let runphreg = 1 ;
      %end;
     

      data _calchazard_ ;
      calchazard = 1 ;
      _sample_ = &bsample ;
      run;

       


        %do intno = %eval(1-&runnc) %to &numint;
         
            %if  (&forrct = 0 and ( &intno = &firstint or &intno = &secondint)) or &forrct = 1 %then %do;
                data hazard&intno ;
                set simulated&intno;
                int = &intno;
                sample = &bsample ;
                keep sample int newtime &outc censor ;
                run;
            %end ;
            %if &forrct = 1 %then %do;
                %if &bsample = 0   %then %do;
                    data hazard_all_&intno ;
                    set hazard&intno ;
                    run;               
                %end;            
                %else %if  (&bootstrap_hazard = 1 OR &chunked = 1 )%then %do; 
                        proc append base = hazard_all_&intno data = hazard&intno ;
                        run;                                        
                %end;
            %end;
       %end;
 
       %if &runphreg = 1 %then %do;
        data both ;
        set hazard&firstint hazard&secondint ;
        run;



        data both ;
        set both ;
        if &outc = 1 then event = 1 ;
        else event = 0 ;
        if int = &firstint then int = 0;
        else int = 1 ;
        run;

 
         ods select none ;
         proc phreg data = both   ;
         ods output ParameterEstimates=_inthr0_  ; 
         model newtime*event(0) =  int / rl  ;
         run;
       ods select all ;

        %if &bsample = 0 %then %do;
            proc sql noprint ;
            select HazardRatio into :sample_hazard from _inthr0_ ;
            quit;
            run;

            %put  hazard = &sample_hazard ;
        %end;

 
        %if &bootstrap_hazard = 1 OR &chunked = 1 %then %do; 

            data _inthr0_ ;
            set _inthr0_ (keep = HazardRatio) ;
            _sample_ = &bsample ;
            run;

            %if &bsample = &sample_start %then %do;
                 data &hazardname ;
                 set _inthr0_;
                 run;
            %end;
            %else %do;
                 data &hazardname ;
                 set &hazardname  _inthr0_ ;
                 run;
            %end;

            proc datasets library = work nolist ;
            delete both   _inthr0_ ;
            quit;

        %end;
     %end;
     %else %do; 
        %let sample_hazard = -1 ;
        %if &bootstrap_hazard = 1 OR &chunked = 1 %then %do; 

            data _inthr0_ ;
            hazardratio = -1;
            _sample_ = &bsample ;
            run;

            %if &bsample = &sample_start %then %do;
                 data &hazardname ;
                 set _inthr0_;
                 run;
            %end;
            %else %do;
                 data &hazardname ;
                 set &hazardname  _inthr0_ ;
                 run;
            %end;
       %end;
     %end;

     proc datasets library = work nolist ;
     delete %do intno =  %eval(1-&runnc) %to &numint ; hazard&intno %end ; ;
     quit;

     data _calchazard_ ;
     calchazard = 0 ;
     _sample_ = &bsample ;
     run;
     
%mend ;

%macro rescaleround2;
   
%if &outctype=bineofu or &outctype=binsurv %then %do;
   
  /* FOR BINARY OUTCOME, MULTIPLY RISKS BY 100 TO GET % AND ROUND OFF TO TWO DECIMAL PLACES.*/
   

   &outputname.A = round(&outputname.A*10000)/100;
   &outputname.A_mean = round(&outputname.A_mean*10000)/100;
   &outputname.A_std  = round(&outputname.A_std*10000)/100;  
   &outputname.A_ulim95 = round(&outputname.A_ulim95*10000)/100;
   &outputname.A_llim95 = round(&outputname.A_llim95*10000)/100;


   
   &outputname.B = round(&outputname.B*10000)/100;
   &outputname.B_mean = round(&outputname.B_mean*10000)/100;
   &outputname.B_std  = round(&outputname.B_std*10000)/100;  
   &outputname.B_ulim95 = round(&outputname.B_ulim95*10000)/100;
   &outputname.B_llim95 = round(&outputname.B_llim95*10000)/100;
   
   RD = round(RD*10000)/100;

   RD_mean = round(RD_mean*10000)/100;
   RD_std  = round(RD_std*10000)/100;

   
   RD_ulim95 = round(RD_ulim95*10000)/100;
   RD_llim95 = round(RD_llim95*10000)/100;
   

   NNT = round(NNT);

   NNT_mean = round(NNT_mean);
   NNT_std  = round(NNT_std);

   

   
      NNT_ulim95 = round(NNT_ulim95);
      NNT_llim95 = round(NNT_llim95);
   
   
%end;
%else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
   /**
   obsp= round(&obsp*100)/100;
   call symput('obspm',obsp);
**/
   s&outc = round(s&outc*100)/100;

   s&outc._mean = round(s&outc._mean*100)/100;
   s&outc._std  = round(s&outc._std*100)/100;
   
   s&outc._ulim95 = round(s&outc._ulim95*100)/100;
   s&outc._llim95 = round(s&outc._llim95*100)/100;
   

   RD = round(RD*100)/100;

   RD_mean = round(RD_mean*100)/100;
   RD_std  = round(RD_std*100)/100;
   
   if int^=0 then do;
      RD_ulim95 = round(RD_ulim95*100)/100;
      RD_llim95 = round(RD_llim95*100)/100;
      end;
   
%end;

/* in both cases do the following: */
   RR = round(RR*100)/100;

      RR_ulim95 = round(RR_ulim95*100)/100;
      RR_llim95 = round(RR_llim95*100)/100;
   
   intervened = round(intervened*10000)/100;
   averinterv = round(averinterv*10000)/100;
   
%mend rescaleround2;

%macro labels2;
   label int2_A= "Interventions under A";
   label int2_B= "Interventions under B";

   label &outputname.A_ulim95 = 'Upper limit 95% CI under A';
   label &outputname.A_llim95 = 'Lower limit 95% CI under A';
   label &outputname.B_ulim95 = 'Upper limit 95% CI under B';
   label &outputname.B_llim95 = 'Lower limit 95% CI under B';
   label RR_ulim95 = 'Upper limit 95% CI';
   label RR_llim95 = 'Lower limit 95% CI';
   label RD_ulim95 = 'Upper limit 95% CI';
   label RD_llim95 = 'Lower limit 95% CI';
   label intervened = '% Intervened On';
   label averinterv = 'Aver % Intervened On';
   label int       = 'Interv.';
   label int2      = 'Description';
   
%if &outctype=binsurv or &outctype=bineofu %then %do;
   label &outputname.A        = 'Risk (%) under A';
   label &outputname.A_std    = 'Bootstrap Risk SE under A';
   label &outputname.A_mean   = 'Bootstrap Risk Mean under A';
   label &outputname.B        = 'Risk (%) under B';
   label &outputname.B_std    = 'Bootstrap Risk SE under B';
   label &outputname.B_mean   = 'Bootstrap Risk Mean under B';
   label rr        = 'Risk ratio';
   label RD        = 'Risk difference';
   label NNT       = '# Needed to Treat';
   label NNT_ulim95 = 'Upper limit 95% CI';
   label NNT_llim95 = 'Lower limit 95% CI';
   %end;

%else %if &outctype=conteofu %then %do;
   label pD        = 'Mean';
   label pD_std    = 'Bootstrap Mean SE';
   label pD_mean   = 'Bootstrap Mean Mean';
   label rr        = 'Ratio of means';
   label RD        = 'Difference of Means';
   %end;
   
%mend labels2;

%macro results2 ;
  
    %let numint = %sysfunc(min(&numinta,&numintb));
    %let runnc = %sysfunc(min(&runnca , &runncb)) ;

    %put (numint, runnc)  &numint &runnc ;
    data fin ;  run ;
    data temp ; run ;
    %if &forrct = 1 %then %do ; 
        data hrtemp ; run ; 
        data hrfin ; run;
    %end;

%if &outctype = binsurv %then %do;
    %let inputname = risk&timepoints ;
    %let outputname = risk ;
%end;
%else %do ;
    %let inputname = s&outc ;
    %let outputname = s&outc ;
%end;
 
     proc sql noprint ;
     select obsp into :obsp_a from _resultsA_ ;
     select obsp into :obsp_b from _resultsB_  ;
     quit;


     %do i= %eval(1-&runnc) %to &numint;
  
        data intA_&i ;
        set &survdata.a ;
        where int = &i ;
        keep &inputname int int2 _sample_ ;
        rename &inputname = &outputname.A int2 = int2_A ;
        run;

        data intB_&i ;
        set &survdata.b ;
        where int = &i ;
        keep &inputname int int2 _sample_ ;
        rename &inputname = &outputname.B int2 = int2_B ;
        run;
 

        data interv&i; 
        merge intA_&i intB_&i;
        by _sample_; 
        if &outputname.A^=0 then rr=&outputname.B / &outputname.A;
        if riskA^=0 then rd=&outputname.B - &outputname.A ;
        if rd^=. and rd^=0 then nnt = 1/rd;       
        run;

     
          %*Appending intervention datasets;
          data fin; 
          set fin interv&i; 
          if _sample_=0;
          if int=. then int=&i;
          run;
      


          %*Calculating bootstrap mean, variance and confidence intervals;
          proc univariate data=interv&i noprint;
          where _sample_ ne 0;
          var &outputname.A &outputname.B rr rd nnt;
          output out = temp&i
          mean = &outputname.A_mean &outputname.B_mean RR_mean RD_mean NNT_mean
          std =  &outputname.A_std  &outputname.B_std RR_std  RD_std  NNT_std
          pctlpre = &outputname.A_ &outputname.B_ RR_     RD_     NNT_
          pctlname = llim95 ulim95  pctlpts = 2.5 97.5;
          run;

          data temp&i;
          set temp&i;
          int = &i;
          run;

          data temp ;
          set temp temp&i;
          run;

          %if &forrct = 1 %then %do;

              data hrfin ;
              set hrfin _inthr_rct&i._ (where = (sample = 0)  keep = sample hazardratio ) ;
              if int = . then int = &i ;
              if hazardratio ne . ;
              drop sample ;
              run;

              proc univariate data = _inthr_rct&i._ noprint ;
              where sample ne 0 ;
              var hazardratio ;
              output out = hrtemp&i 
              mean = rh_mean 
              std = hr_std 
              pctlpre = hr_  
              pctlname = llim95 ulim95  pctlpts = 2.5 97.5;
              run;

              data hrtemp&i;
              set hrtemp&i ;
              int = &i ;
              run;

              data hrtemp ;
              set hrtemp hrtemp&i;
              run;

          %end;
   %end;
         

     



     %if &forrct = 1 %then %do;
         data hrfin ;
         merge hrfin hrtemp ;
         by int ;
         if int ne . ;
         label hr_ulim95 = 'Upper limit 95% CI';
         label hr_llim95 = 'Lower limit 95% CI';
         hazardratio = round(hazardratio,0.01);
         HR_ulim95 = round(HR_ulim95,0.01 )  ;
         HR_llim95 = round(HR_llim95, 0.01)  ;
         run;
     %end;
     %*Cleaning up results to print nicely; 


     
     data fin ;
     merge fin temp ;
     by int ;
     if int ne . ;
     run;

  

     data finfin ;
     set fin;      
     
     %rescaleround2; /* RESCALE AND ROUND OFF THE OUTPUT */
     %labels2;       /* LABEL THE OUTPUT */
     run;
 

     %*Outputting results dataset;
     %if %bquote(&resultsdata)^= %then %do;
          %if &printlogstats = 1 %then %put ;
          %if &printlogstats = 1 %then %put  Outputting results to &resultsdata;
          data &resultsdata;
          set finfin;
          run;

          %if &forrct = 1 %then %do;
            data myresults_hr ;
            set hrfin ;
            run;
          %end;
     %end;
     %if &printlogstats = 1 %then %put ;

     %*Printing results;

     %if &outctype=binsurv or &outctype=bineofu %then %do;    
          title4 "PREDICTED RISK UNDER SEVERAL INTERVENTIONS UNDER TWO ARMS :"; title5 " A (&arm = 0)  and B (&arm = 1 )";
     %end;
     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
          title4 "PREDICTED MEAN &outc UNDER SEVERAL INTERVENTIONS UNDER TWO ARMS : A (&arm = 0) and B (&arm = 1)";
     %end;

 
     proc print data=finfin noobs label double;
     var int  int2_A int2_B ;
     run;


      %if &outctype=binsurv or &outctype=bineofu %then %do;
          title6 "Observed risk= %sysevalf(&obsp_a) for arm A and %sysevalf(&obsp_b) for arm B.  ";
     %end;      
     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
          title6 "Observed mean= %sysevalf(&obsp_a) for arm A and %sysevalf(&obsp_b) for arm B. ";
     %end;
     title7 "Data= &data ";
     title8 "Number of bootstrap samples= &nsamples";
     title9 "Reference is using arm A";


     proc print data=finfin noobs label double; 
     var int &outputname.A  &outputname.A_mean &outputname.A_llim95 &outputname.A_ulim95 &outputname.B &outputname.B_mean  &outputname.B_llim95 &outputname.B_ulim95  /* intervened averinterv */; 
     run;

     proc print data=finfin noobs label double; 
     var int &outputname.A    &outputname.B   rr rr_llim95 rr_ulim95  /* intervened averinterv */; 
     run;

     proc print data=finfin noobs label double;    
     var int &outputname.A &outputname.B rd rd_llim95 rd_ulim95 nnt nnt_llim95 nnt_ulim95;
     run;


     %if &forrct = 1 %then %do;
                
        title4 "PREDICTED HAZARD RATIO UNDER SEVERAL INTERVENTIONS UNDER TWO ARMS :  A (&arm = 0)  and B (&arm = 1 )";
        title5 "Data= &data ";
        title6 "Number of bootstrap samples= &nsamples";
        title7 "Reference is using arm A";
        proc print data = hrfin noobs label double ;
        var int hazardratio hr_llim95 hr_ulim95 ;
        run;

     %end;
     title;
 
     %* Deleting no longer needed datasets;
      
     proc datasets library=work nolist; 
     delete   _ref_ fin finfin dataa datab  temp       
         %if &forrct = 1 %then %do ;
               _inthr_  hrfin hrtemp 
               %do i = %eval(1-&runnc) %to &numint ; _inthr_rct&i._ hrtemp&i hazard_all_&i._a hazard_all_&i._b %end;
         %end;
         %do i = %eval(1-&runnc) %to &numint ; temp&i inta_&i intb_&i %end;
         
         ;             
     quit; 


%mend ;
