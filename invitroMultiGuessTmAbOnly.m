%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Begin invitroMultiTmAbOnly.m

%Code for generating calibrated parameters from the trastuzumab single drug
%dose data using multiple initial guesses for the calibration function.

%Methods described in

%2019 Scientific Reports 
%Experimentally-driven mathematical modeling to improve combination 
%targeted and cytotoxic therapy for HER2+ breast cancer

%in the subsection "Parameter Calibration" in the "Methods" of the
%manuscript.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%This file generates calibrated parameters for trastuzumab treated data and
%saves the results as .mat files.

%Specifically this file shows how to do a sequential calibration over
%certain intervals of the time course data for particular parameter groups
%in the model.

%The results are visualized in a figure at the end.

%Files required: data file, optimization function file, and model
%definition file

%Angela M. Jarrett (ajarret@utexas.edu)
%The University of Texas at Austin
% https://cco.oden.utexas.edu/
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% INITIALIZING WORKSPACE
clear all
close all

matname = 'TmAb24hr10ugmL.mat';
if exist(matname,'file') == 0
   button = questdlg(['Thank you for your interest. This file will not' ...
                      ' run without a corresponding .mat data file. If' ...
                      ' an investigator would like to work with our ' ...
                      'treatment data, please contact us to establish ' ...
                      'a formal collaboration. Otherwise, please ' ...
                      'modify the code to load an appropriate data ' ...
                      'file.'],'Disclaimer','Continue','Continue');
else

    %Define values for cell size to approximate number of cells
    cellradius = 0.00125; %BT474
    areapercell = pi*(cellradius^2); %cm^2
    areaofdish = 0.32; %cm^2
    CellNumConvert = areaofdish/areapercell;

    %Define trastuzumab parameters for drug concentrations
    TmAbmolweight = 145531.5; %g/mol
    AvoNum = 6.02*10^23; %particles/mol Avogadro's number
    TmAbweight = (TmAbmolweight/AvoNum)*10^6; %ug
    Adosevolume = 0.003125; %mL

    %Define the number of intial random guesses for the calibration
    multiguesssize = 100;

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %CALIBRATE GROWTH RATE AND CARRYING CAPACITY
    
    %For the first 24 hours, cells are free to grow, so the inheret growth
    %rate and carrying capacity is calibrated from the data from the data
    %points from the first day of the in vitro results
    
    %Final time, 24 hours
    tf = 1;
    %Free trastuzumab and paclitaxel drugs are zero
    Af = 0;  %trastuzumab
    Pf = 0;  %paclitaxel
    
    %Flag to pass into the calibration function file to indicate
    %calibration of only the two inherent growth parameters (see
    %invitro_optifun.m)
    flag = 0;  
    
    load(matname);
    %params= [k, theta]
    %Intial values for the parameters
    params = [0.76,0.67];
    %Indicates which of the parameters are free for calibration, for the
    %first 24 hours there are only two
    which  = [1,  1];
    free = find(which == 1);
    
    %Lower and upper bounds for the calibration of each parameter defined
    %from the control data sets
    lower  = [0.36,  0.44];
    upper  = [1.12,  0.91];
    
    %The intial condition for the confluence of cells, and the two zeros
    %following are the drug concentrations (none for first 24 hours). Use
    %intial time point two after intial seeding density
    init = [means(2),0,0];
    %Define the start of the simulation (not all data sets start exactly at
    %zero)
    starts = datatimes(2);
    %Find the final data point closest to the defined final simulation time
    %and redefine the vectors to include only those points for the
    %calibration
    a = max(find(datatimes<tf));
    means = means(2:a);
    conf95 = conf95(2:a);
    datatimes = datatimes(2:a);
    
    %Create a matrix of random initial guesses for the calibration
    for nn = 1:multiguesssize
        %Use a uniform distribution for both parameters
        %r = a + (b-a).*rand(N,1);
        RandNums(nn,:) = lower+(upper-lower).*rand(1,size(params,2));
    end
    
    %Using each of the initial guesses, calibrate the model to the data
    for yy = 1:multiguesssize
        %Define  initial parameter guess
        params = RandNums(yy,:);
        
        %Send in data and parameter values to fmincon, which calls the
        %function invitro_optifun where the confidence of the data points
        %is used to weight the calibration
        [control,FVAL,EXITFLAG,OUTPUT] = fmincon(@(p) invitro_optifun(p,init,tf,Af,Pf,free,params,means,conf95',datatimes,flag,0,starts),params(free),[],[],[],[],lower(free),upper(free),[]);
        
        %Save the resulting calibrated parameters and record the error for
        %the fit of the simulation to the data
        FirstDayParams(yy,1:size(params,2)) = control;
        errors(yy) = FVAL;
    end


    %Save the calibrated parameters and their corresponding errors
    save(['FirstDayParams' matname(1:end-4) '.mat'],'FirstDayParams')
    save(['FirstDayParamsErrors' matname(1:end-4) '.mat'],'errors')

    %Determine the parameters that resulted in the least error between the data
    %and the model and save in a separate file
    Bests = zeros(size(params,2),1);
    z = find(errors==min(errors));
    Bests = FirstDayParams(z(1),:);
    save(['FirstDayParamsBests' matname(1:end-4) '.mat'],'Bests');
    
    %Note that the growth and carrying capacity parameters will be fixed
    %for the second calibration over the remaining time course.
    
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% CALIBRATE DRUG DOSE PARAMETERS
    
    %For the remaining time course data of the treated in vitro results,
    %calibrate the parameters associated with the drug dose.
    
    %Flag to pass into the calibration function file to indicate
    %calibration of the drug associated parameters (see invitro_optifun.m)
    flag = 1;
    
    %Redefine our final time as 4 days
    tf = 4;
    
    %Parameters: Note the S parameters are synergy parameters that can be
    %explored with multiple drug doses, S1 is the synergy parameter
    %presented in detail in the paper cited above
    %params = [   k,  S1,  nA, dP, theta,  S2,   bA, S3, aP, gP, tP, tA, S4, S5, S6]
    params  = [0.58,   1, 0.8,  0,  0.78,   1, 7.97,  1,  0,  0,  0,  2,  1,  1,  1];

    %Indicates which of the parameters are free for calibration, for the
    %trastuzumab there are 3 parameters
    which  = [   0,   0,   1,  0,     0,   0,    1,  0,  0,  0,   0,   1,  0,  0,  0];
    free = find(which == 1);
    
    %General lower and upper bounds for the calibration of each parameter
    lower  = [0,0,0,0,0,0,5,0,0,0,0,1,0];
    upper  = [0,0,5,0,0,2,10,5,0,0,0,3,0];

    %Create a matrix of random initial guesses for the calibration for only
    %the freely calibrated parameters
    RandNums = [];%clearing variable for new set
    for nn = 1:multiguesssize
        %Use a uniform distribution for both parameters
        %r = a + (b-a).*rand(N,1);
        RandNums(nn,:) = params;
        RandNums(nn,free) = lower(free)+(upper(free)-lower(free)).*rand(1,size(free,2));
    end

    %Reload file and reintialize variables
    load(matname);
    
    %Find the final data point closest to the defined final simulation time
    %and redefine the vectors to include only those points for the
    %calibration
    b = max(find(datatimes<tf));
    
    %Here we want to only calibrate over the treated time interval (not the
    %first 24 hours, so define our data over that second time period
    means = means((a+2):b)*CellNumConvert;
    conf95 = conf95((a+2):b)*CellNumConvert;
    datatimes = datatimes((a+2):b);
    init = [means(1)*CellNumConvert,0,0];
    starts = datatimes(1);
    
    %Define the dosage given from the name of the .mat file
    if(size(matname,2) > 20) %ug/mL
        Dosages = str2num(matname(9:11));
    else
        Dosages = str2num(matname(9:10));
    end
    %Format is [dose,duration] where the units are [0 or 1 day]
    Af = Dosages'*Adosevolume/TmAbweight; %number of antibodies
    Af = cat(2,Af,ones(size(Dosages,2),1));
    %Free paclitaxel is zero
    Pf = Af;
    Pf(:,1) = Pf(:,1)*0;
    
    %Using each of the initial guesses, calibrate the model to the data
    Guesses = zeros(multiguesssize,size(free,2));
    for yy = 1:multiguesssize
        params = RandNums(yy,:);
        [control,FVAL,EXITFLAG,OUTPUT] = fmincon(@(p) invitro_optifun(p,init,tf,Af,Pf,free,params,means,conf95',datatimes,flag,Bests,starts),params(free),[],[],[],[],lower(free),upper(free),[]);
        Guesses(yy,:) = control;
        Gerrors(yy) = FVAL;
    end

    %Save the calibrated parameters and their corresponding errors    
    save(['Guesses' matname(1:end-4) '.mat'],'Guesses')
    save(['GErrors' matname(1:end-4) '.mat'],'Gerrors')

    %Determine the parameters that resulted in the least error between the data
    %and the model and save in a separate file
    GBests = zeros(size(params,2));
    z = find(Gerrors==min(Gerrors));
    GBests = Guesses(z(1),:);
    save(['GBests' matname(1:end-4) '.mat'],'GBests');


%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% RESULTS

    %Visulize the results for the best calibration
    %(Similar to figure 5 in the manuscript)
    figure;
    %Reload the data
    load(matname);

    %Define the growth rate and carrying capacity parameters from the intial
    %calibration over the first 24 hours
    params([1,5]) = Bests(1:end);
    %Define the separate time intervals for the simulation
    a = max(find(datatimes<1));
    c = max(find(datatimes<tf));
    %Solve the ODE for the simulation
    [t1,yop1] = ode45(@rhsinvitro, [datatimes(2) datatimes(a)], [means(2)*CellNumConvert,init(2:end)], [], GBests, Af, Pf, free, params);
    [t2,yop2] = ode45(@rhsinvitro, [datatimes(a+2) tf], [means(a+2)*CellNumConvert,init(2:end)], [], GBests, Af, Pf, free, params);
    %80% confluence level for reference converted to cell number
    conf80 = 0.8*ones(size([t1',t2']))*CellNumConvert;
    
    %Plotting commands
    plot([t1',t2'],[yop1(:,1)',yop2(:,1)'],'m',[t1',t2'],conf80,'r--','linewidth',3)
    xlabel('Time in days')
    ylabel('Cells')
    hold
    errorbar(datatimes(2:c),means(2:c)'*CellNumConvert,conf95(2:c)'*CellNumConvert,'Linewidth',3)
    legend(['Drug Model Fitted, Dose = ' num2str(Dosages) ' ug/mL'],'80% Confluence','in vitro data','location','northwest')
    set(gca,'FontSize',20,'FontName','Times New Roman')

    %Create a vector of points from the simulation for the corresponding
    %data times
    b = interp1([t1',t2'],[yop1(:,1)',yop2(:,1)']/CellNumConvert,datatimes(2:c),'pchip');

    %The simulation can be compared to the data means at the specific times
    %of data collection (datatimes) using the array b using correlation
    %coeffecients and/or calculating the errors.
    
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%end of file