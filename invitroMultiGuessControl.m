%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Begin invitroMultiGuessControl.m

%Code for generating calibrated parameters from the control data using
%multiple initial guesses for the calibration function

%Methods described in

%2019 Scientific Reports 
%Experimentally-driven mathematical modeling to improve combination 
%targeted and cytotoxic therapy for HER2+ breast cancer

%in the subsection "Parameter Calibration" in the "Methods" of the
%manuscript.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%This file generates calibrated parameters for control data and saves the
%results as .mat files.

%The results are visualized in a figure at the end.

%Please see the additional file for multi-guess calibration for the single
%drug dose of trastuzumab for an example of the sequential calibration
%found in our paper.

%Files required: control data file, optimization function file, and model
%definition file

%Angela M. Jarrett (ajarret@utexas.edu)
%The University of Texas at Austin
% https://cco.oden.utexas.edu/
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% INITIALIZING WORKSPACE
clear all
close all

%Can run in batch if multiple control sets are available
numberofControls = 1;

%Define the final time for the simulation
tf = 4;

%Free drug concentrations are zero for control data
Af = 0;  %trastuzumab
Pf = 0;  %paclitaxel

%Parameters: 
%params = [k, theta]
%Intial values for the parameters
params = [0.58,0.78];

%Define the number of intial random guesses for the calibration
multiguesssize = 100;

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% CALIBRATION

%If there are multiple control sets, the following loop can be modified to
%calibrate the model to each in batch
for ee = 1:numberofControls
    
    %Calibrate controls 1 at a time
    load('ControlSet.mat'); %modify if multiple controls are available
    %Data file provides data times, means of the data set at each time, and
    %the corresponing 95% confidence intervals and standard deviations per
    %time point
    
    %Indicates which of the parameters are free for calibration, for the
    %control there are only two
    which  = [1,  1];
    free = find(which == 1);
    
    %General lower and upper bounds for the calibration of each parameter
    lower  = [0,  0];
    upper  = [2,  1];
    
    %The intial condition for the confluence of cells, and the two zeros
    %following are the drug concentrations (none for controls)
    init = [means(1),0,0];
    %Define the start of the simulation (not all data sets start exactly at
    %zero)
    starts = datatimes(1);
    %Find the final data point closest to the defined final simulation time
    %and redefine the vectors to include only those points for the
    %calibration
    a(ee) = max(find(datatimes<tf));
    means = means(1:a(ee));
    conf95 = conf95(1:a(ee));
    datatimes = datatimes(1:a(ee));
    
    %Create a matrix of random initial guesses for the calibration
    for nn = 1:multiguesssize
        %Use a uniform distribution for both parameters
        %r = a + (b-a).*rand(N,1);
        RandNums(nn,:) = lower+(upper-lower).*rand(1,size(params,2));
    end
    
    %Using each of the initial guesses, calibrate the model to the data
    for yy = 1:multiguesssize
        %Define initial parameter guess
        params = RandNums(yy,:);
        %Send in data and parameter values to fmincon, which calls the
        %function invitro_optifun where the confidence of the data points
        %is used to weight the calibration
        [control,FVAL,EXITFLAG,OUTPUT] = fmincon(@(p) invitro_optifun(p,init,tf,Af,Pf,free,params,means,conf95',datatimes,0,0,starts),params(free),[],[],[],[],lower(free),upper(free),[]);
        
        %Save the resulting calibrated parameters and record the error for
        %the fit of the simulation to the data
        Guesses(yy,1:size(params,2),ee) = control;
        errors(yy,ee) = FVAL;
    end
end

%Save the calibrated parameters and their corresponding errors
save('GuessesControls.mat','Guesses')
save('ErrorsControls.mat','errors')

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% RESULTS

%Visualize the results for the best calibration for each control set
%(Similar to figure 5 in the manuscript)
for sets = 1:numberofControls
    %Open a figure window
    figure;
    
    %Define holding vectors for the calibrated parameters of the control
    %set
    ks = Guesses(:,1,sets);
    thetas = Guesses(:,2,sets);
    %Find the calibration set that resulted in the lowest error between the
    %model simulation and the data
    minimum = min(errors(:,sets));
    %Remove any sets that have errors greater than the minimum
    z = find(errors(:,sets)>ceil(minimum));
    c = unique(sort(z));
    ks(c) = [];
    thetas(c) = [];
    
    %For one control set this is not needed, but for looping through
    %multiple sets the data will need to be reloaded
%     load('ControlSet.mat'); %modify for multiple control sets
    
    %Simulate the simple logisitic growth model using the mean for each
    %calibrated parameter that resulted in the minimum error
    [t,yop] = ode45(@(t,y) mean(ks)*y*(1-y/mean(thetas)), [datatimes(1) tf], means(1));
    %80% confluence level for reference
    conf80 = 0.8*ones(size(t));
    
    %Plotting commands
    plot(t,yop(:,1),'m',t,conf80,'r--','linewidth',3)
    axis([0 tf 0 1])
    xlabel('Time in days')
    ylabel('Confluence (fractionated)')
    hold
    errorbar(datatimes,means',conf95','Linewidth',3)
    legend('Fitted model simulation','80% confluence','Control data','location','southeast')
    legend boxoff
    set(gca,'FontSize',20,'FontName','Times New Roman')
    hold
    
    %Create a vector of points from the simulation for the corresponding
    %data times
    b = interp1(t,yop(:,1),datatimes,'pchip');
    
    %The simulation can be compared to the data means at the specific times
    %of data collection (datatimes) using the array b using correlation
    %coeffecients and/or calculating the errors.

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%end of file