%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Begin CalibVerifControl.m

%Code for assessing the calibration method for its ability to recover
%parameters for growth and carrying capacity using only the first 24 hours
%of data

%Methods described in

%2019 Scientific Reports 
%Experimentally-driven mathematical modeling to improve combination 
%targeted and cytotoxic therapy for HER2+ breast cancer

%in the subsection "Parameter Calibration" in the "Methods" of the
%manuscript.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%This file generates in silico simulations of the logisitic growth model
%using randomly sampled parameter values for growth rate (k) and carrying
%capacity (theta).

%The script then runs the calibration method only considering the first 24
%hours of data points to recover the random parameters for all in silico
%data sets.

%All of the in silico data, randomly selected parameter values, and
%resulting calibrations and their corresponding errors are saved as text or
%.mat files.

%Finally, the results are visualized in a figure.

%To perform calibration verifcation over a different time period (such as
%over the whole time period or data sets where drugs are administered after
%an intial period of time) change the intial and final times for the
%calibration section of the code below along with any other changes that
%would be required such as defining additional parameters.

%Files required: control data file, optimization function file, and model
%definition file

%Angela M. Jarrett (ajarret@utexas.edu)
%The University of Texas at Austin
% https://cco.oden.utexas.edu/
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% INITIALIZING WORKSPACE
clear all
close all

%Can run in batch if multiple control sets are available, only using
%control data to define experimentally relevant data times and weights on
%the data for the calibration method
load('ControlSet.mat');%modify if multiple controls are available

%Define the final time for the simulation
tf = 4;

%Free drug concentrations are zero for control data
Af = 0;  %trastuzumab
Pf = 0;  %paclitaxel

%Number of parameters fitted for the control data (k and theta) 
Nparameters = 2;

%The sample size here is the number of in silico simulations we want to
%test the calibration method against
SampleSize = 100;

%Define the number of intial random guesses for the calibration
multiguesssize = 100;

%The intial condition for the confluence of cells, and the two zeros
%following are the drug concentrations (none for controls)
initials = [means(1),0,0];
%Define the start of the simulation (not all data sets start exactly at
%zero)
starts = datatimes(1);
%Find the final data point closest to the defined final simulation time
%and redefine the vectors to include only those points for the
%calibration
a = max(find(datatimes<tf));
conf95 = conf95(1:a);
datatimes = datatimes(1:a);

%Upper and lower bounds for the calibration of each parameter
lower  = [0,  0];
upper  = [2,  1]; 

%Creating a matrix of random numbers to create test parameter sets (using
%uniform distributions)
RandNums = lower'+(upper'-lower').*rand(Nparameters,SampleSize);  
%Saving test parameters
delete('RandNums.txt');
save('RandNums.txt','RandNums','-ascii');

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% IN SILICO DATA GENERATION

%Define a holding matrix for the resulting simulations from the test
%parameters
TrainingPointsT = zeros(SampleSize,size(datatimes,1));
%Generate in silico data using each of the test parameter sets for the
%simple logistic equation to simulate control data
for j = 1:SampleSize
    clear x t
    %Define random parameter values for simulation
    k = RandNums(1,j);
    theta = RandNums(2,j);
    %Solve ODE
    [t,x] = ode45(@(t,y) k*y*(1-y/theta), [starts tf], initials(1));
    %Record resulting in silico results for the datatimes for the control
    %data reference
    TrainingPointsT(j,:) = interp1(t,x(:,1),datatimes,'pchip');
end
%Save the in silico data
save TrainingPointsT.mat TrainingPointsT -v7.3;


%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% CALIBRATION

%Calibrate the model parameters for each of the in silico generated data

%Holding matrix for each of the parameters recovered for each of the in
%silico data sets
Optimals = zeros(Nparameters,SampleSize);

%Redefine the final time to one day to determine how well the calibration
%method is able to recover the parameters of the in silico data sets using
%only the data points from the first 24 hours.
tf = 1;
initials = [means(1),0,0];
starts = datatimes(1);
a = max(find(datatimes<tf));
conf95 = conf95(1:a);
datatimes = datatimes(1:a);

%Loop through all of the in silico data sets for calibration
for ee = 1:SampleSize
    
    %Reintialize the mulitple first guesses for each iteration of the loop
    clear MultiGuessNums

    %Indicates which of the parameters are free for calibration, for the
    %control there are only two    
    which  = [1,  1]; 
    free = find(which == 1);
    
    %Create a matrix of random initial guesses for the calibration
    MultiGuessNums = zeros(multiguesssize,Nparameters);
    for nn = 1:multiguesssize
        MultiGuessNums(nn,:) = lower(free)+(upper(free)-lower(free)).*rand(1,sum(which));
    end
    
    %Define holding matrices for all the resulting calibrated parameters
    %and the corresponding errors of the resulting simulation compared to
    %the in silico data
    Guesses = zeros(multiguesssize,sum(which));
    errors = zeros(multiguesssize,1);
    %Using each of the initial guesses, calibrate the model to the in
    %silico data set
    for yy = 1:multiguesssize
        clear control FVAL params
        %Define our initial parameter guess
        params = MultiGuessNums(yy,:);
        %Send in data and parameter values to fmincon, which calls the
        %function invitro_optifun where the confidence of the data points
        %is used to weight the calibration
        [control,FVAL,EXITFLAG,OUTPUT] = fmincon(@(p) invitro_optifun(p,initials,tf,Af,Pf,free,params,TrainingPointsT(ee,1:a)',conf95',datatimes,0,0,starts),params(free),[],[],[],[],lower(free),upper(free),[]);
        
        %Save the resulting calibrated parameters and record the error for
        %the fit of the simulation to the data
        Guesses(yy,:) = control;
        errors(yy) = FVAL;
    end
    %Find the result with the minimum error
    minimum = min(errors);
    z = find(errors==minimum);
    Optimals(:,ee) = Guesses(z(1),:)';    
end

%The resulting parameters from the calibrations (Optimals) can be compared
%to the random numbers generated for the in silico data (RandNums) to
%determine the parameter recovery errors

%Save the calibrated parameters
save('OptimalsCalibVerifControl.mat','Optimals')

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% RESULTS

%Visulaize the results for the best calibration for each of the in silico
%data sets

%Reload and reinitialize values for full 4 days
tf = 4;
load('ControlSet.mat');
initials = [means(1),0,0];
starts = datatimes(1);
a = max(find(datatimes<tf));
conf95 = conf95(1:a);
datatimes = datatimes(1:a);

%Open a figure window
figure;
hold;
%Loop through for all the in in silico data sets
for sets = 1:SampleSize
    clear y t
    %Define parameter values from calibration result
    k = Optimals(1,sets);
    theta = Optimals(2,sets);
    %Simulate the simple logisitic growth model using the parameters that 
    %resulted in the minimum error between the model and in silico data
    [t,yop] = ode45(@(t,y) k*y*(1-y/theta), [datatimes(1) tf], initials(1));
    
    %Plotting commands
    %in silico data
    errorbar(datatimes,TrainingPointsT(sets,:),conf95','Linewidth',3)
    %80% confluence level for reference
    conf80 = 0.8*ones(size(t));
    plot(t,yop(:,1),'m','linewidth',1)
    plot(t,conf80,'r--','linewidth',3)
    axis([0 tf 0 1])
    xlabel('Time in days')
    ylabel('Confluence (fractionated)')
    
    legend('in silico data','Fitted model simulation','80% confluence','location','southeast')
    legend boxoff
    set(gca,'FontSize',20,'FontName','Times New Roman')
    
    %Create a vector of points from the simulation for the corresponding
    %data times
    b = interp1(t,yop(:,1),datatimes,'pchip');
    
    %the predicted points b can be compared to the corresponding points of the
    %in silico data ("training points") for correlation coeffecients and
    %errors
    
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%end of file