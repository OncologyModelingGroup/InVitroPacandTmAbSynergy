%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Begin invitro_optifun.m

%Function file for called by the optimization/minimization solver for
%calibrating the parameters of the model

%Methods described in

%2019 Scientific Reports 
%Experimentally-driven mathematical modeling to improve combination 
%targeted and cytotoxic therapy for HER2+ breast cancer

%in the subsection "Parameter Calibration" in the "Methods" of the
%manuscript.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%This file is called to evalute the model for the test parameters against
%the data and return a weighted error evaluation based on the confidences
%of the data at each time point

%The function receives all of the input parameters required to run the
%model as well as the corresponding data to compare the simulation.
%Additional parameter indicate to the function which parameters are free
%for optimization

%The function returns the tested parameter set and the resulting errors for
%the resulting simulation compared to the data to the solver that called it
%for calibrating the parameters.

%Angela M. Jarrett (ajarret@utexas.edu)
%The University of Texas at Austin
% https://cco.oden.utexas.edu/
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


function [out,p] = invitro_optifun(p,init,tf,Af,Pf,which,params,datapoints,dataweights,datatimes,flag,Bests,starts)

    %This flag indicates we are only calibrating a control data set or the
    %first 24 hours of a data set where drug was delivered after the first
    %day--therefore, only the growth rate and carrying capacity are
    %considered.
    
    %The p variable are the test parameter values given by the optimization
    %function
    
    %Here the calibration is over the first 24 hours or for a control data
    %set
    if flag == 0
        k = p(1);
        theta = p(2);
    end

    %Making call to ODE solver for solutions while passing in the 
    %control parameter p
    
    %Here the calibration is over the first 24 hours or for a control data
    %set
    if flag == 0
        [T,X] = ode45(@(t,y) k*y*(1-y/theta), [starts tf], init(1));
    
    %Here the calibration is for a post drug dose period
    elseif flag == 1
        %Defining the growth rate and carrying capacity parameters
        %determined from an earlier calibration of the first 24 hours
        params([1,5]) = Bests(:,nn);
        clear t x
        [T,X] = ode45(@rhsinvitro,[starts(nn) tf],init,[],p,Af(nn,:),Pf(nn,:),which,params);
    end
    %Checking for undefined values
    X(isnan(X)) = 0;
    
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Comparing ODE solution to data
    %Determing corresponding time points for each model output to compare  
    %back to the data. 
    OutputPoints = interp1(T,X(:,1),datatimes,'pchip');
    %Calculate the difference between the simulation output points and the
    %data points
    diff = OutputPoints - datapoints;
    %Divide by the data confidence to weigh the error of each point
    %comparison
    diff = diff./dataweights';
    %caculate the norm
    err = norm(diff);
    
%Sending out the norm of the error between the model and the data points 
out = err;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%end of file