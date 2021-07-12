
clear all

% ############## Parameters ################
CommunicationTime = 5; %timestep in minute for realtime communication between VPP platform and batteries
ForecastTimeStep = 1/2; %Time Step in hour (between 1/2 and 12)
Days_Considered = 2;
Forecast_Inaccuracy = 0; %0; %1 if we ake into account a 25% max error in the forecasts, 0 if we consider perfect forecast
Safety_Power_for_Critical_Commitment = 3/4;  %this value corresponds to the percetnage of total batteries power of an export commitment that requires to have all batteries with energy.
%An example is : if we have 10 batteries at max power of 10kW, then if we commit for energy that corresponds to 100kW during 30 minutes (50kWh), then we need all batteries to have enough energy so they all contribute to the required power.
% It is important to consider that power as if we only consider energy, 8 batteries might have enough energy, but they cannot export such power --> the commitment will not be met. 
Safety_Coefficient = 1; %this represents the percentage of the big export that we can commit to without doubt.
stepRT = 2;  %Step for MPC change in recommendation
Type_of_Selling_Price = 1;  %1 for Market Day-ahead price, 2 for no selling tariff
Type_of_Buying_Price = 1; %1 for Agile tarif, 2 for ToU with 2 prices (Peak and low price), 3 for flat tariff
load('All_Households_Data.mat');
% We remove house 44 as the owner has been removed from the project trial
Porigin(:,44)=[];
Dorigin(:,44)=[];

%We add the VPP (PV and Wind)
%Scaling of Wind corresponding to the approximation given by the owner from Findhorn
for k=1:size(Wind,1)
   if Wind(k,1)>5
       Wind(k,1)=5;
   end
end
PV_VPP =  Porigin(:,size(Porigin,2)-14);

Porigin = [Porigin, PV_VPP*12, Wind/max(Wind)*max(PV_VPP)*18];
Dorigin = [Dorigin, zeros(size(Dorigin,1),2)];
Dorigin= Dorigin(:,23:24);
Porigin= Porigin(:,23:24);
epsilon = 0.1;


N_Bat = size(Dorigin,2);
Efficiency = 0.9;
EfficiencyC = Efficiency;
EfficiencyD = EfficiencyC;


ForecastTimeStep = 60*ForecastTimeStep; %30; %forecast for every 30 minutes
MarketTimeStep = 30; %Time Step in minutes %60*Timestep; %30; %market bid for every 30 minutes
Timestep = MarketTimeStep/60;

Window = 3*24;% Optimization for one day
Length = Window/Timestep; %number of data points

    NdataConsidered = 24/Timestep*Days_Considered; 
    Start =0*24/Timestep;
kWhCost = 200; %Battery cost /kWh
Capacity_Bat = 10300;
PowerRating_Bat = 5300;
MBCRange =   N_Bat*Capacity_Bat; 
Pbatmax = N_Bat*PowerRating_Bat*Timestep;% Limit of battery power
WHCRange = MBCRange'; % Range of Battery capacity in Watt-Hour (WH)

MaxBcap = MBCRange ;% Maximum battery capcity in WHr [Variable Parameter]
MinBcap = 0*0.001*MaxBcap;% Minimum battery capacity at 80% DoD
IBC = MaxBcap*0.51; % Initial Battery Capacity (IBC)of 80%MacBcap
RealTimeDataTimeStep = 1; %1 min data time step for RT


t = (Timestamp(1):minutes(RealTimeDataTimeStep):Timestamp(size(Timestamp,1)))';
t2 = (Timestamp(1):minutes(ForecastTimeStep):Timestamp(size(Timestamp,1)))';
t3 = (Timestamp(1):minutes(MarketTimeStep):Timestamp(size(Timestamp,1)))';  

%% Price Definition           
BPorigin = ones(size(Timestamp));
SPorigin = zeros(size(Timestamp));

j =1;

switch Type_of_Selling_Price 
    case 1 %Selling Price = Day ahead Auction Price from Nordpool
         for i = 1:size(Timestamp,1)
            SPorigin(i) = NordPoolSellingPrice(j);
           if  Timestamp(i)== t3(j+1)
               j=j+1;
           end

        end
%                     SP = interp1(t,SPorigin,t3);
        SP3 = interp1(t,SPorigin,t3);             

    case 2 %No selling Price
%                      SP = interp1(t,SPorigin,t3);
        SP3 = interp1(t,SPorigin,t3);                                 
end
j = 1;
switch Type_of_Buying_Price 
    case 1 %Buying Price = Agile ToU (Octopus)
         for i = 1:size(Timestamp,1)
            BPorigin(i) = BuyingPriceOctopus30min(j);
           if  Timestamp(i)== t3(j+1)
               j=j+1;
           end
        end
%                     BP = interp1(t,BPorigin,t3);
        BP3 = interp1(t,BPorigin,t3); 
    case 2
        HHP1 = 7.5;
        HHP2 = 22.5;
        T1 = t(1:HHP1*60+1);% Hours of  LHBP tariff 
        T2 = t(HHP1*60+2:HHP2*60+1);% Hours of HHBP tariff
        T3 = t(HHP2*60+2:1440);% Hours of LHBP tariff
        LHBP1 = Low_Buying_Price*ones(length(T1),1);% Low buying price for all the hours
        HHBP1 = Peak_Buying_Price*ones(length(T2),1);% High buying price for all the hours
        LHBP2 = Low_Buying_Price*ones(length(T3),1);% Low buying price for all the hours
        Daybuyprice = [LHBP1;HHBP1;LHBP2];
        BPorigin = repmat(Daybuyprice,30+1,1);
        BPorigin = [BPorigin; BPorigin(size(BPorigin,1))];
%                     BP = interp1(t,BPorigin,t3);
        BP3 = interp1(t,BPorigin,t3); 

    case 3 %Flat Buying Price
        BPorigin = Flat_Buying_Price*BPorigin;
%                     BP = interp1(t,BPorigin,t3);
        BP3 = interp1(t,BPorigin,t3); 
end


SP_Factor = 3;  %(SP = BP/SP_Factor) --> should be changed if we want to use elexon prices for example, in which case SP = import from elexon data

BPorigin = zeros(size(Timestamp));
SPorigin = zeros(size(Timestamp));

j =1;
for ki = 1:size(Timestamp,1)
    BPorigin(ki) = BuyingPriceOctopus30min(j);
    SPorigin(ki) = BuyingPriceOctopus30min(j)/SP_Factor;
   if  Timestamp(ki)== t3(j+1)
       j=j+1;
   end

end




SumDorigin = sum(Dorigin,2);
SumPorigin = sum(Porigin,2);

Ndata =  size(Porigin,1);
j = 1;
Pforecast = zeros(size(t2));
Dforecast = zeros(size(t2));


for i =1:ForecastTimeStep/RealTimeDataTimeStep:size(Porigin,1)-1
    Pforecast(j)= sum(SumPorigin(i:i+ForecastTimeStep/RealTimeDataTimeStep-1))*(1+Forecast_Inaccuracy*(0.25-rand(1)/2));
    Dforecast(j)= sum(SumDorigin(i:i+ForecastTimeStep/RealTimeDataTimeStep-1))*(1+Forecast_Inaccuracy*(0.25-rand(1)/2));
    j = j+1;
end

Pbid = zeros(size(t3));
Dbid = zeros(size(t3));

j = 1;
Pbid(1) = Pforecast(j);
Dbid(1) = Dforecast(j);
for i = 2:size(Pbid,1)
    if i== (j)*ForecastTimeStep/MarketTimeStep+1
        Pbid(i) = Pforecast(j+1)*MarketTimeStep/ForecastTimeStep;
        Dbid(i) = Dforecast(j+1)*MarketTimeStep/ForecastTimeStep;
        j = j+1;
    else
     Pbid(i) =  Pbid(i-1);
     Dbid(i) =  Dbid(i-1);
    end     
end

Pref = Pbid;
Dref = Dbid;
t2 = t2(Start*MarketTimeStep/ForecastTimeStep+1:Start*MarketTimeStep/ForecastTimeStep+NdataConsidered*MarketTimeStep/ForecastTimeStep+1);
t3 = t3(Start+1:Start+NdataConsidered+1);
Pforecast = Pforecast(Start*MarketTimeStep/ForecastTimeStep+1:Start*MarketTimeStep/ForecastTimeStep+NdataConsidered*MarketTimeStep/ForecastTimeStep+1);
Dforecast = Dforecast(Start*MarketTimeStep/ForecastTimeStep+1:Start*MarketTimeStep/ForecastTimeStep+NdataConsidered*MarketTimeStep/ForecastTimeStep+1);
Pbid = Pbid(Start+1:Start+NdataConsidered+1);
Dbid= Dbid(Start+1:Start+NdataConsidered+1);
P = Pbid;
D = Dbid;

BP = interp1(t,BPorigin,t3);
SP = interp1(t,SPorigin,t3);
BP3 = interp1(t,BPorigin,t3);
SP3 = interp1(t,SPorigin,t3);

SumPorigin = SumPorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);
Porigin = Porigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1,:);
Dorigin = Dorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1,:);
SumDorigin= SumDorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);
t = t(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);
BPorigin = BPorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);
SPorigin = SPorigin(Start*MarketTimeStep/RealTimeDataTimeStep+1:Start*MarketTimeStep/RealTimeDataTimeStep+NdataConsidered*MarketTimeStep/RealTimeDataTimeStep+1);

Porigininit = Porigin;
Dorigininit = Dorigin;

Pmax = 3*max(max(P),max(D));

    SoC = zeros(length(P),1);
    Energysold = zeros(length(P),1);
    Energybought = zeros(length(P),1);
    FromGrid = zeros(length(P),1);
    ToGrid = zeros(length(P),1);
    SoCmax = MBCRange ; %kWh 
    SoCmin = MinBcap;
    SoCInit = IBC; %kWh
    Powerbat = zeros(length(P),1);
    ellapsedtime = zeros(100,1);
    indice = 0;
    for i = 1:Length:length(P)
            indice = indice +1;
            tic
            LengthOptim = min(i+Length,length(D))-i+1 ; % Because the size of P might not be a multiple of the Length used for the optimization
            
%             t=cputime;
            Demand = D(i:i+LengthOptim-1);
            Production = P(i:i+LengthOptim-1);
            BPrice = BP(i:i+LengthOptim-1);
            SPrice = SP(i:i+LengthOptim-1);
            f = [BPrice;-SPrice;zeros(LengthOptim,1); zeros(LengthOptim,1);zeros(LengthOptim,1);zeros(LengthOptim,1);zeros(LengthOptim,1); zeros(LengthOptim,1) ];
            
            Index_binary = size(f,1)-2*LengthOptim+1:size(f,1); %[121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144];
            Aeq = [eye(LengthOptim) zeros(LengthOptim) EfficiencyD*eye(LengthOptim) -eye(LengthOptim)/EfficiencyC -eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %Pgrin in = Demand - Production + Pbat_recharge/efficiency -Pbat_Discharge*Efficiency +X1   || X1 is here to make t possible to have Pgrid IN !=0 when Pgrid out =0
                   zeros(LengthOptim) eye(LengthOptim) -EfficiencyD*eye(LengthOptim) eye(LengthOptim)/EfficiencyC zeros(LengthOptim) -eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim)]; %Pgrin out = - Demand + Production - Pbat_recharge/efficiency + Pbat_Discharge*Efficiency +X2 || X2 is here to make t possible to have Pgrid OUT !=0 when Pgrid in =0
            beq = [ Demand - Production ; Production - Demand ];

            A = [eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) -Pmax*eye(LengthOptim) zeros(LengthOptim); %makes sure Pgrid in = 0 when Pgrid out != 0
                 zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) Pmax*eye(LengthOptim) zeros(LengthOptim);  %makes sure Pgrid out = 0 when Pgrid in != 0
                 zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) eye(LengthOptim) -Pmax*eye(LengthOptim) zeros(LengthOptim); %makes sure X2  = 0 when X1 != 0
                 zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim) Pmax*eye(LengthOptim) zeros(LengthOptim); %makes sure X1  = 0 when X2 != 0
                 zeros(LengthOptim) zeros(LengthOptim) -tril(ones(LengthOptim)) tril(ones(LengthOptim)) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %makes sure SoC <= SoC max
                 zeros(LengthOptim) zeros(LengthOptim) tril(ones(LengthOptim)) -tril(ones(LengthOptim)) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %makes sure SoC >= SoCmin
                 zeros(LengthOptim) zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) -Pbatmax*eye(LengthOptim) ; %makes sure Pbatdischarge  = 0 when Pcharge != 0
                 zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) Pbatmax*eye(LengthOptim)]; %makes sure Pbat_charge  = 0 when Pdischarge != 0

            b = [zeros(LengthOptim,1); Pmax*ones(LengthOptim,1); zeros(LengthOptim,1); Pmax*ones(LengthOptim,1); (SoCmax-SoCInit)*ones(LengthOptim,1); (SoCInit - SoCmin)*ones(LengthOptim,1); zeros(LengthOptim,1); Pbatmax*ones(LengthOptim,1)]; 
            lb = [zeros(LengthOptim,1); zeros(LengthOptim,1); zeros(LengthOptim,1); zeros(LengthOptim,1); zeros(LengthOptim,1);  zeros(LengthOptim,1);  zeros(LengthOptim,1); zeros(LengthOptim,1)]  ;
            ub = [Pmax*ones(LengthOptim,1); Pmax*ones(LengthOptim,1); Pbatmax*ones(LengthOptim,1);  Pbatmax*ones(LengthOptim,1); Pmax*ones(LengthOptim,1);  Pmax*ones(LengthOptim,1);  ones(LengthOptim,1);  ones(LengthOptim,1) ]  ;

            %if we wanted to change the battery power for discharge for the night (to
            %follow a simulation where price was very high for selling and see if the simulation is able to find a good behaviour)
%             ub(LengthOptim*2+21:LengthOptim*2+24)=3;
            % ub(LengthOptim*3+21:LengthOptim*3+24)=3;
%             x0=[];
%              options=optimoptions('intlinprog','Display','off');
%             options.Display = 'off';
            [x,fval1,exitflag1,output1]  = intlinprog(f, Index_binary, A,b, Aeq, beq, lb, ub);
            elapsedTime = toc ;
%             ellapsedtime(indice) = cputime - t;
            ellapsedtime(indice) = elapsedTime;
            
            SoC(i:i+LengthOptim-1) = -tril(ones(LengthOptim))*x(LengthOptim*2+1:LengthOptim*3)+tril(ones(LengthOptim))*x(LengthOptim*3+1:LengthOptim*4)+SoCInit*ones(LengthOptim,1);
            PgridIN = x(1:LengthOptim);
            PgridOUT = x(LengthOptim+1:LengthOptim*2);
            Pbat = x(LengthOptim*2+1:LengthOptim*3)-x(LengthOptim*3+1:LengthOptim*4);
            X1 = x(LengthOptim*4+1:LengthOptim*5);
            X2 = x(LengthOptim*5+1:LengthOptim*6);
            Alpha = x(LengthOptim*6+1:LengthOptim*7);
            Beta = x(LengthOptim*7+1:LengthOptim*8);
            
            SoCInit = SoC(i+LengthOptim-1);
            for j = 0:LengthOptim-1

            Energybought(i+j) = PgridIN(1+j)*BPrice(j+1);
            Energysold(i+j) =PgridOUT(1+j)*SPrice(j+1);          
            FromGrid(i+j) = max(0,(D(i+j)-P(i+j)-max(0,Pbat(j+1))*EfficiencyD+max(0,-Pbat(j+1)/EfficiencyC)));
            ToGrid(i+j) = max(0, (P(i+j)-D(i+j)+max(0,Pbat(j+1)*EfficiencyD)-max(0,-Pbat(j+1)/EfficiencyC)));
            end
            Powerbat (i:i+LengthOptim-1) = Pbat;

    end

  

   
%% Market Bid

% Determine which mode we are on at the start(0.5 for charging at power < difference between P and D,
%1 for charging with Power = the %difference between production and demand, 2 for charging with other values or -1 for discharging mode with Power = the 
%difference between production and demand, -2 for discharging with other values,- 0.5 for discharging at power < difference between P and D,)
%and 0 if not changing.
% we put specification values back to normal
% MBCRange =   6300; % 13500; %Range of different battery capacity
% Pbatmax = N_Bat*3300*MarketTimeStep/60;% max battery Power = 3.3kW% 3*max(P); %W  We consider there is no limit of battery power
    SoCmax = MBCRange ; %kWh 
    SoCmin = MinBcap;
Pbatmax = N_Bat*PowerRating_Bat*Timestep;% 2*3300*Timestep max battery Power = 3.3kW% 3*max(P); %W  We consider there is no limit of battery power

    
      %we will now determine which mode the battery is in.
    % We list 4 modes:
    %1. Battery follows P-D
    %    1.a: SoC increases (mode 11)
    %    1.b: SoC decreases (mode 12)
    % 2. Battery must ensure a SoC at the end of the time interval
    %    2.a: SoC must be equal or higher than the one obtained during the (mode 21)
    %    optimization at the end of the considered time interval --> use
    %    the algo to reach the SoC
    %    2.b: SoC must be equal or lower than the one obtained during the  (mode 22)
    %    optimization at the end of the considered time interval --> use
    %    the algo to reach the SoC
    
    
    % Note that sometimes, the optimization does not follow a logic choice
    % and stops recharging to recharge later as it gives the same result in
    % terms of revenue. The algorithm does identify this
    
    
   %Determine the time intervals (everytime there is a new price or that the battery does not  
  % first, we determine the SoC targets (every time a SoC evolution changes)  
   k=0;
InitialMode =0;
while InitialMode==0
    k=k+1;
    if SoC(k) - SoC(k+1) <-epsilon
        InitialMode =1; %1 for Charging Mode 

    elseif SoC(k) - SoC(k+1) > epsilon
       InitialMode =-1; %-1 for discharging Mode 

    end
end


    j=1;
SoCTargets =  zeros(size(t3));
SoCTargetsType =  zeros(size(t3)); %0 if target does not correspond to a change of price, 1 otherwise
PreviousMode = InitialMode;
last_index = 1;
ChargeMode = zeros(size(P));
ChargeMode(1) = InitialMode;
CurrentMode = InitialMode;
for k = 1:1:size(P,1)-1
    if SoC(k) - SoC(k+1) <-epsilon
        CurrentMode =1; %1 for Charging Mode 

    elseif SoC(k) - SoC(k+1) > epsilon
       CurrentMode =-1; %-1 for discharging Mode 

    end
    if CurrentMode*PreviousMode <0  || abs(BP(k)-BP(k+1))>0 ||  abs(SP(k+1)-SP(k))>0 %5|| (t3(k).Hour== t3(startfrr).Hour) || t3(k).Hour == t3(endfrr).Hour
       SoCTargets(last_index:k)=SoC(k);
        last_index = k+1;
    end
    PreviousMode = CurrentMode;
    ChargeMode(k+1)=PreviousMode;
     
end
    
    
 



% Determine the times of change of time interval:
Counter = 0; %counter that indicates the number of time interval 
for k = 1:1:size(P,1)-1
    
    if abs(SoCTargets(k)- SoCTargets(k+1))>epsilon || abs(BP(k)-BP(k+1))>0 ||  abs(SP(k+1)-SP(k))>0 %|| (t3(k).Hour== t3(startfrr).Hour) || t3(k).Hour == t3(endfrr).Hour
       Counter = Counter+1;
    end
end

TimeInterval = datetime(zeros(Counter+2,1),0,0,0,0,0);
IndexTimeInterval = ones(Counter+2,1);
j = 1;
TimeInterval(j)=t3(1);
IndexTimeInterval(j)=1;
j = j+1;
for k = 2:1:size(P,1)-1
    if abs(SoCTargets(k)- SoCTargets(k-1))>epsilon || abs(BP(k-1)-BP(k))>0 ||  abs(SP(k-1)-SP(k))>0 %|| (t3(k).Hour== t3(startfrr).Hour) || t3(k).Hour == t3(endfrr).Hour
       TimeInterval(j)=t3(k);
       IndexTimeInterval(j)=k;
       j = j+1;
    elseif abs(SoCTargets(k)- SoCTargets(k-1))<epsilon/10  && k>3%if there was no change in the target during this timeslot, we assume that the target is to be met for the next timeslot
        TimeInterval(j-1)=t3(k);
        IndexTimeInterval(j-1)=k;
    end
end

TimeInterval(j)=t3(k+1);
IndexTimeInterval(j)=k+1;
% Determine the Mode within each time interval
j = 1;
CurrentModeVector = zeros(size(P));
for k = 1:size(P,1)-1
    if t3(k) > TimeInterval(j) || k==1
        %Notice that P and D are energy vectors. not power... Hence we do
        %not need to include the time step. However, Efficiency must be
        %added.

        if (ChargeMode(k) ==1  && abs(abs(sum(P(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))-sum(D(IndexTimeInterval(j)+1:IndexTimeInterval(j+1))))*Efficiency - ...
                   abs(SoCTargets(IndexTimeInterval(j+1))-SoCTargets(IndexTimeInterval(j))))<= epsilon) || (ChargeMode(k) ==1 && abs(SoCTargets(IndexTimeInterval(j+1)) - SoCmax)<epsilon...
                   &&  sum(FromGrid(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))<epsilon    )% || ( ChargeMode(k) ==1 && t3(k).Hour>= t3(startfrr).Hour && t3(k).Hour <= t3(endfrr).Hour)%abs(SoCTargets(IndexTimeInterval(j+1)-1) - SoCTargets(IndexTimeInterval(j+1)) )<epsilon) %charging mode && charge by P -D (mode 11) or if the battery was charged by P until it was full
        CurrentModeVector(k) = 11;
        elseif (ChargeMode(k) ==-1 && abs(abs(sum(P(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))-sum(D(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))/Efficiency - ...
                   (SoCTargets(IndexTimeInterval(j+1))-SoCTargets(IndexTimeInterval(j)))))<= epsilon) || (ChargeMode(k) ==-1 &&  abs(SoCTargets(IndexTimeInterval(j+1)) - SoCmin)<epsilon...
                          &&  sum(ToGrid(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))<epsilon ) %|| (ChargeMode(k) ==-1 && t3(k).Hour>= t3(startfrr).Hour && t3(k).Hour <= t3(endfrr).Hour) %abs(SoCTargets(IndexTimeInterval(j+1)-1) - SoCTargets(IndexTimeInterval(j+1)) )<epsilon) %discharge mode
        CurrentModeVector(k) = 12;    
        elseif ChargeMode(k) ==1  && abs(abs(sum(P(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))-sum(D(IndexTimeInterval(j)+1:IndexTimeInterval(j+1))))*Efficiency - ...
                   abs(SoCTargets(IndexTimeInterval(j+1))-SoCTargets(IndexTimeInterval(j))))> epsilon  %charging mode && charge by P -D (mode 11)
        CurrentModeVector(k) = 21;
        elseif ChargeMode(k) ==-1 && abs(abs(sum(P(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))-sum(D(IndexTimeInterval(j)+1:IndexTimeInterval(j+1)))/Efficiency - ...
                   (SoCTargets(IndexTimeInterval(j+1))-SoCTargets(IndexTimeInterval(j)))))> epsilon  %discharge mode
        CurrentModeVector(k) = 22;    
        end
        j = j+1;        
    else
       CurrentModeVector(k) = CurrentModeVector(k-1); 
    end
end

% We now compute the bid curves
% t3, Pbid, Dbid are the vectors to be considered

%We determine SoCTargetBids, same as SoCTarget, but with the time scale of
%the Market
SoCTargetBid = zeros(size(t3));
ChargeModeBid = zeros(size(t3));
CurrentModeVectorbid = zeros(size(t3));

SoCTargetBid(1) = SoCTargets(1);
ChargeModeBid(1) = ChargeMode(1);
CurrentModeVectorbid(1) = CurrentModeVector(1);
LengthTimeIntervalBid = zeros(size(t3));
j=1;
previousindex = 1;
for k = 1:size(t3,1)-1
if k>=3*48+11.5*2
    tempre=0;
end
    if t3(k) == t3(IndexTimeInterval(j))
        SoCTargetBid(k) = SoCTargets(IndexTimeInterval(j+1));
        ChargeModeBid(k) = ChargeMode(IndexTimeInterval(j+1)); 
        LengthTimeIntervalBid(previousindex:k-1)=k-previousindex;
        CurrentModeVectorbid(k) = CurrentModeVector(IndexTimeInterval(j+1));
        previousindex = k;
        j = j+1;
    else
        SoCTargetBid(k) = SoCTargetBid(k-1);
        ChargeModeBid(k) = ChargeModeBid(k-1);
        CurrentModeVectorbid(k) = CurrentModeVectorbid(k-1);
    end
end
SoCTargetBid(k+1)= SoCTargets(size(SoCTargets,1));


%#################################################################################################
% ##################   when more than
% one battery
SoCTargetBid = [SoCTargetBid(1);SoCTargetBid];

LengthTimeIntervalBid(previousindex:k+1) = k-previousindex;
%We now compute the forecast of the SoC for the battery 
SoCbid = zeros(size(t3));
Ebatbid = zeros(size(t3));
j=1;

SoCbid(1) = SoC(1);
%Initialization
       RealisedOutputprevious = 0;   
SoCprevious = SoC(1);
timeinterval = LengthTimeIntervalBid(1);
LastSoCtarget = SoCprevious;
NextSoCtarget = SoCTargetBid(1);
        if CurrentModeVectorbid(1)==21 || CurrentModeVectorbid(1)==22
            soctargettype = 1; %0 if SoC target is not mandatory, 1 if it is madatory
        else 
            soctargettype = 0;
        end
          StartTime = 0; % time (indice) at which we start this market timeslot 
        Et=0;  %energy that is exported from the battery at the time t (j)
        slopeprevious = 0; %slope of the battery power at the previous time step k
       time1previous = 0;  %time1 determines the the time at which we need to use the battery at its maximum power in order to meet the target
       RealisedOutput = 0;  
       NeededEnergy = SoCTargetBid(1)  - SoCprevious; 
CurrentModeVectorbid = [CurrentModeVectorbid(1);CurrentModeVectorbid];
CurrentModeVector = [CurrentModeVector;CurrentModeVector(size(CurrentModeVector,1))];
IndexTimeInterval = [IndexTimeInterval; IndexTimeInterval(size(IndexTimeInterval,1))];
for k = 1:size(t3,1)
    
    if t3(k) == t3(IndexTimeInterval(j)) %if we are in a new timestep
        StartTime = k-1; % time (indice) at which we start this market timeslot 
        Et=0;  %energy that is exported from the battery at the time t (j)
        slopeprevious = 0; %slope of the battery power at the previous time step k
       time1previous = 0;  %time1 determines the the time at which we need to use the battery at its maximum power in order to meet the target
       RealisedOutput = 0;
       RealisedOutputprevious = 0;        
        timeinterval = LengthTimeIntervalBid(k);
        LastSoCtarget = SoCprevious;
        NextSoCtarget = SoCTargetBid(k);
         NeededEnergy = NextSoCtarget - LastSoCtarget; %SoCTargetBid(k)  - SoCTargetBid(k-1); 
     
        if CurrentModeVector(IndexTimeInterval(j))==21 || CurrentModeVector(IndexTimeInterval(j))==22% ||  (t3(k).Hour>= t3(startfrr).Hour && t3(k).Hour <= t3(endfrr).Hour)
           soctargettype = 1; %0 if SoC target is not mandatory, 1 if it is madatory
        else 
            soctargettype = 0;
        end
        j = j+1;

    end
    CurrentTime = k-1-StartTime; %the time from which we recompute the energy that we still need
     time1 = 2*(sign(NeededEnergy)*Pbatmax*Efficiency*(timeinterval-CurrentTime/2+1/2)-Et*(CurrentTime/2+1/2)+RealisedOutput-sign(NeededEnergy)*min(abs(NeededEnergy),Pbatmax*Efficiency*timeinterval))/(sign(NeededEnergy)*Pbatmax*Efficiency-Et+0.000000001); %determines the the time at which we need to use the battery at its maximum power
    if time1 > timeinterval %if we do not need to be at full power within this time interval
        slope = 2*(sign(NeededEnergy)*min(abs(NeededEnergy),Pbatmax*Efficiency*timeinterval)-RealisedOutput - Et*(timeinterval-CurrentTime))/((timeinterval-CurrentTime)*(timeinterval-CurrentTime+1));
    else
        slope = (sign(NeededEnergy)*Pbatmax*Efficiency-Et)^2/(sign(NeededEnergy)*Pbatmax*Efficiency*(2*timeinterval+1-2*CurrentTime)-Et+2*(RealisedOutput-sign(NeededEnergy)*min(abs(NeededEnergy),Pbatmax*Efficiency*timeinterval))+0.000001);
    end 
    if abs(sign(NeededEnergy)*Pbatmax*Efficiency-Et ) < epsilon  %if we already reached the maximum power at previous step while we were still ramping (happens when time1 is an integer
        time1 = time1previous; %we keep the same value as previous iteration, as it will still be the maximum power for the battery
        slope=slopeprevious;
    end

    if k >=9.5*2
        tempb=0;
    end
    Ebatbid(k) =slope+Et;% Ebatbid(k)+slope; %The battery power is equal to what should be produced to reach the bid minus what is already produced by the PV
    Et = Ebatbid(k);
    s = sign(Pbid(k)-Dbid(k));
    if sign(NeededEnergy)*s>0
       Ebatbid(k) = sign(Et)*max(abs(Ebatbid(k)),abs( Pbid(k)-Dbid(k))*((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency) );
    end
    Ebatbid(k)= min(Pbatmax, max(-Pbatmax,Ebatbid(k)));
    if soctargettype ==0
        Ebatbid(k) = max(-Pbatmax,min(Pbatmax,(Pbid(k)-Dbid(k))*((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency))) ;
        
    end
    
    if k-StartTime - time1 >= -epsilon  && time1 >0 % if we have already reached the time at which the battery should provide all its energy, then we must set its power to Pmax
       Ebatbid(k)=sign(NeededEnergy)*Pbatmax;
    end   
    SoCbid(k) = SoCprevious+Ebatbid(k);
    [Ebatbid(k),SoCbid(k)] = coercsoc(SoCbid(k),SoCmax,SoCmin,Ebatbid(k),1,LastSoCtarget, NextSoCtarget,soctargettype,Pbatmax);  %1 because we consider energy and not power.
    SoCbid(k) = SoCprevious+Ebatbid(k);
    RealisedOutput = RealisedOutputprevious + Ebatbid(k);
    
%     if RealisedOutput > NeededEnergy %if we already exported too much energy compared to our bid, we need to reduce the amount of eergy exported
    if sign( NeededEnergy - RealisedOutput)*sign(NeededEnergy) <0  %if there is more energy provided than required %if we already exported too much energy compared to our bid, we need to reduce the amount of eergy exported
       Ebatbid(k) = Ebatbid(k) + NeededEnergy - RealisedOutput ;%max(-Pbatmax, min(Pbatmax,NeededEnergy-RealisedOutput)); %so we charge the battery with the surplus of energy. making sure the battery's max power is below Pmax
    end    
    SoCbid(k) = SoCprevious+Ebatbid(k);
    RealisedOutput = RealisedOutputprevious + Ebatbid(k);
    RealisedOutputprevious = RealisedOutput;
    time1previous = time1;
    slopeprevious = slope;    
    SoCprevious = SoCbid(k);
end

Extra_Power_Needed = zeros(size(t3));
To_Grid_bid = zeros(size(t3));
From_Grid_bid = zeros(size(t3));

for k = 1:size(t3,1)

    s = -sign(Ebatbid(k));
    a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);
    
    
    From_Grid_bid(k) = max(0,Dbid(k)-Pbid(k)+ Ebatbid(k)*a);
    To_Grid_bid(k) = max(0,Pbid(k)-Dbid(k)- Ebatbid(k)*a);
%Now we detect when there is a need for big power
%(= when the requested power = MaxPower*Percentage_Power_for_Critical_Commitment)
 
    if -Ebatbid(k)>= Pbatmax*Safety_Power_for_Critical_Commitment
        Extra_Power_Needed(k)=1;
        To_Grid_bid(k) = To_Grid_bid(k)*Safety_Coefficient;
    end       
    
end






%% Market bids computation



% We now compute the real Battery operation curves
% t, Porigin, Dorigin are the vectors to be considered


Pbatmax = PowerRating_Bat*RealTimeDataTimeStep/60*1;% max battery Power = 3.3kW% 3*max(P); %W  We consider there is no limit of battery power
SoCmax= SoCmax/N_Bat;
SoCmin = SoCmin/N_Bat;

%We determine SoCTargetreals, same as SoCTarget, but with the time scale of
%the Market
SoCTargetBid = SoCTargetBid/N_Bat;  %We determine the SoC target for all batteries (assumed identical here)
SoCbidscaled = SoCbid/N_Bat;
SoCTargetreal = zeros(size(t));
ChargeModereal = zeros(size(t));
CurrentModeVectorreal = zeros(size(t));
SoCTargets = SoCTargets/N_Bat;
SoCTargetreal(1) = SoCTargetBid(1);
ChargeModereal(1) = ChargeMode(1);
CurrentModeVectorreal(1) = CurrentModeVector(1);
LengthTimeIntervalreal = zeros(size(t));
j=1;
previousindex = 1;
SoCTargetrealprevious = SoCTargetreal(1);
ChargeModerealprevious = ChargeModereal(1);
CurrentModeVectorrealprevious = CurrentModeVectorreal(1);

for k = 1:size(t,1)-1
    if k >=24*60+11.5*60%8*60+1
        tempb=0;
    end
    if t(k) == t3(IndexTimeInterval(j))
%         SoCTargetreal(k) = SoCTargetBid(IndexTimeInterval(j+1));
        SoCTargetreal(k) = SoCbidscaled(IndexTimeInterval(j+1));
        ChargeModereal(k) = ChargeMode(IndexTimeInterval(j+1)); 
        LengthTimeIntervalreal(previousindex:k-1)=k-previousindex;
        CurrentModeVectorreal(k) = CurrentModeVector(IndexTimeInterval(j+1));
        previousindex = k;
        SoCTargetrealprevious = SoCTargetreal(k);
        ChargeModerealprevious = ChargeModereal(k);
        CurrentModeVectorrealprevious = CurrentModeVectorreal(k);
        j = j+1;
    else
        SoCTargetreal(k) = SoCTargetrealprevious; 
        ChargeModereal(k) = ChargeModerealprevious; %ChargeModereal(k-1);
        CurrentModeVectorreal(k) = CurrentModeVectorrealprevious;%CurrentModeVectorreal(k-1);
    end
end
SoCTargetreal(k+1)= SoCTargetBid(size(SoCTargetBid,1));

LengthTimeIntervalreal(previousindex:k+1) = k-previousindex;
%We now compute the forecast of the SoC for the battery 
SoCreal = zeros(size(t,1),N_Bat);
Ebatreal = zeros(size(t,1),N_Bat);


%% Real Time Operation

%Initialization for the SoC timeframe
RealisedOutputprevious = zeros(1,N_Bat);   
SoCprevious = SoC(1)*ones(1,N_Bat)/N_Bat;
timeinterval = LengthTimeIntervalreal(1);
LastSoCtarget = SoCprevious;
NextSoCtarget = SoCTargetreal(1);
if CurrentModeVectorreal(1)==21 || CurrentModeVectorreal(1)==22
    soctargettype = 1; %0 if SoC target is not mandatory, 1 if it is madatory
else 
    soctargettype = 0;
end
StartTime = 0; % time (indice) at which we start this market timeslot 
Et=0;  %energy that is exported from the battery at the time t (j)
slopeprevious = zeros(1,N_Bat); %slope of the battery power at the previous time step k
time1previous = zeros(1,N_Bat);  %time1 determines the the time at which we need to use the battery at its maximum power in order to meet the target
RealisedOutput = zeros(1,N_Bat);  
NeededEnergy = SoCTargetBid(1)*ones(1,N_Bat)  - SoCprevious; 
Required_Export_bid = max(0,To_Grid_bid-From_Grid_bid);



save('Target_MPC_2Batteries2a.mat','Required_Export_bid','CurrentModeVectorreal','SoCTargetreal');
save('TestMPC_Data2a.mat','Porigin','Dorigin','N_Bat','Days_Considered','Timestep','RealTimeDataTimeStep','PowerRating_Bat','epsilon','Timestep','IBC','IBC',...
    'MBCRange','MinBcap','MarketTimeStep','t','t2','t3','BPorigin','SPorigin','Pbatmax','Efficiency','stepRT');
%% Real Time MPC
load('Target_MPC_2Batteries2a.mat');
load('TestMPC_Data2a.mat');

Window = 30;% Optimization for 30 min
Length = 30;
NdataConsidered = 24/Timestep*Days_Considered; %round((size(Timestamp,1)-1)/60/24); %computation for the whole week
% Length = Window/Timestep; %number of data points

%%optimization constraints
Pmax = 3*max(max(Porigin),max(Dorigin));
% BPrice = BPorigin(1:min(Length,size(BP,1))); %0.1* [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1.9 1.9 1.9 1.9 1 1 1 1 1]'; %0.9*ones(size(Demand));
% SPrice = SPorigin(1:min(Length,size(BP,1))); %0.1*[1 1 1 1 1 1 1 1 1 0.1 0.1 0.1 1 1 1 1 1 1 1 1 1 10000.9 10000.9 1.0]'; %0.01*ones(size(Demand));
P = Porigin;
D = Dorigin;
 D(1:30,:)=[]; 
 D = [D;zeros(30,N_Bat)];
 P(1:30,:)=[]; 
 P = [P;zeros(30,N_Bat)];
 
    SoC = zeros(length(P),N_Bat);
%     Energysold = zeros(length(P),1);
%     Energybought = zeros(length(P),1);
%     FromGrid = zeros(length(P),1);
%     ToGrid = zeros(length(P),1);
%     SoC(1) = IBC;
    SoCmax = MBCRange/N_Bat ; %kWh 
    SoCmin = MinBcap/N_Bat;
    Powerbat = zeros(length(P),N_Bat);
    Pbatmax = PowerRating_Bat*RealTimeDataTimeStep/60*1;% max battery Power = 3.3kW% 3*max(P); %W  We consider there is no limit of battery power
    Energysold = zeros(length(P),1);
    Energybought = zeros(length(P),1);
    FromGrid = zeros(length(P),1);
    ToGrid = zeros(length(P),1);
    ComputingTime = zeros(length(P),N_Bat);
    EfficiencyD=Efficiency;
    EfficiencyC=Efficiency;
    indiceExportTarget = 1;
    SoCinit = IBC*ones(N_Bat,1)/N_Bat; %kWh
Export_Target = zeros(N_Bat,1);
for i = 1:30:size(Porigin,1)-31
SoCTarget = SoCTargetreal(i+1);
% ##############  Little trick to be changed ##############################
%############  if the SoCTarget is the same in the next timeinterval, we
%do not want the battery to reach that target in this time interval, but we
%want it to take 2 time interval to meet the SoC target. We should actually
%optimize not  through only one interval, but through all the interval
%where the SoC target is unchanged... but in flat price tariff, that would
%lead to very big interval, so in this trick, we prefer to set
%sub-SoCTarget equal to half the SoC interval needed

Numb_Interval = 1;
while  abs(SoCTargetreal(i+1+Numb_Interval*30)-SoCTarget)<epsilon
   Numb_Interval=Numb_Interval+1;
end
if Numb_Interval>1
    SoCTarget = SoCTarget-(SoCTarget-mean(SoCinit))/Numb_Interval; %##### mean(SoCinit)  To be changed (just do not want to waste time do it right now
end


Already_Realised_Export_Target=zeros(N_Bat,1);
Export_Target_ref = max(0,Required_Export_bid(indiceExportTarget+1)-sum(Already_Realised_Export_Target));

if CurrentModeVectorreal(i+1)==21 || CurrentModeVectorreal(i+1)==22 
    soctargettype = ones(N_Bat,1); %0 if SoC target is not mandatory, 1 if it is madatory
else 
    soctargettype = zeros(N_Bat,1);
end




for k = 0:stepRT:30-1
    weight_Split_Export_Effort = zeros(N_Bat,1);
        for d = 1:N_Bat
            if Already_Realised_Export_Target(d)>=0 || SoCinit(d)-SoCTarget>0%if the battery is not importing electricity or it has too much energy
            weight_Split_Export_Effort(d) = SoCinit(d); %We assume the battery that has the most energy takes the most of the effort to achieve the export target
            end
        end
        weight_Split_Export_Effort=weight_Split_Export_Effort/(0.000001+sum(weight_Split_Export_Effort));
    weight_Split_Import_Effort = zeros(N_Bat,1);
        for d = 1:N_Bat
            weight_Split_Import_Effort(d) = SoCmax-SoCinit(d); %We assume the battery that has the most energy left takes the most of the effort to achieve the export target
        end
        weight_Split_Import_Effort=weight_Split_Import_Effort/sum(weight_Split_Import_Effort);        
    for d = 1:N_Bat
       Export_Target(d) = (Export_Target_ref - sum(Already_Realised_Export_Target));%*weight_Split_Export_Effort(d);
       if Export_Target(d)>=0
            Export_Target(d) = Export_Target(d) *weight_Split_Export_Effort(d);
        else
            Export_Target(d) = Export_Target(d) *weight_Split_Import_Effort(d);
       end
    end 
    
    
    
    
    for d = 1:N_Bat
        LengthOptim = min(i+Length,length(D))-i-(k);  % Because the size of P might not be a multiple of the Length used for the optimization
        if LengthOptim ==0
            LengthOptim=1;
        end
        if i ==1
            Demand_Previous = D(i+k,d);
            Production_Previous = P(i+k,d);
        else
            Demand_Previous = D(i+k-1,d);
            Production_Previous = P(i+k-1,d);            
        end
      
        if LengthOptim > MarketTimeStep/2
            Demand = [Demand_Previous; mean(D(i+k:i+MarketTimeStep/2,d))*ones(LengthOptim-MarketTimeStep/2-1,1); mean(D(i+MarketTimeStep/2+1:i+MarketTimeStep-1,d))*ones(MarketTimeStep/2,1)];
            Production = [Production_Previous; mean(P(i+k:i+MarketTimeStep/2,d))*ones(LengthOptim-MarketTimeStep/2-1,1); mean(P(i+MarketTimeStep/2+1:i+MarketTimeStep-1,d))*ones(MarketTimeStep/2,1)];
        else 
            Demand = [Demand_Previous; mean(D(i+MarketTimeStep/2:i+MarketTimeStep-1,d))*ones(LengthOptim-1,1)];
            Production = [Production_Previous; mean(P(i+MarketTimeStep/2:i+MarketTimeStep-1,d))*ones(LengthOptim-1,1)];

        end
                BPrice = BPorigin(i+k:i+k+LengthOptim-1);
                SPrice = SPorigin(i+k:i+k+LengthOptim-1);

    Export_Flexibility_Down = 0;
    Export_Flexibility_Up = 0;
    SoCTarget_Flexibilityadd = 0;
        if Export_Target(d)<1
            weightMarket = 0;
            SoCTarget_Flexibility =0;
            if sum(Production)-sum(Demand)> max(0,(SoCTarget-SoCinit(d)))
                SoCTarget_Flexibilityadd = -( sum(Production)-sum(Demand)-(SoCTarget-SoCinit(d)));
            end

        else 
            SoCTarget_Flexibility =0;
             weightMarket = 1;  
             soctargettype(d)= 0; %zeros(N_Bat,1);
             if (SoCinit(d)-SoCmin)*Efficiency+sum(Production)-sum(Demand)-Export_Target(d)<50  %######### Safety cofficient
                 Export_Target(d) =max(0,(SoCinit(d)-SoCmin)*Efficiency- (sum(Demand)-sum(Production)))*0.9; %safety coefficient
                 if Export_Target(d)==0
                     weightMarket = 0;
                    SoCTarget_Flexibility =0;                
                 end
             elseif (SoCmax-SoCinit(d))<sum(Production)-sum(Demand)-Export_Target(d)
                 Export_Flexibility_Up = (sum(Production)-sum(Demand)-Export_Target(d)-(SoCmax-SoCinit(d)))*1.2;
             end
        end

    if SoCmax-SoCinit(d)<sum(Production)-sum(Demand)-Export_Target(d) %if the battery cannot take all the extra energy, we relax the Price in the objective function, otherwise simulation times become too big
        Export_Target(d) = Export_Target(d)-(SoCmax-SoCinit(d)-(sum(Production)-sum(Demand)));
        if Export_Target(d) >0
            weightMarket = 1;
            Export_Flexibility_Up = 200;%some relaxation. Not an issue as the wight = -1 will tend the objective function to decrease the negative export
            Export_Flexibility_Down = max(0,sum(Production)-sum(Demand) - (SoCmax-SoCinit(d)));  %some relaxation. Not an issue as the wight = -1 will tend the objective function to decrease the negative export
        end
    end
    if Pbatmax*LengthOptim*Efficiency < sum(Demand)-sum(Production)+Export_Target(d)   %if the battery does not have enough power/time to give everything, we need to relax the constraints
        Export_Target(d) =  - (sum(Demand)-sum(Production)-Pbatmax*LengthOptim*Efficiency);
         if Export_Target(d)<0
          Export_Flexibility_Down = max(100,(sum(Demand)-sum(Production)+Export_Target(d)-Pbatmax*LengthOptim*Efficiency)*1.5);%200;  %some relaxation. Not an issue as the wight = -1 will tend the objective function to decrease the negative export      
         else
          Export_Flexibility_Up = max(100,(sum(Demand)-sum(Production)+Export_Target(d)-Pbatmax*LengthOptim*Efficiency)*1.5); %200;  %some relaxation. Not an issue as the wight = -1 will tend the objective function to decrease the negative export 
         end
    end
    
    if Pbatmax*LengthOptim*Efficiency < sum(Production)-sum(Demand)-Export_Target(d)   %if the battery does not have enough power/time to receive everything, we need to relax the constraints
        Export_Target(d) =  (sum(Production)-sum(Demand)-Pbatmax*LengthOptim*Efficiency);
         if Export_Target(d)<0
          Export_Flexibility_Down = max(100,(sum(Demand)-sum(Production)+Export_Target(d)-Pbatmax*LengthOptim*Efficiency)*1.5);%200;  %some relaxation. Not an issue as the wight = -1 will tend the objective function to decrease the negative export      
         else
          Export_Flexibility_Up = max(100,(sum(Demand)-sum(Production)+Export_Target(d)-Pbatmax*LengthOptim*Efficiency)*1.5); %200;  %some relaxation. Not an issue as the wight = -1 will tend the objective function to decrease the negative export 
         end
    end
    
    
    if Export_Target(d)<-epsilon
        weightMarket = -1;  
%         soctargettype(d) = 0;
        Export_Target(d) = Export_Target(d)+ min(0,sum(Production)-sum(Demand));
        Export_Flexibility_Down = 200;  %some relaxation for solving the optimisation. Not an issue as the wight = -1 will tend the objective function to decrease the negative export
    elseif Export_Target(d)<0  %we do not want to have issues with values close to 0
        Export_Target(d)=0;
    end
        Pmax = 2*max(max(max(Demand),max(Production)),Pbatmax);
% if soctargettype(d) ==1 && (sign(SoCTarget-SoCinit(d))*(SoCTarget-SoCinit(d))>Pbatmax*LengthOptim)
%     %If the battery could never reach the target, we have to relax the
%     %optimization
% %     SoCTarget_Flexibilityadd = SoCTarget - (SoCinit(d)+sign(SoCTarget-SoCinit(d))*Pbatmax*LengthOptim);
% end
    %% Fifth Implementation with MILP: integrating P-D and export target
    tempscpu = cputime;

    H = -diag(ones(LengthOptim,1),0)+diag(ones(LengthOptim-1,1),1);
    H1 = [1 zeros(1,LengthOptim-1);H(1:LengthOptim-1,:)]; %used to minimize the difference between 2 charging powers, so the battery does not go back to 0 spiriously
    
%         Pgrid In ;  Pgrid Out ; Pbat discharge ;     Pbat charge    ;   SoC
%         at time 30 min ;  SoCtarget;  x = abs(SoC at time 30min - SoC target);
%         Export Target ;  Y = abs(Export Target - Total Export); Xi_1; Xi_2; a;
%         b; c = count of number of time steps where the battery power is greater
%         than during previous timestep
    f = [BPrice*1;1*(1-weightMarket)*SPrice;zeros(LengthOptim,1); zeros(LengthOptim,1);0; 0; soctargettype(d)*(1+weightMarket)/2*(1-weightMarket)/SoCmax*1 ;0; weightMarket*1/(Pbatmax*LengthOptim*2); zeros(LengthOptim,1);zeros(LengthOptim,1);zeros(LengthOptim,1); zeros(LengthOptim,1) ; 1/100*zeros(LengthOptim,1)]; %60 is arbitrary. it corresponds to 2 intervals, but the bigger, the faster the optim
%     f = [BPrice*1;1*(1-weightMarket)*SPrice;zeros(LengthOptim,1); zeros(LengthOptim,1);0; 0; soctargettype(d)*(1+weightMarket)*(1-weightMarket) ;0; weightMarket; zeros(LengthOptim,1);zeros(LengthOptim,1);zeros(LengthOptim,1); zeros(LengthOptim,1) ; 1/100*zeros(LengthOptim,1)]; %60 is arbitrary. it corresponds to 2 intervals, but the bigger, the faster the optim

    Index_binary = size(f,1)-3*LengthOptim+1:size(f,1); %[121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144];
    Aeq = [eye(LengthOptim) zeros(LengthOptim) EfficiencyD*eye(LengthOptim) -eye(LengthOptim)/EfficiencyC zeros(LengthOptim,3) zeros(LengthOptim,2) -eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %Pgrin in = Demand - Production + Pbat_recharge/efficiency -Pbat_Discharge*Efficiency +X1   || X1 is here to make t possible to have Pgrid IN !=0 when Pgrid out =0
           zeros(LengthOptim) eye(LengthOptim) -EfficiencyD*eye(LengthOptim) eye(LengthOptim)/EfficiencyC zeros(LengthOptim,3) zeros(LengthOptim,2) zeros(LengthOptim) -eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %Pgrin out = - Demand + Production - Pbat_recharge/efficiency + Pbat_Discharge*Efficiency +X2 || X2 is here to make t possible to have Pgrid OUT !=0 when Pgrid in =0
           (1+weightMarket)/2*soctargettype(d)*[zeros(1,LengthOptim) zeros(1,LengthOptim) ones(1,LengthOptim) -ones(1,LengthOptim) 1 0 0 0 0 zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim)]; % SoC at time 30 min = sum of all previous batterie energy/power
           zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) 1 -1 1 0 0 zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim)]; % defines x =  SoCtarget - SoC_at_time_30min 
        %   abs(weightMarket)*[-ones(1,LengthOptim) ones(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) 0 0 0 -1 1 zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim)]];% defines Y = Export - Export_Target

       beq = [ Demand-Production ; Production-Demand ; (1+weightMarket)/2*soctargettype(d)*SoCinit(d);0];%0;0];
    A1 = -diag(ones(LengthOptim,1),0)+diag(ones(LengthOptim-1,1),1);
    A = [eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim,3) zeros(LengthOptim,2) zeros(LengthOptim) zeros(LengthOptim) -Pmax*eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %makes sure Pgrid in = 0 when Pgrid out != 0
         zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim,3)  zeros(LengthOptim,2) zeros(LengthOptim) zeros(LengthOptim) Pmax*eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim);  %makes sure Pgrid out = 0 when Pgrid in != 0
         zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim,3)  zeros(LengthOptim,2) zeros(LengthOptim) eye(LengthOptim) -Pmax*eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %makes sure X2  = 0 when X1 != 0
         zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim,3)  zeros(LengthOptim,2) eye(LengthOptim) zeros(LengthOptim) Pmax*eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %makes sure X1  = 0 when X2 != 0
         zeros(LengthOptim) zeros(LengthOptim) -tril(ones(LengthOptim)) tril(ones(LengthOptim)) zeros(LengthOptim,3)  zeros(LengthOptim,2) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %makes sure SoC <= SoC max
         zeros(LengthOptim) zeros(LengthOptim) tril(ones(LengthOptim)) -tril(ones(LengthOptim)) zeros(LengthOptim,3)  zeros(LengthOptim,2) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim); %makes sure SoC >= SoCmin
         zeros(LengthOptim) zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim) zeros(LengthOptim,3)  zeros(LengthOptim,2) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) -Pbatmax*eye(LengthOptim)  zeros(LengthOptim); %makes sure Pbatdischarge  = 0 when Pcharge != 0
         zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) eye(LengthOptim) zeros(LengthOptim,3)  zeros(LengthOptim,2) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) Pbatmax*eye(LengthOptim) zeros(LengthOptim); %makes sure Pbat_charge  = 0 when Pdischarge != 0
        abs(weightMarket)*[-ones(1,LengthOptim) ones(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) 0 0 0 -1 -1 zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim)];% defines Y = Export - Export_Target
        abs(weightMarket)*[ones(1,LengthOptim) -ones(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) 0 0 0 1 -1 zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim) zeros(1,LengthOptim)];% defines Y = Export - Export_Target
         0*abs(weightMarket)*[zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) [1 zeros(1,LengthOptim-1);A1(1:LengthOptim-1,:)] zeros(LengthOptim,3) zeros(LengthOptim,2)  zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) -eye(LengthOptim); % a>=PbatCharging(t)-PbatCharging(t-1)
         zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) 1/Pmax*[1 zeros(1,LengthOptim-1);A1(1:LengthOptim-1,:)] zeros(LengthOptim,3) zeros(LengthOptim,2)  zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) zeros(LengthOptim) 1/Pmax*eye(LengthOptim)]]; % a/M<=1+1/M*(PbatCharging(t)-PbatCharging(t-1))


    b = [zeros(LengthOptim,1); Pmax*ones(LengthOptim,1); zeros(LengthOptim,1); Pmax*ones(LengthOptim,1); (SoCmax-SoCinit(d))*ones(LengthOptim,1); (SoCinit(d) - SoCmin)*ones(LengthOptim,1); zeros(LengthOptim,1); Pbatmax*ones(LengthOptim,1); 0;0; zeros(LengthOptim,1); ones(LengthOptim,1)]; 
    lb = [zeros(LengthOptim,1); zeros(LengthOptim,1); zeros(LengthOptim,1); zeros(LengthOptim,1); SoCmin;                      ...
         SoCTarget; (SoCTarget_Flexibilityadd);                      Export_Target(d);  0; zeros(LengthOptim,1);  zeros(LengthOptim,1);  zeros(LengthOptim,1); zeros(LengthOptim,1); zeros(LengthOptim,1)]  ;
    ub = [Pmax*ones(LengthOptim,1); Pmax*ones(LengthOptim,1); Pbatmax*ones(LengthOptim,1);  Pbatmax*ones(LengthOptim,1); SoCmax;  ...
         SoCTarget; SoCmax;    Export_Target(d); Pmax;       Pmax*ones(LengthOptim,1);  Pmax*ones(LengthOptim,1);  ones(LengthOptim,1);  ones(LengthOptim,1) ;  Pbatmax*ones(LengthOptim,1)]  ;

    [x,fval1,exitflag1,output1]  = intlinprog(f, Index_binary, A,b, Aeq, beq, lb, ub);
%      x3 = -x(2*LengthOptim+1:3*LengthOptim) +   x(3*LengthOptim+1:4*LengthOptim)   
%      SoCfinal = x(4*LengthOptim+1)
%      Final_Export = sum(-x(1:LengthOptim)+x(LengthOptim+1:2*LengthOptim))
%     fval1
%     SoC = -tril(ones(LengthOptim))*x(LengthOptim*2+1:LengthOptim*3)+tril(ones(LengthOptim))*x(LengthOptim*3+1:LengthOptim*4)+SoCinit(d)*ones(LengthOptim,1);
%     PgridIN = x(1:LengthOptim);
%     PgridOUT = x(LengthOptim+1:LengthOptim*2);
%     Pbat = x(LengthOptim*2+1:LengthOptim*3)-x(LengthOptim*3+1:LengthOptim*4);
% 
        ComputingTime(i+k,d) =  cputime - tempscpu;

            PgridIN = zeros(stepRT,1);
            PgridOUT = zeros(stepRT,1);
            Pbat = x(LengthOptim*2+1:LengthOptim*3)-x(LengthOptim*3+1:LengthOptim*4);
          if   abs(Pbat(1) - (Demand_Previous-Production_Previous))<epsilon
             for indx = 1:stepRT
                 Pbat(indx) = D(i+k-1+indx,d)-P(i+k-1+indx,d);
             end
          end
            X1 = x(LengthOptim*4+1:LengthOptim*5);
            X2 = x(LengthOptim*5+1:LengthOptim*6);
            Alpha = x(LengthOptim*6+1:LengthOptim*7);
            Beta = x(LengthOptim*7+1:LengthOptim*8);
            SoC_previous = SoCinit(d);

            for indx = 1:stepRT
                SoC(i+k+indx-1,d) = SoC_previous - (Pbat(indx));%-tril(ones(LengthOptim))*x(LengthOptim*2+1:LengthOptim*3)+tril(ones(LengthOptim))*x(LengthOptim*3+1:LengthOptim*4)+SoCinit(d)*ones(LengthOptim,1);
                s = sign(Pbat(indx));
                a = ((1-s/(abs(s)+0.0001))/2/Efficiency + (1+s/(abs(s)+0.0001))/2*Efficiency);
                PgridOUT(indx) =  max(0,max(0,P(i+k-1+indx,d)-D(i+k-1+indx,d)+ Pbat(indx)*a)); %what we export at this time step
                PgridIN(indx) = max(0,max(0,D(i+k-1+indx,d)-P(i+k-1+indx,d)- Pbat(indx)*a));
%                 Energybought(i+k-1+indx) = PgridIN(indx)*BPrice(indx);
%                 Energysold(i+k-1+indx) = PgridIN(indx)*BPrice(indx);
                FromGrid(i+k-1+indx) = FromGrid(i+k-1+indx)+ PgridIN(indx); %max(0,max(0,P(i+k-1+indx)-D(i+k-1+indx)- Pbat(indx)*a)); %what we export at this time step
                ToGrid(i+k-1+indx) = ToGrid(i+k-1+indx)+ PgridOUT(indx);%max(0,max(0,D(i+k-1+indx)-P(i+k-1+indx)+ Pbat(indx)*a));
                Powerbat (i+k-1+indx,d) = Pbat(indx);
                SoC_previous = SoC(i+k+indx-1,d);
                Already_Realised_Export_Target(d) = Already_Realised_Export_Target(d)+PgridOUT(indx)-PgridIN(indx) ;
            end  
            SoCinit(d) = SoC_previous;

    end
end


indiceExportTarget = indiceExportTarget+1;

end


%% Display
close all
ToGridTotal = ToGrid;
FromGridTotal = -FromGrid;
Export_Market_real = zeros(size(t3));
Import_Market_real = zeros(size(t3));
% Export_Market_real2 = zeros(size(t3));

j = 1;
for kindx = 1:size(t3,1)-1
        Export_Market_real(kindx+1) = max(0,sum(ToGridTotal((kindx-1)*MarketTimeStep+1:(kindx)*MarketTimeStep)) + sum(FromGridTotal((kindx-1)*MarketTimeStep+1:(kindx)*MarketTimeStep)));
        Import_Market_real(kindx+1) = min(0,sum(ToGridTotal((kindx-1)*MarketTimeStep+1:(kindx)*MarketTimeStep)) + sum(FromGridTotal((kindx-1)*MarketTimeStep+1:(kindx)*MarketTimeStep)));
%         Export_Market_real2(k) = sum(RealisedOutputMarketVector((k-1)*MarketTimeStep+1:(k)*MarketTimeStep)); 
end 



% close all
figure4 = figure('Name','Timeseries','Color',[1 1 1]);
plot(t,SoC(:,1),t,SoC(:,2),t,(Porigin(:,1)-Dorigin(:,1))*30,t,(Porigin(:,2)-Dorigin(:,2))*30,t,Powerbat(:,1)*30,t,Powerbat(:,2)*30,t3,Required_Export_bid(1:size(t3,1))*0.9,...
    t3, Export_Market_real ,t,SoCTargetreal(1:size(t,1)), t3, Import_Market_real ,t,CurrentModeVectorreal*500, 'LineWidth',2)
legend('SoC 1','SoC 2','P-D 1','P-D 2','Pbat 1','Pbat 2','Export Target','Export Real','SoC Target','Imported Power','Mode')

mean(ComputingTime(:,1))
mean(ComputingTime(:,2))


%% Functions
function [pbat,soc] = coercsoc(soc,socmax,socmin,pbat,timestep,lastsoc, nextsoc,soctargetstype,pbatmax)
bool1 = 1;
bool2 = 1;
 
    if (lastsoc - nextsoc < 0.001) && soc>nextsoc && soctargetstype==1  %if we already charged too much
        bool1=0;
        pbat = min(pbatmax, max(-pbatmax,pbat - (soc-nextsoc)*timestep)); %pbat + (soc-nextsoc)*timestep;
%         soc = nextsoc;
    elseif    (lastsoc - nextsoc >0.001) && soc< nextsoc && soctargetstype==1  %if we already discharged too much
        bool2 = 0;
        pbat =  min(pbatmax, max(-pbatmax, pbat-(soc - nextsoc)*timestep)); %pbat +(soc - nextsoc)*timestep; 
%         pbat =  min(pbatmax, max(-pbatmax,pbat +(soc - nextsoc)*timestep)); %pbat +(soc - nextsoc)*timestep; 

%         soc = nextsoc;
    end

    if soc > socmax && bool1
        pbat = min(pbatmax, max(-pbatmax,pbat - (soc-socmax)*timestep)); %pbat + (soc-socmax)*timestep;
%         soc = socmax;
    elseif soc<socmin && bool2
        pbat = min(pbatmax, max(-pbatmax,pbat -(soc - socmin)*timestep)); %pbat +(soc - socmin)*timestep;
%         soc = socmin;
    end
    

end


