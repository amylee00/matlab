%% 1. Import data

Data.fileName='X:\investments\Chika\Quantitative\World Growth\Economic_Cycle_rawdata.xlsx';

Sheetname='EMGDP';
%  'G10GDP'
%  'EuroGDP'
%  'EMGDP'
[IndicesRawData,Indicestxt] = xlsread(Data.fileName,Sheetname);
DateStrings=Indicestxt(3:end,1);

%% 2. Declare all variables

MAwindow_size=4; %Moving average of 4 quarters
Cyclewindow_size=16; % 4*4 years as a cycle
HalfCycle=Cyclewindow_size/2; %4 quarters per year, over 2 years
Results=zeros(16,size(IndicesRawData,2));
DatesOutput=cell(3,size(IndicesRawData,2));
OutputContents=cell(15,1);
Assumptions=zeros(2,1);
XYplane=zeros(2,size(IndicesRawData,2));

% Row 1: Current Date
% Row 2: Current Smoothed GDP growth rate (smoothed over the MAwindow_size)
% Row 3: Stochastic over Cyclewindow_size
% Row 4: Last date with Min GDP growth %
% Row 5: Last Min GDP growth % (date as of Row 3)
% Row 6: Last date with Max GDP growth %
% Row 7: Last Max GDP growth % (date as of Row 5)
% Row 8: economic trend (1 for expansion, -1 for contraction)
% Row 9: Time into the cycle (expansion or contraction), over the HalfCycle*2
% Row 10: Time from Min
% Row 11: Time from Max
% Row 12: Adjusted Time from Min (max value =8)
% Row 13: Time from Max (max value =8)
% Row 14: Cyclicality - (max-min)/stdev
% Row 15: x plane value
% Row 16: y plane value

%% 3. Calculate stochastic, decide economic trend, store into results
            
for i = 1:size(IndicesRawData,2)
            
            tempdate=DateStrings(1:end,1);
            tempraw=IndicesRawData(:,i);
            
            tempdate(isnan(tempraw(:,:)))=[];
            tempraw(isnan(tempraw(:,:)))=[];
            
            CurDate=tempdate(end,1);     %%need to be different date for those deleted dates       
            DatesOutput(1,i)=CurDate;   
            [r,c]=size(tempdate);
                        
            tempGDPSmoothed=tsmovavg(tempraw,'s',MAwindow_size,1);
            
            Results(2,i)=tempGDPSmoothed(end,1);
            
            tempmin=nanmin(tempGDPSmoothed(end-Cyclewindow_size+1:end,1));
            tempmax=nanmax(tempGDPSmoothed(end-Cyclewindow_size+1:end,1)); 

            MinDate_indx=find(tempGDPSmoothed(end-Cyclewindow_size+1:end,1)==tempmin,1,'last');           
            MinDate=(tempdate(r-Cyclewindow_size+MinDate_indx,:));
            
            MaxDate_indx=find(tempGDPSmoothed(end-Cyclewindow_size+1:end,1)==tempmax,1,'last');           
            MaxDate=(tempdate(r-Cyclewindow_size+MaxDate_indx,:));         
            
            tempStoc=(tempGDPSmoothed(end)-tempmin)/(tempmax-tempmin);
            Results(3,i)=tempStoc;
            DatesOutput(2,i)=MinDate;
            Results(5,i)=tempmin;
            DatesOutput(3,i)=MaxDate;
            Results(7,i)=tempmax;           
          
            QuartersFromMin = Cyclewindow_size - MinDate_indx;
            QuartersFromMax = Cyclewindow_size - MaxDate_indx;
            %DaysFromMax = datenum(CurDate,'dd/mm/yy') - datenum(MaxDate,'dd/mm/yy');            
                        
            if QuartersFromMax > QuartersFromMin
            Sign=1;
            else
            Sign=-1;      
            end            
            Results(8,i)=Sign;
            
            if Sign==1          
            QuartersIntoCycle=QuartersFromMax/HalfCycle;
            else        
            QuartersIntoCycle=QuartersFromMin/HalfCycle;
            end
            Results(9,i)=QuartersIntoCycle;
            
            Results(10,i)=QuartersFromMin;
            Results(11,i)=QuartersFromMax;
           
            
            if QuartersFromMin<HalfCycle
            Results(12,i)=QuartersFromMin;
            else
            Results(12,i)=HalfCycle;  
            end
            
            if QuartersFromMax<HalfCycle
            Results(13,i)=QuartersFromMax;
            else
            Results(13,i)=HalfCycle;
            end
            
            tempstd=nanstd(tempGDPSmoothed,0,1);
            Results(14,i)=(tempmax-tempmin)/tempstd;
            
            if Sign==1          
            x=-(1-Results(12,i)/HalfCycle);
            else        
            x= Results(13,i)/HalfCycle;
            end
            
            Results(15,i)=x;
            Results(16,i)=tempStoc;
            XYplane(1,i)=x;
            XYplane(2,i)=tempStoc;
            
end
            
%% 4. Output into Excel

Data.OutputFileName='X:\investments\Chika\Quantitative\World Growth\Economic_Cycle_output.xlsx';
Banner=Indicestxt(2,3:end);
Results=num2cell(Results);
Results(1,:)=DatesOutput(1,:);
Results(4,:)=DatesOutput(2,:);
Results(6,:)=DatesOutput(3,:);
OutputContents={'Country';'Latest Quarter';'Latest Smoothed GDP Growth(%)';'Stochastic';'Last Smoothed Min GDP Date';'Last Smoothed Min GDP Growth(%)';'Last Smoothed Max GDP Date';'Last Smoothed Max GDP Growth(%)';'Economic Cycle(1:Expansion,-1:Contraction)';'Time Into half cycle(=2yrs)';'Quarters since last Min';'Quarters since last Max';'Adjusted Quarters since last Min';'Adjusted Quarters since last Max';'(max-min)/stdev - cyclicality';'x';'y'};
Output=[Banner;Results];
TotalOutput=[OutputContents,Output];
Assumptions(1,1)=MAwindow_size;
Assumptions(2,1)=Cyclewindow_size;
Assumptions=num2cell(Assumptions);

x=XYplane(1,:);
y=XYplane(2,:);
XYplane=num2cell(XYplane);
XYplaneOutput=[Banner;XYplane];

AssumptionsContents={'Moving Average Window(quarters)';'Cycle length(quarters)'};
AssumptionsOutput=[AssumptionsContents,Assumptions];

if strcmp(Sheetname,'G10GDP')==1
xlswrite(Data.OutputFileName,TotalOutput,'G10','A1');
xlswrite(Data.OutputFileName,AssumptionsOutput,'G10','A26');
elseif strcmp(Sheetname,'EuroGDP')==1
xlswrite(Data.OutputFileName,TotalOutput,'Euro','A1');
xlswrite(Data.OutputFileName,AssumptionsOutput,'Euro','A26');
elseif strcmp(Sheetname,'EMGDP')==1
xlswrite(Data.OutputFileName,TotalOutput,'EM','A1');
xlswrite(Data.OutputFileName,AssumptionsOutput,'EM','A26');
else  
end

%% 5. Plot x, y and Sine Curve 

% plot(x,y,'o');
% labels=Indicestxt(1,3:end);
% h = labelpoints(x,y,labels); 
% hold on
% 
% %sine curve (not a fitting, actual sine curve)
% x1=0:0.05:1;    %this is x used for calculating y value
% x2=-1:0.05:0;   %from x = -1 to x=0 (left half)
% x3=1:-0.05:0;   %from x = 1 to x=0 (right half)
% A=1; % amptitude
% f=pi/2; % frequency
% y1=A*sin(f*x1);
% y2=A*sin(f*x3);
% plot(x2,y1);
% plot(x1,y2);
% hold on
 
% [param]=sine_fit(x,y);

%This part is for plotting fitting curve into the x y data 
% yu = max(y);
% yl = min(y);    
% yr = (yu-yl);                               % Range of ‘y’
% yz = y-yu+(yr/2);
% zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
% per = 32*mean(diff(zx));                     % Estimate period - need to be able to automatically calculate this 
% ym = mean(y);        
% fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);    % Function to fit
% fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
% s = fminsearch(fcn, [yr;  per;  -1;  ym]);
% xp = linspace(min(x),max(x));
% plot(xp,fit(s,xp), 'r');
% grid;
% hold on;
