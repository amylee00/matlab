
%% 1. Import data & mapping

%import raw index & ratio data and the list of reference index
Data.fileName='X:\investments\Chika\Valuation\Matlab Project\Input_rawdata_new';

Data.Sheetname='data';
[Data.Indices, Data.Indicesraw,Data.Indicestxt] = xlsread(Data.fileName,Data.Sheetname);
Data.DateStrings=Data.Indicestxt(6:end,1);

Data.Sheetname='list';
[Data.list, Data.listraw,Data.listtxt] = xlsread(Data.fileName,Data.Sheetname);

tic
%% 2. Declare variables

%weekly data
Stat.MA_Wsize = 52/4*3; %in weeks - moving average of 9 months(default)
Stat.Stoc_Wsize= Stat.MA_Wsize*2; %Stochastic calculated on 1.5 year (default)
Stat.Cycle_Length= 52*4; % 4 years - used for calculating the cycle average & Z-score
Stat.Results_Length= 52*2; % 2 years - used for rebasing the recent 2 years of indices data
Stat.Threshold=0.30; % above and below 4 year average 0.3 is default
Stat.Period=52; % how many period in 1 year - for calculating annualised return
Stat.OW=1.5; %overweight when the signal is buy - this is for performance calculation
Stat.UW=-0.5; %underweight when the signal is sell - this is for performance calculation
Stat.Results=20; %show top 20 results

%this is for the look up later
Data.Index_combined = strcat(Data.Indicestxt(1,:),Data.Indicestxt(3,:)); 
Data.Index_combined = Data.Index_combined(1,2:end); %this is for looking for the referece index location
Data.Index_combined3 = strcat(Data.Indicestxt(3,:),Data.Indicestxt(4,:),Data.Indicestxt(1,:),Data.Indicestxt(2,:));
Data.Index_combined3 = Data.Index_combined3(1,2:end); %this is for looking up the unique z-scores location

%This is for output
Output.OutputContents={'Indices','Description','Type','Valuation Type','Column #','Info Ratio','Signal count','Latest regime','Latest signal(1=BUY,-1=SELL)','4-yr z-score','Ref Index','Weighted Zscore'};
Stat.table_title={'Ref Index';'Sum of Weighted Z-score';'Sum of Trade Sign'};
Output.Titles_Unique={'Index - Price/Valuation-Valuation method','Frequency seen in separate dataset'};
Output.Titles_Zscore={'Z-score'};

%% 3. Calculation of moving average (MA) of indices and valuations

%Calculate MA of the raw data
Data.Smoothed=tsmovavg(Data.Indices,'s',Stat.MA_Wsize,1);
%Calculate cycle length average
Data.cycle_avg=tsmovavg(Data.Indices,'s',Stat.Cycle_Length,1);

%% 4. Getting the stochastic
%result of this section to be used in chapter 7. for getting the location in xy plane to get the
%current regime

Output.Stoc=zeros(size(Data.Smoothed));
    
    for j = Stat.Stoc_Wsize+Stat.MA_Wsize:size(Data.Smoothed,1)
        Temp.tempmin=nanmin(Data.Smoothed(j-Stat.Stoc_Wsize+1:j,:));
        Temp.tempmax=nanmax(Data.Smoothed(j-Stat.Stoc_Wsize+1:j,:));
        Output.Stoc(j,:)=(Data.Smoothed(j,:)-Temp.tempmin)./(Temp.tempmax-Temp.tempmin);
    end

%% 5. Calculating 4-yr z-score

%this is for the level
for k = Stat.Cycle_Length:size(Data.Indices,1) %run loop over time, not over different columns
        Results.Zscore(k,:) = (Data.Indices(k,:)- nanmean(Data.Indices(k-Stat.Cycle_Length+1:k,:)))./nanstd(Data.Indices(k-Stat.Cycle_Length+1:k,:));
end

%Calculate return for z-score
Data.Indices_chg=Data.Indices(2:end,:)./Data.Indices(1:end-1,:)-1; %calculate the return of the ref index
%this is for calculating z-score on the change
for kk = Stat.Cycle_Length:size(Data.Indices_chg,1) %run loop over time, not over different columns
        Results.Zscore_chg(kk,:) = (Data.Indices_chg(kk,:)- nanmean(Data.Indices_chg(kk-Stat.Cycle_Length+1:kk,:)))./nanstd(Data.Indices_chg(kk-Stat.Cycle_Length+1:kk,:));
end

%% 6. Detecting bull or bear phase
%result of this section to be used in getting the regime - chapter 7.
Output.Sign = zeros(size(Data.Indices,1),size(Data.Indices,2));
    
for l = 1:size(Data.Indices,2)
    for m = 2:size(Data.Indices,1)
        if (floor((Data.Indices(m,l)/Data.cycle_avg(m,l)-1)*100)/100)>Stat.Threshold
            Output.Sign(m,l)=1;  %larger than threshold - get bull sign
        elseif (ceil((Data.Indices(m,l)/Data.cycle_avg(m,l)-1)*100)/100)<Stat.Threshold*-1
            Output.Sign(m,l)=-1;  %smaller than threshold - get bera sign
        else
            Output.Sign(m,l)=Output.Sign(m-1,l); %no bull nor bear sign - keep previus value
        end
    end
end

%% 7. Getting the regime
%result of this section to be used in getting the OW/UW trading decisions in chapter 8.

%Below multiplication gives unique values depending on where you are in the four regimes
Signal.Signals=Output.Stoc.*Output.Sign; %Stochastic value multiplied by bull/bear sign
Output.Regime=Signal.Signals;

Output.Regime(Output.Regime~=1 & Output.Regime>0)=3;
Output.Regime(Output.Regime==1)=2;
Output.Regime(Output.Regime<0)=1;
Output.Regime(Output.Regime==0)=4;

%Restore zeros to rows where no stochastic exists
Output.Regime(1:Stat.Stoc_Wsize+Stat.MA_Wsize-1,:)=0;

%% 8. Allocate OW/UW to different regimes
%this is for trading decisions rules

Signal.trades=zeros(size(Output.Regime,1),size(Output.Regime,2));

for p = 1:size(Output.Regime,2)
    for q = 1:size(Output.Regime,1)
        if or(Output.Regime(q,p)== 1,Output.Regime(q,p)== 2)
        Signal.trades(q,p)=Stat.OW;  % to buy at regime 1 & 2 - this worked the best among other combinations
        elseif or(Output.Regime(q,p)== 3,Output.Regime(q,p)== 4)
        Signal.trades(q,p)=Stat.UW;  % to sell at regime 3 & 4 - this worked the best among other combinations
        else
        Signal.trades(q,p)=0;
        end
    end
end

%% 9. Convert trading signal to continuous signal, count turn of signals
%so we only trade when there is clear and longer term change in direction

Signal.trades_cont=zeros(size(Signal.trades)); %continuous trades signal
Output.Excess_return=zeros(size(Signal.trades));

for r=1:size(Signal.trades,2)
    counter = 0;
        for s=2:size(Signal.trades,1)
            if and(Signal.trades(s,r)~=0,Signal.trades(s,r)==Signal.trades(s-1,r))
               Signal.trades_cont(s,r) = Signal.trades(s-1,r);

            elseif Signal.trades(s,r)~=Signal.trades(s-1,r)
               Signal.trades_cont(s,r) = Signal.trades(s,r);
               counter = counter +1; % this counts turn of continuous signal
               Output.Results(r,3)=counter; % counter out put for the results
            else
            end
        end        
        Output.Results(r,1)=r; % output of column number from the raw data for the results
        Output.Results(r,4)=Output.Regime(end,r); % latest regime out put for the results

        if or(Output.Results(r,4)==1,Output.Results(r,4)==2) % latest trading signal out put for the results
               Output.Results(r,5)=1;
        else
               Output.Results(r,5)=-1;
        end    
        Output.Results(r,6)=Results.Zscore(end,r); % z-score out put for the results
end

toc
tic
 
%% 10. Calculate the returns and IR, sort by IR, output results

% store results in 3D, defining the size of the results
Output.SummaryTable=num2cell(zeros(Stat.Results,(Stat.Results/2+1),size(Data.listtxt,1)));
Stat.table=num2cell(zeros(3,size(Data.listtxt,1)));

%loop for the number of reference indices - the loop is over chapter 10. and 11.
for t = 1:size(Data.listtxt,1)   % loop for the output so we don't need to manually change reference index

 [rr,cc]=find(strcmp(Data.Index_combined,Data.listtxt(t,2))); %find the address of the reference index
 Return.Index_key=Data.Indices(2:end,cc)./Data.Indices(1:end-1,cc)-1; %calculate the return of the ref index
 Return.Index_return=[0;Return.Index_key]; %adding zero to make the lengh equal with number of trade signal
 Return.Index_return_matrix=repmat(Return.Index_return,1,size(Signal.trades,2)); %prepare matrix of the same to calculate trading sign adjusted return
 Return.Index_return_matrix2=repmat(Return.Index_key,1,size(Signal.trades,2));
 Export.Z_score(:,t)=num2cell(Results.Zscore(:,cc));
 
 %Delaying index return by one period to signal, multiplying to reference
 %index to get the strategy return
 Return.Strategy_return = Signal.trades_cont(1:end-1,:).*Return.Index_return_matrix(2:end,:);
 
 %Calculating information ratio
 Return.Relative_return=Return.Strategy_return-Return.Index_return_matrix2; %Strategy return (from above) minus reference index return, for information ratio calculation
 Return.Relative_annualised=(1+nansum(log(1+Return.Relative_return)))/(size(Return.Relative_return,1)/Stat.Period);%calculating annualised relative reutrn
 Return.Vol_annualised=nanstd(Return.Relative_return).*sqrt(Stat.Period);%calculating annualised standard deviation of the relative return
 Return.InfoRatio_annualised=Return.Relative_annualised./Return.Vol_annualised; %calculating so called sharp ratio
 Output.Results(:,2)=Return.InfoRatio_annualised'; %output of IR for the results
 
 %Get details of underlying signaling indices from the raw data
 Output.Results_name=[Data.Indicestxt(3:4,Output.Results(:,1)+1);Data.Indicestxt(1:2,Output.Results(:,1)+1)]';
 Results.Results_combined=[Output.Results_name,num2cell(Output.Results)]; %add index descriptions to the Output
 
 %Sort by excess return, extract top 10%, sort rows in descending order
 Results.Results_combined=sortrows(Results.Results_combined,6);
 Results.Results_combined(1:floor(size(Results.Results_combined,1))*(size(Data.Indicestxt,2)-21)/(size(Data.Indicestxt,2)-1),:)=[];
 Results.Results_combined=sortrows(Results.Results_combined,-6);
 
 %Prepare matrix of reference name for output
 Tab=char(Data.listtxt(t,1));
 Output.RefIndex=repmat({Tab},[1 Stat.Results])';
 Results.Results_combined=[Results.Results_combined,Output.RefIndex];
 
 %Separate top 10 list for the output 'Top10'
 Output.RefIndex_top10=repmat({Tab},[1 Stat.Results/2+2]);
 Results.Results_combined_top10=Results.Results_combined(1:Stat.Results/2,:);
 Stat.Top10_IRsum=sum(cell2mat(Results.Results_combined_top10(:,6))); %sum IR ratio
 Stat.Top10_IRweightedZ=cell2mat(Results.Results_combined_top10(:,6))./Stat.Top10_IRsum.*cell2mat(Results.Results_combined_top10(:,10));
 Results.Results_combined_top10=[Results.Results_combined_top10,num2cell(Stat.Top10_IRweightedZ)];
 Stat.SumofWeightedZ=sum(Stat.Top10_IRweightedZ);
 Stat.table(2,t)=num2cell(Stat.SumofWeightedZ); %fill stat table with the sum of the weighted z-score
 Stat.table(1,t)=cellstr(Tab); %fill stat table with the reference index name
 Stat.SumofTradeSign=sum(cell2mat(Results.Results_combined_top10(:,9))); %Calculate sum of trades sign for output
 Stat.table(3,t)=num2cell(Stat.SumofTradeSign); %fill stat table with sum of the trade sign
 
toc
tic

%% 10. To record results in 3D, then stack the output by data type, prepare result for output

Output.SummaryTable(:,:,t)=Results.Results_combined; %store result into 3D matrix

    if t==1
        Export.Stack=Results.Results_combined;
        Export.Stack_top10=[Output.RefIndex_top10;Output.OutputContents;Results.Results_combined_top10];
    else
        Export.Stack=[Export.Stack;Results.Results_combined];
        Export.Stack_top10=[Export.Stack_top10;Output.RefIndex_top10;Output.OutputContents;Results.Results_combined_top10];
    end

end

%preparing for outpout, adding titles
Export.end=repmat(cellstr('xxx end xxx'),1,size(Export.Stack,2));
Export.Stack=[Output.OutputContents(:,1:end-1);Export.Stack;Export.end];
Stat.table=[Stat.table_title,Stat.table];
%Export.Z_score_title=repmat({Output.Titles_Zscore},1,size(Data.listtxt,1));
Export.Date=[0;0;Data.DateStrings];
Export.Z_score=[repmat(cellstr('Z-score'),1,size(Export.Z_score,2));Data.listtxt(:,1)';Export.Z_score];
Export.Z_score=[Export.Date,Export.Z_score];
toc

%% 11. Find unique indices which appears across different reference indices as signal - prepare data for output

Data.Index_combined2 = strcat(Export.Stack(:,1),Export.Stack(:,2),Export.Stack(:,3),Export.Stack(:,4)); %combining the strings to make unique index names
Data.Unique=unique(Data.Index_combined2,'stable'); %looking for unique values

for x=1:size(Data.Unique,1)
    Stat.Unique_value=Data.Unique(x,1); %find unique value
    Stat.Unique_value_matrix=repmat(Stat.Unique_value,size(Data.Index_combined2,1),1); %make matrix of the single unique value
    Stat.Quantity=sum(strcmp(Data.Index_combined2,Stat.Unique_value_matrix),1); % compare the single unique value matrix with the original data
    Stat.Unique(x,1)=num2cell(Stat.Quantity); %record it if it exists
end

Results.Unique=[Data.Unique,Stat.Unique]; %combining data label with the frequency
Results.Unique=sortrows(Results.Unique,-2); %sorting by the frequency
Stat.HighestFreq=cell2mat(Results.Unique(1,2)); %pick the value of highest frequency
Results.Unique((cell2mat(Results.Unique(:,2))<(Stat.HighestFreq-1)),:)=[]; %delete those other than the highest and highest minus 1
Export.Results_Unique_summary=Results.Unique;
Results.Unique=[Output.Titles_Unique;Results.Unique;repmat(cellstr('xxx end xxx'),1,2)]; %preparing for output, adding titles
Output.Unique_Zscore=zeros(size(Results.Zscore,1),size(Results.Unique,1)-2);
Output.Unique_Zscore_chg=zeros(size(Results.Zscore_chg,1),size(Results.Unique,1)-2);

%getting the unique index price z-scores
for u = 2:size(Results.Unique,1)-1
    [rrr,ccc]=find(strcmp(Data.Index_combined3(1,:),Results.Unique(u,1)));
    Output.Unique_Zscore(:,u-1)= Results.Zscore(:,ccc);
    Output.Unique_Zscore_chg(:,u-1)= Results.Zscore_chg(:,ccc);
    Output.Unique_Zscore_title(1,u-1)=Data.Index_combined3(1,ccc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Output.Unique_Zscore_chg(isnan(Output.Unique_Zscore_chg))=0;
[coeff,score,latent,tsquared,explained,mu]=pca(Output.Unique_Zscore_chg,'Algorithm','svd','NumComponents',1, 'Centered',true);

Output.score=zeros(size(Output.Unique_Zscore_chg,1),1);

for a=1:size(Output.Unique_Zscore_chg,1)
        if abs(Output.Unique_Zscore_chg(a,:))>0
        [coeff,score,latent,tsquared,explained,mu]=pca(Output.Unique_Zscore_chg(1:a,:),'NumComponents',1);
        Output.score(a,1)=score(a,1);
        else
        Output.score(a,1)=0;
        end    
end

%prepare plot output for PCA
Start_line=211;%start line of the data in the graph - where there are all data
Output.score1=Output.score(Start_line:end,1);
Output.score1(:,2)=1+Output.score1/100;
Output.score1(1,2)=1;
Output.score1(:,3)=cumprod(Output.score1(:,2));
Output.score1=Output.score1(:,3);
Export.score1=[{'PC1','PC1'};[{'Dates','Factor Index(PC1)'};[Data.DateStrings(Start_line+1:end,1),num2cell(Output.score1)]]];

%this is for charting in Matlab
startDate=datetime(Data.DateStrings(Start_line+1),'InputFormat','dd/MM/yyyy'); %just get the first date where every data exists
endDate=datetime(Data.DateStrings(end,:),'InputFormat','dd/MM/yyyy');
x=linspace(startDate,endDate,size(Output.score1(:,:),1));
plot(x,Output.score1(:,:),'DisplayName','Output.score1','LineWidth',2);
xlabel('Time');
ylabel('PC1 Factor Index');
h = legend('PC1 Factor Index','Location','northeast');
title({'Principal Component Analysis of normalised values','- using Z-scores of unique indices with high IR among S&P500, FTSE All Share, MSCI EMU, TOPIX Index and their respective sector indices -','(raw data: Index Prices and their valuation metrices)'});
annotation('textbox','String',{'Current Value as of',char(Data.DateStrings(end,1)),(':'),+Output.score1(end,1),'(PC1)'});

%prepare output for Zscore
Export.Unique_Zscore=[repmat(cellstr('Z-score'),1,size(Output.Unique_Zscore,2));Output.Unique_Zscore_title;num2cell(Output.Unique_Zscore)];
Export.Unique_Zscore=[Export.Date,Export.Unique_Zscore];
Export.Unique_Zscore_chg=[repmat(cellstr('Z-score of chg'),1,size(Output.Unique_Zscore_chg,2));Output.Unique_Zscore_title;num2cell(Output.Unique_Zscore_chg)];
Export.Unique_Zscore_chg=[Export.Date(2:end),Export.Unique_Zscore_chg];
Stat.Component_titles=[cellstr('PC1'),cellstr('PC2'),cellstr('PC3')];
Export.Stat.coeff=[repmat(cellstr('Coeff of PCA'),1,size(Output.Unique_Zscore_chg,2));Output.Unique_Zscore_title;num2cell(coeff')];

%check if the top10 results contain unique index, if so record 1
Data.Index_combined2 = strcat(Export.Stack_top10(:,1),Export.Stack_top10(:,2),Export.Stack_top10(:,3),Export.Stack_top10(:,4));
for y=1:size(Data.Index_combined2,1);
    Data.IsFrequent=ismember(Data.Index_combined2(y,1),Results.Unique(:,1));
    Data.Frequency_column(y,1)=Data.IsFrequent;
end

%adding frequency column to the output
Export.Stack_top10=[Export.Stack_top10,num2cell(Data.Frequency_column)];
Export.Stack_top10=[Export.Stack_top10;repmat(cellstr('xxx end xxx'),1,size(Export.Stack_top10,2))];

%Exporting the summary for output - arrange table in the way I want
Results.Unique_summary=num2cell(zeros(size(Data.Index_combined2,1),size(Results.Unique,1)-2));
for z=2:size(Results.Unique,1)-1
    [rrrr,cccc]=find(strcmp(Results.Unique(z,1),Data.Index_combined2(:,1)));
    if isempty(rrrr)
    Results.Unique_summary(:,z-1)={0};
    else
        for zz=1:size(rrrr)
        Results.Unique_summary(rrrr(zz),z-1)=Export.Stack_top10(rrrr(zz,1),6);
        end
    end
end

Results.Unique_summary2=num2cell(zeros(size(Data.Index_combined2,1),size(Results.Unique,1)-2));
for z=2:size(Results.Unique,1)-1
    [rrrrr,ccccc]=find(strcmp(Results.Unique(z,1),Data.Index_combined2(:,1)));
    if isempty(rrrrr)
    Results.Unique_summary2(:,z-1)={0};
    else
        for zzz=1:size(rrrrr)
        Results.Unique_summary2(rrrrr(zzz),z-1)=Export.Stack_top10(rrrrr(zzz,1),8);
        end
    end
end

Export.Summary_table=zeros(size(Results.Unique,1)-2,size(Data.listraw,1)*2);
for w=1:size(Results.Unique,1)-2
    Results.Unique_summary_1=sum(cell2mat(reshape(Results.Unique_summary(:,w),[12,5])));
    Results.Unique_summary2_1=sum(cell2mat(reshape(Results.Unique_summary2(:,w),[12,5])));
    Results.Unique_summary3_1=[Results.Unique_summary_1;Results.Unique_summary2_1];
    Results.Unique_summary3_1=reshape(Results.Unique_summary3_1,[1,10]);
    Export.Summary_table(w,:)=Results.Unique_summary3_1;
end

%Tidying up the Output for Export
Export.Summary_table2=[Results.Unique(2:end-1,1),num2cell(Export.Summary_table)];
Export.Summary_table_titles1=[{''},Data.listraw(1,1),{''},Data.listraw(2,1),{''},Data.listraw(3,1),{''},Data.listraw(4,1),{''},Data.listraw(5,1),{''}];
Export.Summary_table_titles2=[{''},repmat({'IR','Regime'},[1,5])];
Export.Summary_table2=[Export.Summary_table_titles2;Export.Summary_table2];
Export.Summary_table2=[Export.Summary_table_titles1;Export.Summary_table2;repmat(cellstr('xxx end xxx'),1,size(Export.Summary_table2,2))];

%% 12. Export data to Excel
tic
        Data.OutputFileName='X:\investments\Chika\Valuation\Matlab Project\Output_valuation_v3.xlsx';
        %%%give location - where to export data in Excel
        xlswrite(Data.OutputFileName,Export.Stack,'Top20','A1'); %export summary table - stacked up table from 3d matrix
        xlswrite(Data.OutputFileName,Results.Unique,'Unique','A2'); % export unique indices and frequency
        xlswrite(Data.OutputFileName,Export.Stack_top10,'Top10','A1'); %export top 10 summary data
        xlswrite(Data.OutputFileName,Stat.table,'Top10','O1'); % export statistics table from top 10 data
        xlswrite(Data.OutputFileName,Export.Z_score,'Unique','D1'); % export statistics table from top 10 data
        xlswrite(Data.OutputFileName,Export.Unique_Zscore,'Unique','K1');% export statistics table from top 10 data
        xlswrite(Data.OutputFileName,Export.Unique_Zscore_chg,'Unique_pca','A1'); % export statistics table from top 10 data
        xlswrite(Data.OutputFileName,Export.Stat.coeff,'Unique_pca','N1'); % export statistics table from top 10 data
        xlswrite(Data.OutputFileName,Export.score1,'Unique_pca','Z1'); % export factor index for pc1
        xlswrite(Data.OutputFileName,Export.Summary_table2,'Summary','A1');% export statistics table from top 10 data

        toc
