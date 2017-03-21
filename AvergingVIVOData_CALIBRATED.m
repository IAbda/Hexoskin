clear all
close all
clc


%%-------------------
%%INPUT

% SPECIFY AVERAGING WINDOW SIZE
windowSize_seconds = 60;

dataDirName = '_ALL_DATA_VIVO_CALIBRATED\'; % specify directory where extracted vivo is located
writeDataTo = strcat('_INDVS_DATA_VIVO_CALIBRATED_',num2str(windowSize_seconds),'s\');
mkdir(writeDataTo);


%%-------------------
ALL_METRICS = {'VENTILATION', 'VT', 'BREATHING_RATE', 'HEART_RATE', 'INSP_VOLUME', 'EXP_VOLUME'};

ALL_CONDITIONS = {'SEATED', 'SUPINE', 'STANDING', '6MIN_BASELINE', '3MST', 'SPIRO1', 'SPIRO2', 'SPIRO3'};



%%
col_time = [1 8 15 22 29 36 43 50];

for kkk = 1:length(ALL_METRICS)
    
    METRIC = ALL_METRICS{kkk}
    
    % DATA COLUMN NUMBER FROM FILES
    switch METRIC
        case 'VENTILATION'
            col_num = [2     9    16    23    30    37    44    51];
        case 'VT'
            col_num = [3    10    17    24    31    38    45    52];
        case 'BREATHING_RATE'
            col_num = [4    11    18    25    32    39    46    53]; % 4=breathing rate seated, 9=breathing rate supine, etc.
        case 'HEART_RATE'
            col_num = [5    12    19    26    33    40    47    54];
        case 'INSP_VOLUME'
            col_num = [6    13    20    27    34    41    48    55];
        case 'EXP_VOLUME'
            col_num = [7    14    21    28    35    42    49    56];
    end
    
    for jj=1:length(ALL_CONDITIONS)
        
        CONDITION = ALL_CONDITIONS{jj}
        
        col_data = col_num(jj);
        
        
        %--------------------------------------------------------------------------
        % Read interpolated data from text file
        %--------------------------------------------------------------------------
        switch CONDITION
            case {'SEATED', 'SUPINE', 'STANDING'}
                fileName_VIVO = {'AG-008_VIVO_CALIBRATED.dat','AJ-017_VIVO_CALIBRATED.dat','CB-013_VIVO_CALIBRATED.dat','DB-012_VIVO_CALIBRATED.dat','EK-016_VIVO_CALIBRATED.dat',...
                    'GC-007_VIVO_CALIBRATED.dat','JI-004_VIVO_CALIBRATED.dat','JL-011_VIVO_CALIBRATED.dat','JM-009_VIVO_CALIBRATED.dat','KJ-010_VIVO_CALIBRATED.dat','LC-003_VIVO_CALIBRATED.dat',...
                    'ML-006_VIVO_CALIBRATED.dat','MR-002_VIVO_CALIBRATED.dat','PB-014_VIVO_CALIBRATED.dat','PM-015_VIVO_CALIBRATED.dat'};
                
                %% DATA FILES FOR 6 MINUTE BASELINE, 3MST, PM-015.dat IS NOT INCLUDED HERE
            case {'6MIN_BASELINE'}
                fileName_VIVO = {'AG-008_VIVO_CALIBRATED.dat','AJ-017_VIVO_CALIBRATED.dat','CB-013_VIVO_CALIBRATED.dat','DB-012_VIVO_CALIBRATED.dat','EK-016_VIVO_CALIBRATED.dat',...
                    'GC-007_VIVO_CALIBRATED.dat','JI-004_VIVO_CALIBRATED.dat','JL-011_VIVO_CALIBRATED.dat','JM-009_VIVO_CALIBRATED.dat','KJ-010_VIVO_CALIBRATED.dat','LC-003_VIVO_CALIBRATED.dat',...
                    'ML-006_VIVO_CALIBRATED.dat','MR-002_VIVO_CALIBRATED.dat','PB-014_VIVO_CALIBRATED.dat'};
            case {'3MST'}
                fileName_VIVO = {'AG-008_VIVO_CALIBRATED.dat','AJ-017_VIVO_CALIBRATED.dat','CB-013_VIVO_CALIBRATED.dat','DB-012_VIVO_CALIBRATED.dat','EK-016_VIVO_CALIBRATED.dat',...
                    'GC-007_VIVO_CALIBRATED.dat','JI-004_VIVO_CALIBRATED.dat','JL-011_VIVO_CALIBRATED.dat','JM-009_VIVO_CALIBRATED.dat','KJ-010_VIVO_CALIBRATED.dat',...
                    'ML-006_VIVO_CALIBRATED.dat','MR-002_VIVO_CALIBRATED.dat','PB-014_VIVO_CALIBRATED.dat'};
                %% DATA FILES FOR SPIRO2, EK-016.DATA IS NOT INCLUDED HERE
                %% DATA FILES FOR SPIRO2, SEVERAL ARE NOT INCLUDED HERE
            case {'SPIRO1','SPIRO2','SPIRO3'}
                fileName_VIVO = {'AG-008_VIVO_CALIBRATED.dat','AJ-017_VIVO_CALIBRATED.dat','CB-013_VIVO_CALIBRATED.dat','DB-012_VIVO_CALIBRATED.dat',...
                    'GC-007_VIVO_CALIBRATED.dat','JL-011_VIVO_CALIBRATED.dat','JM-009_VIVO_CALIBRATED.dat','KJ-010_VIVO_CALIBRATED.dat',...
                    'PB-014_VIVO_CALIBRATED.dat','PM-015_VIVO_CALIBRATED.dat'};
        end
        
        
        for i=1:length(fileName_VIVO)
            fileName_VIVO{i}
            
            
            %--------------------------------------------------------------------------
            % Read data files
            %--------------------------------------------------------------------------
            
            
            M_VIVO = dlmread(strcat(dataDirName,fileName_VIVO{i}));
            
            Y_VIVO{:,i}=M_VIVO(:,col_data);
            indneg1=Y_VIVO{:,i}<=-999999999;     % filter -99999999 data from file
            indneg2=find(~indneg1);
            Y_VIVO{:,i}=Y_VIVO{:,i}(indneg2);
            
            time_VIVO{:,i} = M_VIVO(:,col_time(jj));
            indneg3=time_VIVO{:,i}<0;
            indneg4=find(~indneg3);
            time_VIVO{:,i}=time_VIVO{:,i}(indneg4);
            
            % COLLECT DATA IN CELL ARRAY
            data2file_VIVO{:,i}=[time_VIVO{:,i} Y_VIVO{:,i}];
            
            
            % MAKE MOVING AVERAGE OF REF DATA FOR EVERY N-SECONDS
            modtt = mod(time_VIVO{:,i},windowSize_seconds);
            [B,IX]=sort(modtt);
            ttstop = sort(time_VIVO{:,i}(IX(2:ceil(time_VIVO{:,i}(end)/windowSize_seconds))));
            if isempty(ttstop)
                Y_VIVO_avg_mov_Y_VIVO{:,i}(1,:)=mean(Y_VIVO{:,i}(1,:));
            else
                indttstop=[0;find(ismember(time_VIVO{:,i}(:),ttstop))];
                indttstop1=indttstop+1;
                for cc=1:length(indttstop)-1
                    Y_VIVO_avg_mov_Y_VIVO{:,i}(cc,:) = mean(Y_VIVO{:,i}(indttstop1(cc):1:indttstop(cc+1)));
                end
                Y_VIVO_avg_mov_Y_VIVO{:,i}(end+1,:)=mean(Y_VIVO{:,i}(indttstop1(end):end));
            end
            
            if strcmp(CONDITION,'SEATED') | strcmp(CONDITION,'SUPINE')
                % AVERAGE ALL DATA FROM 60S TO 240S AND EXTRACT STANDARD ERROR
                % AVERAGE ALL DATA FROM 180S TO 240S AND EXTRACT STANDARD ERROR
                [c index60] = min(abs(time_VIVO{:,i}-60));
                [c index180] = min(abs(time_VIVO{:,i}-180));
                [c index240] = min(abs(time_VIVO{:,i}-240));
                % mean
                Y_VIVO_avg_all_60_240{i}(1,:)=mean(Y_VIVO{:,i}(index60:index240));
                % standard deviation
                Y_VIVO_avg_all_60_240{i}(2,:)=std(Y_VIVO{:,i}(index60:index240));
                % standard error
                Y_VIVO_avg_all_60_240{i}(3,:)=std(Y_VIVO{:,i}(index60:index240))/sqrt(length(Y_VIVO{:,i}(index60:index240)));
                % coefficient of variation
                Y_VIVO_avg_all_60_240{i}(4,:)=100*std(Y_VIVO{:,i}(index60:index240))/mean(Y_VIVO{:,i}(index60:index240));
                % mean
                Y_VIVO_avg_all_60_240{i}(5,:)=mean(Y_VIVO{:,i}(index180:index240));
                % standard deviation
                Y_VIVO_avg_all_60_240{i}(6,:)=std(Y_VIVO{:,i}(index180:index240));
                % standard error
                Y_VIVO_avg_all_60_240{i}(7,:)=std(Y_VIVO{:,i}(index180:index240))/sqrt(length(Y_VIVO{:,i}(index180:index240)));
                % coefficient of variation
                Y_VIVO_avg_all_60_240{i}(8,:)=100*std(Y_VIVO{:,i}(index180:index240))/mean(Y_VIVO{:,i}(index180:index240));
                
            end
            if strcmp(CONDITION,'STANDING') | strcmp(CONDITION,'3MST')
                % AVERAGE ALL DATA FROM 60S TO 240S AND EXTRACT STANDARD ERROR
                [c index60] = min(abs(time_VIVO{:,i}-60));
                [c index120] = min(abs(time_VIVO{:,i}-120));
                % mean
                Y_VIVO_avg_all_60_120{i}(1,:)=mean(Y_VIVO{:,i}(index60:index120));
                % standard deviation
                Y_VIVO_avg_all_60_120{i}(2,:)=std(Y_VIVO{:,i}(index60:index120));
                % standard error
                Y_VIVO_avg_all_60_120{i}(3,:)=std(Y_VIVO{:,i}(index60:index120))/sqrt(length(Y_VIVO{:,i}(index60:index120)));
                % coefficient of variation
                Y_VIVO_avg_all_60_120{i}(4,:)=100*std(Y_VIVO{:,i}(index60:index120))/mean(Y_VIVO{:,i}(index60:index120));
                
            end
            
            
            %             leftover=mod(length(Y_VIVO{:,i}),windowSize_seconds);
            %             Y_VIVO_reshape=reshape(Y_VIVO{:,i}(1:end-leftover),windowSize_seconds,[]);
            %             if strcmp(CONDITION, 'SEATED')  | strcmp(CONDITION, 'SUPINE')
            %                 pairedVIVO(:,i)=reshape(Y_VIVO_reshape(:,2:4),size(Y_VIVO_reshape(:,2:4),1)*size(Y_VIVO_reshape(:,2:4),2),1);
            %             end
            %             if strcmp(CONDITION, 'STANDING')
            %                 pairedVIVO(:,i)=reshape(Y_VIVO_reshape(:,2),size(Y_VIVO_reshape(:,2),1)*size(Y_VIVO_reshape(:,2),2),1);
            %             end
            %             if strcmp(CONDITION, '3MST')
            %                 pairedVIVO(:,i)=reshape(Y_VIVO_reshape(:,2),size(Y_VIVO_reshape(:,2),1)*size(Y_VIVO_reshape(:,2),2),1);
            %             end
            
            %clc
        end
        
        
        
        %         if strcmp(CONDITION, 'SEATED')  | strcmp(CONDITION, 'SUPINE') | strcmp(CONDITION, 'STANDING') | strcmp(CONDITION, '3MST')
        %             % COMPUTE THE COV
        %             for i=1:length(fileName_ref)
        %                 cov_VIVO(i) = std(pairedVIVO(:,i))/mean(pairedVIVO(:,i));
        %             end
        %             fprintf('cov_VIVO\n');
        %             fprintf('%f\n',cov_VIVO);
        %
        %             % COMPUTE THE INTER-CLASS CORRELATION COEFFICIENT
        %             %   from the data matrix X. X is a n-by-w matrix with colums corresponding
        %             %   to persons and rows corresponding to measues, respectively.
        %             ICC_VIVO = ICC(pairedVIVO','1-k');
        %             fprintf('ICC_VIVO: %f\n',ICC_VIVO);
        %             fprintf('\n');
        %         end
        
        
        
        
        
        
        
        
        %% WRITE ALL 60-240 and 180-240 AVERAGED Y_REF DATA TO FILE
        if strcmp(CONDITION,'SEATED') | strcmp(CONDITION,'SUPINE')
            
            intl = cellfun( @size, Y_VIVO_avg_all_60_240, 'uni',false) ;
            for i=1:length(intl)
                rownum(i) = intl{i}(1) ;
                colnum(i) = intl{i}(2) ;
            end
            nRow = max(rownum) ;
            nCol = sum(colnum) ;
            % fill out a cell array
            content_Y_VIVO_avg_all_60_240 = cell( nRow, nCol ) ;
            for cId = 1 : length(Y_VIVO_avg_all_60_240)
                ne = size( Y_VIVO_avg_all_60_240{cId} ) ;
                content_Y_VIVO_avg_all_60_240(1:ne(1),(cId-1)*ne(2)+1:ne(2)*cId) = num2cell(Y_VIVO_avg_all_60_240{cId}(1:ne(1),1:ne(2))) ;
            end
            % convert into a cell array of printed (to string) content.
            contentStr_Y_VIVO_avg_all_60_240 = cellfun( @(x)sprintf('%f', x), content_Y_VIVO_avg_all_60_240,'UniformOutput', false ) ;
            % Then we get the max string length per column. Note that we could skip this
            %operation if we wanted to enforce a fixed length format
            strlen = max( cellfun( @length, contentStr_Y_VIVO_avg_all_60_240 )) ;
            
            % FILL OUT EMPTY CELLS WITH -999999
            iidEMP = cellfun( @isempty, contentStr_Y_VIVO_avg_all_60_240 );
            contentStr_Y_VIVO_avg_all_60_240(iidEMP>0) = {'-999999'};
            
            % WRITE DATA TO FILE
            
            fileID = fopen(strcat(writeDataTo,'_All_Indivs_',CONDITION,'_',METRIC,'_AVG_60-240_VIVO_CALIBRATED.dat'),'w');
            format = sprintf( '%%%ds    ', strlen ) ;
            for rId = 1 : size( contentStr_Y_VIVO_avg_all_60_240, 1 )
                fprintf( fileID, format, contentStr_Y_VIVO_avg_all_60_240{rId,:} ) ;
                fprintf( fileID, '\r\n' ) ;
            end
            fclose( fileID ) ;
        end
        
        
        
        
        
        %% WRITE ALL 60-120s AVERAGED Y_VIVO DATA TO FILE
        if strcmp(CONDITION,'STANDING') | strcmp(CONDITION,'3MST')
            
            intl = cellfun( @size, Y_VIVO_avg_all_60_120, 'uni',false) ;
            for i=1:length(intl)
                rownum(i) = intl{i}(1) ;
                colnum(i) = intl{i}(2) ;
            end
            nRow = max(rownum) ;
            nCol = sum(colnum) ;
            % fill out a cell array
            content_Y_VIVO_avg_all_60_120 = cell( nRow, nCol ) ;
            for cId = 1 : length(Y_VIVO_avg_all_60_120)
                ne = size( Y_VIVO_avg_all_60_120{cId} ) ;
                content_Y_VIVO_avg_all_60_120(1:ne(1),(cId-1)*ne(2)+1:ne(2)*cId) = num2cell(Y_VIVO_avg_all_60_120{cId}(1:ne(1),1:ne(2))) ;
            end
            % convert into a cell array of printed (to string) content.
            contentStr_Y_VIVO_avg_all_60_120 = cellfun( @(x)sprintf('%f', x), content_Y_VIVO_avg_all_60_120,'UniformOutput', false ) ;
            % Then we get the max string length per column. Note that we could skip this
            %operation if we wanted to enforce a fixed length format
            strlen = max( cellfun( @length, contentStr_Y_VIVO_avg_all_60_120 )) ;
            
            % FILL OUT EMPTY CELLS WITH -999999
            iidEMP = cellfun( @isempty, contentStr_Y_VIVO_avg_all_60_120 );
            contentStr_Y_VIVO_avg_all_60_120(iidEMP>0) = {'-999999'};
            
            % WRITE DATA TO FILE
            
            fileID = fopen(strcat(writeDataTo,'_All_Indivs_',CONDITION,'_',METRIC,'_AVG_60-120_VIVO_CALIBRATED.dat'),'w');
            format = sprintf( '%%%ds    ', strlen ) ;
            for rId = 1 : size( contentStr_Y_VIVO_avg_all_60_120, 1 )
                fprintf( fileID, format, contentStr_Y_VIVO_avg_all_60_120{rId,:} ) ;
                fprintf( fileID, '\r\n' ) ;
            end
            fclose( fileID ) ;
        end
        
        
        
        
        %% WRITE ALL INDIVIDUALS AVERAGED Y_VIVO DATA TO FILE
        % PREPARE data2file TO BE WRITTEN INTO AN EXTERNAL *.DAT FILE
        % Data varies in length from block to block (seated vs suppine vs etc...)
        % The idea is to write the whole data sets (all blocks) into one huge cell
        % array.
        
        intl = cellfun( @size, Y_VIVO_avg_mov_Y_VIVO, 'uni',false) ;
        for i=1:length(intl)
            rownum(i) = intl{i}(1) ;
            colnum(i) = intl{i}(2) ;
        end
        nRow = max(rownum) ;
        nCol = sum(colnum) ;
        % fill out a cell array
        content_Y_VIVO_avg_mov_Y_VIVO = cell( nRow, nCol ) ;
        for cId = 1 : length(Y_VIVO_avg_mov_Y_VIVO)
            ne = size( Y_VIVO_avg_mov_Y_VIVO{cId} ) ;
            content_Y_VIVO_avg_mov_Y_VIVO(1:ne(1),(cId-1)*ne(2)+1:ne(2)*cId) = num2cell(Y_VIVO_avg_mov_Y_VIVO{cId}(1:ne(1),1:ne(2))) ;
        end
        % convert into a cell array of printed (to string) content.
        contentStr_Y_VIVO_avg_mov_Y_VIVO = cellfun( @(x)sprintf('%f', x), content_Y_VIVO_avg_mov_Y_VIVO,'UniformOutput', false ) ;
        % Then we get the max string length per column. Note that we could skip this
        %operation if we wanted to enforce a fixed length format
        strlen = max( cellfun( @length, contentStr_Y_VIVO_avg_mov_Y_VIVO )) ;
        
        
        % FILL OUT EMPTY CELLS WITH -999999
        iidEMP = cellfun( @isempty, contentStr_Y_VIVO_avg_mov_Y_VIVO );
        contentStr_Y_VIVO_avg_mov_Y_VIVO(iidEMP>0) = {'-999999'};
        
        % WRITE DATA TO FILE
        
        fileID = fopen(strcat(writeDataTo,'_All_Indivs_',CONDITION,'_',METRIC,'_',num2str(windowSize_seconds),'sAVG_VIVO_CALIBRATED.dat'),'w');
        format = sprintf( '%%%ds    ', strlen ) ;
        for rId = 1 : size( contentStr_Y_VIVO_avg_mov_Y_VIVO, 1 )
            fprintf( fileID, format, contentStr_Y_VIVO_avg_mov_Y_VIVO{rId,:} ) ;
            fprintf( fileID, '\r\n' ) ;
        end
        fclose( fileID ) ;
        
        
        
        
        %% WRITE ALL INDIVIDUALS Y_VIVO DATA TO FILE
        % PREPARE data2file TO BE WRITTEN INTO AN EXTERNAL *.DAT FILE
        % Data varies in length from block to block (seated vs suppine vs etc...)
        % The idea is to write the whole data sets (all blocks) into one huge cell
        % array.
        
        intl = cellfun( @size, data2file_VIVO, 'uni',false) ;
        for i=1:length(intl)
            rownum(i) = intl{i}(1) ;
            colnum(i) = intl{i}(2) ;
        end
        nRow = max(rownum) ;
        nCol = sum(colnum) ;
        % fill out a cell array
        content_VIVO = cell( nRow, nCol ) ;
        for cId = 1 : length(data2file_VIVO)
            ne = size( data2file_VIVO{cId} ) ;
            content_VIVO(1:ne(1),(cId-1)*ne(2)+1:ne(2)*cId) = num2cell(data2file_VIVO{cId}(1:ne(1),1:ne(2))) ;
        end
        % convert into a cell array of printed (to string) content.
        contentStr = cellfun( @(x)sprintf('%f', x), content_VIVO,'UniformOutput', false ) ;
        % Then we get the max string length per column. Note that we could skip this
        %operation if we wanted to enforce a fixed length format
        strlen = max( cellfun( @length, contentStr )) ;
        
        % FILL OUT EMPTY CELLS WITH -999999
        iidEMP = cellfun( @isempty, contentStr );
        contentStr(iidEMP>0) = {'-999999'};
        
        % WRITE DATA TO FILE
        fileID = fopen(strcat(writeDataTo,'_All_Indivs_',CONDITION,'_',METRIC,'_VIVO_CALIBRATED.dat'),'w');
        format = sprintf( '%%%ds    ', strlen ) ;
        for rId = 1 : size( contentStr, 1 )
            fprintf( fileID, format, contentStr{rId,:} ) ;
            fprintf( fileID, '\r\n' ) ;
        end
        fclose( fileID ) ;
        
        
        
        %%
        
        clearvars -except writeDataTo dataDirName kkk jj METRIC col_time col_num windowSize_seconds ALL_METRICS ALL_CONDITIONS
        
    end
end
