clear all
close all
clc

%% SPECIFY EXCEL FILE 
filename = '1- VivoSense Master Sheet -6March2017.xlsx';

dataDirName = '_ALL_DATA_VIVO_CALIBRATED\'; % directory name where output files will be dumped
mkdir(dataDirName);


%% FILE NAME WHERE THE PROCESSED DATA WILL BE WRITTEN
FILENAME = {'AG-008_VIVO_CALIBRATED.dat','AJ-017_VIVO_CALIBRATED.dat','CB-013_VIVO_CALIBRATED.dat','DB-012_VIVO_CALIBRATED.dat','EK-016_VIVO_CALIBRATED.dat',...
    'GC-007_VIVO_CALIBRATED.dat','JI-004_VIVO_CALIBRATED.dat','JL-011_VIVO_CALIBRATED.dat','JM-009_VIVO_CALIBRATED.dat','KJ-010_VIVO_CALIBRATED.dat','LC-003_VIVO_CALIBRATED.dat',...
    'ML-006_VIVO_CALIBRATED.dat','MR-002_VIVO_CALIBRATED.dat','PB-014_VIVO_CALIBRATED.dat','PM-015_VIVO_CALIBRATED.dat'};


for sheet = 1:15
    sheet
    [num,txt,raw] = xlsread(filename,sheet);
    
    
    %% EXTRACT TIME AND Y-DATA FROM EXCEL SHEET
    
    % FIND WHERE THE TIME COLUMN IS LOCATED IN THE EXCEL SHEET
    TF=strcmp(raw, raw(2,2)); % RAW(2,2)='Time (sec)'
    [row,col] = find(TF);
    colvivo_time = [col' col(end)+7]'; % PUT ALL TIME STAMPS IN AN ARRAY
    
    % FIND WHERE THE Y DATA IS LOCATED IN THE EXCEL SHEET
    % Y data is located in between time stamps
    n=8;
    for i=1:length(colvivo_time)-1
        colvivo_y((i-1)*n+1:n*i) = linspace(colvivo_time(i),colvivo_time(i+1),n);
    end
    % remove time stamps from the data columns
    colvivo_y = setdiff(colvivo_y,colvivo_time);
    
    % extract data and time stamps
    tvivo=raw(:,colvivo_time(1:end-1));
    yvivo = raw(:,colvivo_y);
    
    % PROCESS TIME AND DATA
    for i=1:size(tvivo,2) % extract data for each block (seated, suppine, etc...)
        i
        % Change time to a running time (running cumulative sum)
        t_tmp = tvivo(:,i);
        fh = @(x) all(isnan(x(:)));
        t_tmp(cellfun(fh, t_tmp)) = [];
        t_tmp = t_tmp(3:end);
        time = cell2mat(t_tmp);
        time = [0; diff(time )*(24*60*60)];
        time =sum(triu(repmat(time',[prod(size(time')) 1])'))';       
        %data
        n=6;
        k=1;
        for j=(i-1)*n+1:n*i
            y_tmp = yvivo(:,j);
            fh = @(x) all(isnan(x(:)));
            y_tmp(cellfun(fh, y_tmp)) = [];
            y_tmp = y_tmp(2:end);
            y(:,k) = cell2mat(y_tmp);
            k=k+1;
        end
        
        % COLLECT DATA IN CELL ARRAY
        data2file{:,i}=[time y];
        
        %PLOT DATA
        hFig = figure;
        plot(time,y,'k-s','LineWidth',1,'MarkerEdgeColor','k',...
            'MarkerFaceColor','k',...
            'MarkerSize',3)
        xlabel('time[s]');
        set(hFig, 'Position', [300 300 1000 400])
        
        clear y;
        clear time;
    end
    
    % PREPARE data2file TO BE WRITTEN INTO AN EXTERNAL *.DAT FILE
    % Data varies in length from block to block (seated vs suppine vs etc...)
    % The idea is to write the whole data sets (all blocks) into one huge cell
    % array.
    
    intl = cellfun( @size, data2file, 'uni',false) ;
    for i=1:length(intl)
        rownum(i) = intl{i}(1) ;
        colnum(i) = intl{i}(2) ;
    end
    nRow = max(rownum) ;
    nCol = sum(colnum) ;
    % fill out a cell array
    content = cell( nRow, nCol ) ;
    for cId = 1 : length(data2file)
        ne = size( data2file{cId} ) ;
        content(1:ne(1),(cId-1)*ne(2)+1:ne(2)*cId) = num2cell(data2file{cId}(1:ne(1),1:ne(2))) ;
    end
    % convert into a cell array of printed (to string) content.
    contentStr = cellfun( @(x)sprintf('%f', x), content,'UniformOutput', false ) ;
    % Then we get the max string length per column. Note that we could skip this
    %operation if we wanted to enforce a fixed length format
    strlen = max( cellfun( @length, contentStr )) ;
    
    % find empty fields and replace them with a 'filler string'
    indempt=find(cellfun(@isempty,contentStr));
    contentStr(indempt)={'-999999999'};

    
    % WRITE DATA TO FILE
    fileID = fopen(strcat(dataDirName,FILENAME{sheet}),'w');
    format = sprintf( '%%%ds    ', strlen ) ;
    for rId = 1 : size( contentStr, 1 )
        fprintf( fileID, format, contentStr{rId,:} ) ;
        fprintf( fileID, '\r\n' ) ;
    end
    fclose( fileID ) ;
    
    clear  data2file
    close all
    
end