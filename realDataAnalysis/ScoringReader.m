function hyp = ScoringReader(path)
%% REMlogic report reader
% this file reads .txt files with sleep macrostructure and microstructure
% annotations exported from REMlogic
% returns the vectors
% hyp___________containing the hypnogram evaluated for each 30 s epochs
% time_tot______containing the starting time in seconds of CAP phases A
% duration______containing the duration of each phase A in seconds
% type_ar_______containing the type of phase A (A1, A2 or A3)

% The format of the header is
% Sleep Stage	Position	Time [hh:mm:ss]	Event	Duration[s]	Location

% Developed by Sara Mariani, MSc (sara1.mariani@mail.polimi.it), under the
% supervision of Professor Anna M. Bianchi and Professor Sergio Cerutti,
% at the Bioengineering Department of Politecnico di Milano, Italy. 

%clear all
close all

start=0;
hyp= []; %hypnogram
h=[];    %hours
m=[];    %minutes
s=[];    %seconds
%[file, path]=uigetfile('*.txt','Select report with annotations');
%path(end)=[];
%cd (path)
ptr = fopen(path,'r'); %opens the .txt files and returns the id

display('Reading file...')

while 1
    tline = fgetl(ptr); 
   
    if ~ischar(tline),
        break,
    end
    if numel(tline)>10
    if tline(1:11)=='Sleep Stage'
       start=1;
       k=0; %hyp counter
       j=0; %CAP counter
    end
    if start==1

            p=strfind(tline,':');
            if numel(p)==0
                p=strfind(tline,'.');
            end
            if tline(p(2)+4:p(2)+5)=='SL'             %sleep stage: write on hyp
                k=k+1;
                if tline(p(2)+11)=='E',   %REM
                    hyp(k,1)=5;
                elseif tline(p(2)+11)=='T',  %MT
                    hyp(k,1)=7;
                    display('MT')
                else
                    hyp(k,1)=str2num(tline(p(2)+11)); %sleep stage 0 1 2 3 4
                end

                h(k)=str2num(tline(p(1)-2:p(1)-1));
                m(k)=str2num(tline(p(1)+1:p(1)+2));
                s(k)=str2num(tline(p(2)+1:p(2)+2));
                if h(k)<10
                    hyp(k,2)=(h(k)+24)*3600+m(k)*60+s(k);
                else
                    hyp(k,2)=h(k)*3600+m(k)*60+s(k);
                end
            elseif (tline(p(2)+4:p(2)+5)=='MC') %CAP A phase: write on time_tot, duration, type
                j=j+1;
                t=strfind(tline,'-');
                type_ar(j)=str2num(tline(t(1)+2));
                duration(j,1)=str2num(tline(t(1)+4:t(1)+5));

                hCAP(j)=str2num(tline(p(1)-2:p(1)-1));
                mCAP(j)=str2num(tline(p(1)+1:p(1)+2));
                sCAP(j)=str2num(tline(p(2)+1:p(2)+2));
                if hCAP(j)<10
                    timevector(j,1)=(hCAP(j)+24)*3600+mCAP(j)*60+sCAP(j);
                else
                    timevector(j,1)=hCAP(j)*3600+mCAP(j)*60+sCAP(j);
                end
            end
    end
    end
end

display('Processing the hypnogram')

hypo=hyp;
hyp2=[];
jump=[];
nj=[];
x=1;
for j=1:length(hyp)-1
    if hyp(j+1,2)-hyp(j,2)>30
        jump(x)=j;
        if hyp(j+1,2)-hyp(j,2)==60
            nj(x)=1;
        elseif hyp(j+1,2)-hyp(j,2)==90
            nj(x)=2;
        elseif hyp(j+1,2)-hyp(j,2)==120
            nj(x)=3;
        end
        x=x+1;
    end
end
for i=1:length(jump)
    if nj(i)==1
        hyp2=[hyp(1:jump(i),:);[hyp(jump(i),1) hyp(jump(i),2)+30]; hyp(jump(i)+1:end,:)]; 
        jump(i:end)=jump(i:end)+1;
    elseif nj(i)==2
        hyp2=[hyp(1:jump(i),:);[hyp(jump(i),1) hyp(jump(i),2)+30];[hyp(jump(i),1) hyp(jump(i),2)+60]; hyp(jump(i)+1:end,:)]; 
        jump(i:end)=jump(i:end)+2;
    else
        hyp2=[hyp(1:jump(i),:);[hyp(jump(i),1) hyp(jump(i),2)+30];[hyp(jump(i),1) hyp(jump(i),2)+60];[hyp(jump(i),1) hyp(jump(i),2)+90];hyp(jump(i)+1:end,:)]; 
        jump(i:end)=jump(i:end)+3;
    end
        hyp=hyp2;
end
di=diff(hyp(:,2));
if di==30*ones(length(hyp)-1,1)
    display('check completed')
else
    display('error')
end
fclose(ptr); 
%cd ..
%display('saving')
start_time.h=h(1);
start_time.m=m(1);
start_time.s=s(1);
timestart=hyp(1,2);
time_tot=timevector-timestart;
hyp(:,2)=hyp(:,2)-hyp(1,2);
%eval (['save micro_str' file(1:end-4) ' time_tot duration type_ar start_time'])
%eval (['save hyp' file(1:end-4) ' hyp'])