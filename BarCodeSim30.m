tic
%{
BarCode sim for Dr. Jongmin Nam
Sean McQuade and Jongmin Nam 1.08.2015 to 5.26.2015
using Sea Urchin Lineage from "Development of Sea Urchins, Ascidians, and
other invertebrate Deuterostomes: Experimental Approaches"
Authors: Charles A. Ettensohn, Gary M. Wessel, Gregory A. Wray
Elsevier academic press
 
% i=cellgroup integrated; 
% j=cell in group; 
% int = number of integrated cells
% exp = number of expressed cells
%}

%All oblique cell cleavages are coded as horizontal cleavages!

%Total_Targets is a cell array with all targets to be compared
%the simulation loops through this cell array
Total_Targets{1} = [51,56];
% Total_Targets{2} = 51;
% Total_Targets{3} = 52;
% Total_Targets{4} = [51,52];

for a = 1:size(Total_Targets,2)
global Number_of_Iterations
%parameters
Number_of_Trials = 1;
Number_of_Iterations = 1000;
percent_trim = 0.01;

%initialize
C{64}=[];          %preallocate memory for Cell array
Integration_Stage = zeros(Number_of_Iterations,1);
Target1 = 51;       %cells in which integration provides transcription
Target2 = 52;
Integrated_Cells = zeros(Number_of_Iterations,1);
Expressed_Cells = zeros(Number_of_Iterations,1);
Ratio_Exp_to_Int = zeros(Number_of_Iterations,1);
A_plot = [0, 0, 0]; %place holder
A_raw = [0, 0, 0]; %place holder
% one_cell_data = [one_cell_data; CRM_raw(i,:)];
% two_cell_data = [two_cell_data; CRM_raw(i,:)];
% four_cell_data = [four_cell_data; CRM_raw(i,:)];
% eight_cell_data = [eight_cell_data; CRM_raw(i,:)];
% sixteen_cell_data = [sixteen_cell_data; CRM_raw(i,:)];
% thirtytwo_cell_data = [thirtytwo_cell_data; CRM_raw(i,:)];
% sixtyfour_cell_data = [four_cell_data; CRM_raw(i,:)];
% one_twentyeight_cell_data = [one_twentyeight_cell_data; CRM_raw(i,:)];
% two_fiftysix_cell_data = [two_fiftysix_cell_data; CRM_raw(i,:)];

%save file name
path = sprintf('Sim26_target_%s.txt',mat2str(Total_Targets{a}));

for b = 1:Number_of_Trials       
for c = 1:Number_of_Iterations     
%{  
The probability of DNA integration occuring in a stage is modeled on 
empiracally obtained profile from Dr. Nam. For a given stage, the probability 
of DNA integration into the ith group is proportional to the number of 
cells in the ith group.
%}

% create integration profile (this must change for oblique clev

h = randi(95148);
if h == 1 %0th stage, yields ~424 integrated cells
    i=1;
    m=1;
elseif h >= 2 && h<=31 %Stage I, yields ~212 integrated cells
    i=2;
    m=2;
elseif h>=32 && h<=958   %Stage II, yields ~106 integrated cells
    i=3;      
    m=3;
elseif h>=959 && h<=9589 %Stage III, yields ~64 integrated cells
    i=4;
    m=4;
elseif h>=9590 && h<=28880  %Stage IV, yields ~32-42 integrated cells
    i=randi([5,8]);
    m=5;
elseif h>=28881 && h<=46066 %Stage V, yields ~16 integrated cells 
    p=randi(6);
    if p==1
        i=10;
    elseif p==2
        i=11;
    elseif p==3 || p==4
        i = 12;
    elseif p==5 || p==6
        i=13;
    end
    m=6;
elseif h>=46067 && h<=58496 %Stage VI, yields ~8-10 integrated cells
    p=randi(14);
    if p==1
        i=9;
    elseif p==2
        i=14;
    elseif p==3
        i = 16;
    elseif p==4
        i=17;
    elseif p==5 || p==6
        i=18;
    elseif p>=7 && p<=10
        i=19;
    elseif p==11 || p==12
        i=20;
    elseif p==13 || p==14
        i=21;
    end
    m=7;
elseif h>=58497 && h<=68321 %Stage VII, yields ~4 integrated cells
    p=randi(26);
    if p==1
        i=22;
    elseif p==2
        i=23;
    elseif p==3 || p==4
        i=24;
    elseif p==5 || p==6
        i=25;
    elseif p>=7 && p<=10
        i=26;
    elseif p>=11 && p<=14
        i=27;
    elseif p>=15 && p<=18
        i=28;
    elseif p>=19 && p<=22
        i=29;
    elseif p>=23 && p<=26
        i=30;
    end
    m=8;
elseif h>=68322 && h<75147  %Stage VIII, yields ~2 integrated cells
    p=randi(53);

    if p == 1
        i = 15;
    elseif p>=2 && p<=5
        i=randi([31,34]);
    elseif p==6 || p==7
        i = 35;
    elseif p==8 || p==9
        i = 36;
    elseif p>=10 && p<=13
        i = 37;
    elseif p>=14 && p<=37
        i=randi([38,40]);
    elseif p>=38 && p<=53
        i=randi([41,44]);
    end
    m=9;
elseif h>=75148 && h<=95148 %Stage IX, yields 1 integrated cell
    i=randi([45,64]);
    m=10;
end

%make vector for hist: Integration_Stage
Integration_Stage(c)=m;

%  Cell in group C{i} chosen uniformly.
if i == 1
    DNA_integrated_cell = 1;
elseif i == 2
    DNA_integrated_cell = randi(2);
elseif i>=3 && i<=11 || i>=14 && i<=17 || i>=31 && i<=34 || i>=45 && i<=46 || ...
        i>=22 && i<=23
    DNA_integrated_cell = randi(4);
elseif i>=12 && i<=13 || i==18 || i>=20 && i<=21 || i>=24 && i<=25 || ... 
        i>=35 && i<=36 || i>=47 && i<=48 || i>=61 
    DNA_integrated_cell = randi(8);
elseif i==19 || i>=26 && i<=30 || i>=41 && i<=44 || i==49 
    DNA_integrated_cell = randi(16);
elseif i>=37 && i<=40 || i>=50 && i<=60
    DNA_integrated_cell = randi(32);
end

%the DNA is integrated in the jth cell in the ith group.
j = DNA_integrated_cell;
    
Cells0 = zeros(1,1); 

C{1}= Cells0;
C{i}(1)=1;
C{2}= cat(2,C{1},C{1});              
C{i}(j)=1;
C{3}= cat(2,C{2},C{2});              
C{i}(j)=1;

C{4}= C{3}; 
C{5}= C{3}; 

C{i}(j)=1;
C{6}= C{4};                                                              %4
C{7}= C{4};                                                              %4  
C{8}= C{5};                                                              %4
C{9}= C{5};                                                              %4

C{i}(j)=1; 
C{10} = C{6};  
C{11} = C{6};  
C{12} = cat(2,C{7},C{7});       
C{13} = cat(2,C{8},C{8});         
C{14} = C{9};  
C{15} = C{9};  

C{i}(j)=1; 
C{16} = C{10};                    
C{17} = C{10};                    
C{18} = cat(2,C{11},C{11});     
C{19} = cat(2,C{12},C{12});      
C{20} = C{13};                  
C{21} = C{13};                   
C{22} = C{14};                        
C{23} = C{14};                        

%116 Total Cells   

C{i}(j)=1;                                                   
C{24}= cat(2,C{16},C{16});
C{25}= cat(2,C{17},C{17});
C{26}= cat(2,C{18},C{18});
C{27}= C{19};
C{28}= C{19};
C{29}=  cat(2,C{20},C{20});
C{30}= cat(2,C{21},C{21});
C{31}= C{22};
C{32}= C{22};
C{33}= C{23};
C{34}= C{23};

%232 total cells

C{i}(j)=1; 
C{35}= C{24};
C{36}= C{24};
C{37}=  cat(2,C{25},C{25});
C{38}= cat(2,C{26},C{26});
C{39}= cat(2,C{27},C{27});
C{40}= cat(2,C{28},C{28});
C{41}= C{29};
C{42}= C{29};
C{43}= C{30};
C{44}= C{30};
C{45}= C{15};
C{46}= C{15};

%424 total cells

C{i}(j)=1; 
C{47}= C{35};
C{48}= C{35};
C{49}= cat(2,C{36},C{36});
C{50}= cat(2,C{37},C{37});
C{51}= C{38};
C{52}= C{38};
C{53}= C{39};
C{54}= C{39};
C{55}= C{40};
C{56}= C{40};
C{57}= cat(2,C{41},C{41});
C{58}= cat(2,C{42},C{42});
C{59}= cat(2,C{43},C{43});
C{60}= cat(2,C{44},C{44}); 
C{61}= cat(2,C{31},C{31});
C{62}= cat(2,C{32},C{32});
C{63}= cat(2,C{33},C{33});
C{64}= cat(2,C{34},C{34});

C{i}(j)=1;

%Count integrated cells within specified "target" groups
Expressed_Cells = 0; 
for d=1:size(Total_Targets{a},2)
Expressed_Cells = Expressed_Cells + nnz(C{Total_Targets{a}(d)});
end

% half1 = round(size(Target1)/2);
% half2 = round(size(Target2)/2);
% %Expressed_Cells(c) = nnz(C{Target1})+ nnz(C{Target2});
% Expressed_Cells(c) = nnz(C{Target1}(half1:end)) + nnz(C{Target2}(1:half2));

%Counts integrated cells in all stage 9 groups(also C{45} and C{46} (SM))
Integrated_Cells  = nnz(C{47}) + nnz(C{48}) + nnz(C{49}) + ...
    nnz(C{50}) + nnz(C{51}) + nnz(C{52}) + nnz(C{53}) + nnz(C{54}) +...
    nnz(C{55}) + nnz(C{56}) + nnz(C{57}) + nnz(C{58}) + nnz(C{59}) +...
    nnz(C{60}) + nnz(C{61}) + nnz(C{62}) + nnz(C{63}) + nnz(C{64}) +...
    nnz(C{45}) + nnz(C{46}); %these last two terms add SM1 and SM2
%log the ratio of this particular iteration in the experiment
Ratio_Exp_to_Int = Expressed_Cells/Integrated_Cells;
A_raw = [A_raw; Expressed_Cells, Integrated_Cells, Ratio_Exp_to_Int];
if Integrated_Cells >= 4
    A_plot = [A_plot; Expressed_Cells, Integrated_Cells, Ratio_Exp_to_Int];
end
%1celldata, 2celldata, 4celldata, etc.
end
%Data is ordered by column 3 (expression ratio), then put into descending
B_raw = sortrows(A_raw,3); %orders array with respect to(WRT) col 3
D_raw = flip(B_raw,1) ;    %puts into descending order WRT col 3
D_raw(end,:) = [];     %removes the last row from array, initial zeros for A.
CRM_raw = D_raw;
avg_exp_raw = sum(CRM_raw(:,1))/sum(CRM_raw(:,2)); %from (eq2)



  
%write information to file
fileID = fopen(path,'a');
fprintf(fileID,'%s%s\t%s\t%i\t%s\t%i\t%s\t%f\n','#Target_',mat2str(Total_Targets{a}),'cDNAall'...
    ,sum(CRM_raw(:,1)),'gDNAall',sum(CRM_raw(:,2)),'c/g_all', avg_exp_raw);
fclose(fileID);
end

one_cell_data{a} = [];
two_cell_data{a} = [];
four_cell_data{a} = [];
eight_cell_data{a} = [];
sixteen_cell_data{a} = [];
thirtytwo_cell_data{a} = [];
fourtytwo_cell_data{a} = [];
sixtyfour_cell_data{a} = [];
one_hundred_six_cell_data{a} = [];
two_twelve_cell_data{a} = [];
four_twentyfour_cell_data{a} = [];

for i = 1:size(CRM_raw,1)
    if CRM_raw(i,2) == 1
        one_cell_data{a} = [one_cell_data{a}; CRM_raw(i,:)];
    elseif CRM_raw(i,2) == 2
        two_cell_data{a} = [two_cell_data{a}; CRM_raw(i,:)];
    elseif CRM_raw(i,2) == 4        
        four_cell_data{a} = [four_cell_data{a}; CRM_raw(i,:)];
    elseif CRM_raw(i,2) == 8
        eight_cell_data{a} = [eight_cell_data{a}; CRM_raw(i,:)];
    elseif CRM_raw(i,2) == 16  
        sixteen_cell_data{a} = [sixteen_cell_data{a}; CRM_raw(i,:)];
    elseif CRM_raw(i,2) == 32
        thirtytwo_cell_data{a} = [thirtytwo_cell_data{a}; CRM_raw(i,:)];
    elseif CRM_raw(i,2) == 42
        fourtytwo_cell_data{a} = [fourtytwo_cell_data{a}; CRM_raw(i,:)];
    elseif CRM_raw(i,2) == 64        
        sixtyfour_cell_data{a} = [sixtyfour_cell_data{a}; CRM_raw(i,:)];
    elseif CRM_raw(i,2) == 106
        one_hundred_six_cell_data{a} = [one_hundred_six_cell_data{a}; CRM_raw(i,:)];
    elseif CRM_raw(i,2) == 212        
        two_twelve_cell_data{a} = [two_twelve_cell_data{a}; CRM_raw(i,:)];
    elseif CRM_raw(i,2) == 424        
        four_twentyfour_cell_data{a} = [four_twentyfour_cell_data{a}; CRM_raw(i,:)];
    end
end

end
linestyle = ['-'];
for a = 1:size(Total_Targets,2)
hold on
%figure 1
if ~isempty(one_cell_data{a})
figure(1)
simple_vector = (1:size(one_cell_data{a},1));
title1 = sprintf('normalized 1 cell expression profile');
plot(simple_vector, one_cell_data{a}(:,3),'Color',[.3 0 .5],...
    'LineWidth',3,'LineStyle',linestyle(a))
xlabel('number of trials');
ylabel({'Rank-Ordered Expression Profile','(cDNA/gDNA)'});
legend('A')
title(title1)
axis auto
end
hold on
%figure 2
if ~isempty(two_cell_data{a})
figure(2)
simple_vector = (1:size(two_cell_data{a},1));
title2 = sprintf('normalized 2 cell expression profile');
plot(simple_vector, two_cell_data{a}(:,3),'Color',[0 .8 0],...
    'LineWidth',3,'LineStyle',linestyle(a))
xlabel('number of trials');
ylabel({'Rank-Ordered Expression Profile','(cDNA/gDNA)'});
legend('A')
title(title2)
axis auto
end
hold on
%figure 3
if ~isempty(four_cell_data{a})
figure(3)
simple_vector = (1:size(four_cell_data{a},1));
title3 = sprintf('normalized 4 cell expression profile');
plot(simple_vector, four_cell_data{a}(:,3),'Color',[1 .5 0],...
    'LineWidth',3,'LineStyle',linestyle(a))
xlabel('number of trials');
ylabel({'Rank-Ordered Expression Profile','(cDNA/gDNA)'});
legend('A')
title(title3)
axis auto
end
hold on
%figure 4
if ~isempty(eight_cell_data{a})
figure(4)
simple_vector = (1:size(eight_cell_data{a},1));
title4 = sprintf('normalized 8 cell expression profile');
plot(simple_vector, eight_cell_data{a}(:,3),'Color',[.8 0 0],...
    'LineWidth',3,'LineStyle',linestyle(a))
xlabel('number of trials');
ylabel({'Rank-Ordered Expression Profile','(cDNA/gDNA)'});
legend('A')
title(title4)
axis auto
end
hold on
%figure 5
if ~isempty(sixteen_cell_data{a})
figure(5)
simple_vector = (1:size(sixteen_cell_data{a},1));
title5 = sprintf('normalized 16 cell expression profile');
plot(simple_vector, sixteen_cell_data{a}(:,3),'Color',[.5 0 0],...
    'LineWidth',3,'LineStyle',linestyle(a))
xlabel('number of trials');
ylabel({'Rank-Ordered Expression Profile','(cDNA/gDNA)'});
legend('A')
title(title5)
axis auto
end
hold on
%figure 6
if ~isempty(thirtytwo_cell_data{a})
figure(6)
simple_vector = (1:size(thirtytwo_cell_data{a},1));
title6 = sprintf('normalized 32 cell expression profile');
plot(simple_vector, thirtytwo_cell_data{a}(:,3),'Color',[0 .3 0],...
    'LineWidth',3,'LineStyle',linestyle(a))
xlabel('number of trials');
ylabel({'Rank-Ordered Expression Profile','(cDNA/gDNA)'});
legend('A')
title(title6)
axis auto
end
hold on

end

toc
