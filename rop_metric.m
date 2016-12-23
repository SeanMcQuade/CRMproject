function [ myanswer ] = rop_metric( ROP_1, ROP_2 )
global Number_of_Iterations
%Computes the distance between two rank order profiles
%   for loop through min amount of Rows between a and b 
%   (though they will have the same number of Rows if from my sim).
running_total = 0; %initialize sum
%running total of differences between two ROP
for i=1:min(size(ROP_1,1),size(ROP_2,1))
    running_total = running_total + abs(ROP_1(i,3) - ROP_2(i,3));

end
Number_of_Iterations
size(running_total)
size(Number_of_Iterations)
myanswer = running_total/Number_of_Iterations;
end
