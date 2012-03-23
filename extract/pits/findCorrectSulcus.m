function sulc = findCorrectSulcus(label1,label2)

labs = sort([label1 label2]);

sulc = 0;

if (labs(1) == 3 && labs(2) == 28)
   sulc = 1;   
elseif (labs(1) == 3 && labs(2) == 18)
   sulc = 2;   
elseif (labs(1) == 18 && labs(2) == 30)
   sulc = 3;   
elseif (labs(1) == 15 && labs(2) == 30)
   sulc = 4;   
elseif (labs(1) == 9 && labs(2) == 15)
   sulc = 5;   
elseif (labs(1) == 24 && labs(2) == 28)
   sulc = 6;   
elseif (labs(1) == 3 && labs(2) == 24)
   sulc = 6;   
elseif (labs(1) == 18 && labs(2) == 24)
   sulc = 6;   
elseif (labs(1) == 22 && labs(2) == 24)
   sulc = 7;   
elseif (labs(1) == 8 && labs(2) == 29)
   sulc = 8;   
elseif (labs(1) == 2 && labs(2) == 28)
   sulc = 9;   
elseif (labs(1) == 2 && labs(2) == 17)
   sulc = 9;   
elseif (labs(1) == 17 && labs(2) == 25)
   sulc = 9;   
elseif (labs(1) == 5 && labs(2) == 25)
   sulc = 10;   
elseif (labs(1) == 5 && labs(2) == 13)
   sulc = 11;   
end















