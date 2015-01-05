function [string] = three_digit_string(num)

temp = num;

hundreds = (temp-mod(temp,100))/100;
temp = temp - hundreds*100;
tens = (temp-mod(temp,10))/10;
ones = temp - tens*10;

string = [num2str(hundreds) num2str(tens) num2str(ones)];


