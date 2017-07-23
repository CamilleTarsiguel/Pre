%%
function printEnergyOnly(energy)
[totEn, D, S, E, P, L, hreg, hlin, hang, hper, hocc, hseg] = ...
    getEnergyValues(energy);

printMessage(2,'\n| ------------- ENERGY  VALUES ----------------| --------------- Label Cost ------------- |||');
printMessage(2,'\n|  Energy|   Data| Smth| DetExc| TrjExc|  Lcost|  hreg|  hlin|  hang|  hper|   hocc|  hseg|||\n');
printMessage(2,'|%8.1f|%7.1f|%5.1f|%7.1f|%7.1f|%7.1f|%6.1f|%6.1f|%6.1f|%6.1f|%7.1f|%6.1f|||\n', ...
    totEn,D,S,E,P,L,hreg,hlin,hang,hper,hocc,hseg); %%% iter output
printMessage(2,  '| ---------------------------------------------| ---------------------------------------- |||\n');