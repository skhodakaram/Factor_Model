function [ Output ] = PlotTube(data, param, ColorTube)
% plot a tube for data given the size "param.interval"
% data: experiments * param.values

values = param.values;
interval = param.interval;
I = prctile(data,interval,1); 

Output.MIN  = I(1,:);
Output.MAX  = I(2,:);
Output.plot = fill([values,fliplr(values)], [Output.MIN fliplr(Output.MAX)], ColorTube);
set(Output.plot,'EdgeColor','None');

Output.AVR  = mean(data,1) % Mean of all of the elements of the columns of 'data'

end

