function [ ha ] = CreateTitleForSubplots( title )
%CREATETITLEFORSUBPLOTS Summary of this function goes here
%   Detailed explanation goes here

    ha = axes();
    ha.Position = [0 0 1 1];
    ha.XLim = [0 1];
    ha.YLim = [0 1];
    ha.Box = 'off';
    ha.Visible = 'off';
    ha.Units = 'normalized';
    ha.Clipping = 'off';

    text(0.5, 1,title ...
        ,'HorizontalAlignment' ,'center','VerticalAlignment', 'top', 'FontSize', 14);
end

