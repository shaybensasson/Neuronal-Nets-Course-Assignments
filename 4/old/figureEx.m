function FigHandle = figureEx(varargin)
%figureMax renders a figure on single screen or on secondary screen if
%exists
%{
run: figureMax('units','normalized','outerposition',[0 0 1 1])
to maximize on the primary screen, as fallback, if no sec exists
%}
%{
original code, before modification found here
%http://www.mathworks.com/matlabcentral/answers/16663-is-it-possible-to-viewing-the-figure-window-on-second-display
%}
    MAXIMIZE = 0;
    indexCustom = find(strncmpi(varargin(1:end), 'Custom', 5),1);
    if ~isempty(indexCustom)
        indexMax = strncmpi(varargin(indexCustom+1), 'Maximize', 8);
        MAXIMIZE = ~isempty(indexMax);
    end
    
    MP = get(0, 'MonitorPositions');
    if size(MP, 1) == 1  % Single monitor
      if (MAXIMIZE)
          FigHandle = figure('units','normalized','outerposition',[0 0 1 1]);
      else
          FigHandle = figure(varargin{:});
      end
        
    else                 % Multiple monitors
      paramVisible = get(0, 'DefaultFigureVisible');


      if (MAXIMIZE)
        FigHandle = figure('Visible', 'off');
        set(FigHandle, 'Units', 'pixels');
        
        left = MP(2, 1);
        bottom = MP(1, 4)-MP(2, 4);
        height = MP(2, 4);
        width = MP(2, 3)-left;
      
        p = [left+2 bottom+5 width-2 height-74];
      else
        FigHandle     = figure(varargin{:}, 'Visible', 'off');
        set(FigHandle, 'Units', 'pixels');
      
        posShift = MP(2, 1:2);
        pos = get(FigHandle, 'Position');
        p = [pos(1:2) + posShift, pos(3:4)];
      end  
      %{
      where 'left' and 'bottom' define the distance from the lower-left corner 
    of the screen to the lower-left corner of the figure window. 
    'width' and 'height' define the dimensions of the window.
      %}
      set(FigHandle, 'Position', p, ...
          'Visible', paramVisible);
    end