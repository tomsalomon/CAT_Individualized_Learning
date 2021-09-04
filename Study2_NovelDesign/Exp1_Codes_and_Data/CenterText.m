% David/Wolfelab's CenterText function code

function [newX, newY] = CenterText (window, message, color, xoffset, yoffset)

if nargin < 2
   error([mfilename ' requires at least two arguments']);
end
if nargin < 3 | isempty(color)
   color = [];
end
if nargin < 4 | isempty(xoffset)
   xoffset = 0;
end
if nargin < 5 | isempty(yoffset)
   yoffset = 0;
end

rectScreen = Screen('Rect', window);
[rectText, rectDummy] = Screen('TextBounds', window, message);
width = RectWidth(rectText);
Screen('DrawText', window, message, ((rectScreen(3)/2)-(width/2))+xoffset, (rectScreen(4)/2)+yoffset, color);
newX = 0;
newY = 0;