function oldAlpha = setFigTransparency(hFig, alpha, fadeDuration, blockingFlag)
% setFigTransparency sets the transparency/opacity of a figure window
%
% Syntax:
%    oldAlpha = setFigTransparency(hFig, alpha, fadeDuration, blockingFlag)
%
% Description:
%    setFigTransparency sets the figure hFig's transparency value.
%    The entire figure window, including all internal menus, toolbars
%    and components, is made transparent according to the alpha value.
%
%    oldAlpha = setFigTransparency(...) returns the old transparency
%    value of the specified figure, prior to its modification.
%
%    This submission is based on an original idea implemented by
%    Malcolm Lidierth in his MUtilities submission:
%    http://www.mathworks.com/matlabcentral/fileexchange/28326-mutilities
%
% Input parameters:  (all parameters are optional)
%
%    hFig -         (default=gcf) Handle(s) of the modified figure(s).
%                   If component handle(s) is/are specified, the containing
%                   figure(s) will be inferred and used.
%
%    alpha -        (default=0.5) Transparency value, between 0.0 (=fully
%                   transparent) and 1.0 (=fully opaque).
%                   Note that all Matlab figure windows are created opaque.
%                   alpha<0 indicates that alpha value should not be modified.
%
%    fadeDuration - (default=0) Number of seconds for fade-in/fade-out effect
%                   Note: default value of 0 means immediately (no fading)
%
%    blockingFlag - (default=true) Whether or not the function should wait for
%                   the fade-in/fade-out effect to complete before returning
%
% Examples:
%    oldAlpha = setFigTransparency(hFig,-1);   % get hFig's current alpha
%    oldAlpha = setFigTransparency(hFig);      % set hFig's alpha to 0.5 (semi-transparent)
%    oldAlpha = setFigTransparency(hFig,0.7);  % set hFig's alpha to 0.7
%    oldAlpha = setFigTransparency([hFig1,hFig2],0.7);  % set transparency for several figures
%    oldAlpha = setFigTransparency(hFig,0.3,1.5,false); % non-blocking fade over 1.5 secs
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% Warning:
%    This code heavily relies on undocumented and unsupported Matlab functionality.
%    It works on Matlab 7.9 (R2009b) and higher, but use at your own risk!
%
% Change log:
%    2011-03-01: First version posted on Matlab's File Exchange: <a href="http://www.mathworks.com/matlabcentral/fileexchange/?term=authorid%3A27420">http://www.mathworks.com/matlabcentral/fileexchange/?term=authorid%3A27420</a>
%    2011-03-02: Prevented setting transparency for docked figures
%    2011-10-14: Fix for R2011b
%
% See also:
%    MUtilities, enableDisableFig, setFigDockGroup, getJFrame (all of them on the File Exchange)

% Programmed by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.2 $  $Date: 2011/10/14 03:46:26 $

  % Require Java 1.6.0_10 to run
  if ~usejava('awt')
      error([mfilename ' requires Java to run']);
  elseif ~ismethod('com.sun.awt.AWTUtilities','setWindowOpacity')
      error([mfilename ' requires a newer Java engine (1.6.0_10 minimum) to run']);
  end

  % Set default input parameters values
  if nargin < 1,  hFig = gcf;           end
  if nargin < 2,  alpha = 0.5;          end
  if nargin < 3,  fadeDuration = 0.0;   end
  if nargin < 4,  blockingFlag = true;  end

  % Check valid scalar values
  if ~isscalar(alpha) || alpha>1
      error([mfilename ' requires a scalar alpha value <= 1']);
  elseif ~isscalar(fadeDuration) || fadeDuration<0
      error([mfilename ' requires a scalar duration value >= 0']);
  end

  % Normalize blockingFlag
  if ischar(blockingFlag)
      blockingFlag = strcmpi(blockingFlag,'on');
  else
      blockingFlag = logical(blockingFlag);
  end

  % Get unique figure handles
  if ~all(ishghandle(hFig))
      error('hFig must be a valid GUI handle or array of handles');
  end
  hFigCell = {};
  for hIdx = 1 : length(hFig)
      if hFig(hIdx) == 0
          hFigCell{hIdx} = 0;
      else
          hFigCell{hIdx} = ancestor(hFig(hIdx),'figure');
      end
  end
  [unigueHandles,uniqueIdx,foldedIdx] = unique(hFig,'first');  %#ok unused
  currentAlpha = zeros(size(hFig));

  % Loop over all requested figures
  for figIdx = 1 : length(hFigCell)

      % Get the root Java frame
      jff = getJFrame(hFigCell{figIdx});

      % If this figure has already been processed - bail out
      if ~ismember(figIdx,uniqueIdx)
          currentAlpha(figIdx) = currentAlpha(foldedIdx(figIdx));
          continue;
      end

      % Get the current frame's state
      currentAlpha(figIdx) = com.sun.awt.AWTUtilities.getWindowOpacity(jff);

      % If same alpha or negative alpha - bail out
      if currentAlpha(figIdx) == alpha  || (alpha < 0),  continue;  end

      % Fade-in/out - automatically compute step size and duration
      deltaAlpha = alpha - currentAlpha(figIdx);
      maxStepAlpha = 0.03;
      if fadeDuration <= 0 || abs(deltaAlpha) < maxStepAlpha
          stepAlpha = deltaAlpha;
          steps = 1;
      else
          steps = fix(abs(deltaAlpha) / maxStepAlpha) + 1;
          stepAlpha = deltaAlpha / steps;
          stepDuration = fadeDuration / (steps-1);
      end

      % If blocking, do the fade effect immediately
      if blockingFlag || steps==1
          for stepIdx = 1 : steps
              newAlpha = currentAlpha(figIdx) + stepAlpha*stepIdx;
              com.sun.awt.AWTUtilities.setWindowOpacity(jff,newAlpha);
              jff.repaint;
              if stepIdx < steps,  pause(stepDuration);  end
          end
      else  % non-blocking: fade in/out asynchronously using a dedicated timer
          start(timer('ExecutionMode','fixedRate', 'Period',0.1, 'TasksToExecute',steps, 'TimerFcn', {@timerFcn,jff,currentAlpha(figIdx),stepAlpha}));
      end
  end

  % Return old alpha value(s), if requested
  if nargout
      oldAlpha = currentAlpha;
  end

end  % setFigTransparency


%% Get the root Java frame (up to 10 tries, to wait for figure to become responsive)
function jframe = getJFrame(hFig)

  % Ensure that hFig is a non-empty handle...
  if isempty(hFig)
      error(['Cannot retrieve the figure handle for handle ' num2str(hFig)]);
  end

  % Check for the desktop handle
  if hFig == 0
      %jframe = com.mathworks.mde.desk.MLDesktop.getInstance.getMainFrame; return;
      error('Setting the Desktop''s transparency is not supported');
  end

  % Check whether the figure is invisible
  if strcmpi(get(hFig,'Visible'),'off')
      error(['Cannot set the transparency for non-visible figure ' num2str(hFig)]);
  end

  % Check whether the figure is docked
  if strcmpi(get(hFig,'WindowStyle'),'docked')
      error(['Cannot set the transparency for docked figure ' num2str(hFig)]);
  end

  % Retrieve the figure window (JFrame) handle
  jframe = [];
  maxTries = 10;
  while maxTries > 0
      try
          % Get the figure's underlying Java frame
          jf = get(handle(hFig),'JavaFrame');

          % Get the Java frame's root frame handle
          %jframe = jf.getFigurePanelContainer.getComponent(0).getRootPane.getParent;
          try
              jframe = jf.fFigureClient.getWindow;  % equivalent to above...
          catch
              jframe = jf.fHG1Client.getWindow;  % equivalent to above...
          end
          if ~isempty(jframe)
              break;
          else
              maxTries = maxTries - 1;
              drawnow; pause(0.1);
          end
      catch
          maxTries = maxTries - 1;
          drawnow; pause(0.1);
      end
  end
  if isempty(jframe)
      error(['Cannot retrieve the Java Frame reference for figure ' num2str(hFig)]);
  end
  
end  % getJFrame


%% Timer function for non-blocking fade-in/fade-out effect
function timerFcn(hTimer,eventData,jFrame,currentAlpha,stepAlpha)  %#ok eventData is unused
  stepIdx = hTimer.TasksExecuted;
  newAlpha = currentAlpha + stepAlpha*stepIdx;
  com.sun.awt.AWTUtilities.setWindowOpacity(jFrame,newAlpha);
  jFrame.repaint;
  if stepIdx == hTimer.TasksToExecute
      stop(hTimer);
      delete(hTimer);
  end
end  % timerFcn