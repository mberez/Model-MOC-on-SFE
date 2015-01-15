
function yok = index(y,ok)
% function yok = index(y,ok)
% Returns y(ok) or y{ok} if y is a cell array
% Useful when you don't want to create a dummy variable to index.
% Example,
%   Yok = index(fft(y),2:20);

  yok = y;
  if (isempty(y))
    return
  end
  if (iscell(y(1)))
    yok = y{ok};
  else
    yok = y(ok);
  end