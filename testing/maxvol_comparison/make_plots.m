function updates = make_plots(probname, startiter)
  figwidth = 410;
  figheight = 280;
  maxvol_tol = [10.0; 4.0; 2.0; 1.1];
  updates_hrt = [];
  kktiter_hrt = [];
  updates_seq = [];
  kktiter_seq = [];
  for i=1:length(maxvol_tol)
      v = maxvol_tol(i);
      %% heuristic
      probname_method = sprintf('%s_heuristic', probname);
      [updates,kktiter] = readlogfile(probname_method, v, startiter);
      if size(updates,1) > size(updates_hrt,1)
          pad = zeros(size(updates,1)-size(updates_hrt,1), size(updates_hrt,2));
          pad(:) = NaN;
          updates_hrt = [updates_hrt; pad];
          kktiter_hrt = [kktiter_hrt; pad];
      else
          pad = zeros(size(updates_hrt,1)-size(updates,1), 1);
          pad(:) = NaN;
          updates = [updates; pad];
          kktiter = [kktiter; pad];
      end
      updates_hrt = [updates_hrt updates];
      kktiter_hrt = [kktiter_hrt kktiter];
      %% sequential
      probname_method = sprintf('%s_sequential', probname);
      [updates,kktiter] = readlogfile(probname_method, v, startiter);
      if size(updates,1) > size(updates_seq,1)
          pad = zeros(size(updates,1)-size(updates_seq,1), size(updates_seq,2));
          pad(:) = NaN;
          updates_seq = [updates_seq; pad];
          kktiter_seq = [kktiter_seq; pad];
      else
          pad = zeros(size(updates_seq,1)-size(updates,1), 1);
          pad(:) = NaN;
          updates = [updates; pad];
          kktiter = [kktiter; pad];
      end
      updates_seq = [updates_seq updates];
      kktiter_seq = [kktiter_seq kktiter];
  end
  %% Plot updates.
  maxiter = max(size(updates_hrt,1), size(updates_seq,1));
  figure('units', 'points', 'position', [0,0,figwidth,figheight]);
  %colors = {[0 0 0]; [0 0 1]; [1 0 0]; [0 1 0]};
  colors = {[0 0 0]; [0.25 0.25 0.25]; [0.5 0.5 0.5]; [0.75 0.75 0.75]};
  h = plot(updates_seq, '-');
  set(h, {'color'}, colors);
  hold on;
  h = plot(updates_hrt, '--');
  set(h, {'color'}, colors);
  xlim([1 maxiter]);
  xtk = (1:5:maxiter)';
  xticks(xtk);
  set(gca, 'xticklabel', num2str(xtk+startiter-1));
  legend('\rho=10.0','\rho=4.0','\rho=2.0','\rho=1.1');
  xlabel('interior point iteration');
  ylabel('basis updates');
  crop_margins();
  %% Plot kktiter.
  figure('units', 'points', 'position', [0,0,figwidth,figheight]);
  h = plot(kktiter_seq, '-');
  set(h, {'color'}, colors);
  hold on;
  h = plot(kktiter_hrt, '--');
  set(h, {'color'}, colors);
  xlim([1 maxiter]);
  xtk = (1:5:maxiter)';
  xticks(xtk);
  set(gca, 'xticklabel', num2str(xtk+startiter-1));
  legend('\rho=10.0','\rho=4.0','\rho=2.0','\rho=1.1');
  xlabel('interior point iteration');
  ylabel('CR iterations');
  crop_margins();
end

function [updates,kktiter] = readlogfile(probname_method, v, startiter)
  filename = sprintf('logs/%s_%05.1f.log', probname_method, v);
  fid = fopen(filename, 'r');
  line_in = [];
  while true
      line_in = fgetl(fid);
      if strcmp(line_in(2:4), sprintf('%3d', startiter))
          break
      end
  end
  updates = [];
  kktiter = [];
  while true
      fields = strsplit(line_in);
      updates = [updates; str2double(fields{11})];
      kktiter = [kktiter; str2double(fields{12})];
      if line_in(5)=='*' break; end
      line_in = fgetl(fid);
  end
  fclose(fid);
end
