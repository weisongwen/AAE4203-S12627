function rinex2csv(obs_file, nav_file, out_dir)
% GPS L1 teaching baseline: RINEX -> five measurement CSVs + epoch metadata.
% Derived from GAO Yixin's MIT-licensed rinex2csv.m; original in legacy/.
% No figures, no Mapping Toolbox, no ground-truth input. All paths explicit.
% Example: rinex2csv('benchmark/...obs','benchmark/...nav','outputs/benchmark_inputs')
if nargin ~= 3
    error('Usage: rinex2csv(obs_file, nav_file, out_dir)');
end
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'matRTKLIB'));
if ~isfile(obs_file) || ~isfile(nav_file)
    error('Observation or navigation file does not exist.');
end
if ~isfolder(out_dir); mkdir(out_dir); end
obs = gt.Gobs(obs_file);
nav = gt.Gnav(nav_file);
obs = obs.selectSat(obs.sys == 1); % GPS only: a single receiver clock state
if obs.nsat < 4; error('Fewer than four GPS satellites in the file.'); end
n = obs.n; m = obs.nsat;
P = NaN(m,n); D = P; I = P; T = P; B = P; E = P;
S = NaN(m,3*n);
xlast = [];
success = false(n,1);
for k = 1:n
    o = obs.selectTime(k);
    sat = gt.Gsat(o,nav);
    sp = [sat.x(:), sat.y(:), sat.z(:)];
    pr = o.L1.P(:); clk = sat.dts(:);
    tgd = nav.getTGD(sat.sat); tgd = tgd(:);
    mask = isfinite(pr) & pr > 1e7 & isfinite(clk) & all(isfinite(sp),2) & sat.svh(:)==0;
    if nnz(mask) < 4; continue; end
    % Coarse observation-only solution, also used to recover after failed epochs.
    if isempty(xlast); x = zeros(4,1); else; x = xlast; end
    good = true;
    for j = 1:15
        delta = sp(mask,:) - x(1:3)';
        range = vecnorm(delta,2,2);
        H = [-delta ./ range, ones(nnz(mask),1)];
        if rank(H)<4; good=false; break; end
        dx = H \ (pr(mask)+clk(mask)-tgd(mask)-range-x(4));
        x = x+dx;
        if norm(dx)<1e-3; break; end
    end
    if ~good || ~all(isfinite(x)) || norm(x(1:3))<6.2e6 || norm(x(1:3))>6.5e6
        continue;
    end
    % Iterate atmospheric corrections at an observation-derived approximate position.
    % Unlike the old script, ionosphere/troposphere are subtracted exactly once.
    converged = false;
    for j = 1:20
        sat.setRcvPos(gt.Gpos(x(1:3)', 'xyz'));
        ion = sat.ionL1(:); trp = sat.trp(:); el = sat.el(:);
        use = mask & isfinite(ion) & isfinite(trp) & el>=15 & o.L1.S(:)>=15;
        if nnz(use)<4; break; end
        H = [-sat.ex(use)', -sat.ey(use)', -sat.ez(use)', ones(nnz(use),1)];
        if rank(H)<4; break; end
        v = pr(use)+clk(use)-ion(use)-trp(use)-tgd(use)-sat.rng(use)'-x(4);
        w = sind(el(use)).^2;
        dx = lscov(H,v,w);
        x = x+dx;
        if norm(dx)<1e-3; converged=true; break; end
    end
    if ~converged; continue; end
    xlast = x;
    sat.setRcvPos(gt.Gpos(x(1:3)', 'xyz'));
    el=sat.el(:); ion=sat.ionL1(:); trp=sat.trp(:);
    use=mask & isfinite(ion) & isfinite(trp) & el>=15 & o.L1.S(:)>=15;
    P(use,k)=pr(use); D(use,k)=clk(use); I(use,k)=ion(use); T(use,k)=trp(use);
    B(use,k)=tgd(use); E(use,k)=el(use);
    block=NaN(m,3); block(use,:)=sp(use,:); S(:,3*k-2:3*k)=block;
    success(k)=true;
    if mod(k,200)==0; fprintf('Converted %d/%d epochs\n',k,n); end
end
writematrix(P,fullfile(out_dir,'pseudoranges_meas.csv'));
writematrix(D,fullfile(out_dir,'satellite_clock_bias.csv'));
writematrix(I,fullfile(out_dir,'ionospheric_delay.csv'));
writematrix(T,fullfile(out_dir,'tropospheric_delay.csv'));
writematrix(S,fullfile(out_dir,'satellite_positions.csv'));
writematrix(B,fullfile(out_dir,'group_delay.csv'));
writematrix(E,fullfile(out_dir,'elevation_deg.csv'));
epochs=table(obs.time.week(:),obs.time.tow(:),'VariableNames',{'gps_week','gps_tow_s'});
writetable(epochs,fullfile(out_dir,'epochs.csv'));
writetable(table(obs.satstr(:),'VariableNames',{'satellite_id'}),fullfile(out_dir,'satellites.csv'));
fid=fopen(fullfile(out_dir,'conversion_summary.txt'),'w');
fprintf(fid,'GPS L1 only; elevation >=15 deg; SNR >=15 dB-Hz.\n');
fprintf(fid,'Input epochs: %d; converted: %d; omitted: %d.\n',n,nnz(success),n-nnz(success));
fprintf(fid,'Atmosphere and TGD are separate arrays; satellite clocks are metres.\n');
fprintf(fid,'Unsolved epochs remain NaN, never replaced by the previous position.\n');
fclose(fid);
fprintf('RINEX conversion complete: %d/%d epochs. Output: %s\n',nnz(success),n,out_dir);
end
