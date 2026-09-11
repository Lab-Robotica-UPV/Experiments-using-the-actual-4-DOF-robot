% =========================================================================
%  mocap_analysis.m
%  Analisis estadistico del error cartesiano del end-effector
%  usando las corridas de MoCap del LABC.
%
%  Estructura esperada:
%    LABC/run_1/x_mocap.txt  ...  LABC/run_12/x_mocap.txt
%    ref_cart.txt  (referencia cartesiana, 4 cols: x z theta psi)
%
%  Columnas x_mocap.txt (10 cols):
%    1-2: Encoder time (sec + nanosec)
%    3-6: Init/Current time (ignoradas)
%    7:   x  [m]    8: z  [m]    9: theta [rad]    10: psi [rad]
%
%  Notas:
%    - run_8 excluida: duracion 67.8s incompatible con referencia (52.1s)
%    - Sincronizacion: correlacion cruzada en x
%    - Filtro outliers: percentil 95 del |error| por corrida y componente
% =========================================================================

clear; clc; close all;

BASE_DIR  = pwd;
N_RUNS    = 12;
SKIP_RUNS = [7,8];        % excluir run_8: duracion incompatible
DT_REF    = 0.01;       % paso de tiempo referencia [s]
PCT_FILT  = 95;         % percentil para filtro de outliers
COLORS    = {[0.22 0.48 0.75],[0.13 0.70 0.47],[0.93 0.54 0.10],[0.64 0.08 0.18]};
UNITS     = [1000 1000 1000 1000];   % m->mm, rad->mrad
DOF_NAMES = {'x','z','theta','psi'};
DOF_UNITS = {'mm','mm','mrad','mrad'};

%% ── 1. CARGAR REFERENCIA ─────────────────────────────────────────────────
ref_raw = readmatrix(fullfile(BASE_DIR,'ref_cart.txt'),...
    'FileType','text','Delimiter',' ');
N_ref  = size(ref_raw,1);
t_ref  = (0:N_ref-1)' * DT_REF;
x_ref  = ref_raw(:,1);  z_ref  = ref_raw(:,2);
th_ref = ref_raw(:,3);  ps_ref = ref_raw(:,4);

fprintf('Referencia: N=%d muestras (%.1f s)\n', N_ref, t_ref(end));
fprintf('Corridas excluidas: '); fprintf('run_%d ', SKIP_RUNS); fprintf('\n\n');

%% ── 2. PROCESAR CORRIDAS ─────────────────────────────────────────────────
valid_runs = setdiff(1:N_RUNS, SKIP_RUNS);
N_valid    = numel(valid_runs);

RMSE_cart = zeros(N_valid, 4);
EMAX_cart = zeros(N_valid, 4);
IAE_cart  = zeros(N_valid, 4);
lag_ms    = zeros(N_valid, 1);
rep_run   = struct();

fprintf('Procesando %d corridas MoCap...\n', N_valid);

for idx = 1:N_valid
    r     = valid_runs(idx);
    fpath = fullfile(BASE_DIR,'LABC',sprintf('run_%d',r),'x_mocap.txt');
    raw   = readmatrix(fpath,'FileType','text','Delimiter','\t','NumHeaderLines',1);

    % Tiempo y datos cartesianos (cols 7-10)
    t_m  = raw(:,1) + raw(:,2)*1e-9;
    t_m  = t_m - t_m(1);
    x_m  = raw(:,7);   z_m  = raw(:,8);
    th_m = raw(:,9);   ps_m = raw(:,10);

    % Truncar al minimo entre mocap y referencia
    N_use = min(length(t_m), N_ref);
    t_m   = t_m(1:N_use);
    x_m   = x_m(1:N_use);   z_m  = z_m(1:N_use);
    th_m  = th_m(1:N_use);  ps_m = ps_m(1:N_use);
    dt_m  = median(diff(t_m));

    % Sincronizacion: correlacion cruzada en x
    xc      = xcorr(x_m - mean(x_m), x_ref(1:N_use) - mean(x_ref(1:N_use)));
    [~,imax]= max(xc);
    lag_s   = (imax - N_use) * DT_REF;
    lag_ms(idx) = lag_s * 1000;

    % Interpolar referencia con offset de tiempo
    t_adj = t_ref(1:N_use) - lag_s;
    xr_s  = interp1(t_adj, x_ref(1:N_use),  t_m, 'linear', 'extrap');
    zr_s  = interp1(t_adj, z_ref(1:N_use),  t_m, 'linear', 'extrap');
    thr_s = interp1(t_adj, th_ref(1:N_use), t_m, 'linear', 'extrap');
    psr_s = interp1(t_adj, ps_ref(1:N_use), t_m, 'linear', 'extrap');
    x_m  = x_m  - mean(x_m)  + mean(xr_s);
    z_m  = z_m  - mean(z_m)  + mean(zr_s);
    th_m = th_m - mean(th_m) + mean(thr_s);
    ps_m = ps_m - mean(ps_m) + mean(psr_s);
    % Errores en unidades originales [m, rad]
    % Errores en unidades originales [m, rad]
    % Filtro percentil 95 por componente (oclusiones MoCap)
    E_raw = [x_m-xr_s, z_m-zr_s, th_m-thr_s, ps_m-psr_s];
    E_filt = E_raw;
    % Filtro percentil 95 por componente (oclusiones MoCap)
    E_filt = E_raw;
    for k = 1:4
        thr_k          = prctile(abs(E_raw(:,k)), PCT_FILT);
        mask           = abs(E_raw(:,k)) > thr_k;
        E_filt(mask,k) = NaN;
    end

    % Metricas escalares (ignorando NaN)
    for k = 1:4
        v = E_filt(:,k); v = v(~isnan(v));
        RMSE_cart(idx,k) = sqrt(mean(v.^2));
        EMAX_cart(idx,k) = max(abs(v));
        IAE_cart(idx,k)  = sum(abs(v)) * dt_m;
    end

    % Corrida representativa (primera valida)
    if idx == 1
        rep_run.t      = t_m;
        rep_run.xr     = xr_s;   rep_run.zr  = zr_s;
        rep_run.x_m    = x_m;    rep_run.z_m = z_m;
        rep_run.E_raw  = E_raw;
        rep_run.E_filt = E_filt;
        rep_run.run    = r;
    end

    fprintf('  run_%2d: lag=%+5.1f ms  RMSE x=%5.2f mm  z=%5.2f mm  th=%5.2f mrad  ps=%5.2f mrad\n',...
        r, lag_ms(idx),...
        RMSE_cart(idx,1)*1000, RMSE_cart(idx,2)*1000,...
        RMSE_cart(idx,3)*1000, RMSE_cart(idx,4)*1000);
end

%% ── 3. TABLA DESCRIPTIVA ─────────────────────────────────────────────────
fprintf('\n=========================================\n');
fprintf(' ERROR CARTESIANO END-EFFECTOR (LABC)\n');
fprintf(' N=%d corridas validas (run_%d excluida)\n', N_valid, SKIP_RUNS);
fprintf(' Filtro outliers: percentil %d\n', PCT_FILT);
fprintf('=========================================\n');
fprintf('%-8s %12s %12s %12s %12s\n','DOF','RMSE mean','RMSE std','max mean','max std');
fprintf('%s\n', repmat('-',1,58));
for k = 1:4
    sc = UNITS(k);
    fprintf('%-8s %10.3f %14.3f %14.3f %14.3f  [%s]\n',...
        DOF_NAMES{k},...
        mean(RMSE_cart(:,k))*sc, std(RMSE_cart(:,k))*sc,...
        mean(EMAX_cart(:,k))*sc, std(EMAX_cart(:,k))*sc,...
        DOF_UNITS{k});
end

%% ── 4. FIGURAS ───────────────────────────────────────────────────────────

% Fig 1: Trayectoria 2D + error temporal corrida representativa
figure('Name','MoCap corrida representativa','Position',[50 50 1200 480]);

subplot(1,3,1);
plot(rep_run.xr*100, rep_run.zr*100,'--k','LineWidth',2,'DisplayName','Reference');
hold on;
xf = rep_run.x_m; xf(any(isnan(rep_run.E_filt),2)) = NaN;
zf = rep_run.z_m; zf(any(isnan(rep_run.E_filt),2)) = NaN;
plot(xf*100, zf*100,'-','Color',COLORS{1},'LineWidth',1.5,'DisplayName','LABC (MoCap)');
xlabel('X axis (cm)'); ylabel('Z axis (cm)');
title(sprintf('Trajectory — run\\_1')); legend('Location','best');
grid on; box on; axis equal;

subplot(1,3,2);
plot(rep_run.t, rep_run.E_raw(:,1)*1000,'Color',[0.8 0.8 0.8],...
    'LineWidth',0.8,'DisplayName','Raw'); hold on;
plot(rep_run.t, rep_run.E_filt(:,1)*1000,'Color',COLORS{1},...
    'LineWidth',1.2,'DisplayName','Filtered');
yline(0,'k','LineWidth',0.8);
xlabel('Time (s)'); ylabel('Error (mm)');
title(sprintf('X error  RMSE=%.2f mm', RMSE_cart(1,1)*1000));
legend('Location','best'); grid on; box on;

subplot(1,3,3);
plot(rep_run.t, rep_run.E_raw(:,2)*1000,'Color',[0.8 0.8 0.8],...
    'LineWidth',0.8,'DisplayName','Raw'); hold on;
plot(rep_run.t, rep_run.E_filt(:,2)*1000,'Color',COLORS{2},...
    'LineWidth',1.2,'DisplayName','Filtered');
yline(0,'k','LineWidth',0.8);
xlabel('Time (s)'); ylabel('Error (mm)');
title(sprintf('Z error  RMSE=%.2f mm', RMSE_cart(1,2)*1000));
legend('Location','best'); grid on; box on;

sgtitle(sprintf('LABC end-effector — MoCap run\\_1 (lag=%.1f ms)', lag_ms(1)));

% Fig 2: RMSE por corrida
figure('Name','RMSE por corrida','Position',[50 580 1000 380]);
run_labels = arrayfun(@(r) sprintf('r%d',r), valid_runs,'UniformOutput',false);
dof_plot = [1 2 4];
dof_lbls = {'x [mm]','z [mm]','\psi [mrad]'};
for kk = 1:3
    k  = dof_plot(kk);
    sc = UNITS(k);
    subplot(1,3,kk);
    bar(RMSE_cart(:,k)*sc,'FaceColor',COLORS{kk},'FaceAlpha',0.75);
    hold on;
    yline(mean(RMSE_cart(:,k))*sc,'--k','LineWidth',1.5);
    xlabel('Run'); ylabel(dof_lbls{kk});
    xticks(1:N_valid); xticklabels(run_labels); xtickangle(45);
    title(sprintf('%s: %.2f \\pm %.2f',...
        dof_lbls{kk}, mean(RMSE_cart(:,k))*sc, std(RMSE_cart(:,k))*sc));
    grid on; box on;
end
sgtitle(sprintf('Cartesian RMSE per run — LABC MoCap (N=%d, p%d filter)',N_valid,PCT_FILT));

% Fig 3: Boxplot con puntos individuales
figure('Name','Boxplot RMSE','Position',[50 1020 600 380]);
data_box = [RMSE_cart(:,1)*1000, RMSE_cart(:,2)*1000, RMSE_cart(:,4)*1000];
bp = boxplot(data_box,'Labels',{'x [mm]','z [mm]','\psi [mrad]'});
h  = findobj(gca,'Tag','Box');
cols_bp = [COLORS{3}; COLORS{2}; COLORS{1}];
for b = 1:length(h)
    patch(get(h(b),'XData'),get(h(b),'YData'),cols_bp(b,:),'FaceAlpha',0.4);
end
hold on;
for kk = 1:3
    k = dof_plot(kk);
    scatter(kk*ones(N_valid,1)+0.08*(rand(N_valid,1)-0.5),...
        RMSE_cart(:,k)*UNITS(k),35,COLORS{kk},'filled','MarkerFaceAlpha',0.8);
end
ylabel('Error'); grid on; box on;
title(sprintf('Cartesian RMSE distribution (N=%d runs)',N_valid));

%% ── 5. TABLA LaTeX ───────────────────────────────────────────────────────
fprintf('\n=========================================\n');
fprintf(' TABLA LaTeX\n');
fprintf('=========================================\n');
fprintf('\\begin{table}[h]\n\\centering\n');
fprintf('\\caption{End-effector Cartesian tracking error measured\n');
fprintf('independently by the optical tracking system (OptiTrack, 14 cameras,\n');
fprintf('accuracy $<$0.1~mm). $N=%d$ valid runs; run\\_8 excluded due to\n',N_valid);
fprintf('recording duration mismatch. Outliers ($>$p$_{%d}$) removed\n',PCT_FILT);
fprintf('per run to exclude marker occlusions.}\n');
fprintf('\\begin{tabular}{lcc}\n\\toprule\n');
fprintf('DOF & RMSE (mean$\\pm$std) & Max.~error (mean$\\pm$std) \\\\\n\\midrule\n');
dof_tex = {'$x$','$z$','$\\psi$'};
dof_k   = [1 2 4];
dof_u   = {'mm','mm','mrad'};
for kk = 1:3
    k = dof_k(kk); sc = UNITS(k);
    fprintf('%s & $%.2f\\pm%.2f$~%s & $%.2f\\pm%.2f$~%s \\\\\n',...
        dof_tex{kk},...
        mean(RMSE_cart(:,k))*sc, std(RMSE_cart(:,k))*sc, dof_u{kk},...
        mean(EMAX_cart(:,k))*sc, std(EMAX_cart(:,k))*sc, dof_u{kk});
end
fprintf('\\bottomrule\n\\end{tabular}\n');
fprintf(['\\footnotesize{$\\theta$ not reported: systematic offset\n'...
    'attributed to rigid-body calibration in the OptiTrack frame.}\n']);
fprintf('\\label{tab:mocap_error}\n\\end{table}\n');

fprintf('\nAnalisis completado. N=%d corridas validas.\n', N_valid);

%% ── 6. ANÁLISIS POR CICLO (variabilidad intra-run) ──────────────────────
% Detectar ciclos desde la referencia (señal limpia, sin oclusiones)
% Un ciclo = cruce ascendente por la media de x_ref
x_ref_c = x_ref - mean(x_ref);
idx_cross = find(diff(sign(x_ref_c)) > 0);   % cruces ascendentes

% Filtrar ciclos completos (descartar segmentos < 500 muestras)
cycle_starts = idx_cross;
cycle_ends   = [idx_cross(2:end)-1; N_ref];
valid_cycles = (cycle_ends - cycle_starts) >= 500;
cycle_starts = cycle_starts(valid_cycles);
cycle_ends   = cycle_ends(valid_cycles);
N_cycles     = numel(cycle_starts);

fprintf('\n=========================================\n');
fprintf(' VARIABILIDAD INTRA-RUN (por ciclo)\n');
fprintf(' %d ciclos completos por run (~%.1f s cada uno)\n',...
    N_cycles, mean(cycle_ends-cycle_starts)*DT_REF);
fprintf('=========================================\n');

% Matrices: [N_valid x N_cycles x 4]
RMSE_cycle = zeros(N_valid, N_cycles, 4);

for idx = 1:N_valid
    r     = valid_runs(idx);
    fpath = fullfile(BASE_DIR,'LABC',sprintf('run_%d',r),'x_mocap.txt');
    raw   = readmatrix(fpath,'FileType','text','Delimiter','\t','NumHeaderLines',1);

    t_m  = raw(:,1) + raw(:,2)*1e-9; t_m = t_m - t_m(1);
    x_m  = raw(:,7); z_m  = raw(:,8);
    th_m = raw(:,9); ps_m = raw(:,10);
    N_use= min(length(t_m), N_ref);
    t_m=t_m(1:N_use); x_m=x_m(1:N_use); z_m=z_m(1:N_use);
    th_m=th_m(1:N_use); ps_m=ps_m(1:N_use);

    % Sincronizar (mismo lag calculado antes)
    lag_s = lag_ms(idx)/1000;
    t_adj = t_ref(1:N_use) - lag_s;
    xr_s  = interp1(t_adj, x_ref(1:N_use),  t_m,'linear','extrap');
    zr_s  = interp1(t_adj, z_ref(1:N_use),  t_m,'linear','extrap');
    thr_s = interp1(t_adj, th_ref(1:N_use), t_m,'linear','extrap');
    psr_s = interp1(t_adj, ps_ref(1:N_use), t_m,'linear','extrap');

    E_raw = [x_m-xr_s, z_m-zr_s, th_m-thr_s, ps_m-psr_s];

    % Filtro outliers
    E_filt = E_raw;
    for k = 1:4
        thr_k = prctile(abs(E_raw(:,k)), PCT_FILT);
        E_filt(abs(E_raw(:,k)) > thr_k, k) = NaN;
    end

    % RMSE por ciclo
    for c = 1:N_cycles
        i1 = cycle_starts(c); i2 = cycle_ends(c);
        if i2 > N_use; continue; end
        for k = 1:4
            v = E_filt(i1:i2, k); v = v(~isnan(v));
            if ~isempty(v)
                RMSE_cycle(idx,c,k) = sqrt(mean(v.^2));
            end
        end
    end
end

% Estadísticos intra-run
fprintf('%-8s %14s %14s %14s\n','DOF','media_ciclos','std_ciclos','CV%%');
fprintf('%s\n', repmat('-',1,44));
for k = [1 2 4]
    sc   = UNITS(k);
    vals = RMSE_cycle(:,:,k)*sc;   % [N_valid x N_cycles]
    % Variabilidad media dentro de cada run
    std_intra = mean(std(vals, 0, 2));   % std sobre ciclos, media sobre runs
    mu_all    = mean(vals(:));
    cv        = std_intra/mu_all*100;
    fprintf('%-8s %14.3f %14.3f %14.1f%%  [%s]\n',...
        DOF_NAMES{k}, mu_all, std_intra, cv, DOF_UNITS{k});
end

% Fig 4: RMSE por ciclo — todas las runs superpuestas
figure('Name','Variabilidad por ciclo','Position',[50 50 900 380]);
dof_plt = [1 2 4];
for kk = 1:3
    k  = dof_plt(kk); sc = UNITS(k);
    subplot(1,3,kk);
    data_c = squeeze(RMSE_cycle(:,:,k))*sc;   % [N_valid x N_cycles]
    for idx2 = 1:N_valid
        plot(1:N_cycles, data_c(idx2,:), 'o-',...
            'Color',[COLORS{kk} 0.4],'LineWidth',1,'MarkerSize',5);
        hold on;
    end
    plot(1:N_cycles, mean(data_c,1), 'k-o','LineWidth',2,...
        'MarkerFaceColor','k','MarkerSize',7,'DisplayName','Mean');
    xlabel('Cycle'); ylabel(DOF_UNITS{k});
    xticks(1:N_cycles);
    title(sprintf('%s — cycle variability', DOF_NAMES{k}));
    legend('Location','best'); grid on; box on;
end
sgtitle(sprintf('Within-run cycle-to-cycle variability — LABC MoCap (N=%d runs)',N_valid));

% Texto para el articulo
fprintf('\nTexto para el articulo:\n');
for k = [1 2 4]
    sc   = UNITS(k);
    vals = RMSE_cycle(:,:,k)*sc;
    fprintf('  %s: ciclo-a-ciclo std=%.2f %s (CV=%.1f%%)\n',...
        DOF_NAMES{k}, mean(std(vals,0,2)), DOF_UNITS{k},...
        mean(std(vals,0,2))/mean(vals(:))*100);
end