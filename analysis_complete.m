% =========================================================================
%  analysis_complete.m
%  Análisis estadístico completo: PID vs CTC vs LABC
%  Fundamento metodológico: Nehmzow (2006), Scientific Methods in
%  Mobile Robotics, Sections 2.6.3, 3.3.5, 3.4.3, 3.4.4
%
%  Estructura de carpetas esperada (en el directorio de trabajo):
%    PID/run_1/ ... PID/run_10/
%    CTC/run_1/ ... CTC/run_10/
%    LABC/run_1/ ... LABC/run_10/
%
%  Cada run_X contiene:
%    joint_position.txt   ref_pose.txt   control_action.txt
%    (formato TSV, 10 columnas: 6 timestamps + 4 articulaciones)
%
%  Flujo del análisis (siguiendo Nehmzow 2006):
%    1. Carga y cálculo de métricas escalares por corrida (anti pseudo-rep.)
%    2. Friedman test global (k=3 controladores, Sec. 3.4.4)
%    3. Wilcoxon pareado post-hoc por pares (Sec. 3.4.3)
%    4. Corrección de Bonferroni (3 pares x 4 articulaciones)
%    5. Effect size e IC 95% bootstrap
%    6. Análisis de potencia post-hoc
%    7. Figuras y tabla LaTeX
% =========================================================================

clear; clc; close all;

%% ── 0. PARÁMETROS ────────────────────────────────────────────────────────
BASE_DIR    = pwd;
N_RUNS      = 12;          % corridas por controlador (mínimo común)
N_JOINTS    = 4;
DT          = 0.01;        % paso de tiempo [s]
ALPHA       = 0.05;
N_PAIRS     = 3;           % PID-CTC, PID-LABC, CTC-LABC
N_COMP      = N_PAIRS * N_JOINTS;          % 12 comparaciones (Bonferroni)
ALPHA_ADJ   = ALPHA / N_COMP;             % 0.0042
CTRLS       = {'PID','CTC','LABC'};
JOINTS      = {'q_{1,3}','q_{2,3}','q_{3,3}','q_{4,2}'};
PAIRS       = {'PID vs CTC','PID vs LABC','CTC vs LABC'};
COLORS      = {[0.22 0.48 0.75],[0.93 0.54 0.10],[0.13 0.70 0.47]};

fprintf('=========================================\n');
fprintf(' ANÁLISIS ESTADÍSTICO COMPLETO\n');
fprintf(' PID vs CTC vs LABC  |  N=%d corridas\n', N_RUNS);
fprintf(' Ref: Nehmzow (2006) Sec. 2.6.3, 3.4.3, 3.4.4\n');
fprintf('=========================================\n');
fprintf('Alpha global=%.2f  |  Bonferroni (%d comp.)  |  alpha_adj=%.4f\n\n',...
    ALPHA, N_COMP, ALPHA_ADJ);

%% ── 1. CARGA Y CÁLCULO DE MÉTRICAS ──────────────────────────────────────
% Métricas: RMSE, MaxError, IAE, TVE, TVC
% Un escalar por corrida por articulación (Nehmzow 2006, Sec. 2.6.3:
% "collapse the observations into one value before analysis")
METRIC_NAMES = {'RMSE','MaxError','IAE','TVE','TVC'};
N_METRICS    = numel(METRIC_NAMES);

% Matrices [N_RUNS x N_JOINTS x N_METRICS] por controlador
for c = 1:3
    eval(sprintf('%s_data = zeros(N_RUNS, N_JOINTS, N_METRICS);', CTRLS{c}));
end

for c = 1:3
    ctrl = CTRLS{c};
    fprintf('Cargando %s...\n', ctrl);
    for r = 1:N_RUNS
        run_dir = fullfile(BASE_DIR, ctrl, sprintf('run_%d', r));
        jp  = readmatrix(fullfile(run_dir,'joint_position.txt'), ...
                         'FileType','text','Delimiter','\t','NumHeaderLines',1);
        rp  = readmatrix(fullfile(run_dir,'ref_pose.txt'),       ...
                         'FileType','text','Delimiter','\t','NumHeaderLines',1);
        ca  = readmatrix(fullfile(run_dir,'control_action.txt'), ...
                         'FileType','text','Delimiter','\t','NumHeaderLines',1);
        n   = min([size(jp,1) size(rp,1) size(ca,1)]);
        q   = jp(1:n, 7:10);
        ref = rp(1:n, 7:10);
        u   = ca(1:n, 7:10);
        e   = q - ref;
        for j = 1:N_JOINTS
            ej = e(:,j); uj = u(:,j);
            metrics = [sqrt(mean(ej.^2));          % RMSE
                       max(abs(ej));               % MaxError
                       sum(abs(ej))*DT;            % IAE
                       sum(abs(uj))*DT;            % TVE
                       sum(abs(diff(uj)))];        % TVC
            eval(sprintf('%s_data(r,j,:) = metrics;', ctrl));
        end
        fprintf('  run_%d OK\n', r);
    end
end
fprintf('\n');

%% ── 2. TABLA DESCRIPTIVA ─────────────────────────────────────────────────
fprintf('=========================================\n');
fprintf(' ESTADÍSTICOS DESCRIPTIVOS (media ± std)\n');
fprintf('=========================================\n');
fprintf('%-10s %-4s %18s %18s %18s %10s %10s\n',...
    'Metrica','Jnt','PID','CTC','LABC','PID-CTC%%','PID-LABC%%');
fprintf('%s\n', repmat('-',1,90));

for m = 1:N_METRICS
    for j = 1:N_JOINTS
        pid_v  = PID_data(:,j,m);
        ctc_v  = CTC_data(:,j,m);
        labc_v = LABC_data(:,j,m);
        pct_ctc  = (mean(pid_v)-mean(ctc_v)) /mean(pid_v)*100;
        pct_labc = (mean(pid_v)-mean(labc_v))/mean(pid_v)*100;
        fprintf('%-10s %-4s %8.5f±%7.5f %8.5f±%7.5f %8.5f±%7.5f %+9.1f%% %+9.1f%%\n',...
            METRIC_NAMES{m}, JOINTS{j},...
            mean(pid_v), std(pid_v),...
            mean(ctc_v), std(ctc_v),...
            mean(labc_v),std(labc_v),...
            pct_ctc, pct_labc);
    end
    fprintf('%s\n', repmat('-',1,90));
end

%% ── 3. FRIEDMAN TEST GLOBAL ──────────────────────────────────────────────
% Nehmzow (2006) Sec. 3.4.4: Kruskal-Wallis / Friedman para k>2 grupos
fprintf('\n=========================================\n');
fprintf(' PASO 1: FRIEDMAN TEST (k=3 controladores)\n');
fprintf(' Pregunta: ¿hay alguna diferencia global?\n');
fprintf('=========================================\n');
fprintf('%-4s %10s %10s %-16s\n','Jnt','Chi2','p-val','Conclusion');
fprintf('%s\n', repmat('-',1,44));

friedman_sig = false(1, N_JOINTS);
for j = 1:N_JOINTS
    X = [PID_data(:,j,1), CTC_data(:,j,1), LABC_data(:,j,1)];
    [p_f, tbl] = friedman(X, 1, 'off');
    chi2 = tbl{2,5};
    friedman_sig(j) = p_f < ALPHA;
    fprintf('%-4s %10.4f %10.4f %-16s\n', JOINTS{j}, chi2, p_f,...
        ternary(p_f<ALPHA,'SI -> post-hoc','no significativo'));
end

%% ── 4. WILCOXON PAREADO POST-HOC ─────────────────────────────────────────
% Nehmzow (2006) Sec. 3.4.3: Wilcoxon test for paired observations
% "This variation can be used when outcomes are paired by some criterion"
fprintf('\n=========================================\n');
fprintf(' PASO 2: WILCOXON PAREADO POST-HOC\n');
fprintf(' alpha_adj=%.4f (Bonferroni: %d pares x %d articulaciones)\n',...
    ALPHA_ADJ, N_PAIRS, N_JOINTS);
fprintf('=========================================\n');
fprintf('%-4s %-16s %8s %8s %8s %-14s\n',...
    'Jnt','Par','W','p-val','r','Significativo?');
fprintf('%s\n', repmat('-',1,62));

pair_A = {PID_data, PID_data, CTC_data};
pair_B = {CTC_data, LABC_data, LABC_data};
results = struct();

for j = 1:N_JOINTS
    for p = 1:N_PAIRS
        A = pair_A{p}; B = pair_B{p};
        d_j = A(:,j,1) - B(:,j,1);   % diferencias RMSE pareadas

        if friedman_sig(j)
            [pval, ~, stats] = signrank(d_j, 0, 'method','exact');
            mu_T    = N_RUNS*(N_RUNS+1)/4;
            sigma_T = sqrt((2*N_RUNS+1)*mu_T/6);
            z       = (stats.signedrank - mu_T) / sigma_T;
            r       = min(abs(z)/sqrt(N_RUNS), 1.0);
            sig = pval < ALPHA_ADJ;
            fprintf('%-4s %-16s %8.1f %8.4f %8.3f %s\n',...
                JOINTS{j}, PAIRS{p}, stats.signedrank, pval, r,...
                ternary(sig,'SI *','no'));
        else
            pval = NaN; r = NaN; sig = false;
            fprintf('%-4s %-16s %8s %8s %8s %s\n',...
                JOINTS{j}, PAIRS{p},'--','--','--','Friedman no sig.');
        end
        results(j,p).diffs = d_j;
        results(j,p).pval  = pval;
        results(j,p).r     = r;
        results(j,p).sig   = sig;
        results(j,p).pair  = PAIRS{p};
    end
end
fprintf('\n* p < %.4f (Bonferroni, %d comparaciones)\n', ALPHA_ADJ, N_COMP);

%% ── 5. EFFECT SIZE E IC 95% BOOTSTRAP ────────────────────────────────────
fprintf('\n=========================================\n');
fprintf(' PASO 3: EFFECT SIZE E IC 95%% BOOTSTRAP\n');
fprintf('=========================================\n');
fprintf('%-4s %-16s %6s %8s %-22s %10s\n',...
    'Jnt','Par','r','magnitud','IC 95%% diferencia','mejora%%');
fprintf('%s\n', repmat('-',1,72));

rng(42);
for j = 1:N_JOINTS
    for p = 1:N_PAIRS
        if ~results(j,p).sig; continue; end
        d_j = results(j,p).diffs;
        r   = results(j,p).r;
        if r < 0.1;     mag = 'trivial';
        elseif r < 0.3; mag = 'pequeno';
        elseif r < 0.5; mag = 'mediano';
        else;           mag = 'grande';
        end
        boot  = bootstrp(5000, @mean, d_j);
        ci_lo = prctile(boot, 2.5);
        ci_hi = prctile(boot, 97.5);
        A_mean = mean(pair_A{p}(:,j,1));
        pct    = mean(d_j)/A_mean*100;
        fprintf('%-4s %-16s %6.3f %8s [%+.5f,%+.5f] %+9.1f%%\n',...
            JOINTS{j}, PAIRS{p}, r, mag, ci_lo, ci_hi, pct);
        results(j,p).ci_lo = ci_lo;
        results(j,p).ci_hi = ci_hi;
        results(j,p).pct   = pct;
    end
end

%% ── 6. ANÁLISIS DE POTENCIA POST-HOC ────────────────────────────────────
fprintf('\n=========================================\n');
fprintf(' PASO 4: POTENCIA POST-HOC (RMSE)\n');
fprintf('=========================================\n');
fprintf('%-4s %-16s %10s %10s %10s\n',...
    'Jnt','Par','Cohen_d','Potencia','N_min 80%%');
fprintf('%s\n', repmat('-',1,54));

for j = 1:N_JOINTS
    for p = 1:N_PAIRS
        if ~results(j,p).sig; continue; end
        d_j     = results(j,p).diffs;
        d_cohen = abs(mean(d_j))/std(d_j);
        ncp     = d_cohen*sqrt(N_RUNS);
        t_crit  = tinv(1-ALPHA_ADJ/2, N_RUNS-1);
        sig_nc  = sqrt(1+ncp^2/(2*(N_RUNS-1)));
        power   = 1-normcdf(t_crit,ncp,sig_nc)+normcdf(-t_crit,ncp,sig_nc);
        n_min   = NaN;
        for nt = 3:100
            ncp_t = d_cohen*sqrt(nt);
            tc_t  = tinv(1-ALPHA_ADJ/2,nt-1);
            sg_t  = sqrt(1+ncp_t^2/(2*(nt-1)));
            pw_t  = 1-normcdf(tc_t,ncp_t,sg_t)+normcdf(-tc_t,ncp_t,sg_t);
            if pw_t>=0.80; n_min=nt; break; end
        end
        fprintf('%-4s %-16s %10.3f %10.3f %10d\n',...
            JOINTS{j}, PAIRS{p}, d_cohen, power, n_min);
    end
end

%% ── 7. FIGURAS ───────────────────────────────────────────────────────────

%% Fig 1: Boxplot RMSE por controlador y articulación
figure('Name','RMSE por controlador','Position',[50 50 1100 450]);
for j = 1:N_JOINTS
    subplot(1,4,j);
    data_box = [PID_data(:,j,1), CTC_data(:,j,1), LABC_data(:,j,1)];
    bp = boxplot(data_box, 'Labels',{'PID','CTC','LABC'});
    h = findobj(gca,'Tag','Box');
    for b = 1:length(h)
        col_idx = length(h)-b+1;
        patch(get(h(b),'XData'),get(h(b),'YData'),COLORS{col_idx},'FaceAlpha',0.4);
    end
    % Marcar pares significativos
    y_max = max(data_box(:))*1.08;
    for p = 1:N_PAIRS
        if results(j,p).sig
          parts = strsplit(PAIRS{p}, ' vs ');
          x1 = find(strcmp({'PID','CTC','LABC'}, strtrim(parts{1})));
          x2 = find(strcmp({'PID','CTC','LABC'}, strtrim(parts{2})));
          y_line = y_max + (p-1)*max(data_box(:))*0.04;
          line([x1 x2],[y_line y_line],'Color','k','LineWidth',1.2);
          text((x1+x2)/2, y_line+max(data_box(:))*0.01,'*',...
            'HorizontalAlignment','center','FontSize',14,'Color','r');
        end
    end
    title(JOINTS{j},'FontSize',11);
    ylabel('RMSE [rad]'); grid on; box on;
end
sgtitle(sprintf('RMSE por controlador — * p < %.4f (Bonferroni)',ALPHA_ADJ),'FontSize',11);

%% Fig 2: Forest plot — diferencias pareadas con IC 95%
figure('Name','Forest plot IC','Position',[50 550 900 500]);
colors_sig = [0.13 0.70 0.47];
colors_ns  = [0.70 0.70 0.70];
y_pos = 0; y_ticks = []; y_labels = {};

for p = 1:N_PAIRS
    for j = 1:N_JOINTS
        y_pos = y_pos + 1;
        if results(j,p).sig && ~isnan(results(j,p).ci_lo)
            col = colors_sig;
            d_m = mean(results(j,p).diffs);
            ci_lo = results(j,p).ci_lo;
            ci_hi = results(j,p).ci_hi;
        else
            col = colors_ns;
            d_m = mean(results(j,p).diffs);
            se  = std(results(j,p).diffs)/sqrt(N_RUNS);
            tc  = tinv(0.975, N_RUNS-1);
            ci_lo = d_m - tc*se;
            ci_hi = d_m + tc*se;
        end
        plot([ci_lo ci_hi],[y_pos y_pos],'-','Color',col,'LineWidth',2); hold on;
        plot(d_m, y_pos,'o','Color',col,'MarkerFaceColor',col,'MarkerSize',7);
        y_ticks(end+1) = y_pos;
        y_labels{end+1} = sprintf('%s %s',JOINTS{j},PAIRS{p});
    end
    y_pos = y_pos + 0.5;
end

xline(0,'--k','LineWidth',1.2);
yticks(y_ticks); yticklabels(y_labels);
xlabel('\Delta RMSE = A - B [rad]  (> 0 \Rightarrow B mejora)');
title('Forest plot: diferencias pareadas RMSE  (IC 95% bootstrap)');
legend([plot(nan,nan,'o-','Color',colors_sig,'MarkerFaceColor',colors_sig)...
        plot(nan,nan,'o-','Color',colors_ns, 'MarkerFaceColor',colors_ns)],...
    {'Significativo','No significativo'},'Location','best');
grid on; box on;

%% Fig 3: Trade-off RMSE vs TVE
figure('Name','Trade-off RMSE vs TVE','Position',[50 1050 700 500]);
hold on;
markers = {'o','s','^'};
for c = 1:3
    for j = 1:N_JOINTS
        rmse_v = eval(sprintf('%s_data(:,j,1)', CTRLS{c}));
        tve_v  = eval(sprintf('%s_data(:,j,4)', CTRLS{c}));
        scatter(tve_v, rmse_v, 50, COLORS{c}, markers{c},...
            'DisplayName', sprintf('%s J%d',CTRLS{c},j-1),...
            'LineWidth',1.5);
    end
end
xlabel('TVE — Esfuerzo de control (V·s)');
ylabel('RMSE (rad)');
title('Trade-off: precisión vs esfuerzo de control');
legend('Location','best','NumColumns',3,'FontSize',7);
grid on; box on;

%% ── 8. TABLA LaTeX ───────────────────────────────────────────────────────
fprintf('\n=========================================\n');
fprintf(' TABLA LaTeX — RMSE (copiar al artículo)\n');
fprintf('=========================================\n');
fprintf('\\begin{table}[h]\n\\centering\n');
fprintf('\\caption{Statistical comparison: PID vs CTC vs LABC ($N=%d$ paired runs)}\n',N_RUNS);
fprintf('\\begin{tabular}{llrrrrrr}\n\\toprule\n');
fprintf('Joint & Controller & RMSE (rad) & MaxErr (rad) & IAE & $W$ & $p$ & $r$ \\\\\n\\midrule\n');

for j = 1:N_JOINTS
    for c = 1:3
        d = eval(sprintf('%s_data(:,j,:)',CTRLS{c}));
        rmse_s = sprintf('%.5f$\\pm$%.5f', mean(d(:,:,1)), std(d(:,:,1)));
        emax_s = sprintf('%.5f$\\pm$%.5f', mean(d(:,:,2)), std(d(:,:,2)));
        iae_s  = sprintf('%.4f$\\pm$%.4f', mean(d(:,:,3)), std(d(:,:,3)));
        if c == 1
            fprintf('%s & %s & %s & %s & %s & & & \\\\\n',...
                JOINTS{j}, CTRLS{c}, rmse_s, emax_s, iae_s);
        elseif c == 2
            p_idx = 1; % PID vs CTC
            sig_m = ternary(results(j,p_idx).sig,'$^*$','');
            W_str = ternary(results(j,p_idx).sig,...
                sprintf('%.1f',results(j,p_idx).diffs'*0+0),'--');
            fprintf('   & %s & %s & %s & %s & & & %s\\\\\n',...
                CTRLS{c}, rmse_s, emax_s, iae_s, sig_m);
        else
            p_idx = 2; % PID vs LABC
            sig_m = ternary(results(j,p_idx).sig,'$^*$','');
            fprintf('   & %s & %s & %s & %s & & & %s\\\\\n',...
                CTRLS{c}, rmse_s, emax_s, iae_s, sig_m);
        end
    end
    if j < N_JOINTS; fprintf('\\midrule\n'); end
end

fprintf('\\bottomrule\n');
fprintf('\\end{tabular}\n');
fprintf(['\\footnotesize{$^*$Significant at $\\alpha_{adj}=%.4f$ ',...
    '(Bonferroni, %d comparisons). Friedman + Wilcoxon signed-rank ',...
    '(exact, Nehmzow 2006).}\n'], ALPHA_ADJ, N_COMP);
fprintf('\\label{tab:statistical_comparison}\n\\end{table}\n');

%% ── FUNCIÓN AUXILIAR ─────────────────────────────────────────────────────
function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end

function [H, pValue, W] = swtest(x, alpha)
    if nargin<2; alpha=0.05; end
    x=real(x(:)); x=sort(x); n=length(x);
    if n<3; H=false; pValue=1; W=1; return; end
    m=norminv(((1:n)'-3/8)/(n+1/4)); m2=m'*m;
    c=m/sqrt(m2); u=1/sqrt(n);
    a_n=polyval([-2.706056,4.434685,-2.071190,-0.147981,0.221157,c(n)],u);
    if n==3
        a=zeros(n,1); a(n)=a_n; a(1)=-a_n;
    else
        a_n2=polyval([-3.582633,5.682633,-1.752461,-0.293762,0.042981,c(n-1)],u);
        a=zeros(n,1); a(n)=a_n; a(n-1)=a_n2; a(1)=-a_n; a(2)=-a_n2;
        phi=m2-2*a_n^2-2*a_n2^2;
        if phi<=0; H=false; pValue=1; W=1; return; end
        a(3:n-2)=m(3:n-2)/sqrt(phi);
    end
    SS=(x-mean(x))'*(x-mean(x));
    if SS<=0; H=false; pValue=1; W=1; return; end
    W=min(real((a'*x)^2/SS),1-eps);
    if n<=11
        gamma=polyval([-2.273,0.459],n);
        y=-log(gamma-log(1-W));
        mu=polyval([-1.5861,-0.31082,-0.083751,0.0038915],log(n));
        sigma=exp(polyval([-0.4803,-0.082676,0.0030302],log(n)));
    else
        y=log(1-W);
        mu=polyval([4.3881,-1.7531,0.20788,-0.006714],log(n));
        sigma=exp(polyval([-3.4280,1.7888,-0.22178,0.0092478],log(n)));
    end
    sigma=real(sigma);
    if sigma<=0||~isfinite(sigma); pValue=0.5; H=false; return; end
    pValue=max(0,min(1,1-normcdf(real((y-mu)/sigma))));
    H=(pValue<alpha);
end
