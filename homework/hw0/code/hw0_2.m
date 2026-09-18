close all
clc
clear

%% Problem 2: Periodic second derivative

Ns = [25 50 100 200 400 800 1600];
a = 0;
b = 2*pi;

hs          = zeros(length(Ns),1);
relInfError = zeros(length(Ns),1);
rel2Error   = zeros(length(Ns),1);

f   = @(x) exp(sin(x));
ddf = @(x) ((cos(x)).^2 - sin(x)).*exp(sin(x));

%% Convergence study

for k = 1:length(Ns)

    N = Ns(k);
    h = (b-a)/N;
    hs(k) = h;

    x = linspace(a,b,N+1);
    x = x(1:end-1);

    u = f(x);
    dduExact  = ddf(x);
    dduApprox = zeros(size(x));

    % Periodic wrap
    dduApprox(1) = ...
        (u(N) - 2*u(1) + u(2))/h^2;

    for i = 2:N-1
        dduApprox(i) = ...
            (u(i-1) - 2*u(i) + u(i+1))/h^2;
    end

    dduApprox(N) = ...
        (u(N-1) - 2*u(N) + u(1))/h^2;

    error = dduExact - dduApprox;

    relInfError(k) = ...
        norm(error,inf) / norm(dduExact,inf);

    rel2Error(k) = ...
        norm(error,2) / norm(dduExact,2);

end
%% Observed convergence rates

relInfOrder = NaN(length(Ns),1);
rel2Order   = NaN(length(Ns),1);

for k = 2:length(Ns)

    ratioH = hs(k-1)/hs(k);

    relInfOrder(k) = ...
        log(relInfError(k-1)/relInfError(k)) / log(ratioH);

    rel2Order(k) = ...
        log(rel2Error(k-1)/rel2Error(k)) / log(ratioH);

end

Results = table( ...
    Ns.', hs, ...
    relInfError, relInfOrder, ...
    rel2Error, rel2Order, ...
    'VariableNames', ...
    {'N','h','RelInfError','RelInfOrder','Rel2Error','Rel2Order'});

disp(Results)

%% Representative mesh

Nplot = 100;
hplot = (b-a)/Nplot;

xPlot = linspace(a,b,Nplot+1);
xPlot = xPlot(1:end-1);

uPlot = f(xPlot);
dduExactPlot  = ddf(xPlot);
dduApproxPlot = zeros(size(xPlot));

dduApproxPlot(1) = ...
    (uPlot(Nplot) - 2*uPlot(1) + uPlot(2))/hplot^2;

for i = 2:Nplot-1
    dduApproxPlot(i) = ...
        (uPlot(i-1) - 2*uPlot(i) + uPlot(i+1))/hplot^2;
end

dduApproxPlot(Nplot) = ...
    (uPlot(Nplot-1) - 2*uPlot(Nplot) + uPlot(1))/hplot^2;


%% Figure 1: Exact and finite-difference second derivative

figure(1)
clf

% Exact second derivative
plot(xPlot, dduExactPlot, '-', ...
    'LineWidth', 1.6, ...
    'DisplayName', 'Exact $f''''(x)$')

hold on

% Show only selected numerical points
markerStepPlot = max(1,round(length(xPlot)/25));
idxPlot = 1:markerStepPlot:length(xPlot);

% Approximate second derivative
plot(xPlot(idxPlot), dduApproxPlot(idxPlot), 'o', ...
    'LineStyle', 'none', ...
    'MarkerSize', 5, ...
    'LineWidth', 1.0, ...
    'DisplayName', 'Finite difference $D_h^2 f$')

xlabel('$x$', ...
    'Interpreter','latex', 'FontSize',18)

ylabel('$f''''(x)$', ...
    'Interpreter','latex', 'FontSize',18)

title('Exact and Finite-Difference Second Derivatives', ...
    'Interpreter','latex', 'FontSize',21)

legend('Location','best', ...
    'Interpreter','latex', 'FontSize',15)

grid on
box on

ax = gca;
ax.FontSize = 18;
ax.TickLabelInterpreter = 'latex';
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

xlim([0 2*pi])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$0$','$\pi/2$','$\pi$', ...
             '$3\pi/2$','$2\pi$'})

set(gcf,'Units','inches')
set(gcf,'Position',[1 1 6.5 4])

exportgraphics(gcf, ...
    'figures/p2_second_derivative_comparison.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');

%% Figure 2: Relative error convergence

figure(2)
clf

loglog(Ns, relInfError, 'o-', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6, ...
    'DisplayName', 'Relative $\|\cdot\|_\infty$ error')

hold on

loglog(Ns, rel2Error, 's-', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6, ...
    'DisplayName', 'Relative $\|\cdot\|_2$ error')

Nref = Ns(3:6);

ref2 = 0.7 * relInfError(3) .* ...
    (Nref/Nref(1)).^(-2);

loglog(Nref, ref2, '--', ...
    'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.5, ...
    'DisplayName', '$O(N^{-2})$')

xlabel('$N$ (grid intervals)', ...
    'Interpreter','latex', 'FontSize',18)

ylabel('Relative error', ...
    'Interpreter','latex', 'FontSize',18)

title('Convergence of Relative Error', ...
    'Interpreter','latex', 'FontSize',21)

legend('Location','best', ...
    'Interpreter','latex', ...
    'FontSize',15)

grid on
box on

ax = gca;
ax.FontSize = 18;
ax.TickLabelInterpreter = 'latex';
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

set(gcf,'Units','inches')
set(gcf,'Position',[1 1 6.5 4])

exportgraphics(gcf, ...
    'figures/p2_convergence.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');