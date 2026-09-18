close all
clc
clear

%% Problem 1: First derivative with mixed-order boundary treatment

% Convergence meshes
Ns = [25 50 100 200 400 800 1600];

a = 0;
b = 2*pi;

hs             = zeros(length(Ns),1);
endpointErrors = zeros(length(Ns),1);
interiorErrors = zeros(length(Ns),1);
relInfError    = zeros(length(Ns),1);
rel2Error      = zeros(length(Ns),1);

% Function and exact derivative
f  = @(x) exp(sin(x));
df = @(x) cos(x).*exp(sin(x));

%% Convergence study

for k = 1:length(Ns)

    N = Ns(k);
    h = (b-a)/N;
    hs(k) = h;

    x = linspace(a,b,N+1);

    u       = f(x);
    duExact = df(x);
    duApprox = zeros(size(x));

    % First-order one-sided differences at endpoints
    duApprox(1) = (u(2)-u(1))/h;
    duApprox(end) = (u(end)-u(end-1))/h;

    % Second-order centered difference in the interior
    for i = 2:N
        duApprox(i) = (u(i+1)-u(i-1))/(2*h);
    end

    error = duExact - duApprox;

    endpointErrors(k) = max(abs(error([1 end])));
    interiorErrors(k) = norm(error(2:N),inf);

    relInfError(k) = norm(error,inf) / norm(duExact,inf);
    rel2Error(k)   = norm(error,2)   / norm(duExact,2);

end

%% Observed convergence rates

endpointOrder = NaN(length(Ns),1);
interiorOrder = NaN(length(Ns),1);
relInfOrder   = NaN(length(Ns),1);
rel2Order     = NaN(length(Ns),1);

for k = 2:length(Ns)

    ratioH = hs(k-1)/hs(k);

    endpointOrder(k) = ...
        log(endpointErrors(k-1)/endpointErrors(k)) / log(ratioH);

    interiorOrder(k) = ...
        log(interiorErrors(k-1)/interiorErrors(k)) / log(ratioH);

    relInfOrder(k) = ...
        log(relInfError(k-1)/relInfError(k)) / log(ratioH);

    rel2Order(k) = ...
        log(rel2Error(k-1)/rel2Error(k)) / log(ratioH);

end

%% Representative mesh for solution/error figures

Nplot = 100;
hplot = (b-a)/Nplot;

xPlot = linspace(a,b,Nplot+1);

uPlot = f(xPlot);
duExactPlot = df(xPlot);
duApproxPlot = zeros(size(xPlot));

duApproxPlot(1) = (uPlot(2)-uPlot(1))/hplot;

for i = 2:Nplot
    duApproxPlot(i) = ...
        (uPlot(i+1)-uPlot(i-1))/(2*hplot);
end

duApproxPlot(end) = ...
    (uPlot(end)-uPlot(end-1))/hplot;

errorPlot = duExactPlot - duApproxPlot;

Results = table( ...
    Ns.', hs, ...
    endpointErrors, endpointOrder, ...
    interiorErrors, interiorOrder, ...
    relInfError, relInfOrder, ...
    rel2Error, rel2Order, ...
    'VariableNames', ...
    {'N','h', ...
    'EndpointError','EndpointOrder', ...
    'InteriorError','InteriorOrder', ...
    'RelInfError','RelInfOrder', ...
    'Rel2Error','Rel2Order'});

disp(Results)


%% Figure 1: Exact and finite-difference derivative

figure(1)
clf

plot(xPlot, duExactPlot, '-', ...
    'LineWidth', 1.6, ...
    'DisplayName', 'Exact $f''(x)$')

hold on

markerStepPlot = max(1,round(length(xPlot)/25));
idxPlot = 1:markerStepPlot:length(xPlot);

plot(xPlot(idxPlot), duApproxPlot(idxPlot), 'o', ...
    'LineStyle', 'none', ...
    'MarkerSize', 4.5, ...
    'LineWidth', 1.0, ...
    'DisplayName', 'Finite difference $D_hf$')

xlabel('$x$', 'Interpreter','latex','FontSize',18)
ylabel('$f''(x)$', 'Interpreter','latex','FontSize',18)

title('Exact and Finite-Difference Derivatives', ...
    'Interpreter','latex','FontSize',21)

legend('Location','best', ...
    'Interpreter','latex','FontSize',15)

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
    'figures/p1_derivative_comparison.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');


%% Figure 2: Point-wise absolute error

figure(2)
clf

hErr = semilogy(xPlot, abs(errorPlot), '-', ...
    'LineWidth', 1.2, ...
    'DisplayName', '$|e_h(x_i)|$');

hold on

hEnd = semilogy([xPlot(1) xPlot(end)], ...
                [abs(errorPlot(1)) abs(errorPlot(end))], ...
                'o', ...
                'LineStyle', 'none', ...
                'MarkerSize', 7, ...
                'LineWidth', 1.2, ...
                'DisplayName', 'First-order endpoints');

xlabel('$x$', 'Interpreter','latex','FontSize',18)
ylabel('$|e_h(x_i)|$', 'Interpreter','latex','FontSize',18)

legend([hErr hEnd], ...
    'Location','best', ...
    'Interpreter','latex','FontSize',15)

title('Pointwise Error', 'Interpreter','latex','FontSize',21)

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
    'figures/p1_pointwise_error.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');



%% Figure 3: Relative-error convergence

figure(3)
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


% Reference rates anchored to the first numerical error
refInf = 0.65 * relInfError(1) .* (Ns/Ns(1)).^(-1);

ref2 = 1.35 * rel2Error(1) .* (Ns/Ns(1)).^(-3/2);

loglog(Ns, refInf, '--', ...
    'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.5, ...
    'DisplayName', '$O(N^{-1})$')

loglog(Ns, ref2, ':', ...
    'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.5, ...
    'DisplayName', '$O(N^{-3/2})$')

xlabel('$N$ (grid intervals)', 'Interpreter','latex','FontSize',18)
ylabel('Relative error', 'Interpreter','latex','FontSize',18)

title('Convergence of Relative Error', ...
    'Interpreter','latex','FontSize',21)

legend('Location','best', 'Interpreter','latex','FontSize',15)

grid on
ax = gca;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

box on

ax = gca;
ax.FontSize = 18;
ax.TickLabelInterpreter = 'latex';
ax.GridAlpha = 0.15;

set(gcf,'Units','inches')
set(gcf,'Position',[1 1 6.5 4])

exportgraphics(gcf, ...
    'figures/p1_convergence.pdf', ...
    'ContentType','vector',...
    'BackgroundColor','white');
