%% Problem 3: Finite-difference solution of a boundary-value ODE

Ns = [25 50 100 200 400 800 1600];
x0 = 0;
xN = 5;

hs          = zeros(length(Ns),1);
relInfError = zeros(length(Ns),1);
rel2Error   = zeros(length(Ns),1);

uExactFun = @(x) sin(x);
g         = @(x) sin(x);
f         = @(x) sin(x).*cos(x);

for k = 1:length(Ns)

    N = Ns(k);
    h = (xN-x0)/N;
    hs(k) = h;

    x = linspace(x0,xN,N+1);
    xInterior = x(2:N).';

    gInterior = g(xInterior);
    fInterior = f(xInterior);

    mainDiag = (h^2 - 2)*ones(N-1,1);
    upDiag   = 1 + (h/2)*gInterior(1:end-1);
    lowDiag  = 1 - (h/2)*gInterior(2:end);

    A = diag(mainDiag) ...
      + diag(upDiag,1) ...
      + diag(lowDiag,-1);

    b = h^2*fInterior;

    % Dirichlet boundary contributions
    b(1) = b(1) ...
        - (1-(h/2)*gInterior(1))*uExactFun(x0);

    b(end) = b(end) ...
        - (1+(h/2)*gInterior(end))*uExactFun(xN);

    uInterior = A\b;

    uNum   = [uExactFun(x0); uInterior; uExactFun(xN)];
    uExact = uExactFun(x).';

    error = uExact - uNum;

    relInfError(k) = norm(error,inf) / norm(uExact,inf);
    rel2Error(k)   = norm(error,2)   / norm(uExact,2);

end

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
    {'N','h', ...
    'RelInfError','RelInfOrder', ...
    'Rel2Error','Rel2Order'});

disp(Results)

%% Representative mesh

Nplot = 100;
hplot = (xN-x0)/Nplot;

xPlot = linspace(x0,xN,Nplot+1);
xPlot = xPlot.';
xInteriorPlot = xPlot(2:Nplot);

gInteriorPlot = sin(xInteriorPlot).';
fInteriorPlot = (sin(xInteriorPlot).*cos(xInteriorPlot)).';

% Construct tridiagonal system
mainDiagPlot = (hplot^2 - 2)*ones(Nplot-1,1);

upDiagPlot = ...
    1 + (hplot/2)*gInteriorPlot(1:end-1);

lowDiagPlot = ...
    1 - (hplot/2)*gInteriorPlot(2:end);

Aplot = diag(mainDiagPlot) ...
    + diag(upDiagPlot,1) ...
    + diag(lowDiagPlot,-1);

% Right-hand side
bplot = hplot^2*fInteriorPlot;

% Dirichlet boundary contributions
bplot(1) = bplot(1) ...
    - (1-(hplot/2)*gInteriorPlot(1))*sin(x0);

bplot(end) = bplot(end) ...
    - (1+(hplot/2)*gInteriorPlot(end))*sin(xN);



% Solve
bplot = bplot(:);
uInteriorPlot = Aplot\bplot;

uNumPlot = [sin(x0); uInteriorPlot; sin(xN)];
uExactPlot = sin(xPlot).';


%% Figure 1: Exact and finite-difference solution

figure(1)
clf

% Exact derivative
plot(xPlot, uExactPlot, '-', ...
    'LineWidth', 1.6, ...
    'DisplayName', 'Exact $u(x)=\sin(x)$')

hold on

% Show only selected numerical points
markerStepPlot = max(1,round(length(xPlot)/25));
idxPlot = 1:markerStepPlot:length(xPlot);

% Numerical Solution
plot(xPlot(idxPlot), uNumPlot(idxPlot), 'o', ...
    'MarkerSize', 4, ...
    'LineWidth', 1.0, ...
    'DisplayName', 'Finite-difference $u_h(x)$')

xlabel('$x$', ...
    'Interpreter','latex','FontSize',18)

ylabel('$u(x)$', ...
    'Interpreter','latex','FontSize',18)

title('Exact and Finite-Difference Solutions', ...
    'Interpreter','latex','FontSize',21)

legend('Location','southwest', ...
    'Interpreter','latex','FontSize',15)

grid on
box on

ax = gca;
ax.FontSize = 18;
ax.TickLabelInterpreter = 'latex';
ax.GridAlpha = 0.15;

xlim([0 5])
xticks([0 1 2 3 4 5])
xticklabels({'$0$','$1$','$2$','$3$','$4$','$5$'})

set(gcf,'Units','inches')
set(gcf,'Position',[1 1 6.5 4])

exportgraphics(gcf, ...
    'figures/p3_solution_comparison.pdf', ...
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

% Second-order reference
Nref = Ns(3:6);

ref2 = 0.7*relInfError(3) .* ...
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
    'Interpreter','latex', 'FontSize',15)

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
    'figures/p3_convergence.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');