Ns = [20 40 80 160];

errors = zeros(length(Ns),1);
hs = zeros(length(Ns),1);

for k = 1:length(Ns)

    N = Ns(k);
    h = 5/N;
    hs(k) = h;
    m = (N-1);
    x = linspace(0,5,N+1);
    y = linspace(0,5,N+1);

    [X,Y] = meshgrid(x,y);

    u = sin(X).*cos(Y);

    lapExact = -2*sin(X).*cos(Y);

    uNum = zeros(N+1,N+1);
        uNum(1,:) = u(1,:);
        uNum(end,:) = u(end,:);
        uNum(:,1) = u(:,1);
        uNum(:,end) = u(:,end);

    A = sparse(m^2,m^2);
    b = zeros(m^2,1);

    for j = 2:N
        for i = 2:N
            p = (i-1) + (j-2)*m;
            b(p) = h^2*lapExact(j,i);
            A(p,p) = -4;
                
            if i > 2
                A(p,p-1) = 1;
            else
                b(p) = b(p) - u(j,1);
            end

            if i < N
                A(p, p+1) = 1;
            else
                b(p) = b(p) - u(j,N+1);
            end

            if j > 2
                A(p,p-m) = 1;
            else
                b(p) = b(p) - u(1,i);
            end

            if j < N
                A(p,p+m) = 1;
            else
                b(p) = b(p) - u(N+1,i);
            end

        end
    end

uInterior = A \ b;
uInteriorGrid = reshape(uInterior, m, m).';

uNum(2:N,2:N) = uInteriorGrid;


    errorInterior = u(2:N,2:N) - uNum(2:N,2:N);
    errors(k) = norm(errorInterior(:),inf);

end

ref = .5*errors(1)*(hs/hs(1)).^2;

orders = NaN(length(Ns),1);

for k = 2:length(Ns)
    orders(k) = log(errors(k-1)/errors(k))/log(2);
end

Results = table(Ns', hs, errors, orders, ...
    'VariableNames', {'N','h','Inf_Error','Order'});

disp(Results)

%% Representative mesh for Dirichlet solution

Nplot = 100;
hplot = 5/Nplot;
mplot = Nplot - 1;

xPlot = linspace(0,5,Nplot+1);
yPlot = linspace(0,5,Nplot+1);

[Xplot,Yplot] = meshgrid(xPlot,yPlot);

uExactPlot = sin(Xplot).*cos(Yplot);
fPlot = -2*sin(Xplot).*cos(Yplot);

% Interior finite-difference system
Aplot = sparse(mplot^2,mplot^2);
bplot = zeros(mplot^2,1);

for j = 2:Nplot
    for i = 2:Nplot

        p = (i-1) + (j-2)*mplot;

        Aplot(p,p) = -4;
        bplot(p) = hplot^2*fPlot(j,i);

        % Left neighbor
        if i > 2
            Aplot(p,p-1) = 1;
        else
            bplot(p) = bplot(p) - uExactPlot(j,1);
        end

        % Right neighbor
        if i < Nplot
            Aplot(p,p+1) = 1;
        else
            bplot(p) = bplot(p) - uExactPlot(j,Nplot+1);
        end

        % Lower neighbor
        if j > 2
            Aplot(p,p-mplot) = 1;
        else
            bplot(p) = bplot(p) - uExactPlot(1,i);
        end

        % Upper neighbor
        if j < Nplot
            Aplot(p,p+mplot) = 1;
        else
            bplot(p) = bplot(p) - uExactPlot(Nplot+1,i);
        end

    end
end

% Solve
uVecPlot = Aplot\bplot;

% Insert boundary values
uNumPlot = zeros(Nplot+1,Nplot+1);

uNumPlot(1,:)   = uExactPlot(1,:);
uNumPlot(end,:) = uExactPlot(end,:);
uNumPlot(:,1)   = uExactPlot(:,1);
uNumPlot(:,end) = uExactPlot(:,end);

uNumPlot(2:Nplot,2:Nplot) = ...
    reshape(uVecPlot,mplot,mplot).';


%% Figure: Representative Dirichlet solution

figure(1)
clf

surf(Xplot,Yplot,uNumPlot, ...
    'EdgeColor','none')

xlabel('$x$', ...
    'Interpreter','latex', 'FontSize',18)

ylabel('$y$', ...
    'Interpreter','latex', 'FontSize',18)

zlabel('$u_{\mathrm{Num}}(x,y)$', ...
    'Interpreter','latex', 'FontSize',18)

title('Poisson Solution with Dirichlet Conditions', ...
    'Interpreter','latex', 'FontSize',21)

set(gcf,'Units','inches')
set(gcf,'Position',[1 1 7 4.5])

%% Diverging blue-white-red colormap

nColors = 256;
nHalf = nColors/2;

% Muted endpoint colors
blue  = [0.15 0.35 0.60];
white = [0.95 0.95 0.95];
red   = [0.70 0.20 0.20];

% Blue -> white
lower = [linspace(blue(1),white(1),nHalf)', ...
         linspace(blue(2),white(2),nHalf)', ...
         linspace(blue(3),white(3),nHalf)'];

% White -> red
upper = [linspace(white(1),red(1),nHalf)', ...
         linspace(white(2),red(2),nHalf)', ...
         linspace(white(3),red(3),nHalf)'];

divMap = [lower; upper];

colormap(divMap)

% Make zero the center of the color scale
cmax = max(abs(uNumPlot(:)));
clim([-cmax cmax])

colormap(divMap)
cmax = max(abs(uNumPlot(:)));
clim([-cmax cmax])

cb = colorbar;
cb.FontSize = 14;

view(45,30)

grid on
box on

ax = gca;
ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';
ax.GridAlpha = 0.15;

set(gcf,'Units','inches')
set(gcf,'Position',[1 1 6.9 4.5])

exportgraphics(gcf, ...
    'figures/p4_dirichlet_solution.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');

%% Convergence of Dirichlet solution

figure(2)
loglog(Ns, errors, 'o-', 'LineWidth', 1.5, 'MarkerSize', 7)

xlabel('$N$ (grid intervals)','Interpreter','latex', 'FontSize',18)
ylabel('$\|u-u_{\mathrm{Num}}\|_\infty$', 'Interpreter', 'latex', 'FontSize',18)
title('Convergence of the 2D Poisson Solver', 'Interpreter', 'latex', 'FontSize',21)

hold on
loglog(Ns, ref, '--', ...
    'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.5)
hold off

legend('Numerical error', '$O(N^{-2})$', ...
    'Interpreter', 'latex', 'Location', 'best', 'FontSize',15)

grid on
box on

exportgraphics(gcf, ...
    'figures/p4_dirichlet_convergence.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');