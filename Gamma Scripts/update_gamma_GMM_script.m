%% Update gamma_GMM_fitLFP.m
%
% Creates a new version that:
%   1. Does not skip an electrode merely because its data do not cross 30 uV.
%   2. Fits K=1 first and retains K=1 statistics if K=2 fails.
%   3. Retains K=2 parameters in the table even when mu2 is rejected.
%   4. Displays mu2 only when the same final validity criteria pass.
%   5. Optionally requires the descriptive t-test to pass.
%   6. Includes one-component-only rows in mu1 summaries.
%
% The original source file is not modified.

clear
clc

sourceFile = "gamma_GMM_fitLFP.m";
outputFile = "gamma_GMM_fitLFP_fitAll_test.m";

if ~isfile(sourceFile)
    error("Source file not found: %s",sourceFile);
end

% Preserve blank lines and the original physical line structure.
lines = readlines(sourceFile,'EmptyLineRule','read');

fprintf('Reading source file:\n  %s\n',sourceFile);

%% ------------------------------------------------------------------------
% 1. Add shared validity settings before the fitting loop
% -------------------------------------------------------------------------

thresholdAnchor = find( ...
    strtrim(lines) == "minMu2_log = log10(minMu2_uV);", ...
    1,'first');

if isempty(thresholdAnchor)
    error('Could not find the minMu2_log definition.');
end

% Do not add the block twice.
alreadyHasEarlyCriteria = any( ...
    strtrim(lines(thresholdAnchor+1: ...
    min(thresholdAnchor+15,numel(lines)))) == ...
    "requireTtestForMu2 = true;");

if ~alreadyHasEarlyCriteria
    settingsBlock = [
        ""
        "% Operational criteria used consistently for plotting and final validity."
        "minWeight2    = 0.01;"
        "minCount2     = 50;"
        "minSeparation = 2;"
        ""
        "% If true, an unperformed or nonsignificant descriptive t-test hides mu2."
        "% Note: this test is descriptive because the fitted boundary defines"
        "% the two groups being compared."
        "requireTtestForMu2 = true;"
    ];

    lines = [
        lines(1:thresholdAnchor)
        settingsBlock
        lines(thresholdAnchor+1:end)
    ];
end

%% ------------------------------------------------------------------------
% 2. Preallocate model-completion and plotting columns
% -------------------------------------------------------------------------

preallocAnchor = find( ...
    strtrim(lines) == ...
    "gofTable.std2Larger = false(nRows,1);", ...
    1,'first');

if isempty(preallocAnchor)
    error('Could not find the std2Larger preallocation.');
end

if ~any(strtrim(lines) == ...
        "gofTable.fit1Completed = false(nRows,1);")

    preallocBlock = [
        "gofTable.fit1Completed = false(nRows,1);"
        "gofTable.fit2Completed = false(nRows,1);"
        "gofTable.plotMu2 = false(nRows,1);"
    ];

    lines = [
        lines(1:preallocAnchor)
        preallocBlock
        lines(preallocAnchor+1:end)
    ];
end

%% ------------------------------------------------------------------------
% 3. Replace threshold-dependent initialization and skipping
% -------------------------------------------------------------------------

initializationStart = find( ...
    contains(lines,"lowStartData  = myData"), ...
    1,'first');

if isempty(initializationStart)
    error(['Could not find the old lowStartData initialization. ' ...
        'Use the original/current gamma_GMM_fitLFP.m as the source.']);
end

componentProportionStart = find( ...
    (1:numel(lines))' > initializationStart & ...
    contains(lines,"start.ComponentProportion ="), ...
    1,'first');

if isempty(componentProportionStart)
    error('Could not find start.ComponentProportion.');
end

initializationEnd = componentProportionStart;

% Include the continuation line containing the two weights.
if initializationEnd < numel(lines) && ...
        contains(lines(initializationEnd+1), ...
        "[1-highStartWeight, highStartWeight]")
    initializationEnd = initializationEnd+1;
end

newInitializationBlock = [
    "    % -------------------------------------------------------------------------"
    "    % Attempt fitting for every electrode containing enough variable data."
    "    %"
    "    % The 30-uV threshold is applied after fitting. It is not a reason to"
    "    % skip an electrode or suppress its fitted parameter estimates."
    "    % -------------------------------------------------------------------------"
    "    dataVariance = var(myData,1);"
    ""
    "    if numel(myData) < 4 || ..."
    "            ~isfinite(dataVariance) || ..."
    "            dataVariance <= 0"
    ""
    "        gofTable.fitStatus(gofRow) = ""skipped"";"
    "        gofTable.skipReason(gofRow) = ..."
    "            ""insufficient_or_constant_data"";"
    ""
    "        gofTable.gmm{gofRow} = [];"
    "        gofTable.myData{gofRow} = myData;"
    "        gofTable.boundaries{gofRow} = [];"
    "        gofTable.originalBounds{gofRow} = [];"
    ""
    "        if isgraphics(tempFig)"
    "            close(tempFig);"
    "        end"
    ""
    "        continue"
    "    end"
    ""
    "    clear start"
    ""
    "    % Initialize both means from the complete distribution."
    "    start.mu = ["
    "        quantile(myData,0.25)"
    "        quantile(myData,0.75)"
    "    ];"
    ""
    "    dataSD = std(myData,1);"
    ""
    "    % Guard against invalid or identical starting means."
    "    if ~all(isfinite(start.mu)) || start.mu(2) <= start.mu(1)"
    "        start.mu = ["
    "            median(myData)-0.25*dataSD"
    "            median(myData)+0.25*dataSD"
    "        ];"
    "    end"
    ""
    "    % Positive, scale-dependent starting variances."
    "    var1Start = max(0.5*dataVariance,1e-4);"
    "    var2Start = max(dataVariance,(1.5^2)*var1Start);"
    ""
    "    start.Sigma = zeros(1,1,2);"
    "    start.Sigma(:,:,1) = var1Start;"
    "    start.Sigma(:,:,2) = var2Start;"
    ""
    "    % The observed fraction above 30 uV informs initialization but is not"
    "    % required in order to attempt the fit."
    "    highStartWeight = mean(myData > minMu2_log);"
    "    highStartWeight = min(max(highStartWeight,0.05),0.40);"
    ""
    "    start.ComponentProportion = ..."
    "        [1-highStartWeight, highStartWeight];"
];

lines = [
    lines(1:initializationStart-1)
    newInitializationBlock
    lines(initializationEnd+1:end)
];

%% ------------------------------------------------------------------------
% 4. Replace the active model-fitting block
% -------------------------------------------------------------------------

% Locate the active, uncommented K=2 call.
fit2Start = [];

for k = 1:numel(lines)
    thisLine = strtrim(lines(k));

    if startsWith(thisLine,"gmm = fitgmdist(myData,2")
        fit2Start = k;
        break
    end
end

if isempty(fit2Start)
    error('Could not find the active K=2 fitgmdist call.');
end

% Replace everything up to the active component-sorting statement.
sortLine = find( ...
    (1:numel(lines))' > fit2Start & ...
    contains(lines,"[mu,idx] = sort(gmm.mu(:),'ascend');"), ...
    1,'first');

if isempty(sortLine)
    error('Could not find the active component-sorting statement.');
end

newFitBlock = [
    "    % -------------------------------------------------------------------------"
    "    % Fit K=1 first so a valid mu1 can be retained if K=2 fails."
    "    % -------------------------------------------------------------------------"
    "    try"
    "        gmm1 = fitgmdist(myData,1, ..."
    "            'Options',options, ..."
    "            'RegularizationValue',regValue);"
    ""
    "        gofTable.fit1Completed(gofRow) = true;"
    ""
    "    catch ME1"
    "        gofTable.fitStatus(gofRow) = ""skipped"";"
    "        gofTable.skipReason(gofRow) = ..."
    "            ""K1_fit_failed: "" + string(ME1.identifier);"
    ""
    "        warning('K=1 fit failed for %s: %s', ..."
    "            axonNames(i),ME1.message);"
    ""
    "        gofTable.myData{gofRow} = myData;"
    "        gofTable.boundaries{gofRow} = [];"
    "        gofTable.originalBounds{gofRow} = [];"
    ""
    "        if isgraphics(tempFig)"
    "            close(tempFig);"
    "        end"
    ""
    "        continue"
    "    end"
    ""
    "    % -------------------------------------------------------------------------"
    "    % Attempt K=2. If it fails, retain and plot the K=1 result."
    "    % -------------------------------------------------------------------------"
    "    try"
    "        gmm = fitgmdist(myData,2, ..."
    "            'Options',options, ..."
    "            'Start',start, ..."
    "            'Replicates',1, ..."
    "            'RegularizationValue',regValue);"
    ""
    "        gofTable.fit2Completed(gofRow) = true;"
    "        gofTable.fitStatus(gofRow) = ""fitted"";"
    "        gofTable.skipReason(gofRow) = """";"
    ""
    "    catch ME2"
    "        gofTable.fitStatus(gofRow) = ""one_component_only"";"
    "        gofTable.skipReason(gofRow) = ..."
    "            ""K2_fit_failed: "" + string(ME2.identifier);"
    ""
    "        warning(['K=2 fit failed for %s. Retaining the ' ..."
    "            'K=1 result: %s'],axonNames(i),ME2.message);"
    ""
    "        oneMu = gmm1.mu(1);"
    "        oneVariance = gmm1.Sigma(1);"
    "        oneSD = sqrt(oneVariance);"
    ""
    "        gofTable.NegativeLogLikelihood(gofRow) = ..."
    "            gmm1.NegativeLogLikelihood;"
    "        gofTable.BIC1(gofRow) = gmm1.BIC;"
    "        gofTable.BIC2(gofRow) = NaN;"
    ""
    "        gofTable.mu1(gofRow) = oneMu;"
    "        gofTable.sigma1(gofRow) = oneVariance;"
    "        gofTable.std1(gofRow) = oneSD;"
    "        gofTable.weights1(gofRow) = 1;"
    "        gofTable.ampUnweighted1(gofRow) = ..."
    "            1/sqrt(2*pi*oneVariance);"
    "        gofTable.ampWeighted1(gofRow) = ..."
    "            gofTable.ampUnweighted1(gofRow);"
    ""
    "        % No K=2 estimate exists, so all component-2 fields stay NaN."
    "        gofTable.mu2(gofRow) = NaN;"
    "        gofTable.sigma2(gofRow) = NaN;"
    "        gofTable.std2(gofRow) = NaN;"
    "        gofTable.weights2(gofRow) = NaN;"
    ""
    "        gofTable.startMu1(gofRow) = start.mu(1);"
    "        gofTable.startMu2(gofRow) = start.mu(2);"
    "        gofTable.startStd1(gofRow) = ..."
    "            sqrt(start.Sigma(:,:,1));"
    "        gofTable.startStd2(gofRow) = ..."
    "            sqrt(start.Sigma(:,:,2));"
    ""
    "        gofTable.gmm{gofRow} = gmm1;"
    "        gofTable.myData{gofRow} = myData;"
    "        gofTable.boundaries{gofRow} = [];"
    "        gofTable.originalBounds{gofRow} = [];"
    ""
    "        % Plot the retained one-component result."
    "        oneFig = figure;"
    "        hold on"
    ""
    "        histogram(myData,log10(binEdges), ..."
    "            'Normalization','pdf', ..."
    "            'FaceColor',[0.55 0.78 1.00], ..."
    "            'EdgeColor',[0.10 0.25 0.40], ..."
    "            'FaceAlpha',0.65, ..."
    "            'DisplayName','Observed Data');"
    ""
    "        onePad = 0.10*range(myData);"
    ""
    "        if ~isfinite(onePad) || onePad <= 0"
    "            onePad = 0.1;"
    "        end"
    ""
    "        oneGrid = linspace( ..."
    "            min(myData)-onePad, ..."
    "            max(myData)+onePad,10000)';"
    ""
    "        onePdf = exp(-0.5*((oneGrid-oneMu).^2)./oneVariance) ./ ..."
    "            sqrt(2*pi*oneVariance);"
    ""
    "        plot(oneGrid,onePdf,'b-', ..."
    "            'LineWidth',3, ..."
    "            'DisplayName',sprintf( ..."
    "            'One-component fit: \mu=%.3f',oneMu));"
    ""
    "        xline(oneMu,'b:', ..."
    "            'LineWidth',2, ..."
    "            'HandleVisibility','off');"
    ""
    "        xlabel('Log_{10} \muV');"
    "        ylabel('Probability Density');"
    ""
    "        title({axonNames(i), ..."
    "            sprintf(['K=2 fit failed; showing K=1 only, ' ..."
    "            '\mu=%.3f, SD=%.3f'],oneMu,oneSD)});"
    ""
    "        legend('show','Location','best');"
    "        box on"
    "        hold off"
    ""
    "        saveas(oneFig, ..."
    "            fullfile(figOutDir, ..."
    "            ""GMM muSort "" + axonNames(i) + "".png""));"
    ""
    "        if isgraphics(oneFig)"
    "            close(oneFig);"
    "        end"
    ""
    "        if isgraphics(tempFig)"
    "            close(tempFig);"
    "        end"
    ""
    "        continue"
    "    end"
    ""
    "    % -------------------------------------------------------------------------"
    "    % Sort successful K=2 components by fitted mean."
    "    % Component 1 = lower-amplitude mode."
    "    % Component 2 = higher-amplitude mode."
    "    % -------------------------------------------------------------------------"
];

lines = [
    lines(1:fit2Start-1)
    newFitBlock
    lines(sortLine:end)
];

%% ------------------------------------------------------------------------
% 5. Replace the t-test-only showMu2 rule
% -------------------------------------------------------------------------

oldShowLine = find( ...
    strtrim(lines) == ...
    "showMu2 = ttestPerformed && (h == 1);", ...
    1,'first');

if isempty(oldShowLine)
    error('Could not find the old showMu2 rule.');
end

newShowBlock = [
    "    % -------------------------------------------------------------------------"
    "    % Apply the same criteria used later for final mu2Valid."
    "    % -------------------------------------------------------------------------"
    "    localThresholdPass = ..."
    "        isfinite(mu(2)) && ..."
    "        mu(2) > minMu2_log;"
    ""
    "    localHasBoundary = ..."
    "        ~isempty(boundaries) && ..."
    "        isfinite(boundaries(1)) && ..."
    "        boundaries(1) > mu(1) && ..."
    "        boundaries(1) < mu(2);"
    ""
    "    localBicPass = ..."
    "        isfinite(gmm1.BIC) && ..."
    "        isfinite(gmm.BIC) && ..."
    "        gmm.BIC < gmm1.BIC;"
    ""
    "    localWeightPass = ..."
    "        isfinite(weights(2)) && ..."
    "        weights(2) >= minWeight2 && ..."
    "        weights(2)*numel(myData) >= minCount2;"
    ""
    "    localSeparationPass = ..."
    "        isfinite(mu(1)) && ..."
    "        isfinite(mu(2)) && ..."
    "        isfinite(stdComponents(1)) && ..."
    "        stdComponents(1) > 0 && ..."
    "        mu(2) > mu(1) + ..."
    "            minSeparation*stdComponents(1);"
    ""
    "    localTtestPass = ..."
    "        ttestPerformed && ..."
    "        isfinite(h) && ..."
    "        h == 1;"
    ""
    "    if requireTtestForMu2"
    "        showMu2 = ..."
    "            localThresholdPass && ..."
    "            localHasBoundary && ..."
    "            localBicPass && ..."
    "            localWeightPass && ..."
    "            localSeparationPass && ..."
    "            localTtestPass;"
    "    else"
    "        showMu2 = ..."
    "            localThresholdPass && ..."
    "            localHasBoundary && ..."
    "            localBicPass && ..."
    "            localWeightPass && ..."
    "            localSeparationPass;"
    "    end"
    ""
    "    gofTable.plotMu2(gofRow) = showMu2;"
];

lines = [
    lines(1:oldShowLine-1)
    newShowBlock
    lines(oldShowLine+1:end)
];

%% ------------------------------------------------------------------------
% 6. Replace the individual-figure title logic
% -------------------------------------------------------------------------

% Locate the main two-component GMM figure, not the earlier K=1 fallback.
gmmFigureLine = find( ...
    strtrim(lines) == "gmmFig = figure;", ...
    1,'first');

if isempty(gmmFigureLine)
    error('Could not locate the main gmmFig figure creation.');
end

% Find the first legend only after the main gmmFig is created.
legendRelative = find( ...
    startsWith( ...
        strtrim(lines(gmmFigureLine+1:end)), ...
        "legend('show','Location','best')"), ...
    1,'first');

if isempty(legendRelative)
    error('Could not locate the main individual GMM figure legend.');
end

legendLine = gmmFigureLine + legendRelative;

% The main plot contains multiple "if showMu2" blocks. The final one before
% this legend is the title block.
lineNumbers = (1:numel(lines))';

showCandidates = find( ...
    lineNumbers > gmmFigureLine & ...
    lineNumbers < legendLine & ...
    strtrim(lines) == "if showMu2");

if isempty(showCandidates)
    error(['Could not locate the title-block "if showMu2" between ' ...
        'gmmFig creation and its legend.']);
end

titleStart = showCandidates(end);
titleEnd = legendLine-1;

newTitleBlock = [
    "    if showMu2"
    "        if ttestPerformed"
    "            pText = sprintf(', descriptive p=%.3e',p);"
    "        else"
    "            pText = '';"
    "        end"
    ""
    "        title({axonNames(i), ..."
    "            sprintf(['\mu_1=%.3f, \mu_2=%.3f, ' ..."
    "            'SD_1=%.3f, SD_2=%.3f, weights=%.2f/%.2f%s'], ..."
    "            mu(1),mu(2), ..."
    "            stdComponents(1),stdComponents(2), ..."
    "            weights(1),weights(2),pText)});"
    "    else"
    "        failedCriteria = strings(0,1);"
    ""
    "        if ~localThresholdPass"
    "            failedCriteria(end+1) = ""30-uV threshold"";"
    "        end"
    ""
    "        if ~localHasBoundary"
    "            failedCriteria(end+1) = ""boundary"";"
    "        end"
    ""
    "        if ~localBicPass"
    "            failedCriteria(end+1) = ""BIC"";"
    "        end"
    ""
    "        if ~localWeightPass"
    "            failedCriteria(end+1) = ""weight/count"";"
    "        end"
    ""
    "        if ~localSeparationPass"
    "            failedCriteria(end+1) = ""separation"";"
    "        end"
    ""
    "        if requireTtestForMu2 && ~localTtestPass"
    "            failedCriteria(end+1) = ""descriptive t-test"";"
    "        end"
    ""
    "        if isempty(failedCriteria)"
    "            failureText = ""unspecified criterion"";"
    "        else"
    "            failureText = strjoin(failedCriteria,"", "");"
    "        end"
    ""
    "        resultText = sprintf( ..."
    "            '\mu_1=%.3f, SD_1=%.3f; \mu_2 hidden: %s', ..."
    "            mu(1),stdComponents(1),char(failureText));"
    ""
    "        title({axonNames(i),resultText});"
    "    end"
    ""
];

lines = [
    lines(1:titleStart-1)
    newTitleBlock
    lines(titleEnd+1:end)
];

%% ------------------------------------------------------------------------
% 7. Remove duplicate late validity threshold definitions
% -------------------------------------------------------------------------

criteriaSection = find( ...
    strtrim(lines) == "%% Per-axon bimodality criteria", ...
    1,'first');

if isempty(criteriaSection)
    error('Could not locate the per-axon bimodality criteria section.');
end

duplicateLines = [
    "minWeight2    = 0.01;"
    "minCount2     = 50;"
    "minSeparation = 2;"
];

for d = 1:numel(duplicateLines)
    candidate = find( ...
        (1:numel(lines))' > criteriaSection & ...
        strtrim(lines) == duplicateLines(d), ...
        1,'first');

    if ~isempty(candidate)
        lines(candidate) = [];
    end
end

%% ------------------------------------------------------------------------
% 8. Add t-test agreement to final mu2Valid
% -------------------------------------------------------------------------

criteriaSection = find( ...
    strtrim(lines) == "%% Per-axon bimodality criteria", ...
    1,'first');

mu2ValidStart = find( ...
    (1:numel(lines))' > criteriaSection & ...
    strtrim(lines) == "gofTable.mu2Valid = ...", ...
    1,'first');

if isempty(mu2ValidStart)
    error('Could not locate the final mu2Valid definition.');
end

mu2ValidEnd = find( ...
    (1:numel(lines))' >= mu2ValidStart & ...
    strtrim(lines) == "sepPass;", ...
    1,'first');

if isempty(mu2ValidEnd)
    error('Could not locate the end of the final mu2Valid definition.');
end

newFinalValidityBlock = [
    "% Descriptive t-test criterion."
    "ttestPass = ..."
    "    gofTable.ttestPerformed & ..."
    "    isfinite(gofTable.ttest_h) & ..."
    "    gofTable.ttest_h == 1;"
    ""
    "gofTable.ttestPass = ttestPass;"
    ""
    "% Use the same rule as the individual figures."
    "if requireTtestForMu2"
    "    gofTable.mu2Valid = ..."
    "        wasFitted & ..."
    "        thresholdPass & ..."
    "        hasBoundary & ..."
    "        bicPass & ..."
    "        wPass & ..."
    "        sepPass & ..."
    "        ttestPass;"
    "else"
    "    gofTable.mu2Valid = ..."
    "        wasFitted & ..."
    "        thresholdPass & ..."
    "        hasBoundary & ..."
    "        bicPass & ..."
    "        wPass & ..."
    "        sepPass;"
    "end"
    ""
    "% Confirm that the figures and final validity decisions agree."
    "plotValidityMismatch = ..."
    "    gofTable.plotMu2 ~= gofTable.mu2Valid;"
    ""
    "fprintf('Plot/final-validity mismatches: %d\n', ..."
    "    sum(plotValidityMismatch));"
    ""
    "if any(plotValidityMismatch)"
    "    warning(['Some plotting decisions differ from mu2Valid. ' ..."
    "        'Inspect rows where plotMu2 ~= mu2Valid.']);"
    "end"
];

lines = [
    lines(1:mu2ValidStart-1)
    newFinalValidityBlock
    lines(mu2ValidEnd+1:end)
];

%% ------------------------------------------------------------------------
% 9. Make final workbook t-test messages more informative
% -------------------------------------------------------------------------

finalTableMarker = find( ...
    strtrim(lines) == "finalXlsxTable = gofTable;", ...
    1,'first');

if isempty(finalTableMarker)
    error('Could not locate finalXlsxTable creation.');
end

finalNotTestedLine = find( ...
    (1:numel(lines))' > finalTableMarker & ...
    strtrim(lines) == ...
    "finalXlsxTable.ttest_p_display(r) = ""not tested"";", ...
    1,'first');

if ~isempty(finalNotTestedLine)
    % The preceding line is:
    % if ~finalXlsxTable.ttestPerformed(r)
    ifLine = finalNotTestedLine-1;

    replacement = [
        "    if ~finalXlsxTable.ttestPerformed(r)"
        "        if finalXlsxTable.fitStatus(r) == ""skipped"""
        "            finalXlsxTable.ttest_p_display(r) = ..."
        "                ""not tested: "" + finalXlsxTable.skipReason(r);"
        ""
        "        elseif finalXlsxTable.fitStatus(r) == ""one_component_only"""
        "            finalXlsxTable.ttest_p_display(r) = ..."
        "                ""not tested: K=1 fit only"";"
        ""
        "        elseif ~finalXlsxTable.hasBoundary(r)"
        "            finalXlsxTable.ttest_p_display(r) = ..."
        "                ""not tested: no between-mean boundary"";"
        ""
        "        else"
        "            finalXlsxTable.ttest_p_display(r) = ..."
        "                ""not tested: insufficient boundary groups"";"
        "        end"
    ];

    lines = [
        lines(1:ifLine-1)
        replacement
        lines(finalNotTestedLine+1:end)
    ];
end

%% ------------------------------------------------------------------------
% 10. Include K=1-only rows in active mu1 summaries
% -------------------------------------------------------------------------

% Replace active validMu1Rows blocks. Commented historical code is ignored.
k = 1;

while k <= numel(lines)-2
    if strtrim(lines(k)) == "validMu1Rows = regRows( ..." && ...
            contains(strtrim(lines(k+1)), ...
            'gofTable.fitStatus(regRows) == "fitted"') && ...
            ~startsWith(strtrim(lines(k)),"%")

        oldSecondLine = k+1;

        replacement = [
            "    validMu1Rows = regRows( ..."
            "        (gofTable.fitStatus(regRows) == ""fitted"" | ..."
            "         gofTable.fitStatus(regRows) == ""one_component_only"") & ..."
            "        isfinite(gofTable.mu1(regRows)));"
        ];

        % Existing block occupies three lines.
        lines = [
            lines(1:k-1)
            replacement
            lines(k+3:end)
        ];

        k = k+numel(replacement);
    else
        k = k+1;
    end
end

% Replace the active bar-chart mu1 count condition.
countStatusLine = find( ...
    contains(lines,'gofTable.fitStatus == "fitted" & ...') & ...
    ~startsWith(strtrim(lines),"%"), ...
    1,'last');

if ~isempty(countStatusLine)
    lines(countStatusLine) = ...
        "        (gofTable.fitStatus == ""fitted"" | ...";

    lines = [
        lines(1:countStatusLine)
        "         gofTable.fitStatus == ""one_component_only"") & ..."
        lines(countStatusLine+1:end)
    ];
end

%% ------------------------------------------------------------------------
% 11. Write the updated test script
% -------------------------------------------------------------------------

writelines(lines,outputFile);

fprintf('\nUpdate completed successfully.\n');
fprintf('Original file was not modified:\n  %s\n',sourceFile);
fprintf('Updated test script created:\n  %s\n\n',outputFile);

fprintf('Expected behavior:\n');
fprintf(['  fitted + valid mu2       -> table has mu1/mu2; ' ...
    'plot shows both\n']);
fprintf(['  fitted + rejected mu2    -> table has mu1/mu2; ' ...
    'plot shows mu1 only\n']);
fprintf(['  K=1 succeeds, K=2 fails  -> table has mu1 only; ' ...
    'plot shows K=1 only\n']);
fprintf(['  unusable/constant data   -> row remains skipped with ' ...
    'missing fit values\n\n']);

fprintf('After the generated script runs, look for:\n');
fprintf('  Plot/final-validity mismatches: 0\n');