function plan = buildfile()
%BUILDFILE Build plan for the Term Premium Modelling project.

% Build plan
plan = buildplan( localfunctions() );
plan.DefaultTasks = "export";

% List Live Script files for conversion
liveScript = fullfile( plan.RootFolder, "TermPremiumModelling.mlx" );

% Set task inputs
plan("md").Inputs = liveScript;

% Set task dependencies
plan("export").Dependencies = "md";

end % buildfile

function mdTask( context )
% Generate Markdown files from the MATLAB Live Scripts.

liveScripts = context.Task.Inputs;

for scriptIdx = 1 : numel( liveScripts )
    currentLiveScript = liveScripts(scriptIdx).Path;
    [currentPath, currentFilename] = fileparts( currentLiveScript );
    exportFile = fullfile( currentPath, currentFilename + ".md" );
    export( currentLiveScript, exportFile, ...
        "Format", "markdown", ...
        "IncludeOutputs", true, ...
        "Run", false, ...
        "AcceptHTML", true, ...
        "RenderLateXOnline", "svg" );
    fprintf( "[+] %s\n", exportFile )
end % for

end % mdTask

function exportTask( ~ )
% Export the project.

prj = currentProject();
mlproj = fullfile( prj.RootFolder, "ACM.mlproj" );
prj.export( mlproj )

end % exportTask