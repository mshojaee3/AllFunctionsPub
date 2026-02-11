%
%% ==========================================================
% Local function: safe cleaning (never crashes if files are locked)
% ==========================================================
function safeCleanRunDir(runDir)

    if ~exist(runDir,'dir')
        mkdir(runDir);
        return;
    end

    fclose('all');
    pause(0.2);

    % Delete common Abaqus lock files first (can block reuse) :contentReference[oaicite:2]{index=2}
    delete(fullfile(runDir,'*.lck'));

    items = dir(runDir);
    for i = 1:numel(items)
        nm = items(i).name;
        if strcmp(nm,'.') || strcmp(nm,'..')
            continue;
        end

        fp = fullfile(runDir,nm);

        try
            if items(i).isdir
                % Try removing folder, but DON'T crash if Windows refuses
                [ok, msg] = rmdir(fp,'s');
                if ~ok
                    warning("Could not remove folder: %s\nReason: %s", fp, msg);
                end
            else
                delete(fp);
            end
        catch ME
            warning("Cleanup failed for: %s\n%s", fp, ME.message);
        end
    end

end
