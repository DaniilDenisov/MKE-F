#requires -Version 7.0
<#
.SYNOPSIS
Build and test the complete release in Docker, without publishing it.
.DESCRIPTION
Run from any directory. Each invocation gets fresh fixtures, logs and a private
Compose project. Independent stages continue after failures; dependent stages
are BLOCKED. PASS requires every planned stage, including cleanup, to pass.
#>
[CmdletBinding()]
param(
    [switch]$PlanOnly,
    [ValidateRange(60, 14400)][int]$StageTimeoutSeconds = 3600
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$runId = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$project = "mkef-release-$runId"
$runDir = Join-Path $root "output/release-runs/$runId"
$fixtures = Join-Path $runDir 'artifacts'
$manifest = Get-Content -LiteralPath (Join-Path $PSScriptRoot 'release-tests.json') -Raw | ConvertFrom-Json
$plan = [System.Collections.Generic.List[object]]::new()
$results = [System.Collections.Generic.List[object]]::new()
$states = @{}
$images = @{}
$builtImages = [System.Collections.Generic.List[string]]::new()
foreach ($name in @('solver-test','browser','solver','section','web')) { $images[$name] = "${project}-${name}:local" }
$composeFile = Join-Path $runDir 'compose.override.json'
$composeArgs = @('compose','-f',(Join-Path $root 'compose.yaml'),'-f',$composeFile,'-p',$project)
$webId = $null
$solverId = $null
$stackAttempted = $false
$activeContainer = $null
$started = [DateTime]::UtcNow
$reportStatus = 'RUNNING'
$currentLog = $null

function Add-Stage($Id, $Kind, $Arguments = @(), $Depends = @(), $Image = '', $Network = '', $Environment = @(), $Output = '') {
    $plan.Add([pscustomobject]@{id=$Id; kind=$Kind; arguments=@($Arguments); depends=@($Depends); image=$Image; network=$Network; environment=@($Environment); output=$Output})
}
function Add-Node($Name, $Depends, $Network = '', $Environment = @(), $Output = '') {
    Add-Stage $Name 'container' @('node',$manifest.entryPoints.$Name) (@('build-browser') + $Depends) 'browser' $Network $Environment $Output
}
function Add-Octave($Id, $Code, $Depends) {
    Add-Stage $Id 'container' @('octave-cli','--no-gui','--quiet','--eval',"setup(); addpath('/app/tests/octave'); $Code") (@('build-solver-test') + $Depends) 'solver-test'
}

Add-Stage 'inventory' 'inventory'
Add-Stage 'runner-checks' 'runner-checks' @() @('inventory')
foreach ($name in @('solver-test','browser','solver','section','web')) {
    $dockerfile = switch ($name) { 'solver-test' {'solver'} 'browser' {'browser-tests'} 'section' {'section-solver'} default {$name} }
    $buildArgs = @('build','-f',"docker/$dockerfile.Dockerfile",'-t',$images[$name])
    if ($name -eq 'solver-test') { $buildArgs += @('--target','test') }
    if ($name -eq 'browser') { $buildArgs += @('--target','release') }
    if ($name -in @('solver','section')) { $buildArgs += @('--target','runtime') }
    Add-Stage "build-$name" 'docker' ($buildArgs + @('.')) @('runner-checks')
}
Add-Stage 'spatial-bundle' 'container' @('node','scripts/build/build-spatial.cjs','--check') @('build-browser') 'browser'
foreach ($name in $manifest.beforeOctave) { Add-Node $name @('inventory') }
Add-Octave 'octave' 'run_octave_tests();' @('inventory')
Add-Octave 'extra-fixtures' 'prepare_release_fixtures();' @('octave')
Add-Octave 'nafems-profiles' "run_nafems_challenge5('quick'); run_nafems_challenge5('full');" @('octave')
Add-Octave 'condition-parity' 'test_spatial_condition_parity();' @('octave','test_spatial_condition_parity')
foreach ($name in $manifest.recovery) {
    $dependencies=@('octave')
    if ($name -eq 'test_spatial_conditions') { $dependencies += 'condition-parity' }
    Add-Node $name $dependencies
}
foreach ($name in $manifest.offline) {
    $dependencies=@('octave','extra-fixtures')
    if ($name -eq 'test_nafems_browser') { $dependencies += 'nafems-profiles' }
    Add-Node $name $dependencies
}
Add-Octave 'downloaded-examples' 'test_loaded_examples();' @('test_examples_browser')
Add-Stage 'mises-offline' 'container' @('node','tests/browser/test_mises_browser.cjs') @('build-browser','octave','extra-fixtures') 'browser'
Add-Stage 'backend' 'container' @('python3','-m','pytest','-q','-p','no:cacheprovider','backend/tests','--junitxml=/app/output/backend.xml') @('build-solver-test') 'solver-test'
Add-Stage 'section-unit' 'container' @('python','-m','unittest','discover','-s','tests/section_solver','-v') @('build-section') 'section'
Add-Stage 'section-api-unit' 'container' @('python','-m','unittest','discover','-s','tests/section_service','-v') @('build-section') 'section'
Add-Stage 'rounded-parity' 'container' @('python','tests/integration/test_rounded_section_parity.py') @('build-section','octave') 'section'
Add-Stage 'section-matrix' 'container' @('python','-m','tests.integration.section_solver_acceptance.matrix') @('build-section') 'section' '' @('SECTION_REPORT_DIR=output/section-acceptance')
Add-Stage 'section-acceptance' 'container' @('python','-m','tests.integration.section_solver_acceptance.verify') @('section-matrix') 'section' '' @('SECTION_REPORT_DIR=output/section-acceptance')
Add-Stage 'stack-up' 'stack-up' @() @('build-solver','build-section','build-web')
Add-Stage 'compose-smoke' 'smoke' @() @('stack-up')
Add-Stage 'reference-hosted' 'container' @('node','tests/browser/test_reference_browser.cjs') @('stack-up','build-browser','octave') 'browser' 'web' @('MKEF_REFERENCE_BASE_URL=http://localhost:8080')
Add-Node 'test_v8_api_browser' @('stack-up') 'web' @() (Join-Path $runDir 'runtime-artifacts')
Add-Node 'test_spatial_real_api' @('stack-up') 'solver'
foreach ($name in @('test_section_inspector_browser','test_mises_browser','test_section_diagnostic_browser','test_rounded_sections_browser')) {
    Add-Stage "reset-$name" 'compose' @('restart','--no-deps','section-solver') @('stack-up')
    Add-Stage "ready-$name" 'compose' @('up','-d','--no-build','--wait','--wait-timeout','180','section-solver') @("reset-$name")
    Add-Node $name @("ready-$name",'octave','extra-fixtures') 'web' @('MKEF_POST_BASE_URL=http://localhost:8080','MKEF_TEST_BACKEND=http://localhost:8080')
}
Add-Stage 'proxy-up' 'proxy' @('up') @('stack-up')
Add-Stage 'section-stop' 'compose' @('stop','section-solver') @('stack-up')
Add-Stage 'proxy-down' 'proxy' @('down') @('section-stop')
Add-Stage 'diagnostic-down' 'container' @('node','tests/browser/test_section_diagnostic_browser.cjs','down') @('section-stop','build-browser') 'browser' 'web'
Add-Stage 'section-recover' 'compose' @('up','-d','--no-build','--wait','--wait-timeout','180','section-solver') @('stack-up')
Add-Stage 'proxy-recovered' 'proxy' @('recovered') @('section-recover')
Add-Stage 'diagnostic-recovered' 'container' @('node','tests/browser/test_section_diagnostic_browser.cjs','recovered') @('section-recover','build-browser') 'browser' 'web'
Add-Stage 'retired-assets' 'compose' @('exec','-T','web','sh','-ec','cd /usr/share/nginx/html && test ! -e shared/section-fem && test ! -e shared/mises.js && test ! -e postprocessor/js/mises-worker.js && test ! -e scripts/section_fem_fvm_reference.cjs') @('stack-up')
Add-Stage 'cleanup' 'cleanup'

function Assert-Inventory {
    $listed = @($manifest.beforeOctave) + @($manifest.recovery) + @($manifest.offline) + @($manifest.live)
    $actual = @(Get-ChildItem -Path (Join-Path $root 'tests'),(Join-Path $root 'scripts') -Recurse -Filter 'test_*.cjs' | ForEach-Object { [IO.Path]::GetRelativePath($root, $_.FullName).Replace('\','/') })
    $entries = @($manifest.entryPoints.PSObject.Properties)
    $paths = @($entries | ForEach-Object Value)
    $difference = @(Compare-Object ($paths | Sort-Object) ($actual | Sort-Object))
    if ($difference.Count -or @(Compare-Object ($listed | Sort-Object) ($entries.Name | Sort-Object)).Count -or
        @($listed | Group-Object | Where-Object Count -gt 1).Count -or
        @($paths | Group-Object | Where-Object Count -gt 1).Count) {
        throw "JavaScript test inventory mismatch. Update scripts/testing/release-tests.json: $($difference | Out-String)"
    }
    foreach ($entry in $entries) {
        if ($entry.Value -cnotmatch '^tests/(browser|node)/test_\w+\.cjs$' -or
            [IO.Path]::GetFileNameWithoutExtension($entry.Value) -cne $entry.Name -or
            !(Test-Path -LiteralPath (Join-Path $root $entry.Value) -PathType Leaf)) { throw "JavaScript test inventory mismatch: invalid entry $($entry.Name)" }
    }
    $scheduled = @($plan | Where-Object { $_.kind -eq 'container' -and $_.arguments[0] -eq 'node' } | ForEach-Object { $_.arguments[1] })
    foreach ($name in $listed) { if ($manifest.entryPoints.$name -cnotin $scheduled) { throw "Unscheduled test: $name" } }
    $runner = Get-Content -LiteralPath (Join-Path $root 'tests/octave/run_verification_tests.m') -Raw
    $octave = @([regex]::Matches($runner, '@(test_\w+)') | ForEach-Object { 'tests/octave/' + $_.Groups[1].Value + '.m' }) + @($manifest.extraOctave)
    $actualM = @(Get-ChildItem -Path (Join-Path $root 'tests'),(Join-Path $root 'scripts') -Recurse -Filter 'test_*.m' | ForEach-Object { [IO.Path]::GetRelativePath($root, $_.FullName).Replace('\','/') })
    if (@(Compare-Object ($octave | Sort-Object) ($actualM | Sort-Object)).Count) { throw 'Octave inventory mismatch: register new tests in run_verification_tests.m or the release plan.' }
    foreach ($file in $manifest.extraOctave) {
        $function=[IO.Path]::GetFileNameWithoutExtension($file)
        $fixtureCode=Get-Content (Join-Path $root 'scripts/testing/prepare_release_fixtures.m') -Raw
        if ($fixtureCode -notmatch ([regex]::Escape($function) + '\(\)') -and !@($plan | Where-Object { ($_.arguments -join ' ') -match ([regex]::Escape($function) + '\(\)') }).Count) { throw "Unscheduled Octave test: $file" }
    }
    $actualPy = @(Get-ChildItem (Join-Path $root 'tests/integration') -Recurse -Filter 'test_*.py' | ForEach-Object { [IO.Path]::GetRelativePath($root, $_.FullName).Replace('\','/') })
    if (@(Compare-Object @($manifest.pythonScripts) $actualPy).Count) { throw 'Python script inventory mismatch.' }
    foreach ($file in $manifest.pythonScripts) {
        if (!@($plan | Where-Object { $file -in $_.arguments }).Count) { throw "Unscheduled Python test: $file" }
    }
    $actualPs = @(Get-ChildItem (Join-Path $root 'scripts/testing') -Filter 'test_*.ps1' | ForEach-Object Name)
    if (@(Compare-Object @($manifest.compositeScripts.PSObject.Properties.Name) $actualPs).Count) { throw 'PowerShell test inventory mismatch: classify the new runner in compositeScripts.' }
    $seen = @{}
    foreach ($stage in $plan) {
        if ($seen.ContainsKey($stage.id)) { throw "Duplicate stage: $($stage.id)" }
        foreach ($dependency in $stage.depends) { if (!$seen.ContainsKey($dependency)) { throw "Dependency must precede $($stage.id): $dependency" } }
        $seen[$stage.id] = $true
    }
    Write-Host "Inventory: $($actual.Count) JavaScript entry points, $($actualM.Count) Octave test files; Python suites use discovery."
}

if ($PlanOnly) {
    Assert-Inventory
    $plan | Select-Object id,kind,@{n='depends';e={$_.depends -join ', '}} | Format-Table -AutoSize
    return
}

New-Item -ItemType Directory -Path $fixtures,(Join-Path $runDir 'logs') -Force | Out-Null
$configuration = @{services=@{
    solver=@{image=$images.solver; environment=@{MKEF_ALLOWED_ORIGINS='http://localhost:8080,http://127.0.0.1:8080'}}
    'section-solver'=@{image=$images.section}; web=@{image=$images.web}
}}
$configuration | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $composeFile -Encoding utf8
$plan | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $runDir 'plan.json') -Encoding utf8
$revision = (& git -C $root rev-parse HEAD 2>$null) -join ''
$dirty = @(& git -C $root status --porcelain=v1 2>$null)

function Save-Report {
    $counts=[ordered]@{}
    foreach ($status in @('PASS','FAIL','BLOCKED','RUNNING','ABORTED')) { $counts[$status]=@($results | Where-Object status -eq $status).Count }
    $report = [ordered]@{schemaVersion=1; status=$reportStatus; startedAt=$started.ToString('o'); updatedAt=[DateTime]::UtcNow.ToString('o'); durationSeconds=[math]::Round(([DateTime]::UtcNow-$started).TotalSeconds,2); revision=$revision; workingTreeChanges=$dirty; project=$project; images=$images; plannedStages=$plan.Count; counts=$counts; stages=@($results.ToArray())}
    $report | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $runDir 'report.json') -Encoding utf8
    $lines = @('# Release acceptance', '', "Status: **$reportStatus**", '', "Stages: $($plan.Count) planned; $($counts.PASS) passed; $($counts.FAIL) failed; $($counts.BLOCKED) blocked.", '', "Revision: $revision", '', '| Stage | Status | Seconds | Log / reason |', '|---|---|---:|---|')
    foreach ($r in $results) {
        $detail = if ($r.log) { "[log]($($r.log))" } else { '' }
        if ($r.error) { $detail += ' ' + ($r.error -replace '[\r\n|]', ' ') }
        $lines += "| $($r.id) | $($r.status) | $($r.seconds) | $detail |"
    }
    $lines | Set-Content -LiteralPath (Join-Path $runDir 'report.md') -Encoding utf8
}

# ArgumentList avoids shell quoting and works with spaces in repository paths.
# Drain both pipes while running, so verbose Octave/Chromium cannot deadlock.
function Invoke-Native([string]$File, [string[]]$Arguments, [string]$InputText) {
    $hasInput=$PSBoundParameters.ContainsKey('InputText')
    $info = [Diagnostics.ProcessStartInfo]::new()
    $info.FileName=$File; $info.WorkingDirectory=$root; $info.UseShellExecute=$false
    $info.CreateNoWindow=$true; $info.RedirectStandardOutput=$true; $info.RedirectStandardError=$true
    $info.StandardOutputEncoding=[Text.UTF8Encoding]::new($false)
    $info.StandardErrorEncoding=[Text.UTF8Encoding]::new($false)
    $info.RedirectStandardInput=$hasInput
    if ($hasInput) { $info.StandardInputEncoding=[Text.UTF8Encoding]::new($false) }
    foreach ($argument in $Arguments) { $info.ArgumentList.Add($argument) }
    $info.Environment['MKEF_PORT']='0'
    $process = [Diagnostics.Process]::new(); $process.StartInfo=$info
    $writer = [IO.StreamWriter]::new($currentLog,$true); $writer.AutoFlush=$true
    $clock = [Diagnostics.Stopwatch]::StartNew()
    $stdout = [Text.StringBuilder]::new()
    $processStarted=$false
    try {
        $writer.WriteLine("COMMAND: $File " + ($Arguments -join ' '))
        [void]$process.Start()
        $processStarted=$true
        if ($hasInput) { $process.StandardInput.Write($InputText); $process.StandardInput.Close() }
        $outTask=$process.StandardOutput.ReadLineAsync(); $errTask=$process.StandardError.ReadLineAsync()
        while ($null -ne $outTask -or $null -ne $errTask -or !$process.HasExited) {
            if ($clock.Elapsed.TotalSeconds -gt $StageTimeoutSeconds) { throw "Command timed out after $StageTimeoutSeconds seconds" }
            $drainedAny=$false
            foreach ($stream in @('out','err')) {
                $task = if ($stream -eq 'out') {$outTask} else {$errTask}
                $drained=0
                while ($null -ne $task -and $task.IsCompleted -and $drained -lt 128) {
                    $drained++; $drainedAny=$true
                    $line=$task.GetAwaiter().GetResult()
                    if ($null -ne $line) {
                        $writer.WriteLine($line); Write-Host $line
                        if ($stream -eq 'out') { [void]$stdout.AppendLine($line) }
                    }
                    if ($stream -eq 'out') { $outTask=if ($null -ne $line) {$process.StandardOutput.ReadLineAsync()} else {$null} }
                    else { $errTask=if ($null -ne $line) {$process.StandardError.ReadLineAsync()} else {$null} }
                    $task = if ($stream -eq 'out') {$outTask} else {$errTask}
                }
            }
            if (!$drainedAny) { Start-Sleep -Milliseconds 10 }
        }
        $process.WaitForExit()
        if ($process.ExitCode -ne 0) { throw "$File exited with code $($process.ExitCode)" }
        return $stdout.ToString().Trim()
    } finally {
        if ($processStarted -and !$process.HasExited) { $process.Kill($true); $process.WaitForExit() }
        $process.Dispose(); $writer.Dispose()
    }
}

function Invoke-Stage($Stage) {
    switch ($Stage.kind) {
        'inventory' { Assert-Inventory }
        'runner-checks' { [void](Invoke-Native (Join-Path $PSHOME 'pwsh') @('-NoProfile','-File',(Join-Path $PSScriptRoot 'check_release_runner.ps1'),'-OutputDirectory',(Join-Path $runDir 'runner-checks'))) }
        'docker' {
            [void](Invoke-Native 'docker' $Stage.arguments)
            if ($Stage.id.StartsWith('build-')) { $builtImages.Add($images[$Stage.id.Substring(6)]) }
        }
        'container' {
            $directory = if ($Stage.output) {$Stage.output} else {$fixtures}
            New-Item -ItemType Directory -Path $directory -Force | Out-Null
            $script:activeContainer = "$project-check"
            # One UID for producers and consumers: browser-created directories
            # must remain writable by Octave and Python on Linux bind mounts.
            $arguments = @('run','--rm','--name',$activeContainer,'--user','10001:10001','-e','HOME=/tmp','-e','PYTHONPATH=/app','--mount',"type=bind,source=$root,target=/app,readonly",'--mount',"type=bind,source=$directory,target=/app/output",'--workdir','/app','--shm-size','1g')
            if ($Stage.network) {
                $containerId=if ($Stage.network -eq 'web') {$webId} else {$solverId}
                $arguments += @('--network',"container:$containerId")
            }
            foreach ($variable in $Stage.environment) { $arguments += @('-e',$variable) }
            $arguments += @('--entrypoint',$Stage.arguments[0],$images[$Stage.image]) + @($Stage.arguments | Select-Object -Skip 1)
            try { [void](Invoke-Native 'docker' $arguments) }
            finally {
                # Only the per-run test container; never a user's container.
                & docker rm -f $activeContainer 2>$null | Out-Null
                $script:activeContainer=$null
            }
        }
        'compose' { [void](Invoke-Native 'docker' ($composeArgs + $Stage.arguments)) }
        'stack-up' {
            $script:stackAttempted=$true
            [void](Invoke-Native 'docker' ($composeArgs + @('up','-d','--no-build','--wait','--wait-timeout','180')))
            $script:webId=Invoke-Native 'docker' ($composeArgs + @('ps','-q','web'))
            $script:solverId=Invoke-Native 'docker' ($composeArgs + @('ps','-q','solver'))
            if (!$webId -or !$solverId) { throw 'Missing test stack container IDs' }
        }
        'smoke' {
            $address=Invoke-Native 'docker' ($composeArgs + @('port','web','8080'))
            # Inspect this project's solver, not an unrelated service on host :8000.
            $ports=(Invoke-Native 'docker' @('inspect',$solverId,'--format','{{json .NetworkSettings.Ports}}')) | ConvertFrom-Json
            foreach ($port in $ports.PSObject.Properties) { if ($port.Value) { throw "Solver has a published host port: $($port.Name)" } }
            [void](Invoke-Native (Join-Path $PSHOME 'pwsh') @('-NoProfile','-File',(Join-Path $PSScriptRoot 'compose_smoke_test.ps1'),'-BaseUri',"http://$address",'-SkipHostPortProbe'))
        }
        'proxy' {
            $probe=Get-Content -LiteralPath (Join-Path $root 'tests/section_service/proxy_probe.py') -Raw
            [void](Invoke-Native 'docker' ($composeArgs + @('exec','-T','solver','python3','-') + $Stage.arguments) -InputText $probe)
        }
        'cleanup' {
            if ($stackAttempted) {
                try { [void](Invoke-Native 'docker' ($composeArgs + @('logs','--no-color'))) }
                finally { [void](Invoke-Native 'docker' ($composeArgs + @('down','--volumes','--remove-orphans'))) }
                $script:stackAttempted=$false
            }
            if ($builtImages.Count) {
                [void](Invoke-Native 'docker' (@('image','rm') + $builtImages.ToArray()))
                $builtImages.Clear()
            }
        }
    }
}

Push-Location $root
try {
    Write-Host "Release acceptance: $runDir"
    foreach ($stage in $plan) {
        $blocked = @($stage.depends | Where-Object {$states[$_] -ne 'PASS'})
        $row=[pscustomobject]@{id=$stage.id; status='RUNNING'; seconds=0; log="logs/$($stage.id).log"; error=$null}
        $results.Add($row)
        $script:currentLog=Join-Path $runDir $row.log
        if ($blocked.Count) {
            $row.status='BLOCKED'; $row.error='Prerequisites did not pass: ' + ($blocked -join ', '); $row.log=$null
        } else {
            New-Item -ItemType File -Path $currentLog -Force | Out-Null
            Save-Report
            Write-Host "[$($results.Count)/$($plan.Count)] $($stage.id)"
            $timer=[Diagnostics.Stopwatch]::StartNew()
            try { Invoke-Stage $stage; $row.status='PASS' }
            catch { $row.status='FAIL'; $row.error=$_.Exception.Message; Add-Content -LiteralPath $currentLog -Value $row.error; Write-Host "FAIL: $($row.error)" }
            $row.seconds=[math]::Round($timer.Elapsed.TotalSeconds,2)
        }
        $states[$stage.id]=$row.status
        Save-Report
    }
    $reportStatus=if (@($results | Where-Object status -ne 'PASS').Count) {'FAIL'} else {'PASS'}
} finally {
    if ($reportStatus -eq 'RUNNING') {
        $reportStatus='ABORTED'
        foreach ($unfinished in @($results | Where-Object status -eq 'RUNNING')) { $unfinished.status='ABORTED'; $unfinished.error='Run interrupted' }
    }
    if ($stackAttempted) {
        try { $currentLog=Join-Path $runDir 'logs/emergency-cleanup.log'; [void](Invoke-Native 'docker' ($composeArgs + @('down','--volumes','--remove-orphans'))) }
        catch { $reportStatus='FAIL'; Write-Warning "Cleanup failed: $_" }
    }
    if ($builtImages.Count) {
        try {
            $currentLog=Join-Path $runDir 'logs/emergency-cleanup.log'
            [void](Invoke-Native 'docker' (@('image','rm') + $builtImages.ToArray()))
            $builtImages.Clear()
        } catch { $reportStatus='FAIL'; Write-Warning "Image cleanup failed: $_" }
    }
    Save-Report
    Pop-Location
}
Write-Host "$reportStatus - $(Join-Path $runDir 'report.md')"
if ($reportStatus -ne 'PASS') { throw "Release acceptance $reportStatus. See $runDir/report.md" }
