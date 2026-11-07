#requires -Version 7.0
# Exercise process handling and inventory checks without invoking Docker.
param([string]$OutputDirectory = (Join-Path $PSScriptRoot '../../output/release-runner-checks'))
$ErrorActionPreference='Stop'
if ([Console]::IsOutputRedirected) { [Console]::OutputEncoding=[Text.UTF8Encoding]::new($false) }
$root=Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
New-Item -ItemType Directory -Path $OutputDirectory -Force | Out-Null
$currentLog=Join-Path $OutputDirectory 'native.log'
$StageTimeoutSeconds=10
$tokens=$null; $errors=$null
$ast=[Management.Automation.Language.Parser]::ParseFile((Join-Path $PSScriptRoot 'test_release.ps1'),[ref]$tokens,[ref]$errors)
if ($errors.Count) { throw ($errors | Out-String) }
foreach ($name in @('Invoke-Native','Assert-Inventory','Save-Report')) {
    $definition=$ast.Find({param($node) $node -is [Management.Automation.Language.FunctionDefinitionAst] -and $node.Name -eq $name},$false)
    . ([scriptblock]::Create($definition.Extent.Text))
}
$pwsh=Join-Path $PSHOME 'pwsh'
$echo=Join-Path $OutputDirectory 'argument probe.ps1'
'param([string]$Value); [Console]::OutputEncoding=[Text.UTF8Encoding]::new($false); [Console]::Out.WriteLine($Value); [Console]::Error.WriteLine("stderr probe")' | Set-Content $echo
$value='space кириллица "quote" $literal; & unchanged'
$result=Invoke-Native $pwsh @('-NoProfile','-File',$echo,'-Value',$value)
if ($result -cne $value -or (Get-Content $currentLog -Raw) -notmatch 'stderr probe') { throw 'Arguments or redirected streams were corrupted' }
$result=Invoke-Native $pwsh @('-NoProfile','-Command','[Console]::InputEncoding=[Text.UTF8Encoding]::new($false); [Console]::OutputEncoding=[Text.UTF8Encoding]::new($false); [Console]::Write([Console]::In.ReadToEnd())') -InputText $value
if ($result -cne $value) { throw 'Standard input was corrupted' }
$failed=$false
try { [void](Invoke-Native $pwsh @('-NoProfile','-Command','exit 7')) } catch { $failed=$_.Exception.Message -match 'code 7' }
if (!$failed) { throw 'Nonzero exit was accepted' }
$failed=$false
try { [void](Invoke-Native (Join-Path $OutputDirectory 'missing-executable') @()) } catch { $failed=$true }
if (!$failed) { throw 'Missing executable was accepted' }
$StageTimeoutSeconds=1
$failed=$false
try { [void](Invoke-Native $pwsh @('-NoProfile','-Command','Start-Sleep -Seconds 30')) } catch { $failed=$_.Exception.Message -match 'timed out' }
if (!$failed) { throw 'Timeout was not enforced' }

# Load the real plan without building images, then exercise each inventory guard.
. (Join-Path $PSScriptRoot 'test_release.ps1') -PlanOnly | Out-Null
$validManifest=$manifest | ConvertTo-Json -Depth 10
$validPlan=@($plan.ToArray())
function Assert-RejectedInventory([scriptblock]$Change, [string]$Expected) {
    $manifest=$validManifest | ConvertFrom-Json
    $plan=@($validPlan)
    & $Change
    $failed=$false
    try { Assert-Inventory } catch { $failed=$_.Exception.Message -match $Expected }
    if (!$failed) { throw "Invalid inventory was accepted: $Expected" }
}
Assert-RejectedInventory { $manifest.offline += 'test_deliberately_missing' } 'inventory mismatch'
Assert-RejectedInventory { $manifest.entryPoints.test_browser='tests/browser/test_deliberately_missing.cjs' } 'inventory mismatch'
Assert-RejectedInventory { $manifest.entryPoints.PSObject.Properties.Remove('test_browser') } 'inventory mismatch'
Assert-RejectedInventory { $manifest.offline += 'test_browser' } 'inventory mismatch'
Assert-RejectedInventory { $manifest.entryPoints.test_editor_adapter=$manifest.entryPoints.test_browser } 'inventory mismatch'
Assert-RejectedInventory { Set-Variable plan @($validPlan | Where-Object id -ne 'test_browser') -Scope 1 } 'Unscheduled test'
Assert-RejectedInventory { $manifest.extraOctave += 'tests/octave/test_deliberately_missing.m' } 'Octave inventory mismatch'
Assert-RejectedInventory { $manifest.pythonScripts += 'tests/integration/test_deliberately_missing.py' } 'Python script inventory mismatch'

# FAIL and BLOCKED must survive JSON/Markdown serialization, including reasons.
$runDir=$OutputDirectory; $started=[DateTime]::UtcNow; $reportStatus='FAIL'
$revision='probe'; $dirty=@(); $project='probe'; $images=@{}
$results=[System.Collections.Generic.List[object]]::new()
$results.Add([pscustomobject]@{id='failure';status='FAIL';seconds=1;log='native.log';error='failure | reason'})
$results.Add([pscustomobject]@{id='dependent';status='BLOCKED';seconds=0;log=$null;error='Prerequisites did not pass: failure'})
$plan=@(1,2)
Save-Report
$report=Get-Content (Join-Path $runDir 'report.json') -Raw | ConvertFrom-Json
if ($report.status -ne 'FAIL' -or $report.stages[1].status -ne 'BLOCKED' -or !$report.stages[0].error) { throw 'Failure report lost status or reason' }
Write-Host 'PASS release runner: argument quoting, stdout/stderr, exit codes, timeout, inventory rejection and failure reports.'
