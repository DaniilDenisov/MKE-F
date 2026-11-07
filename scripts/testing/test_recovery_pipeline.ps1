param(
    [string]$SolverImage = 'mkef-solver-test:latest',
    [string]$BrowserImage = 'mkef-browser-tests:latest'
)

$ErrorActionPreference = 'Stop'
$repositoryRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$runId = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$runDirectory = Join-Path $repositoryRoot "output/recovery-runs/$runId"
New-Item -ItemType Directory -Path $runDirectory -Force | Out-Null
$mounts = @('--rm', '--mount', "type=bind,source=$repositoryRoot,target=/app,readonly",
    '--mount', "type=bind,source=$runDirectory,target=/app/output", '--workdir', '/app')

function Invoke-ContainerCheck([string]$Image, [string]$EntryPoint, [string[]]$CheckArguments) {
    & docker run @mounts --entrypoint $EntryPoint $Image @CheckArguments
    if ($LASTEXITCODE -ne 0) {
        throw "Docker check failed ($LASTEXITCODE): $EntryPoint $CheckArguments. Artifacts: $runDirectory"
    }
}

Write-Host "Recovery checks, fresh output: $runDirectory"
Invoke-ContainerCheck $BrowserImage 'node' @('tests/node/test_spatial_condition_parity.cjs')
Invoke-ContainerCheck $SolverImage 'octave-cli' @('--no-gui', '--quiet', '--eval',
    "setup(); run_octave_tests(); test_spatial_condition_parity();")
Invoke-ContainerCheck $BrowserImage 'node' @('scripts/testing/run_recovery_checks.cjs')
Invoke-ContainerCheck $BrowserImage 'node' @('tests/browser/test_browser.cjs')
Invoke-ContainerCheck $SolverImage 'python3' @('-m','pytest','-q','-p','no:cacheprovider','backend/tests')
Write-Host "PASS: recovery checks in Docker. Artifacts: $runDirectory"
