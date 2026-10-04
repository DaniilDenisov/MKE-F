param(
    [string]$BaseUri = 'http://127.0.0.1:8080'
)

$ErrorActionPreference = 'Stop'
$repositoryRoot = Split-Path -Parent $PSScriptRoot

$health = Invoke-RestMethod -Uri "$BaseUri/api/v1/health" -Method Get
if ($health.status -ne 'ready') {
    throw 'Solver health check did not report ready.'
}
Write-Host "Solver ready: $($health.octaveVersion)"

foreach ($asset in @('/preprocessor/', '/preprocessor/js/case-format.js', '/postprocessor/js/data-model.js')) {
    $response = Invoke-WebRequest -Uri "$BaseUri$asset" -Method Get
    if (($response.Headers['Cache-Control'] -join ',') -notmatch 'no-cache') {
        throw "Static asset $asset must require cache revalidation after rebuilds."
    }
}
Write-Host 'PASS static assets require cache revalidation'

$cases = @(
    @{ File = 'CasePreprocessorStatic.txt'; Type = 'static'; Version = 1 },
    @{ File = 'CasePreprocessorModal.txt'; Type = 'modal'; Version = 1 },
    @{ File = 'CasePreprocessorTransient.txt'; Type = 'transient'; Version = 1 },
    @{ File = 'CaseUniformFrame.txt'; Type = 'static'; Version = 2 }
)

foreach ($case in $cases) {
    $path = Join-Path $repositoryRoot "examples/cases/$($case.File)"
    $body = @{
        name = "Compose $($case.Type)"
        caseText = [System.IO.File]::ReadAllText($path)
    } | ConvertTo-Json -Compress
    $job = Invoke-RestMethod -Uri "$BaseUri/api/v1/jobs" -Method Post -ContentType 'application/json' -Body $body
    $deadline = [DateTime]::UtcNow.AddMinutes(2)
    do {
        Start-Sleep -Milliseconds 250
        $job = Invoke-RestMethod -Uri "$BaseUri/api/v1/jobs/$($job.id)" -Method Get
        if ([DateTime]::UtcNow -gt $deadline) {
            throw "Timed out waiting for $($case.Type) job."
        }
    } while ($job.status -in @('queued', 'running'))
    if ($job.status -ne 'succeeded') {
        throw "$($case.Type) job ended as $($job.status): $($job.error.message)"
    }
    $result = Invoke-RestMethod -Uri "$BaseUri/api/v1/jobs/$($job.id)/result" -Method Get
    if ($result.format -ne 'mkef-postprocessor' -or $result.version -ne $case.Version -or $result.analysis.type -ne $case.Type) {
        throw "Unexpected $($case.Type) result schema."
    }
    Write-Host "PASS $($case.File) (schema v$($case.Version))"
}

$solverExposed = $false
try {
    Invoke-WebRequest -Uri 'http://127.0.0.1:8000/api/v1/health' -TimeoutSec 1 | Out-Null
    $solverExposed = $true
} catch {
    $solverExposed = $false
}
if ($solverExposed) {
    throw 'The private solver service is unexpectedly reachable on host port 8000.'
}
Write-Host 'PASS solver port is private'
