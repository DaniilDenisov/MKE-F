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

$cases = @(
    @{ File = 'CasePreprocessorStatic.txt'; Type = 'static' },
    @{ File = 'CasePreprocessorModal.txt'; Type = 'modal' },
    @{ File = 'CasePreprocessorTransient.txt'; Type = 'transient' }
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
    if ($result.format -ne 'mkef-postprocessor' -or $result.version -ne 1 -or $result.analysis.type -ne $case.Type) {
        throw "Unexpected $($case.Type) result schema."
    }
    Write-Host "PASS $($case.Type)"
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
