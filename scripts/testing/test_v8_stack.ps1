param(
    [string]$SolverImage = 'mkef-v8-solver:verification',
    [string]$WebImage = 'mkef-v8-web:verification',
    [string]$BrowserImage = 'mkef-browser-tests:latest'
)
$ErrorActionPreference = 'Stop'
$repositoryRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$project = 'mkef-v8-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$runDirectory = Join-Path $repositoryRoot "output/$project"
New-Item -ItemType Directory -Path $runDirectory | Out-Null
function Invoke-Docker([string[]]$DockerArguments) {
    & docker @DockerArguments
    if ($LASTEXITCODE -ne 0) { throw "Docker verification failed: $DockerArguments" }
}
Invoke-Docker @('build','-f',"$repositoryRoot/docker/solver.Dockerfile",'--target','runtime','-t',$SolverImage,$repositoryRoot)
Invoke-Docker @('build','-f',"$repositoryRoot/docker/web.Dockerfile",'-t',$WebImage,$repositoryRoot)
# Isolated stack, no host ports and no production volumes/containers.
$configuration = @{
    services = @{
        solver = @{ image=$SolverImage; tmpfs=@('/jobs:uid=10001,gid=10001,mode=0700','/tmp:mode=1777') }
        web = @{ image=$WebImage; depends_on=@('solver') }
    }
}
$composeFile = Join-Path $runDirectory 'compose.json'
$configuration | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $composeFile -Encoding utf8
try {
    Invoke-Docker @('compose','-f',$composeFile,'-p',$project,'up','-d')
    $webContainer = & docker compose -f $composeFile -p $project ps -q web
    if ($LASTEXITCODE -ne 0 -or !$webContainer) { throw 'No verification web container.' }
    # Share web's network namespace: nginx deliberately accepts localhost only.
    Invoke-Docker @('run','--rm','--network',"container:$webContainer",
        '--mount',"type=bind,source=$repositoryRoot,target=/app,readonly",
        '--mount',"type=bind,source=$runDirectory,target=/app/output",
        '--workdir','/app','--entrypoint','node',$BrowserImage,'tests/browser/test_v8_api_browser.cjs')
    Write-Host "PASS: rebuilt v8 runtime images, API and browser. Artifacts: $runDirectory"
} finally {
    & docker compose -f $composeFile -p $project logs --no-color | Set-Content -LiteralPath (Join-Path $runDirectory 'stack.log')
    & docker compose -f $composeFile -p $project down --volumes
    if ($LASTEXITCODE -ne 0) { throw "Could not remove verification stack $project" }
}
