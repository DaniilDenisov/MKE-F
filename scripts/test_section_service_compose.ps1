# Creates a private test stack on a randomly allocated loopback port.
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$project = 'mkef-section-test-' + [guid]::NewGuid().ToString('N').Substring(0,10)
$previousPort = $env:MKEF_PORT
Push-Location $root
try {
    $env:MKEF_PORT = '0'
    docker compose -p $project up -d --build --wait
    if ($LASTEXITCODE) { throw 'Test stack startup failed' }
    docker build -f docker/browser-tests.Dockerfile -t mkef-browser-tests .
    if ($LASTEXITCODE) { throw 'Browser image build failed' }
    $webId = docker compose -p $project ps -q web
    function Probe([string]$state) {
        Get-Content -Raw -Encoding utf8 tests/section_service/proxy_probe.py |
            docker compose -p $project exec -T solver python3 - $state
        if ($LASTEXITCODE) { throw "Proxy probe failed ($state)" }
        docker run --rm --network "container:$webId" --mount "type=bind,source=$root,target=/app,readonly" -w /app mkef-browser-tests node scripts/test_section_diagnostic_browser.cjs $state
        if ($LASTEXITCODE) { throw "Diagnostic browser failed ($state)" }
    }
    Probe 'up'
    docker compose -p $project stop section-solver
    if ($LASTEXITCODE) { throw 'Could not stop section solver' }
    Probe 'down'
    docker compose -p $project up -d --wait section-solver
    if ($LASTEXITCODE) { throw 'Could not restart section solver' }
    Probe 'recovered'
    docker compose -p $project exec -T web sh -ec 'cd /usr/share/nginx/html && test ! -e shared/section-fem && test ! -e shared/mises.js && test ! -e postprocessor/js/mises-worker.js && test ! -e scripts/section_fem_fvm_reference.cjs'
    if ($LASTEXITCODE) { throw 'Web image contains retired solver assets' }
} finally {
    docker compose -p $project down --rmi local
    $env:MKEF_PORT = $previousPort
    Pop-Location
}
