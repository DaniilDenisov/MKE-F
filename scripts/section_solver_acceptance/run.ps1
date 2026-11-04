$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$dir = 'output/section-solver/' + [guid]::NewGuid().ToString('N')
Push-Location $root
try {
    docker build --target numeric -f docker/section-solver.Dockerfile -t mkef-section-test .
    if ($LASTEXITCODE) { throw 'Numeric image build failed' }
    $common = @('--rm','--mount',"type=bind,source=$root,target=/app",'-w','/app','-e',"SECTION_REPORT_DIR=$dir")
    function Run-InDocker([string[]]$arguments) {
        docker run @common @arguments
        if ($LASTEXITCODE) { throw "Failed: $arguments" }
    }
    Run-InDocker @('mkef-section-test','python','-m','unittest','discover','-s','tests/section_solver','-v')
    Run-InDocker @('mkef-section-test','python','-m','scripts.section_solver_acceptance.matrix')
    Run-InDocker @('mkef-section-test','python','-m','scripts.section_solver_acceptance.verify')
    Write-Output "$root/$dir/report.json"
} finally {
    Pop-Location
}
