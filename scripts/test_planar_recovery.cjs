const fs=require('node:fs'),assert=require('node:assert/strict'),path=require('node:path');
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const context=require('./recovery-test-context.cjs'),root=path.resolve(__dirname,'..'),directory=root+'/output/planar-recovery';
const read=name=>JSON.parse(fs.readFileSync(directory+'/'+name+'.json'));
const M=context(),S=M.spatialResults,manifest=read('manifest');assert.equal(manifest.length,23);
const plain=x=>JSON.parse(JSON.stringify(x));
function close(actual,expected,label){
  assert.equal(actual.length,expected.length,label);
  const scale=Math.max(1,...expected.map(Math.abs));
  actual.forEach((v,i)=>assert(Number.isFinite(v)&&Math.abs(v-expected[i])<=1e-8+2e-9*scale,`${label}[${i}]: ${v} != ${expected[i]}`));
}
let frames=0;
for(const row of manifest){
  const references=read(row.name+'-reference');assert.equal(references.length,row.elements);
  let baseline;
  for(const mode of ['inline','table','auto']){
    const raw=read(row.name+'-'+mode),before=JSON.stringify(raw),ds=M.validateDataset(raw),a=raw.analysis;
    if(!baseline)baseline=ds;else assert.deepEqual(a,baseline.raw.analysis);
    for(const [i,e] of raw.model.elements.entries()){
      const count=a.type==='modal'?a.eigenvalues.length:a.type==='transient'?a.time.length:1;
      for(let frame=0;frame<count;frame++){
        const u=S.vector(ds,'displacements',frame),acc=S.vector(ds,'accelerations',frame),r=S.recover(ds,e,u,acc);
        if(!r){assert(row.name.endsWith('partial')||row.name.endsWith('acceleration'));continue;}
        const ref=references[i];
        if(a.type==='static'){
          close(r.end,ref.end,'static end');assert.deepEqual(plain(r.coeff),ref.coeff);
          for(const [j,q]of ['N','V','M'].entries())for(const xi of [0,.25,.5,.75,1]){
            const expected=ref.coeff[j].reduceRight((s,v)=>s*xi+v,0);
            close([M.staticResults.elementDiagram(ds,e,q,xi)],[expected],q);
          }
          if(e.type===113)for(const q of ['N','V','M'])assert(M.staticResults.diagramSamples(ds,e,q,11).every(p=>Number.isFinite(p.value)));
          if(e.sectionId){
            const stress=M.sectionResults.get(ds,e).stress;
            assert.deepEqual(plain(stress.minimum),ds.resultsByElementId.get(e.id).sectionStress.minimum);
            assert.equal(M.normalStressResults.get(ds,e).status,'ok');
          }
        }else{
          close(r.u,ref.local.map(row=>row[frame]),'local');
          close(r.end,ref.elastic.map(row=>row[frame]),'elastic');
          if(r.inertia)close(r.inertia,ref.inertia.map(row=>row[frame]),'inertia');
          const baseElement=baseline.raw.model.elements[i];
          assert.deepEqual(plain(r),plain(S.recover(baseline,baseElement,S.vector(baseline,'displacements',frame),S.vector(baseline,'accelerations',frame))));
          if(e.type===112)assert(r.coeff.slice(1).every(row=>row.every(v=>v===0)));
        }
        frames++;
      }
      if(raw.model.recovery.version===2)assert.strictEqual(S.operators(e,8,raw.model),raw.model.recoveryOperators[e.recoveryOperatorId-1]);
    }
    assert.equal(JSON.stringify(raw),before,'Reader mutated result');
  }
}
// Exported operators are authoritative, not a hint to recompute from properties.
const original=read('CaseMPCModal-inline'),base=M.validateDataset(structuredClone(original));
const changed=structuredClone(original);
for(const e of changed.model.elements)for(const key of ['K','M'])e.recoveryOperators[key]=e.recoveryOperators[key].map(row=>row.map(v=>v*2));
const altered=M.validateDataset(changed),e=altered.raw.model.elements[0],b=base.raw.model.elements[0];
const x=S.recover(altered,e,S.vector(altered,'displacements'),S.vector(altered,'accelerations'));
const y=S.recover(base,b,S.vector(base,'displacements'),S.vector(base,'accelerations'));
close(x.end,y.end.map(v=>2*v),'authoritative K');close(x.inertia,y.inertia.map(v=>2*v),'authoritative M');
const invalid=[
  x=>delete x.model.recovery,x=>x.metadata.exporter.contractVersion=3,
  x=>x.model.recovery.stiffnessFormulation='euler-bernoulli-saint-venant',
  x=>x.model.dofLabels.reverse(),x=>delete x.model.numberOfDOFs,
  x=>x.model.elements[0].length*=2,x=>x.model.elements[0].localFromSystem[0][0]=.7,
  x=>x.model.elements[0].globalDOFs.reverse(),x=>x.model.elements[0].globalDOFs[0]=999,
  x=>delete x.model.elements[0].recoveryOperators,
  x=>x.model.elements[0].recoveryOperators.K.pop(),x=>x.model.elements[0].recoveryOperators.M[0][0]=NaN,
  x=>x.model.elements[0].recoveryOperators.K[0][1]=1e30,
  x=>delete x.analysis.eigenvalues,x=>x.analysis.eigenvalues[0]*=2
];
for(const mutate of invalid){const bad=structuredClone(original);mutate(bad);assert.throws(()=>M.validateDataset(bad));}
const pooled=read('CaseMPCModal-table');
for(const mutate of [x=>x.model.elements[0].recoveryOperatorId=999,x=>x.model.recoveryOperators[0].id=2,x=>x.model.elements[0].recoveryOperators=original.model.elements[0].recoveryOperators]){
  const bad=structuredClone(pooled);mutate(bad);assert.throws(()=>M.validateDataset(bad));
}
const stat=read('CaseReleaseStatic-inline');
for(const mutate of [x=>delete x.analysis.elementResults[0].sectionForceCoefficients,x=>x.analysis.elementResults[0].sectionForceCoefficients[1].pop(),x=>x.model.elements[0].recoveryOperators=original.model.elements[0].recoveryOperators]){
  const bad=structuredClone(stat);mutate(bad);assert.throws(()=>M.validateDataset(bad));
}
const release=read('inclined-modal-inline');release.model.elements[0].globalDOFs[2]=release.model.elements[0].globalDOFs[0];
assert.throws(()=>M.validateDataset(release));
console.log(`Planar recovery: 23 cases, ${frames} element frames, exact storage parity, solver references and 22 malformed contracts PASS`);

(async()=>{
  const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try{
    const page=await browser.newPage({viewport:{width:1500,height:1100}}),errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    await page.goto('file://'+root+'/postprocessor/index.html');
    await page.evaluate(()=>{
      const mount=MKEFPost.WebGLRenderer.prototype.mount;
      MKEFPost.WebGLRenderer.prototype.mount=function(...args){window.resultRenderer=this;return mount.apply(this,args);};
    });
    for(const name of ['CaseReleaseStatic','CaseLinearFrame','inclined-static','CaseSectionBeam','CaseMPCModal','inclined-single','inclined-transient','inclined-partial','inclined-acceleration','truss-modal','truss-transient']){
      const source=read(name+'-auto');
      await page.evaluate(()=>window.previousData=window.resultRenderer?.dataset);
      await page.locator('#file-input').setInputFiles(directory+'/'+name+'-auto.json');
      await page.waitForFunction(()=>window.resultRenderer?.dataset&&resultRenderer.dataset!==window.previousData);
      assert.equal(await page.locator('#error-panel').isVisible(),false,name);
      await page.locator('#planar-result-objects select').selectOption('element');
      await page.locator('#planar-result-rows tr').first().click();
      if(source.analysis.type==='static'){
        await page.locator('#static-result').selectOption('M');
        assert(await page.evaluate(()=>{resultRenderer.flush();return resultRenderer.frame.diagrams.positions.length>0;}));
        if(name==='CaseSectionBeam'){
          assert(await page.locator('#section-result-panel').isVisible());
          await page.locator('#section-position-number').fill('.3');await page.locator('#section-position-number').dispatchEvent('input');
          assert.match(await page.locator('#section-result-info').innerText(),/σmin=/);
          // Inspect the actual request builder; the section solver itself is unchanged.
          const body=await page.evaluate(async()=>{
            const originalFetch=window.fetch;let captured;
            window.fetch=async(url,options)=>{
              if(String(url).endsWith('/member')){captured=JSON.parse(options.body);return new Response(JSON.stringify({id:'test'}));}
              return new Response(JSON.stringify({state:'completed',result:{status:'ok',maximum:{value:1,xi:.5}}}));
            };
            try{const button=document.querySelector('#member-run');button.disabled=false;await button.onclick();return captured;}
            finally{window.fetch=originalFetch;}
          });
          assert.equal(source.model.elements[0].sectionAngle,0);
          const coeff=source.analysis.elementResults[0].sectionForceCoefficients;
          assert.deepEqual(body.coefficients,[coeff[0],coeff[1],[0,0,0,0],[0],[0,0,0,0],coeff[2]]);
        }
      }
      if(source.analysis.type==='modal'){await page.locator('#modal-phase').fill('-0.5');await page.locator('#modal-phase').dispatchEvent('input');}
      if(source.analysis.type==='transient'){
        await page.locator('#time-index').fill('1');await page.locator('#time-index').dispatchEvent('input');
        if(name==='inclined-partial'||name==='inclined-acceleration'){
          assert(await page.locator('#show-deformed').isDisabled());
          assert.match(await page.locator('#notice-panel').innerText(),/charts only/);
          assert(await page.locator('#history-chart .chart-series').count());
        }
      }
      if(['CaseSectionBeam','inclined-transient'].includes(name)){
        const waiting=page.waitForEvent('download');await page.locator('#download-json').click();
        const file=directory+'/'+name+'-roundtrip.json';await(await waiting).saveAs(file);assert.deepEqual(JSON.parse(fs.readFileSync(file)),source);
        const exporting=page.waitForEvent('download');await page.locator('#export-svg').click();
        const svg=directory+'/'+name+'.svg';await(await exporting).saveAs(svg);assert(fs.readFileSync(svg,'utf8').includes('data-layer="deformed"'));
      }
    }
    assert.deepEqual(errors,[]);
    console.log('Planar browser: static diagrams, sections, modes, transient/partial data, JSON round-trip and SVG with solver recovery data PASS');
  }finally{await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
