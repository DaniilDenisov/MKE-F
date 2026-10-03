(function () {
  'use strict';
  var tests = [];
  function test(name, operation) { tests.push({ name: name, operation: operation }); }
  function assert(condition, message) { if (!condition) throw new Error(message || 'Assertion failed.'); }
  function response(status, payload) {
    return { ok: status >= 200 && status < 300, status: status, headers: { get: function () { return 'application/json'; } }, json: function () { return Promise.resolve(payload); } };
  }

  test('validates canonical UUID job IDs', function () {
    assert(MKEFApi.isValidJobId('d9428888-122b-4d7d-9ca4-eefd2f155174'));
    assert(!MKEFApi.isValidJobId('../result.json'));
  });

  test('submits only name and serialized case text', async function () {
    var request = null, client = MKEFApi.createClient({ fetch: function (url, options) { request = { url: url, options: options }; return Promise.resolve(response(202, { id: 'd9428888-122b-4d7d-9ca4-eefd2f155174', status: 'queued' })); } });
    var job = await client.createJob('Frame', 'analysis\nstatic\n');
    assert(job.status === 'queued');
    assert(request.url === '/api/v1/jobs');
    assert(JSON.stringify(JSON.parse(request.options.body)) === JSON.stringify({ name: 'Frame', caseText: 'analysis\nstatic\n' }));
  });

  test('maps structured API failures to readable errors', async function () {
    var client = MKEFApi.createClient({ fetch: function () { return Promise.resolve(response(409, { detail: { status: 'failed', error: { message: 'bad model' } } })); } });
    try { await client.getResult('d9428888-122b-4d7d-9ca4-eefd2f155174'); }
    catch (error) { assert(error.status === 409 && error.message.indexOf('failed') >= 0); return; }
    throw new Error('Expected an API error.');
  });

  test('rejects invalid IDs before making a request', async function () {
    var called = false, client = MKEFApi.createClient({ fetch: function () { called = true; } });
    try { await client.cancelJob('../../bad'); }
    catch (error) { assert(error.message.indexOf('Invalid') >= 0 && !called); return; }
    throw new Error('Expected invalid ID rejection.');
  });

  document.addEventListener('DOMContentLoaded', async function () {
    var output = document.getElementById('test-output'), lines = [], passed = 0;
    for (var index = 0; index < tests.length; index += 1) {
      var item = tests[index];
      try { await item.operation(); passed += 1; lines.push('PASS  ' + item.name); }
      catch (error) { lines.push('FAIL  ' + item.name + '\n      ' + error.message); }
    }
    output.textContent = lines.join('\n');
    document.body.setAttribute('data-test-status', passed === tests.length ? 'passed' : 'failed');
  });
}());
