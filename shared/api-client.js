(function (global) {
  'use strict';

  var UUID_PATTERN = /^[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/i;

  function ApiError(message, status, detail) {
    this.name = 'ApiError';
    this.message = message;
    this.status = status || 0;
    this.detail = detail;
    if (Error.captureStackTrace) Error.captureStackTrace(this, ApiError);
  }
  ApiError.prototype = Object.create(Error.prototype);
  ApiError.prototype.constructor = ApiError;

  function detailMessage(payload, fallback) {
    if (!payload || payload.detail === undefined) return fallback;
    if (typeof payload.detail === 'string') return payload.detail;
    if (payload.detail && typeof payload.detail.message === 'string') return payload.detail.message;
    try { return JSON.stringify(payload.detail); } catch (_) { return fallback; }
  }

  function createClient(options) {
    options = options || {};
    var base = (options.baseUrl || '/api/v1').replace(/\/$/, '');
    var fetcher = options.fetch || (global.fetch && global.fetch.bind(global));

    async function request(path, requestOptions) {
      if (!fetcher) throw new ApiError('This browser does not provide the Fetch API.', 0, null);
      var response;
      try { response = await fetcher(base + path, requestOptions || {}); }
      catch (error) { throw new ApiError('Could not reach the local solver: ' + error.message, 0, null); }
      var contentType = response.headers && response.headers.get ? response.headers.get('content-type') || '' : '';
      var payload = null;
      if (contentType.indexOf('application/json') >= 0) {
        try { payload = await response.json(); } catch (_) { payload = null; }
      } else {
        try { payload = await response.text(); } catch (_) { payload = null; }
      }
      if (!response.ok) throw new ApiError(detailMessage(payload, 'Solver request failed with status ' + response.status + '.'), response.status, payload);
      return payload;
    }

    return {
      health: function () { return request('/health'); },
      createJob: function (name, caseText) {
        return request('/jobs', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ name: name, caseText: caseText }) });
      },
      getJob: function (id) {
        if (!UUID_PATTERN.test(id)) return Promise.reject(new ApiError('Invalid solver job ID.', 0, null));
        return request('/jobs/' + encodeURIComponent(id));
      },
      cancelJob: function (id) {
        if (!UUID_PATTERN.test(id)) return Promise.reject(new ApiError('Invalid solver job ID.', 0, null));
        return request('/jobs/' + encodeURIComponent(id) + '/cancel', { method: 'POST' });
      },
      getResult: function (id) {
        if (!UUID_PATTERN.test(id)) return Promise.reject(new ApiError('Invalid solver job ID.', 0, null));
        return request('/jobs/' + encodeURIComponent(id) + '/result');
      }
    };
  }

  function isHosted() { return global.location && (global.location.protocol === 'http:' || global.location.protocol === 'https:'); }
  function isValidJobId(value) { return typeof value === 'string' && UUID_PATTERN.test(value); }
  function safeFilename(value, fallback) {
    var cleaned = String(value || '').replace(/[^A-Za-z0-9._ -]/g, '_').replace(/^[ .]+|[ .]+$/g, '');
    return cleaned || fallback || 'MKE-F-result';
  }
  function downloadJson(data, name) {
    var blob = new Blob([JSON.stringify(data)], { type: 'application/json;charset=utf-8' });
    var url = URL.createObjectURL(blob), link = document.createElement('a');
    link.href = url; link.download = safeFilename(name, 'MKE-F-result') + '.json'; link.click();
    setTimeout(function () { URL.revokeObjectURL(url); }, 0);
  }

  var client = createClient();
  global.MKEFApi = {
    ApiError: ApiError,
    createClient: createClient,
    isHosted: isHosted,
    isValidJobId: isValidJobId,
    safeFilename: safeFilename,
    downloadJson: downloadJson,
    health: client.health,
    createJob: client.createJob,
    getJob: client.getJob,
    cancelJob: client.cancelJob,
    getResult: client.getResult
  };
}(window));
