// Legacy chapter URLs retain their section anchors in the Russian edition.
(() => {
  const link = document.getElementById('redirect-target');
  if (!link) return;
  const target = new URL(link.href);
  target.search = window.location.search;
  target.hash = window.location.hash;
  link.href = target.href;
  window.location.replace(target.href);
})();
