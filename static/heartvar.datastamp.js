    (function(){
      'use strict';

      var FALLBACK = 'Data is updated monthly.';

      function formatAU(iso){
        var d = new Date(iso);
        if (isNaN(d.getTime())) return '';
        var pad = function(n){ return (n < 10 ? '0' : '') + n; };
        return pad(d.getDate()) + '/' + pad(d.getMonth() + 1) + '/' + d.getFullYear();
      }

      function render(stamp){
        var el = document.querySelector('.src-update-note');
        if (!el) return;
        var iso = stamp && stamp.last_refresh_utc;
        var date = iso ? formatAU(iso) : '';
        el.textContent = date ? (FALLBACK + ' Last update: ' + date + '.') : FALLBACK;
      }

      function load(){
        fetch('/api/data-status', { credentials: 'same-origin' })
          .then(function(r){ return r.ok ? r.json() : null; })

          .catch(function(){ return null; })
          .then(render);
      }

      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', load);
      } else {
        load();
      }
    })();
