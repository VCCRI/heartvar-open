    (function(){
      'use strict';

      var MS_ICON = '<svg viewBox="0 0 20 20" width="17" height="17" aria-hidden="true">'
        + '<rect x="1" y="1" width="8.2" height="8.2" fill="#F25022"/>'
        + '<rect x="10.8" y="1" width="8.2" height="8.2" fill="#7FBA00"/>'
        + '<rect x="1" y="10.8" width="8.2" height="8.2" fill="#00A4EF"/>'
        + '<rect x="10.8" y="10.8" width="8.2" height="8.2" fill="#FFB900"/></svg>';

      var GOOGLE_ICON = '<svg viewBox="0 0 48 48" width="17" height="17" aria-hidden="true">'
        + '<path fill="#EA4335" d="M24 9.5c3.5 0 6.6 1.2 9 3.5l6.7-6.7C35.5 2.5 30.1 0 24 0 14.6 0 6.4 5.4 2.5 13.2l7.8 6.1C12.2 13.2 17.6 9.5 24 9.5z"/>'
        + '<path fill="#4285F4" d="M46.1 24.5c0-1.6-.1-2.8-.4-4.1H24v9.3h12.6c-.3 2.1-1.6 5.2-4.7 7.3l7.6 5.9c4.5-4.2 6.6-10.3 6.6-18.4z"/>'
        + '<path fill="#FBBC05" d="M10.3 28.7a14.6 14.6 0 0 1 0-9.4l-7.8-6.1a24 24 0 0 0 0 21.6l7.8-6.1z"/>'
        + '<path fill="#34A853" d="M24 48c6.1 0 11.3-2 15.5-5.5l-7.6-5.9c-2 1.4-4.7 2.4-7.9 2.4-6.4 0-11.8-3.7-13.7-9.3l-7.8 6.1C6.4 42.6 14.6 48 24 48z"/></svg>';

      var PROVIDER_ICONS = { microsoft: MS_ICON, google: GOOGLE_ICON };

      var PROVIDER_PATHS = { microsoft: '/login', google: '/login/google' };

      var SNAPSHOT_KEY = 'hvAuthFormSnapshot';
      var INTENT_KEY = 'hvAuthAiIntent';

      var FIELD_IDS = [
        'hvl-gene', 'hvl-variant', 'hvl-hpo', 'hvl-zygosity', 'hvl-inheritance',
        'hvl-sex', 'hvl-trio', 'hvl-denovo', 'hvl-famhx',
        'hvl-seg-affected-carriers', 'hvl-in-trans',
      ];

      var state = { ready: false, authenticated: false, label: '', gate: false, providers: [], is_admin: false };
      window.hvAuthState = state;

      function esc(s){
        return String(s == null ? '' : s).replace(/[&<>"]/g, function(c){
          return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c];
        });
      }

      function refresh(){
        return fetch('/api/auth/status', { credentials: 'same-origin' })
          .then(function(r){ return r.ok ? r.json() : null; })
          .catch(function(){ return null; })
          .then(function(s){
            s = s || {};
            state.gate = !!s.ai_requires_signin;
            state.providers = s.providers || [];
            state.authenticated = !!s.authenticated;
            state.label = s.label || '';
            state.is_admin = !!s.is_admin;
            state.ready = true;
            render();
            return state;
          });
      }
      window.hvAuthRefresh = refresh;

      function blocked(){
        return state.ready && state.gate && !state.authenticated;
      }
      window.hvAuthBlocked = blocked;

      function chipHTML(){

        if (!state.gate && !state.authenticated) return '';
        if (state.authenticated) {
          return '<button type="button" class="hv-auth-chip" '
            + 'title="Signed in as ' + esc(state.label) + '" '
            + 'onclick="hvAuthSignOut()">Sign out</button>';
        }
        return '<button type="button" class="hv-auth-chip hv-auth-chip__signin" '
          + 'onclick="hvAuthOpen()">Sign in</button>';
      }

      function render(){
        ['hv-auth-slot-landing', 'hv-auth-slot-app'].forEach(function(id){
          var host = document.getElementById(id);
          if (host) host.innerHTML = chipHTML();
        });

        var adminTab = document.getElementById('nt-admin');
        if (adminTab) adminTab.style.display = state.is_admin ? '' : 'none';

        var adminSlotLanding = document.getElementById('hv-admin-slot-landing');
        if (adminSlotLanding) {
          adminSlotLanding.innerHTML = state.is_admin
            ? '<a onclick="document.body.classList.remove(\'hv-landing-on\'); showPage(\'admin\');">Admin</a>'
            : '';
        }
      }

      window.hvlAiPillChange = function(cb){
        var pill = cb.closest('.hvl-ai-pill');
        if (cb.checked && blocked()) {
          cb.checked = false;
          if (pill) pill.classList.remove('is-on');
          openModal({ intent: true });
          return;
        }
        if (pill) pill.classList.toggle('is-on', cb.checked);
      };

      function ensureModal(){
        var el = document.getElementById('hv-auth-modal');
        if (el) return el;
        el = document.createElement('div');
        el.id = 'hv-auth-modal';
        el.className = 'hv-auth-modal';
        el.setAttribute('role', 'dialog');
        el.setAttribute('aria-modal', 'true');
        el.setAttribute('aria-labelledby', 'hv-auth-modal-title');

        el.innerHTML =
          '<div class="hv-auth-modal__backdrop" onclick="hvAuthClose()"></div>'
          + '<div class="hv-auth-modal__card">'
          +   '<button type="button" class="hv-auth-modal__x" aria-label="Close" '
          +     'onclick="hvAuthClose()">&times;</button>'
          +   '<h2 id="hv-auth-modal-title" class="hv-auth-modal__title">'
          +     'Sign in to use AI interpretation</h2>'
          +   '<p class="hv-auth-modal__why">Sign-in lets us meter AI usage '
          +     'fairly. Everything else on HeartVar works without it.</p>'
          +   '<div class="hv-auth-modal__providers" id="hv-auth-provider-list"></div>'
          +   '<button type="button" class="hv-auth-modal__skip" onclick="hvAuthClose()">'
          +     'Continue without AI</button>'
          + '</div>';
        document.body.appendChild(el);
        return el;
      }

      function openModal(opts){
        var el = ensureModal();
        var list = el.querySelector('#hv-auth-provider-list');
        if (!state.providers.length) {
          list.innerHTML = '<p class="hv-auth-modal__empty">Sign-in is not '
            + 'configured on this server yet.</p>';
        } else {
          list.innerHTML = state.providers.map(function(p){
            return '<button type="button" class="hv-auth-prov" '
              + 'onclick="hvAuthGo(\'' + esc(p.key) + '\')">'
              + (PROVIDER_ICONS[p.key] || '')
              + '<span>Continue with ' + esc(p.label) + '</span></button>';
          }).join('');
        }

        try {
          if (opts && opts.intent) sessionStorage.setItem(INTENT_KEY, '1');
        } catch (e) {}
        el.classList.add('is-open');
        document.body.classList.add('hv-auth-modal-open');

        var card = el.querySelector('.hv-auth-modal__card');
        if (card) { card.setAttribute('tabindex', '-1'); card.focus(); }
      }
      window.hvAuthOpen = function(){ openModal({ intent: false }); };

      window.hvAuthClose = function(){
        var el = document.getElementById('hv-auth-modal');
        if (el) el.classList.remove('is-open');
        document.body.classList.remove('hv-auth-modal-open');
        try { sessionStorage.removeItem(INTENT_KEY); } catch (e) {}
      };

      document.addEventListener('keydown', function(e){
        if (e.key !== 'Escape') return;
        var el = document.getElementById('hv-auth-modal');
        if (el && el.classList.contains('is-open')) window.hvAuthClose();
      });

      function snapshotForm(){
        var data = {};
        FIELD_IDS.forEach(function(id){
          var el = document.getElementById(id);
          if (el) data[id] = el.value;
        });
        var buildRow = document.getElementById('build-toggle-row');
        if (buildRow) data.__build = buildRow.dataset.build || 'GRCh38';
        try { sessionStorage.setItem(SNAPSHOT_KEY, JSON.stringify(data)); } catch (e) {}
      }

      function restoreForm(){
        var raw = null;
        try { raw = sessionStorage.getItem(SNAPSHOT_KEY); } catch (e) {}
        if (!raw) return false;
        try { sessionStorage.removeItem(SNAPSHOT_KEY); } catch (e) {}
        var data;
        try { data = JSON.parse(raw); } catch (e) { return false; }
        if (!data || typeof data !== 'object') return false;
        var restoredAny = false;
        FIELD_IDS.forEach(function(id){
          var el = document.getElementById(id);
          if (el && typeof data[id] === 'string' && data[id] !== '') {
            el.value = data[id];
            restoredAny = true;
          }
        });

        if (typeof window.hvlSyncBuildToggle === 'function') window.hvlSyncBuildToggle();

        if (data.__build === 'GRCh37' && typeof window.hvlSetBuild === 'function') {
          var seg = document.querySelector('.hvl-build-opt[data-build="GRCh37"]');
          if (seg) window.hvlSetBuild(seg);
        }
        if (typeof window.hvlSyncDropdowns === 'function') window.hvlSyncDropdowns();
        return restoredAny;
      }

      window.hvAuthGo = function(providerKey){
        snapshotForm();
        var path = PROVIDER_PATHS[providerKey] || '/login';
        var next = window.location.pathname + window.location.search;
        window.location.href = path + '?next=' + encodeURIComponent(next);
      };

      window.hvAuthSignOut = function(){
        var cb = document.getElementById('hvl-ai-enable');
        if (cb) cb.checked = false;
        var pill = document.querySelector('.hvl-ai-pill');
        if (pill) pill.classList.remove('is-on');
        window.location.href = '/logout';
      };

      function notify(msg){
        if (typeof window.showToast === 'function') { window.showToast(msg); return; }
        var host = document.getElementById('toast');
        if (!host) return;
        host.textContent = msg;
        host.classList.add('show');
        setTimeout(function(){ host.classList.remove('show'); }, 4200);
      }
      window.hvAuthNotify = notify;

      window.hvAuthHandle401 = function(){
        state.authenticated = false;
        state.ready = true;
        render();
        openModal({ intent: true });
      };

      function boot(){
        var hadSnapshot = restoreForm();
        refresh().then(function(){
          var wanted = false;
          try { wanted = sessionStorage.getItem(INTENT_KEY) === '1'; } catch (e) {}
          if (!state.authenticated) return;
          if (wanted) {

            try { sessionStorage.removeItem(INTENT_KEY); } catch (e) {}
            var cb = document.getElementById('hvl-ai-enable');
            if (cb && !cb.checked) {
              cb.checked = true;
              var pill = document.querySelector('.hvl-ai-pill');
              if (pill) pill.classList.add('is-on');
            }
            notify('Signed in as ' + state.label + ' — AI interpretation enabled.');
          } else if (hadSnapshot) {
            notify('Signed in as ' + state.label + '.');
          }
        });
      }

      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', boot);
      } else {
        boot();
      }
    })();

    (function(){
      'use strict';

      function fmtTs(ts){
        if (!ts) return '\u2014';
        try {
          var d = new Date(ts);
          if (isNaN(d)) return ts;
          return d.toLocaleString(undefined, {
            year:'numeric', month:'short', day:'numeric',
            hour:'2-digit', minute:'2-digit', second:'2-digit',
            hour12: false
          });
        } catch(e){ return ts; }
      }

      function fmtSize(bytes){
        if (bytes == null) return '';
        if (bytes < 1024) return bytes + ' B';
        if (bytes < 1048576) return (bytes / 1024).toFixed(1) + ' KB';
        return (bytes / 1048576).toFixed(1) + ' MB';
      }

      function renderEntry(category, entry){
        var ts   = fmtTs(entry.ts);
        var size = fmtSize(entry.size);
        var name = entry.name || '';
        function q(v){ return JSON.stringify(v).replace(/"/g, '&quot;'); }
        return '<li class="admin-log-item">'
          + '<button class="admin-log-item__btn" onclick="hvAdminOpenLog(' + q(category) + ',' + q(name) + ',' + q(ts) + ')">'
          + '<span class="admin-log-item__dot"></span>'
          + '<span class="admin-log-item__body">'
          + '<span class="admin-log-item__ts">' + ts + '</span>'
          + '<span class="admin-log-item__name">' + name + '</span>'
          + '</span>'
          + (size ? '<span class="admin-log-item__size">' + size + '</span>' : '')
          + '<span class="admin-log-item__chevron">&#8250;</span>'
          + '</button></li>';
      }

      function setList(id, category, entries, emptyMsg){
        var el = document.getElementById(id);
        if (!el) return;
        el.classList.remove('admin-log-list--loading');
        if (!entries || !entries.length){
          el.innerHTML = '<p class="admin-empty">' + emptyMsg + '</p>';
          return;
        }
        el.innerHTML = '<ul class="admin-log-ul">'
          + entries.map(function(e){ return renderEntry(category, e); }).join('')
          + '</ul>';
      }

      window.hvAdminLoadLogs = function(){
        var dbEl = document.getElementById('admin-db-logs');
        var dpEl = document.getElementById('admin-deploy-logs');
        if (dbEl) { dbEl.classList.add('admin-log-list--loading'); dbEl.innerHTML = '<span class="admin-spinner"></span> Loading\u2026'; }
        if (dpEl) { dpEl.classList.add('admin-log-list--loading'); dpEl.innerHTML = '<span class="admin-spinner"></span> Loading\u2026'; }
        var iv = document.getElementById('admin-index-view');
        var cv = document.getElementById('admin-content-view');
        if (iv) iv.style.display = '';
        if (cv) cv.style.display = 'none';
        fetch('/api/admin/log-index', { credentials: 'same-origin' })
          .then(function(r){
            if (r.status === 403) throw new Error('forbidden');
            if (!r.ok) throw new Error('http ' + r.status);
            return r.json();
          })
          .then(function(data){
            setList('admin-db-logs', 'db_builder', data.db_builder, 'No database build records found.');
            setList('admin-deploy-logs', 'deployments', data.deployments, 'No deployment records found.');
          })
          .catch(function(err){
            var msg = err.message === 'forbidden'
              ? 'Access denied \u2014 admin credentials required.'
              : 'Could not load logs: ' + err.message;
            var errHtml = '<p class="admin-empty admin-empty--err">' + msg + '</p>';
            if (dbEl) { dbEl.innerHTML = errHtml; dbEl.classList.remove('admin-log-list--loading'); }
            if (dpEl) { dpEl.innerHTML = errHtml; dpEl.classList.remove('admin-log-list--loading'); }
          });
      };

      window.hvAdminOpenLog = function(category, name, title){
        var iv      = document.getElementById('admin-index-view');
        var cv      = document.getElementById('admin-content-view');
        var titleEl = document.getElementById('admin-content-title');
        var bodyEl  = document.getElementById('admin-content-body');
        if (iv) iv.style.display = 'none';
        if (cv) cv.style.display = '';
        if (titleEl) titleEl.textContent = title || name;
        if (bodyEl) { bodyEl.textContent = 'Loading\u2026'; bodyEl.classList.add('admin-log-content--loading'); }
        fetch('/api/admin/log-content?category=' + encodeURIComponent(category) + '&name=' + encodeURIComponent(name),
          { credentials: 'same-origin' })
          .then(function(r){
            if (!r.ok) throw new Error('http ' + r.status);
            return r.json();
          })
          .then(function(data){
            if (bodyEl) { bodyEl.textContent = data.content || '(empty)'; bodyEl.classList.remove('admin-log-content--loading'); }
          })
          .catch(function(err){
            if (bodyEl) { bodyEl.textContent = 'Could not load log: ' + err.message; bodyEl.classList.remove('admin-log-content--loading'); }
          });
      };

      window.hvAdminBackToIndex = function(){
        var iv = document.getElementById('admin-index-view');
        var cv = document.getElementById('admin-content-view');
        if (iv) iv.style.display = '';
        if (cv) cv.style.display = 'none';
      };

      window.hvAdminLoadDbStatus = function(){
        var el = document.getElementById('admin-db-status');
        if (el) el.innerHTML = '<span class="admin-spinner"></span> Loading\u2026';
        fetch('/api/admin/db-status', { credentials: 'same-origin' })
          .then(function(r){
            if (r.status === 403) throw new Error('forbidden');
            if (!r.ok) throw new Error('http ' + r.status);
            return r.json();
          })
          .then(function(data){
            if (!el) return;
            var cadenceNote = {
              monthly:   'Rebuilt monthly by build_all.sh. Run the data builder to create it.',
              versioned: 'Rebuilt on upstream version bumps. Run the builder when a new release is available.',
              static:    'One-time build (static upstream data). Run the builder once to create it.',
              image:     'Ships inside the Docker image. Rebuild the app image to regenerate it.',
            };
            var rows = (data.sources || []).map(function(s){
              var updCell = s.last_updated
                ? '<span class="admin-ds-ts">' + fmtTs(s.last_updated) + '</span>'
                : '<span class="admin-ds-missing">not built</span>';

              var statusCell;
              if (s.present && !s.failed_in_last_run) {
                statusCell = '<span class="admin-ds-ok" title="Artifact found">&#10003;</span>';
              } else {
                var failTitle = s.failed_in_last_run
                  ? 'Failed in last build run \u2014 click for details'
                  : 'Artifact missing \u2014 click for details';
                statusCell = '<button class="admin-ds-fail" title="' + failTitle + '" '
                  + 'onclick="hvAdminToggleDbDetail(\'' + s.name + '\')">&#10007;</button>';
              }

              var mainRow = '<tr>'
                + '<td class="admin-ds-name">' + s.name.replace(/_/g, '_<wbr>') + '</td>'
                + '<td>' + updCell + '</td>'
                + '<td class="admin-ds-status-cell">' + statusCell + '</td>'
                + '</tr>';

              var detailRow = '';
              if (!s.present || s.failed_in_last_run) {
                var note = cadenceNote[s.cadence] || 'Run the data builder to create it.';
                var errorBlock = '';
                if (s.failed_in_last_run && s.error_snippet) {
                  errorBlock = '<pre class="admin-ds-error-pre">' + s.error_snippet.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;') + '</pre>';
                } else if (s.failed_in_last_run) {
                  errorBlock = '<p class="admin-ds-error-note">Failed in the last build run (no snippet available \u2014 check the Database Builder log).</p>';
                }
                var heading = s.failed_in_last_run
                  ? '<strong>Failed in last build run:</strong> <code>' + (s.artifact_path || s.name) + '</code>'
                  : '<strong>Artifact not found:</strong> <code>' + (s.artifact_path || s.name) + '</code>';
                detailRow = '<tr id="admin-ds-detail-' + s.name + '" class="admin-ds-detail-row" style="display:none">'
                  + '<td colspan="3" class="admin-ds-detail-cell">'
                  + heading + '<br>'
                  + '<strong>Cadence:</strong> ' + (s.cadence || 'monthly') + '&nbsp;&mdash;&nbsp;' + note
                  + errorBlock
                  + '</td></tr>';
              }
              return mainRow + detailRow;
            }).join('');

            el.innerHTML = '<table class="admin-ds-table"><thead><tr>'
              + '<th>Source</th><th>Last Updated</th><th>Status</th>'
              + '</tr></thead><tbody>'
              + (rows || '<tr><td colspan="3" class="admin-ds-empty">No build data found.</td></tr>')
              + '</tbody></table>';
          })
          .catch(function(err){
            var msg = err.message === 'forbidden'
              ? 'Access denied \u2014 admin credentials required.'
              : 'Could not load data status: ' + err.message;
            if (el) el.innerHTML = '<p class="admin-empty admin-empty--err">' + msg + '</p>';
          });
      };

      window.hvAdminToggleDbDetail = function(name){
        var row = document.getElementById('admin-ds-detail-' + name);
        if (row) row.style.display = row.style.display === 'none' ? '' : 'none';
      };
    })();
