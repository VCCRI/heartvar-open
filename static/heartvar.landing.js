    (function(){
      'use strict';

      window.hvlTogglePanel = function(btn, id){
        var panel = document.getElementById(id);
        if (!panel) return;
        var open = panel.classList.toggle('open');
        btn.setAttribute('aria-expanded', open ? 'true' : 'false');
        if (open){
          panel.style.maxHeight = panel.scrollHeight + 'px';
          panel.addEventListener('transitionend', function once(){
            if (panel.classList.contains('open')) {
              panel.style.maxHeight = 'none';

              panel.style.overflow = 'visible';
            }
            panel.removeEventListener('transitionend', once);
          });
        } else {
          panel.style.overflow = 'hidden';
          panel.style.maxHeight = panel.scrollHeight + 'px';
          requestAnimationFrame(function(){
            requestAnimationFrame(function(){ panel.style.maxHeight = '0px'; });
          });
        }
      };

      window.hvlInitDropdowns = function(){
        var selects = document.querySelectorAll('#hv-landing select.hvl-ctl');
        Array.prototype.forEach.call(selects, function(sel){
          if (sel.closest('.hvl-dd')) return;
          var wrap = document.createElement('div');
          wrap.className = 'hvl-dd';
          sel.parentNode.insertBefore(wrap, sel);
          wrap.appendChild(sel);
          var btn = document.createElement('button');
          btn.type = 'button'; btn.className = 'hvl-dd__btn';
          btn.setAttribute('aria-haspopup', 'listbox'); btn.setAttribute('aria-expanded', 'false');
          var menu = document.createElement('div');
          menu.className = 'hvl-dd__menu'; menu.setAttribute('role', 'listbox');
          Array.prototype.forEach.call(sel.options, function(opt){
            var o = document.createElement('div');
            o.className = 'hvl-dd__opt'; o.setAttribute('role', 'option');
            o.setAttribute('data-value', opt.value); o.textContent = opt.textContent;
            o.addEventListener('click', function(){
              sel.value = opt.value;
              sel.dispatchEvent(new Event('change', { bubbles: true }));
              wrap.classList.remove('open'); btn.setAttribute('aria-expanded', 'false');
            });
            menu.appendChild(o);
          });
          wrap.appendChild(btn); wrap.appendChild(menu);
          btn.addEventListener('click', function(e){
            e.stopPropagation();
            var isOpen = wrap.classList.contains('open');
            document.querySelectorAll('#hv-landing .hvl-dd.open').forEach(function(w){
              w.classList.remove('open'); var b = w.querySelector('.hvl-dd__btn'); if (b) b.setAttribute('aria-expanded', 'false');
            });
            if (!isOpen){ wrap.classList.add('open'); btn.setAttribute('aria-expanded', 'true'); }
          });
          sel._hvlSync = function(){
            var opt = sel.options[sel.selectedIndex];
            btn.textContent = opt ? opt.textContent : '';
            Array.prototype.forEach.call(menu.children, function(o){
              o.setAttribute('aria-selected', o.getAttribute('data-value') === sel.value ? 'true' : 'false');
            });
          };
          sel.addEventListener('change', sel._hvlSync);
          sel._hvlSync();
        });
        if (!window.__hvlDDCloser){
          window.__hvlDDCloser = true;
          document.addEventListener('click', function(){
            document.querySelectorAll('#hv-landing .hvl-dd.open').forEach(function(w){
              w.classList.remove('open'); var b = w.querySelector('.hvl-dd__btn'); if (b) b.setAttribute('aria-expanded', 'false');
            });
          });
          document.addEventListener('keydown', function(e){
            if (e.key === 'Escape') document.querySelectorAll('#hv-landing .hvl-dd.open').forEach(function(w){ w.classList.remove('open'); });
          });
        }
      };
      window.hvlSyncDropdowns = function(){
        document.querySelectorAll('#hv-landing select.hvl-ctl').forEach(function(sel){ if (sel._hvlSync) sel._hvlSync(); });
      };

      if (document.readyState !== 'loading') hvlInitDropdowns();
      else document.addEventListener('DOMContentLoaded', hvlInitDropdowns);

      window.enterAppFromLanding = function(){
        try {

          document.body.classList.remove('hv-landing-on');
          window.__hvRan = true;
          if (typeof showPage === 'function') showPage('curate');

          if (typeof runCuration === 'function') runCuration();
        } catch (e) {

          console.error('[HeartVar] enterAppFromLanding failed', e);
          document.body.classList.remove('hv-landing-on');
        }
      };

      window.returnToLanding = function(opts){
        if (!(opts && opts.silent) && typeof hvSyncUrlToHome === 'function') hvSyncUrlToHome();
        document.body.classList.add('hv-landing-on');
        document.body.classList.remove('hv-results-full', 'hv-page-about', 'hv-page-contact');

        var ab = document.getElementById('page-about');    if (ab) ab.style.display = 'none';
        var ct = document.getElementById('page-contact');  if (ct) ct.style.display = 'none';

        var l = document.getElementById('hv-landing');
        if (l) l.scrollTop = 0;

        if (typeof fadeInPage === 'function') fadeInPage(l ? l.querySelector('.hvl-main') : null);
        window.scrollTo(0, 0);
      };

      window.HVL_EXAMPLES = [
        { gene:'MYH7', variant:'c.1988G>A', hpo:'Hypertrophic cardiomyopathy',
          inheritance:'AD', zygosity:'het', sex:'male', trio:'trio', denovo:'inherited_affected',
          famhx:'Father and paternal uncle with HCM; paternal grandfather sudden cardiac death age 48.' },
        { gene:'JAG1', variant:'c.703C>T', hpo:'Tetralogy of Fallot, peripheral pulmonary stenosis',
          inheritance:'AD', zygosity:'het', sex:'male', trio:'trio', denovo:'inherited_affected',
          famhx:'Father has Alagille syndrome with chronic cholestasis; paternal grandfather had pulmonary stenosis repaired in childhood.' },
        { gene:'CHD7', variant:'c.5058del', hpo:'Atrioventricular septal defect, coloboma',
          inheritance:'AD', zygosity:'het', sex:'female', trio:'trio', denovo:'confirmed',
          famhx:'Parents unaffected; no family history of CHD, hearing loss, or coloboma.' },
        { gene:'RAF1', variant:'c.770C>T', hpo:'Hypertrophic cardiomyopathy, pulmonary stenosis',
          inheritance:'AD', zygosity:'het', sex:'female', trio:'trio', denovo:'confirmed',
          famhx:'Parents unaffected; no family history of cardiomyopathy or congenital heart disease.' },
      ];
      window.__hvlExIdx = 0;

      var HVL_SEG_IDS = ['hvl-seg-affected-carriers','hvl-in-trans'];
      window.hvlTryExample = function(){
        var set = function(id,v){ var e=document.getElementById(id); if(e){ e.value=v; e.dispatchEvent(new Event('input',{bubbles:true})); e.dispatchEvent(new Event('change',{bubbles:true})); } };
        var ex = window.HVL_EXAMPLES[window.__hvlExIdx % window.HVL_EXAMPLES.length];
        window.__hvlExIdx++;
        set('hvl-gene', ex.gene); set('hvl-variant', ex.variant);
        set('hvl-hpo', ex.hpo);
        set('hvl-inheritance', ex.inheritance); set('hvl-zygosity', ex.zygosity); set('hvl-sex', ex.sex);
        set('hvl-trio', ex.trio || ''); set('hvl-denovo', ex.denovo || 'unconfirmed');
        set('hvl-famhx', ex.famhx);

        HVL_SEG_IDS.forEach(function(id){ var e=document.getElementById(id); if(e){ if(e.tagName==='SELECT'){ e.selectedIndex=0; } else { e.value=''; } } });

        ['hvl-panel-clinical'].forEach(function(pid){
          var panel = document.getElementById(pid);
          var btn = document.querySelector('.hvl-toggle[aria-controls="' + pid + '"]');
          if (panel && btn && !panel.classList.contains('open')) hvlTogglePanel(btn, pid);
        });

        if (typeof hvlSyncDropdowns === 'function') hvlSyncDropdowns();
      };

      window.hvlClearAll = function(){
        ['hvl-gene','hvl-variant','hvl-hpo','hvl-famhx','hvl-zygosity','hvl-inheritance','hvl-sex','hvl-trio','hvl-denovo'].concat(HVL_SEG_IDS).forEach(function(id){
          var e=document.getElementById(id); if(e){ if(e.tagName==='SELECT'){ e.selectedIndex=0; } else { e.value=''; } }
        });

        if (typeof hvlSyncDropdowns === 'function') hvlSyncDropdowns();

        if (typeof hvlResetBuildToggle === 'function') hvlResetBuildToggle();
        var g=document.getElementById('hvl-gene'); if(g) g.focus();
      };

      window.newVariant = function(){
        window.__hvRan = false;

        ['hvl-gene','hvl-variant','hvl-hpo','hvl-famhx','hvl-zygosity','hvl-inheritance','hvl-sex','hvl-trio','hvl-denovo'].concat(HVL_SEG_IDS).forEach(function(id){
          var e=document.getElementById(id); if(e){ if(e.tagName==='SELECT'){ e.selectedIndex=0; } else { e.value=''; } }
        });

        if (typeof hvlSyncDropdowns === 'function') hvlSyncDropdowns();

        if (typeof hvlResetBuildToggle === 'function') hvlResetBuildToggle();

        ['hvl-panel-clinical'].forEach(function(pid){
          var panel = document.getElementById(pid);
          var btn = document.querySelector('.hvl-toggle[aria-controls="' + pid + '"]');
          if (panel && panel.classList.contains('open') && btn && typeof hvlTogglePanel === 'function') hvlTogglePanel(btn, pid);
        });
        if (typeof returnToLanding === 'function') returnToLanding();
      };
    })();
