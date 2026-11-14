// Light and dark theme for the series' Plotly figures, shared by the R and Python scripts.
// Each script draws its figure in the light theme and then calls statgenTheme(gd, minWidth).
// Inside the blog the figure sits in an iframe on the same site, so it follows the site's own
// theme toggle (data-theme on the parent page); opened on its own, it follows the OS setting.
//  - text in the default ink switches ink; text with its own colour (a highlighted label) keeps it
//  - marker rings drawn in the light surface colour switch to the dark surface colour
//  - traces with meta "seq" (heatmaps) get a sequential scale that starts at the surface colour
//  - below minWidth pixels the figure keeps minWidth and scrolls sideways
function statgenTheme(gd, minWidth) {
  var LIGHT = {ink: '#0b0b0b', ink2: '#52514e', line: '#c9c8c2', surface: '#fcfcfb',
               grid: 'rgba(128,128,128,0.16)', seq: [[0, '#fcfcfb'], [0.15, '#cde2fb'], [0.5, '#6da7ec'], [1, '#104281']]};
  var DARK = {ink: '#f0efec', ink2: '#c3c2b7', line: '#4a4a47', surface: '#1a1a19',
              grid: 'rgba(160,160,160,0.16)', seq: [[0, '#1a1a19'], [0.15, '#184f95'], [0.5, '#3987e5'], [1, '#cde2fb']]};
  var INKS = [LIGHT.ink, DARK.ink, LIGHT.ink2, DARK.ink2];
  function isDark() {
    try { var t = window.parent.document.documentElement.dataset.theme; if (t) return t === 'dark'; } catch (err) {}
    return !!(window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches);
  }
  function paint() {
    var t = isDark() ? DARK : LIGHT;
    document.body.style.background = 'transparent';
    document.documentElement.style.colorScheme = isDark() ? 'dark' : 'light';
    var upd = {'font.color': t.ink, 'legend.font.color': t.ink};
    Object.keys(gd.layout).forEach(function (k) {
      if (/^[xy]axis\d*$/.test(k)) {
        upd[k + '.gridcolor'] = t.grid; upd[k + '.linecolor'] = t.line; upd[k + '.zerolinecolor'] = t.line;
      }
    });
    upd.annotations = (gd.layout.annotations || []).map(function (a) {
      var c = a.font && a.font.color;
      if (!c || INKS.indexOf(c) >= 0) {
        var second = c === LIGHT.ink2 || c === DARK.ink2;
        a.font = Object.assign({}, a.font, {color: second ? t.ink2 : t.ink});
      }
      if (a.arrowcolor && INKS.indexOf(a.arrowcolor) >= 0) a.arrowcolor = t.ink2;
      return a;
    });
    upd.shapes = (gd.layout.shapes || []).map(function (s) {
      if (s.line && [LIGHT.line, DARK.line].indexOf(s.line.color) >= 0) s.line = Object.assign({}, s.line, {color: t.line});
      return s;
    });
    Plotly.relayout(gd, upd);
    gd.data.forEach(function (tr, i) {
      var m = tr.marker && tr.marker.line && tr.marker.line.color;
      if (m === LIGHT.surface || m === DARK.surface) Plotly.restyle(gd, {'marker.line.color': t.surface}, [i]);
      if (tr.meta === 'seq') Plotly.restyle(gd, {colorscale: [t.seq]}, [i]);
    });
  }
  if (minWidth && window.innerWidth < minWidth) {
    document.body.style.overflowX = 'auto'; document.body.style.margin = '0';
    Plotly.relayout(gd, {width: minWidth});
  }
  paint();
  if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', paint);
  try {
    new MutationObserver(paint).observe(window.parent.document.documentElement, {attributes: true, attributeFilter: ['data-theme']});
  } catch (err) {}
}
