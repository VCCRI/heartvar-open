"""Tests for the addressable page routes — GET /about and GET /contact.

HeartVar is a single-page app: About and Contact were sibling <div>s that
``showPage()`` toggled with ``display``, so the URL never changed. The cost was
paid by the user, not the code — the browser Back button left the site instead
of returning to the curator, a refresh on About dropped you back at the landing,
and there was no way to send a colleague to the Data & privacy section.

Both paths now serve the same ``index.html`` and the client routes off
``location.pathname``. The properties worth pinning:

  * ``/about`` and ``/contact`` serve the SAME document as ``/`` — one shell, so
    a page can never drift from the app it is mounted in;
  * they are declared as EXPLICIT paths, not a catch-all: a greedy
    ``/{path:path}`` would shadow ``/api/*``, ``/login``, ``/health`` and turn
    every API typo into a 200 with an HTML body. The unknown-path and
    still-routes-correctly tests below are what stop that regression;
  * they are no-store, like ``/`` — a deploy must not leave a browser holding a
    stale shell (the same reason _RevalidateStaticFiles exists);
  * ``?next=`` on /login only accepts a site-relative path. Making ``/about`` a
    real route means the sign-in return path is now a URL a user can influence,
    so the open-redirect guard is tested here beside the routes that motivated
    it.
"""
from __future__ import annotations

import pytest
from starlette.testclient import TestClient

from backend.app import PROJECT_ROOT, _safe_next, app

PAGE_PATHS = ("/about", "/contact")


@pytest.fixture
def client():
    with TestClient(app) as c:
        yield c


@pytest.mark.parametrize("path", PAGE_PATHS)
def test_page_path_serves_the_app_shell(client, path):
    r = client.get(path)
    assert r.status_code == 200
    assert "text/html" in r.headers["content-type"]
    assert 'id="page-about"' in r.text
    assert 'id="page-contact"' in r.text


@pytest.mark.parametrize("path", PAGE_PATHS)
def test_page_path_is_byte_identical_to_root(client, path):
    """One shell for every route — no second copy to drift out of sync."""
    assert client.get(path).content == client.get("/").content


@pytest.mark.parametrize("path", ("/",) + PAGE_PATHS)
def test_shell_is_never_cached(client, path):
    cc = client.get(path).headers.get("cache-control", "")
    assert "no-store" in cc


@pytest.mark.parametrize("path", ("/nonsense", "/about/privacy", "/curate"))
def test_unknown_paths_still_404(client, path):
    """No catch-all. Only the declared pages are addressable; section anchors
    stay client-side fragments (/about#privacy), which never reach the server."""
    assert client.get(path).status_code == 404


def test_api_and_infra_paths_are_not_shadowed(client):
    """A greedy page route would answer these with the HTML shell."""
    for path in ("/api/data-status", "/health", "/healthz"):
        r = client.get(path)
        assert r.status_code == 200, path
        assert "text/html" not in r.headers.get("content-type", ""), path
    r = client.get("/api/definitely-not-a-route")
    assert r.status_code == 404
    assert "text/html" not in r.headers.get("content-type", "")


def test_index_exists_at_project_root():
    """The routes are only as good as the file they serve."""
    assert (PROJECT_ROOT / "index.html").exists()


@pytest.mark.parametrize("raw", ["/", "/about", "/contact", "/about?x=1"])
def test_safe_next_allows_site_relative_paths(raw):
    assert _safe_next(raw) == raw


@pytest.mark.parametrize("raw", [
    "https://evil.example/phish",
    "//evil.example/phish",
    "/\\evil.example",
    "http:/evil.example",
    "about",
    "",
    None,
])
def test_safe_next_rejects_offsite_targets(raw):
    assert _safe_next(raw) == "/"


def test_a_browser_gets_a_styled_404_page():
    c = TestClient(app)
    r = c.get("/no-such-page", headers={"accept": "text/html,application/xhtml+xml"})
    assert r.status_code == 404
    assert "text/html" in r.headers["content-type"]
    assert "We couldn't find that page" in r.text
    assert 'href="/"' in r.text


def test_an_api_path_still_gets_json_even_from_a_browser():
    """🔴 The important one. If this ever returns HTML, an /api/ typo starts
    looking like an empty result instead of a failure."""
    c = TestClient(app)
    r = c.get("/api/no-such-endpoint",
              headers={"accept": "text/html,application/xhtml+xml"})
    assert r.status_code == 404
    assert "application/json" in r.headers["content-type"]
    assert r.json() == {"detail": "Not Found"}


def test_a_non_browser_client_still_gets_json():
    c = TestClient(app)
    r = c.get("/no-such-page", headers={"accept": "*/*"})
    assert r.status_code == 404
    assert "application/json" in r.headers["content-type"]


def test_non_404_errors_are_not_rewritten():
    """The handler intercepts EVERY HTTPException, so anything that is not a
    404 must fall through to FastAPI's own handler untouched.

    This used to POST /api/contact for its non-404.  The contact form was
    removed on 2026-09-08 (VCCRI IT declined to provide an SMTP relay), so the
    403 an unauthenticated caller gets from ``_require_admin`` stands in: same
    handler path, and it needs no credentials and no request body."""
    c = TestClient(app)
    r = c.get("/api/admin/log-index", headers={"accept": "text/html"})
    assert r.status_code == 403
    assert "application/json" in r.headers["content-type"]


def test_the_real_pages_are_unaffected():
    c = TestClient(app)
    for path in ("/", "/about", "/contact"):
        r = c.get(path, headers={"accept": "text/html"})
        assert r.status_code == 200, path
        assert "<html" in r.text.lower(), path
        assert "We couldn't find that page" not in r.text, path
