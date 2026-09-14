"""HTTPS downloads that trust what Windows trusts.

Python checks certificates with its own bundled OpenSSL, which is stricter
than Windows and knows nothing about the certificates a company network
installs. On a PC whose network (or antivirus) inspects HTTPS -- common in
hospitals -- every download in this app failed with an error Windows itself
accepts, e.g.

    CERTIFICATE_VERIFY_FAILED: Missing Authority Key Identifier

truststore hands the check to Windows instead, so whatever Windows trusts
(including a company's own certificate) works here too -- the same approach
pip takes, and the same trust store the Windows installer already uses for
its SpliceAI download. Without truststore this falls back to Python's own
default, i.e. the previous behaviour.
"""
import ssl
import urllib.request

try:
    import truststore
except ImportError:  # falls back to Python's own certificate list
    truststore = None

DEFAULT_TIMEOUT = 60

CERTIFICATE_ADVICE = (
    "The secure connection to this site couldn't be verified, even using Windows' own "
    "certificate list. That usually means this network (or an antivirus program) inspects "
    "secure connections, and its certificate isn't installed on this PC.\n\n"
    "You don't have to download this: copy the files from the USB drive's reference-data "
    "folder and pick them with Browse... (or drag them onto the window). Otherwise ask your "
    "IT department to install the network's certificate on this PC."
)

_context = None


def ssl_context():
    """Certificate checking done by Windows when truststore is available."""
    global _context
    if _context is None:
        _context = truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT) if truststore else ssl.create_default_context()
    return _context


def urlopen(url, timeout=DEFAULT_TIMEOUT):
    """urllib.request.urlopen with Windows-backed certificate checking. `url`
    may be a string or a urllib Request (java_download sends a User-Agent)."""
    return urllib.request.urlopen(url, timeout=timeout, context=ssl_context())


def is_certificate_error(exc):
    return isinstance(exc, ssl.SSLCertVerificationError) or "CERTIFICATE_VERIFY_FAILED" in str(exc)


def download_error_message(exc, context=None):
    """What to show the user for a failed download -- with what to do about
    it when the certificate check is what failed. context: what was being
    downloaded, shown first (e.g. "Could not list <url>.")."""
    if is_certificate_error(exc):
        detail = f"{CERTIFICATE_ADVICE}\n\nTechnical detail: {exc}"
    else:
        detail = f"Download failed: {exc}"
    return f"{context}\n\n{detail}" if context else detail
