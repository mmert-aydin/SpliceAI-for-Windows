"""Downloads: certificates checked by Windows, and a clear message when that
check still fails (a network that inspects HTTPS).

    .venv\Scripts\python tests\test_downloads.py
"""
import ssl
import sys
import urllib.error
import urllib.request

import helpers
from helpers import check, skip

helpers.watchdog(180)
helpers.setup(need_spliceai=False)

from spliceai_gui import java_download, reference_download  # noqa: E402
from spliceai_pipeline import mane, net  # noqa: E402

check("HTTPS is verified through Windows (truststore), not Python's own list",
      type(net.ssl_context()).__module__.startswith("truststore"), type(net.ssl_context()).__module__)

opened = []
real_urlopen = net.urlopen
net.urlopen = lambda url, timeout=None: opened.append(url) or real_urlopen(url, timeout=timeout)
try:
    try:
        reference_download.download_file("https://example.invalid/x", str(helpers.CACHE / "unused.bin"))
    except reference_download.DownloadError:
        pass
    check("the FASTA/SnpEff download path uses it by default", opened == ["https://example.invalid/x"], str(opened))
finally:
    net.urlopen = real_urlopen

check("the MANE and Java download paths use it by default",
      "net.urlopen" in open(mane.__file__, encoding="utf-8").read()
      and "net.urlopen" in open(java_download.__file__, encoding="utf-8").read())

# The error the user hit, as urllib delivers it.
wrapped = urllib.error.URLError(ssl.SSLCertVerificationError(
    "[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: Missing Authority Key Identifier (_ssl.c:1018)"))
check("that error is recognised as a certificate problem", net.is_certificate_error(wrapped))
message = net.download_error_message(wrapped)
check("...and the message says what to do instead of showing only the error",
      "copy the files from the USB drive" in message and "Technical detail" in message, message.splitlines()[0])
check("an ordinary failure keeps its plain message",
      net.download_error_message(OSError("timed out")) == "Download failed: timed out")


def reachable(url):
    try:
        with net.urlopen(urllib.request.Request(url, method="HEAD")) as response:
            return response.status
    except Exception:
        return None


status = reachable(mane.MANE_CURRENT_DIR_URL)
if status is None and reachable("https://example.com") is None:
    skip("live: the download sites answer", "no internet connection")
else:
    check("live: NCBI answers over a Windows-verified connection", status == 200, f"HTTP {status}")
    check("live: UCSC answers over a Windows-verified connection",
          reachable(reference_download.REFERENCE_URLS["hg19"]) == 200)

helpers.finish()
