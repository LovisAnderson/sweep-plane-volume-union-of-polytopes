#!/usr/bin/env python3
"""Local browser viewer. Run after cargo build --release; no Python dependencies."""
import argparse
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MAX_BODY = 8 * 1024 * 1024


class Handler(BaseHTTPRequestHandler):
    def reply(self, status, body, content_type="application/json"):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        paths = {"/": ("index.html", "text/html; charset=utf-8"),
                 "/app.js": ("app.js", "text/javascript; charset=utf-8"),
                 "/math.js": ("math.js", "text/javascript; charset=utf-8")}
        if self.path not in paths:
            self.reply(404, b'{"error":"Not found"}')
            return
        name, mime = paths[self.path]
        self.reply(200, (ROOT / "viewer" / name).read_bytes(), mime)

    def do_POST(self):
        try:
            # Only the local UI may submit work; no cross-origin or DNS-rebound requests.
            host = self.headers.get("Host", "")
            if host not in {f"127.0.0.1:{self.server.server_port}", f"localhost:{self.server.server_port}"}:
                raise ValueError("Invalid host")
            if self.headers.get("Origin", f"http://{host}") != f"http://{host}":
                raise ValueError("Invalid origin")
            if self.path != "/api/compute":
                self.reply(404, b'{"error":"Not found"}')
                return
            size = int(self.headers.get("Content-Length", "0"))
            if not 0 < size <= MAX_BODY:
                raise ValueError("Input must be between 1 byte and 8 MiB")
            payload = json.loads(self.rfile.read(size))
            files = payload.get("files")
            direction = payload.get("direction")
            if not isinstance(files, list) or not files or not all(isinstance(f, str) for f in files):
                raise ValueError("Choose at least one input file")
            if not isinstance(direction, str) or len(direction.split(",")) != 2:
                raise ValueError("Enter two direction coordinates, e.g. 1,2")
            with tempfile.TemporaryDirectory(prefix="nefvol-viewer-") as tmp:
                paths = []
                for i, content in enumerate(files):
                    path = Path(tmp) / f"input-{i}"
                    path.write_text(content)
                    paths.append(str(path))
                result = subprocess.run([str(self.server.binary), "view-data", "--input", *paths,
                                         f"--direction={direction}", "--threads", "2"],
                                        capture_output=True, text=True, timeout=120)
            if result.returncode:
                raise ValueError(result.stderr.strip() or "Computation failed")
            self.reply(200, result.stdout.encode())
        except (ValueError, TypeError, AttributeError, OSError, subprocess.TimeoutExpired) as error:
            self.reply(400, json.dumps({"error": str(error)}).encode())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/nefvol")
    args = parser.parse_args()
    if not args.binary.is_file():
        parser.error("Binary not found; run cargo build --release first or supply --binary")
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    server.binary = args.binary.resolve()
    print(f"Sweep viewer: http://127.0.0.1:{server.server_port} (Ctrl+C to stop)", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
