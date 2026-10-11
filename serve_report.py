import http.server
import socketserver
import webbrowser
import threading
import os
import sys

PORT = 8000
DIRECTORY = "ldmt-experiment-report"


def open_browser():
    webbrowser.open(f"http://localhost:{PORT}")


def main():
    os.chdir(DIRECTORY)
    handler = http.server.SimpleHTTPRequestHandler

    # 启动后延迟 0.8 秒再打开浏览器，确保服务已就绪
    threading.Timer(0.8, open_browser).start()

    with socketserver.TCPServer(("", PORT), handler) as httpd:
        print(f"Serving at http://localhost:{PORT}")
        print("Press Ctrl+C to stop.")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nServer stopped.")


if __name__ == "__main__":
    main()