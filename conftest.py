"""Ensure the project root is importable from anywhere + pytest config."""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def pytest_configure(config):
    config.addinivalue_line("markers", "smoke: lightweight end-to-end run")
