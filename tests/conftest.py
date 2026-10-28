import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from v1lb.config import raw_root  # noqa: E402


def pytest_collection_modifyitems(config, items):
    if raw_root().exists():
        return
    skip = pytest.mark.skip(reason=f"raw data not available at {raw_root()}")
    for item in items:
        if "raw" in item.keywords:
            item.add_marker(skip)
