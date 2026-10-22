"""Check or format JSON schema files used by the project."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

SCHEMA_DIRECTORY = Path(__file__).parents[1] / "src/fairmd/lipids/schema_validation/schema"


def format_json(path: Path) -> str:
    """Return a consistently formatted representation of a JSON file."""
    with path.open(encoding="utf-8") as stream:
        document = json.load(stream)
    return json.dumps(document, indent=2, ensure_ascii=False) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when a JSON schema is invalid or not formatted",
    )
    args = parser.parse_args()

    failed = False
    for path in sorted(SCHEMA_DIRECTORY.glob("*.json")):
        try:
            formatted = format_json(path)
        except (OSError, json.JSONDecodeError) as error:
            print(f"{path}: {error}")
            failed = True
            continue

        current = path.read_text(encoding="utf-8")
        if current != formatted:
            failed = True
            print(f"{path}: not formatted")
            if not args.check:
                path.write_text(formatted, encoding="utf-8")

    return int(failed and args.check)


if __name__ == "__main__":
    raise SystemExit(main())
