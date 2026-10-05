#!/usr/bin/env python3
"""Rebuild the certificates stored with the A340984 representatives.

Copied unmodified from certificate.py in github.com/meghanto/OEIS_A340984,
so that enum_prime.py can de-duplicate on the catalogue's own D4 certificate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def neighbors(rects):
    nb = {}
    for i, r in enumerate(rects):
        right = sorted([j for j, s in enumerate(rects)
                        if s[0] == r[2] and s[1] < r[3] and s[3] > r[1]],
                       key=lambda j: rects[j][1])
        bottom = sorted([j for j, s in enumerate(rects)
                         if s[1] == r[3] and s[0] < r[2] and s[2] > r[0]],
                        key=lambda j: rects[j][0])
        left = sorted([j for j, s in enumerate(rects)
                       if s[2] == r[0] and s[1] < r[3] and s[3] > r[1]],
                      key=lambda j: rects[j][1])
        top = sorted([j for j, s in enumerate(rects)
                      if s[3] == r[1] and s[0] < r[2] and s[2] > r[0]],
                     key=lambda j: rects[j][0])
        nb[i] = (right, bottom, left, top)
    return nb


def cert(rects, W, H):
    nb = neighbors(rects)
    try:
        tl = next(i for i, r in enumerate(rects) if r[0] == 0 and r[1] == 0)
    except StopIteration:
        raise ValueError("No top-left rectangle starting at (0,0) found.")
    ids = {tl: 0}
    order = [tl]
    qi = 0
    while qi < len(order):
        cur = order[qi]
        qi += 1
        for side in nb[cur]:
            for j in side:
                if j not in ids:
                    ids[j] = len(order)
                    order.append(j)
    if len(order) != len(rects):
        raise ValueError(f"Disconnected tiling: traversal reached {len(order)} of {len(rects)} rectangles.")
    out = []
    for i in order:
        out.append(tuple(tuple(ids[j] for j in side) for side in nb[i]))
    return tuple(out)


def transform(rects, W, H, t):
    out = []
    for (x1, y1, x2, y2) in rects:
        if t == 0:
            nx1, ny1, nx2, ny2, WW, HH = x1, y1, x2, y2, W, H
        elif t == 1:  # rot90 cw: (x,y)->(H-y, x)
            nx1, ny1, nx2, ny2, WW, HH = H - y2, x1, H - y1, x2, H, W
        elif t == 2:  # rot180
            nx1, ny1, nx2, ny2, WW, HH = W - x2, H - y2, W - x1, H - y1, W, H
        elif t == 3:  # rot270
            nx1, ny1, nx2, ny2, WW, HH = y1, W - x2, y2, W - x1, H, W
        elif t == 4:  # mirror horizontal (flip x)
            nx1, ny1, nx2, ny2, WW, HH = W - x2, y1, W - x1, y2, W, H
        elif t == 5:  # mirror vertical (flip y)
            nx1, ny1, nx2, ny2, WW, HH = x1, H - y2, x2, H - y1, W, H
        elif t == 6:  # transpose
            nx1, ny1, nx2, ny2, WW, HH = y1, x1, y2, x2, H, W
        else:  # anti-transpose
            nx1, ny1, nx2, ny2, WW, HH = H - y2, W - x2, H - y1, W - x1, H, W
        out.append((nx1, ny1, nx2, ny2))
    return out, WW, HH


def reading_order(rectangles):
    return tuple(sorted(rectangles, key=lambda r: (r[1], r[0], r[3], r[2])))


def exact_cover(rectangles, width: int, height: int) -> bool:
    if width <= 0 or height <= 0:
        return False
    cells = [[0] * width for _ in range(height)]
    for r in rectangles:
        if len(r) != 4:
            return False
        x1, y1, x2, y2 = r
        if not (0 <= x1 < x2 <= width and 0 <= y1 < y2 <= height):
            return False
        for y in range(y1, y2):
            for x in range(x1, x2):
                cells[y][x] += 1
    return all(value == 1 for row in cells for value in row)


def canonical_d4_certificate(rectangles, width: int, height: int):
    certificates = []
    for action in range(8):
        moved, moved_width, moved_height = transform(
            rectangles, width, height, action
        )
        ordered = reading_order(moved)
        certificates.append(cert(list(ordered), moved_width, moved_height))
    return min(certificates)


def jsonable_certificate(value):
    return [[[int(neighbor) for neighbor in side] for side in room]
            for room in value]


def run_jsonl_verification(jsonl_path: str, index: int | None, class_id: str | None) -> int:
    path = Path(jsonl_path)
    if not path.is_file():
        raise FileNotFoundError(f"JSONL file not found: {jsonl_path}")

    target_record = None
    seen_count = 0
    with path.open(encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as e:
                raise ValueError(f"Invalid JSON at line {line_no}: {e}")

            seen_count += 1
            if index is not None and seen_count == index:
                target_record = record
                break
            if class_id is not None and record.get("class_id") == class_id:
                target_record = record
                break

    if target_record is None:
        selector = f"class_id '{class_id}'" if class_id is not None else f"index {index}"
        raise ValueError(f"Representative not found for {selector}")

    width = target_record.get("width")
    height = target_record.get("height")
    rectangles_raw = target_record.get("rectangles")
    stored_cert = target_record.get("d4_certificate")
    rec_class_id = target_record.get("class_id")

    if width is None or height is None or rectangles_raw is None:
        raise ValueError("Record is missing required fields (width, height, or rectangles)")

    rectangles = [tuple(r) for r in rectangles_raw]

    if not exact_cover(rectangles, width, height):
        raise ValueError(f"Record {rec_class_id} does not represent an exact cover or has invalid bounds.")

    recomputed = canonical_d4_certificate(rectangles, width, height)
    recomputed_json = jsonable_certificate(recomputed)

    matches_stored = (recomputed_json == stored_cert)

    result = {
        "class_id": rec_class_id,
        "matches_stored": matches_stored,
        "d4_certificate": recomputed_json,
    }
    print(json.dumps(result, separators=(",", ":")))

    return 0 if matches_stored else 1


def run_raw_verification(raw_path: str, line_index: int) -> int:
    path = Path(raw_path)
    if not path.is_file():
        raise FileNotFoundError(f"Raw source file not found: {raw_path}")

    target_line = None
    with path.open(encoding="utf-8") as f:
        for current_line_no, line in enumerate(f, 1):
            if current_line_no == line_index:
                target_line = line
                break

    if target_line is None:
        raise ValueError(f"Line {line_index} not found in raw file {raw_path}")

    fields = target_line.split()
    if not fields:
        raise ValueError(f"Line {line_index} is empty or only whitespace.")

    marked_crossed = (fields[0] == "X")
    if marked_crossed:
        fields = fields[1:]

    if len(fields) < 2:
        raise ValueError("Malformed raw line: missing width or height.")

    try:
        values = [int(v) for v in fields]
    except ValueError as e:
        raise ValueError(f"Malformed raw line: non-integer field: {e}")

    width, height = values[0], values[1]
    coord_fields = values[2:]

    if len(coord_fields) % 4 != 0:
        raise ValueError(f"Malformed raw line: coordinate values count ({len(coord_fields)}) must be a multiple of 4.")

    rectangles = [tuple(coord_fields[i:i+4]) for i in range(0, len(coord_fields), 4)]

    if not exact_cover(rectangles, width, height):
        raise ValueError(f"Line {line_index} does not represent an exact cover or has invalid bounds.")

    sorted_rectangles = list(reading_order(rectangles))
    oriented = cert(sorted_rectangles, width, height)
    oriented_json = jsonable_certificate(oriented)

    canonical = canonical_d4_certificate(rectangles, width, height)
    canonical_json = jsonable_certificate(canonical)

    result = {
        "width": width,
        "height": height,
        "rectangle_count": len(rectangles),
        "crossed": marked_crossed,
        "oriented_certificate": oriented_json,
        "d4_certificate": canonical_json,
    }
    print(json.dumps(result, separators=(",", ":")))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="A340984 Minimal Certificate Verifier")
    parser.add_argument("--jsonl", help="Path to JSONL representative catalogue file")
    parser.add_argument("--index", type=int, help="One-based record index to select from --jsonl")
    parser.add_argument("--class-id", help="Exact class_id to select from --jsonl")
    parser.add_argument("--raw", help="Path to raw source file")
    parser.add_argument("--line", type=int, help="One-based line number to select from --raw")

    args = parser.parse_args()

    if args.jsonl:
        if args.raw or args.line:
            print("Error: --jsonl and --raw/--line options are mutually exclusive.", file=sys.stderr)
            return 1
        if (args.index is None) == (args.class_id is None):
            print("Error: Specify exactly one of --index or --class-id with --jsonl.", file=sys.stderr)
            return 1
        if args.index is not None and args.index <= 0:
            print("Error: --index must be a positive, one-based integer.", file=sys.stderr)
            return 1

        try:
            return run_jsonl_verification(args.jsonl, args.index, args.class_id)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    elif args.raw:
        if args.index or args.class_id:
            print("Error: --raw and --index/--class-id options are mutually exclusive.", file=sys.stderr)
            return 1
        if args.line is None:
            print("Error: --line must be specified with --raw.", file=sys.stderr)
            return 1
        if args.line <= 0:
            print("Error: --line must be a positive, one-based integer.", file=sys.stderr)
            return 1

        try:
            return run_raw_verification(args.raw, args.line)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
    else:
        print("Error: Must specify either --jsonl or --raw mode.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
