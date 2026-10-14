#!/usr/bin/env python3
"""Print an explicitly hypothetical two-species timing-competition scenario.

Run after editable install, or use:
    PYTHONPATH=src python scripts/simulate_temporal_supply_competition.py
"""
from __future__ import annotations

import json

from butterfly_specialization_ecology.temporal_supply_competition import evaluate


def main() -> None:
    example = {"A": (0.75, 0.25), "B": (0.25, 0.75)}
    print(json.dumps(evaluate(example, 0.20), indent=2))


if __name__ == "__main__":
    main()
