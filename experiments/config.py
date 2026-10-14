"""Single configuration entry point (v3.2, plan 5.3).

The V2.0 plan flagged an inconsistency: ``configs/*.yaml`` existed but was stale
and read by *nothing*, so a user editing YAML saw no effect.  Plan 5.3 requires
picking ONE source of truth.  The choice made here is:

    **YAML is the authoritative experiment definition; CLI flags override it,
    and every override is validated and echoed.**

That combination keeps the "just run it" convenience of CLI flags while giving a
single, reviewable, version-controlled place to record a reported experiment.
Anything not expressible in YAML is rejected rather than silently ignored.

Usage
-----
    from experiments.config import load_config
    cfg = load_config("configs/main_experiment.yaml", overrides={"tau": 25})

The loader returns a plain dict of *resolved* values plus the list of overrides
that were applied, so a run can record exactly what differed from the file.
"""
import os

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
CONFIG_DIR = os.path.join(PROJ, "configs")

#: every key a config file may set.  Unknown keys are an ERROR (plan 5.3: no
#: silent divergence between documentation and behaviour).
ALLOWED_KEYS = {
    "name", "description",
    "problems", "taus", "seeds", "algorithms", "algos",
    "N", "n", "n_gen", "n_gen_total", "n_t", "tau",
    "out", "phase", "code_version", "experiment_group",
    "algo_kwargs", "limit_seeds", "n_points", "change_mode",
    # documentation-only keys that several configs carry; they are validated
    # for presence/shape but never silently reinterpreted
    "authoritative_problems", "variants", "studies", "variants_note",
}


class ConfigError(ValueError):
    """Raised for an unknown key, a bad type, or an unreadable file."""


def _load_yaml(path):
    if not os.path.isabs(path):
        candidate = os.path.join(PROJ, path)
        if not os.path.exists(candidate):
            candidate = os.path.join(CONFIG_DIR, os.path.basename(path))
        path = candidate
    if not os.path.exists(path):
        raise ConfigError(f"config file not found: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ConfigError(f"config root must be a mapping, got {type(data).__name__}")
    unknown = sorted(set(data) - ALLOWED_KEYS)
    if unknown:
        raise ConfigError(
            f"unknown key(s) in {os.path.basename(path)}: {unknown}. "
            f"Allowed: {sorted(ALLOWED_KEYS)}")
    return data


def _parse_list(value):
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    return [v.strip() for v in str(value).split(",") if v.strip()]


def _parse_seeds(value):
    """Accept [0,1,2], "0-29", "0,1,2" or an int."""
    if value is None:
        return None
    if isinstance(value, int):
        return list(range(value))
    if isinstance(value, (list, tuple)):
        return [int(v) for v in value]
    s = str(value)
    if "-" in s:
        lo, hi = s.split("-", 1)
        return list(range(int(lo), int(hi) + 1))
    return [int(v) for v in s.split(",") if v.strip()]


def available_configs():
    if not os.path.isdir(CONFIG_DIR):
        return []
    return sorted(f for f in os.listdir(CONFIG_DIR) if f.endswith((".yaml", ".yml")))


def load_config(path, overrides=None):
    """Load a YAML experiment definition, applying and validating overrides.

    Returns ``(resolved, applied)`` where ``applied`` maps each overridden key to
    ``(old, new)``.  Overrides whose value is ``None`` are ignored, so an
    argparse default of ``None`` means "use the file".
    """
    data = _load_yaml(path)
    overrides = {k: v for k, v in (overrides or {}).items() if v is not None}
    unknown = sorted(set(overrides) - ALLOWED_KEYS)
    if unknown:
        raise ConfigError(f"unknown override key(s): {unknown}")

    applied = {}
    for k, v in overrides.items():
        old = data.get(k)
        if old != v:
            applied[k] = (old, v)
        data[k] = v

    # normalise
    data["problems"] = _parse_list(data.get("problems"))
    data["algorithms"] = _parse_list(
        data.get("algorithms", data.get("algos")))
    data["taus"] = ([int(t) for t in _parse_list(data.get("taus"))]
                    if data.get("taus") is not None
                    else (None if data.get("tau") is None else [int(data["tau"])]))
    data["seeds"] = _parse_seeds(data.get("seeds"))
    for alias, canonical in (("n", "N"), ("n_gen_total", "n_gen")):
        if data.get(canonical) is None and data.get(alias) is not None:
            data[canonical] = data[alias]
    for key in ("N", "n_gen", "n_t"):
        if data.get(key) is not None:
            data[key] = int(data[key])
    if data.get("limit_seeds") is not None:
        data["limit_seeds"] = int(data["limit_seeds"])
    if data.get("algo_kwargs") is None:
        data["algo_kwargs"] = {}
    if not isinstance(data["algo_kwargs"], dict):
        raise ConfigError("algo_kwargs must be a mapping")
    data["_path"] = os.path.abspath(
        path if os.path.isabs(path) else os.path.join(PROJ, path))
    data["_overrides"] = applied
    return data, applied


def describe(resolved):
    """One-line human summary used to echo the effective configuration."""
    return (f"problems={resolved.get('problems')} taus={resolved.get('taus')} "
            f"seeds={len(resolved.get('seeds') or [])} "
            f"algos={resolved.get('algorithms')} N={resolved.get('N')} "
            f"n_gen={resolved.get('n_gen')} n_t={resolved.get('n_t')} "
            f"group={resolved.get('experiment_group')}")
