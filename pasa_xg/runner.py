"""Run directory, manifest, stage checkpoints and external commands.

A stage writes ``<run>/<stage_dir>/stage.done.json`` when it succeeds. That file
records a digest of the stage's inputs and settings; the stage is skipped while the
digest still matches and rerun otherwise, so a run resumes where it stopped and a
stale result is never reused.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shlex
import subprocess
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence

log = logging.getLogger("pasa_xg")

SMALL_FILE = 64 << 20  # hash whole files below this size
SAMPLE = 4 << 20       # otherwise hash head, middle and tail samples


def file_digest(path: str) -> str:
    """Content digest; large files are sampled (size + three 4 MB windows)."""
    h = hashlib.sha256()
    size = os.path.getsize(path)
    h.update(str(size).encode())
    with open(path, "rb") as fh:
        if size <= SMALL_FILE:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        else:
            for off in (0, size // 2, max(0, size - SAMPLE)):
                fh.seek(off)
                h.update(fh.read(SAMPLE))
    return h.hexdigest()[:20]


def digest(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:20]


@dataclass
class Run:
    out_dir: str
    settings: Dict[str, Any]
    force_from: Optional[str] = None
    manifest: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.out_dir = os.path.abspath(self.out_dir)
        os.makedirs(self.out_dir, exist_ok=True)
        self.manifest_path = os.path.join(self.out_dir, "manifest.json")
        if os.path.exists(self.manifest_path):
            with open(self.manifest_path) as fh:
                self.manifest = json.load(fh)
        self.manifest.setdefault("stages", {})
        self.manifest["settings"] = self.settings
        self._forcing = False
        self._file_digests: Dict[str, str] = {}

    def path(self, *parts: str) -> str:
        return os.path.join(self.out_dir, *parts)

    def stage_dir(self, name: str) -> str:
        d = self.path(name)
        os.makedirs(d, exist_ok=True)
        return d

    def fdigest(self, path: str) -> str:
        key = os.path.abspath(path)
        if key not in self._file_digests:
            self._file_digests[key] = file_digest(key)
        return self._file_digests[key]

    def save_manifest(self) -> None:
        tmp = self.manifest_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(self.manifest, fh, indent=2, sort_keys=True)
        os.replace(tmp, self.manifest_path)

    def run_stage(self, name: str, dirname: str, inputs: Dict[str, Any],
                  fn: Callable[[str], Dict[str, Any]]) -> Dict[str, Any]:
        """Run `fn(stage_dir)` unless a matching checkpoint exists.

        `inputs` must name everything the stage's result depends on: file paths
        (digested by content) and settings. `fn` returns the stage's outputs
        (paths relative to the stage dir, plus small summary values)."""
        if self.force_from == name:
            self._forcing = True
        sdir = self.stage_dir(dirname)
        done_path = os.path.join(sdir, "stage.done.json")
        key = digest({k: (self.fdigest(v) if isinstance(v, str) and os.path.isfile(v) else v)
                      for k, v in inputs.items()})
        if not self._forcing and os.path.exists(done_path):
            with open(done_path) as fh:
                done = json.load(fh)
            if done.get("key") == key:
                log.info("stage %s: up to date, skipped", name)
                self.manifest["stages"][name] = {"status": "skipped (up to date)", **done}
                self.save_manifest()
                return {k: (os.path.join(sdir, v) if isinstance(v, str) and k.endswith("_path") else v)
                        for k, v in done["outputs"].items()}
        if os.path.exists(done_path):
            os.remove(done_path)
        self._forcing = True  # everything downstream of a rerun stage reruns too
        log.info("stage %s: running", name)
        t0 = time.time()
        self.manifest["stages"][name] = {"status": "running"}
        self.save_manifest()
        outputs = fn(sdir)
        rel = {k: (os.path.relpath(v, sdir) if isinstance(v, str) and k.endswith("_path") else v)
               for k, v in outputs.items()}
        done = {"key": key, "outputs": rel, "seconds": round(time.time() - t0, 1),
                "finished": time.strftime("%Y-%m-%d %H:%M:%S")}
        with open(done_path, "w") as fh:
            json.dump(done, fh, indent=2, sort_keys=True)
        self.manifest["stages"][name] = {"status": "done", **done}
        self.save_manifest()
        return outputs


def run_cmd(cmd: Sequence[str] | str, log_path: Optional[str] = None, cwd: Optional[str] = None,
            stdout_path: Optional[str] = None) -> None:
    """Run an external command; its stderr (and stdout unless redirected) go to log_path."""
    shell = isinstance(cmd, str)
    text = cmd if shell else " ".join(shlex.quote(str(c)) for c in cmd)
    log.info("cmd: %s", text)
    logf = open(log_path, "a") if log_path else None
    outf = open(stdout_path, "w") if stdout_path else logf
    try:
        if logf:
            logf.write(f"\n# {time.strftime('%Y-%m-%d %H:%M:%S')}  {text}\n")
            logf.flush()
        rc = subprocess.call(text if shell else [str(c) for c in cmd], shell=shell, cwd=cwd,
                             stdout=outf, stderr=logf,
                             executable="/bin/bash" if shell else None)
    finally:
        if stdout_path and outf:
            outf.close()
        if logf:
            logf.close()
    if rc != 0:
        raise RuntimeError(f"command failed (exit {rc}): {text}" + (f"\nsee {log_path}" if log_path else ""))


def tool_version(cmd: List[str]) -> str:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        text = (out.stdout or out.stderr).strip().splitlines()
        return text[0] if text else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unavailable"
