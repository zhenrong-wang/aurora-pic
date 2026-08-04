#!/usr/bin/env python3
"""Bounded synthetic regression for the eduPIC stage runner."""

from __future__ import annotations

import json
import hashlib
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts/run_edupic_stage.py"
ACK = "I_UNDERSTAND_THIS_IS_A_BOUNDED_EDUPIC_STAGE"


def checkpoint(path: Path, cycle: int, electrons: int, ions: int) -> None:
    values = [float(cycle) / 13.56e6, float(cycle), float(electrons)]
    values += [0.0] * (4 * electrons)
    values += [float(ions)] + [0.0] * (4 * ions)
    path.write_bytes(struct.pack(f"={len(values)}d", *values))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="aurorapic_edupic_", dir=ROOT / "tmp") as tmp:
        work = Path(tmp)
        source = work / "source"
        source.mkdir()
        checkpoint(source / "picdata.bin", 1, 3, 4)
        (source / "conv.dat").write_text("       1         3         4\n", encoding="utf-8")
        fake = work / "fake_edupic.py"
        fake.write_text(
            "#!/usr/bin/env python3\n"
            "import pathlib, struct, sys\n"
            "p=pathlib.Path('picdata.bin'); d=p.read_bytes(); cycle=int(struct.unpack_from('=d',d,8)[0]); n=int(sys.argv[1]); e=3+n; i=4+2*n\n"
            "v=[(cycle+n)/13.56e6,float(cycle+n),float(e)]+[0.0]*(4*e)+[float(i)]+[0.0]*(4*i); p.write_bytes(struct.pack(f'={len(v)}d',*v))\n"
            "with open('conv.dat','a') as f:\n"
            "  for c in range(cycle+1,cycle+n+1): f.write(f'{c} {3+c-1} {4+2*(c-1)}\\n')\n",
            encoding="utf-8")
        fake.chmod(0o755)
        binary_hash = hashlib.sha256(fake.read_bytes()).hexdigest()
        input_hash = hashlib.sha256((source / "picdata.bin").read_bytes()).hexdigest()
        common = [sys.executable, str(RUNNER), str(fake), str(source),
                  str(work / "stage"), "--cycles", "2",
                  "--expected-binary-sha256", binary_hash,
                  "--expected-input-sha256", input_hash]
        denied = subprocess.run(common, text=True, capture_output=True)
        require(denied.returncode == 2 and ACK in denied.stderr,
                "runner bypassed acknowledgement")
        completed = subprocess.run(
            [*common, "--acknowledge-cost", ACK], text=True, capture_output=True)
        require(completed.returncode == 0, completed.stderr)
        report = json.loads(completed.stdout)
        require(report["completed"] and report["stage"]["end_cycle"] == 3 and
                report["final_state"]["electrons"] == 5 and
                report["final_state"]["ions"] == 8 and
                len(report["new_cycle_population"]) == 2,
                "runner produced an invalid stage report")
        require(hashlib.sha256((source / "picdata.bin").read_bytes()).hexdigest()
                == input_hash,
                "runner modified its input checkpoint")
    print("eduPIC stage regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
