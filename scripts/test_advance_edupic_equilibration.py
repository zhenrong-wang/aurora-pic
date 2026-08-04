#!/usr/bin/env python3
"""Synthetic regression for adaptive eduPIC equilibration advancement."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
ADVANCE = ROOT / "scripts/advance_edupic_equilibration.py"
ACK = "I_UNDERSTAND_THIS_ADVANCES_BOUNDED_EDUPIC_EQUILIBRATION"


def checkpoint(path: Path, cycle: int, electrons: int, ions: int) -> None:
    values = [cycle / 13.56e6, float(cycle), float(electrons)]
    values += [0.0] * (4 * electrons)
    values += [float(ions)] + [0.0] * (4 * ions)
    path.write_bytes(struct.pack(f"={len(values)}d", *values))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="aurorapic_edupic_advance_",
                                     dir=ROOT / "tmp") as tmp:
        work = Path(tmp)
        source = work / "source"
        source.mkdir()
        checkpoint(source / "picdata.bin", 1, 3, 4)
        (source / "conv.dat").write_text("1 3 4\n", encoding="utf-8")
        fake = work / "fake_edupic.py"
        fake.write_text(
            "#!/usr/bin/env python3\n"
            "import pathlib,struct,sys\n"
            "p=pathlib.Path('picdata.bin');d=p.read_bytes();c=int(struct.unpack_from('=d',d,8)[0]);n=int(sys.argv[1]);e=3+n+c-1;i=4+2*(n+c-1)\n"
            "v=[(c+n)/13.56e6,float(c+n),float(e)]+[0.0]*(4*e)+[float(i)]+[0.0]*(4*i);p.write_bytes(struct.pack(f'={len(v)}d',*v))\n"
            "with open('conv.dat','a') as f:\n"
            "  for x in range(c+1,c+n+1):f.write(f'{x} {3+x-1} {4+2*(x-1)}\\n')\n",
            encoding="utf-8")
        fake.chmod(0o755)
        binary_hash = hashlib.sha256(fake.read_bytes()).hexdigest()
        input_hash = hashlib.sha256((source / "picdata.bin").read_bytes()).hexdigest()
        common = [sys.executable, str(ADVANCE), str(fake), str(source),
                  str(work / "campaign"), "--expected-binary-sha256", binary_hash,
                  "--expected-input-sha256", input_hash, "--target-cycle", "5",
                  "--max-wall-seconds", "10", "--stage-timeout-seconds", "5",
                  "--max-stage-cycles", "2",
                  "--max-stage-initial-particle-steps", "100000"]
        denied = subprocess.run(common, text=True, capture_output=True)
        require(denied.returncode == 2 and ACK in denied.stderr,
                "coordinator bypassed acknowledgement")
        completed = subprocess.run([*common, "--acknowledge-cost", ACK],
                                   text=True, capture_output=True)
        require(completed.returncode == 0, completed.stderr)
        report = json.loads(completed.stdout)
        require(report["target_reached"] and report["latest_state"]["cycles"] == 5
                and len(report["stages"]) == 3 and
                report["stages"][0]["end_cycle"] == 3,
                "coordinator did not complete exact adaptive stages")
        require(hashlib.sha256((source / "picdata.bin").read_bytes()).hexdigest()
                == input_hash, "coordinator modified its input checkpoint")
    print("eduPIC adaptive advancement regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
