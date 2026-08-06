#!/usr/bin/env python3
"""Bounded synthetic regression for the eduPIC measurement-stage runner."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts/run_edupic_measurement_stage.py"
ACK = "I_UNDERSTAND_THIS_IS_A_BOUNDED_EDUPIC_MEASUREMENT_STAGE"


def checkpoint(path: Path, cycle: int, electrons: int, ions: int) -> None:
    values = [float(cycle) / 13.56e6, float(cycle), float(electrons)]
    values += [0.0] * (4 * electrons)
    values += [float(ions)] + [0.0] * (4 * ions)
    path.write_bytes(struct.pack(f"={len(values)}d", *values))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fake_source(valid_diagnostics: bool) -> str:
    diagnostic_code = ""
    if valid_diagnostics:
        diagnostic_code = (
            "density=''.join(f'{i*0.025/399:.12e} 1.0 1.0\\n' for i in range(400)); pathlib.Path('density.dat').write_text(density)\n"
            "eepf=''.join(f'{(i+0.5)*0.05:.8e} {1/(2000*0.05*math.sqrt((i+0.5)*0.05)):.8e}\\n' for i in range(2000)); pathlib.Path('eepf.dat').write_text(eepf)\n"
            "ifed=''.join(f'{(i+0.5):.8e} 5.0e-3 5.0e-3\\n' for i in range(200)); pathlib.Path('ifed.dat').write_text(ifed)\n"
            "row=' '.join(['0']*200)+'\\n'\n"
            "for name in ('pot_xt.dat','efield_xt.dat','ne_xt.dat','ni_xt.dat','je_xt.dat','ji_xt.dat','powere_xt.dat','poweri_xt.dat','meanee_xt.dat','meanei_xt.dat','ioniz_xt.dat'): pathlib.Path(name).write_text(row*400)\n"
            "info=f'''# of simulation cycles in this run    = {n}\\nElectron density @ center             = 1.0e15 [m^{{-3}}]\\nPlasma frequency @ center             = 1.0e9 [rad/s]\\nDebye length @ center                 = 1.0e-4 [m]\\nElectron collision frequency          = 1.0e7 [1/s]\\nIon collision frequency               = 1.0e6 [1/s]\\nIon flux at powered electrode         = 1.0e18 [m^{{-2}} s^{{-1}}]\\nIon flux at grounded electrode        = 1.0e18 [m^{{-2}} s^{{-1}}]\\nMean ion energy at powered electrode  = 10 [eV]\\nMean ion energy at grounded electrode = 11 [eV]\\nElectron flux at powered electrode    = 1.0e18 [m^{{-2}} s^{{-1}}]\\nElectron flux at grounded electrode   = 1.0e18 [m^{{-2}} s^{{-1}}]\\nElectron power density (average)      = 2.0 [W m^{{-3}}]\\nIon power density (average)           = 3.0 [W m^{{-3}}]\\nTotal power density(average)          = 5.0 [W m^{{-3}}]\\n'''\n"
            "pathlib.Path('info.txt').write_text(info)\n")
    return (
        "#!/usr/bin/env python3\n"
        "import math, pathlib, struct, sys\n"
        "assert len(sys.argv)==3 and sys.argv[2]=='m'\n"
        "p=pathlib.Path('picdata.bin'); d=p.read_bytes(); cycle=int(struct.unpack_from('=d',d,8)[0]); n=int(sys.argv[1]); e=3+cycle-1+n; i=4+2*(cycle-1+n)\n"
        "v=[(cycle+n)/13.56e6,float(cycle+n),float(e)]+[0.0]*(4*e)+[float(i)]+[0.0]*(4*i); p.write_bytes(struct.pack(f'={len(v)}d',*v))\n"
        "with open('conv.dat','a') as f:\n"
        "  for c in range(cycle+1,cycle+n+1): f.write(f'{c} {3+c-1} {4+2*(c-1)}\\n')\n"
        + diagnostic_code)


def run_case(work: Path, valid_diagnostics: bool, name: str) -> subprocess.CompletedProcess:
    source = work / f"source-{name}"
    source.mkdir()
    checkpoint(source / "picdata.bin", 1, 3, 4)
    (source / "conv.dat").write_text("1 3 4\n", encoding="utf-8")
    fake = work / f"fake-{name}.py"
    fake.write_text(fake_source(valid_diagnostics), encoding="utf-8")
    fake.chmod(0o755)
    binary_hash = hashlib.sha256(fake.read_bytes()).hexdigest()
    input_hash = hashlib.sha256((source / "picdata.bin").read_bytes()).hexdigest()
    command = [
        sys.executable, str(RUNNER), str(fake), str(source), str(work / f"stage-{name}"),
        "--cycles", "2", "--expected-binary-sha256", binary_hash,
        "--expected-input-sha256", input_hash,
        "--max-initial-particle-steps", "1000000",
        "--acknowledge-cost", ACK,
    ]
    completed = subprocess.run(command, text=True, capture_output=True)
    if valid_diagnostics:
        require(hashlib.sha256((source / "picdata.bin").read_bytes()).hexdigest()
                == input_hash, "runner modified its input checkpoint")
    return completed


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
            prefix="aurorapic_edupic_measurement_", dir=ROOT / "tmp") as tmp:
        work = Path(tmp)
        completed = run_case(work, True, "valid")
        require(completed.returncode == 0, completed.stderr)
        report = json.loads(completed.stdout)
        require(
            report["completed"] and report["stage"]["measurement_mode"] and
            report["stage"]["end_cycle"] == 3 and
            report["reported_observables"]["measurement_cycles"] == 2.0 and
            math.isclose(
                report["outputs"]["diagnostics"]["eepf.dat"]
                ["weighted_normalization"], 1.0, rel_tol=2e-5),
            "runner produced an invalid measurement report")
        rejected = run_case(work, False, "missing")
        require(rejected.returncode == 2 and "density.dat" in rejected.stderr,
                "runner accepted missing native diagnostics")
    print("eduPIC measurement stage regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
