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

from advance_edupic_equilibration import host_guard_violations


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
    policy = {"maximum_load_per_cpu": 0.5,
              "minimum_available_memory_mib": 2048.0,
              "maximum_swap_io_pages_per_stage": 0}
    require(host_guard_violations({
        "load_one_minute_per_cpu": 0.75,
        "available_memory_mib": 1024.0,
        "swap_io_pages_since_previous_check": 2,
    }, policy) == ["host_load_above_maximum",
                   "host_available_memory_below_minimum",
                   "host_swap_activity_above_maximum"],
            "host guard did not identify all pressure signals")
    require(not host_guard_violations({
        "load_one_minute_per_cpu": 0.25,
        "available_memory_mib": 4096.0,
        "swap_io_pages_since_previous_check": 0,
    }, policy), "host guard rejected a healthy synthetic sample")
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

        manifest_path = work / "campaign" / "campaign-report.json"
        interrupted = json.loads(manifest_path.read_text(encoding="utf-8"))
        interrupted["stages"].pop()
        interrupted["latest_state"] = json.loads((
            work / "campaign" / "stage-000003-000004" / "stage-report.json"
        ).read_text(encoding="utf-8"))["final_state"]
        interrupted["completed"] = False
        interrupted["target_reached"] = False
        interrupted["stop_reason"] = "stage_completed"
        interrupted.pop("total_wall_seconds", None)
        manifest_path.write_text(json.dumps(interrupted), encoding="utf-8")
        resumed = subprocess.run([
            *common, "--resume-existing", "--acknowledge-cost", ACK],
            text=True, capture_output=True)
        require(resumed.returncode == 0, resumed.stderr)
        resumed_report = json.loads(resumed.stdout)
        require(resumed_report["target_reached"] and
                resumed_report["recovered_unrecorded_stages"] == 1 and
                resumed_report["stages"][-1].get(
                    "recovered_after_coordinator_interruption") is True and
                resumed_report["total_stage_wall_seconds"] > 0.0,
                "coordinator did not reconcile a completed unrecorded stage")

        limited = list(common)
        limited[4] = str(work / "limited-campaign")
        limited[limited.index("--max-wall-seconds") + 1] = "5"
        stopped = subprocess.run([*limited, "--acknowledge-cost", ACK],
                                 text=True, capture_output=True)
        require(stopped.returncode == 0, stopped.stderr)
        stopped_report = json.loads(stopped.stdout)
        require(stopped_report["completed"] and not stopped_report["target_reached"]
                and not stopped_report["stages"] and
                stopped_report["stop_reason"] == "overall_wall_time_exhausted",
                "coordinator did not stop cleanly at its wall-time reserve")

        if Path("/proc/meminfo").is_file():
            guarded = list(common)
            guarded[4] = str(work / "guarded-campaign")
            guarded.extend(["--min-available-memory-mib", "1000000000000"])
            stopped = subprocess.run([*guarded, "--acknowledge-cost", ACK],
                                     text=True, capture_output=True)
            require(stopped.returncode == 0, stopped.stderr)
            guarded_report = json.loads(stopped.stdout)
            require(guarded_report["completed"] and
                    not guarded_report["target_reached"] and
                    not guarded_report["stages"] and
                    guarded_report["stop_reason"] ==
                        "host_available_memory_below_minimum" and
                    guarded_report["host_health_checks"][0]["phase"] ==
                        "before_stage",
                    "coordinator did not stop cleanly at its memory guard")
    print("eduPIC adaptive advancement regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
