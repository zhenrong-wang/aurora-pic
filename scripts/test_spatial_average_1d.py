#!/usr/bin/env python3
"""Conservative CLI regression for restart-safe 1D spatial averaging."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run(binary: Path, config: Path) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    return subprocess.run(
        [str(binary), str(config)],
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def deck(output: Path, steps: int, restart: Path | None = None,
         interval: int = 1, start_step: int = 1,
         end_step: int = 4, reset: bool = False,
         phase_bins: int = 0) -> str:
    restart_line = "" if restart is None else f"restart_path = {restart}\n"
    reset_line = (
        "spatial_average_reset_on_restart = true\n" if reset else ""
    )
    phase_line = (
        f"spatial_average_phase_bins = {phase_bins}\n"
        if phase_bins else ""
    )
    return (
        "config_version = 1\n"
        "units = normalized\n"
        "dimension = 1\n"
        "velocity_dimensions = 1\n"
        "nx = 5\n"
        "length = 1\n"
        "dt = 0.01\n"
        f"steps = {steps}\n"
        "output_interval = 2\n"
        f"output_dir = {output}\n"
        f"{restart_line}"
        "boundary = periodic\n"
        "mode = transient\n"
        "runtime_backend = serial\n"
        "runtime_threads = 1\n"
        "checkpoint_output = true\n"
        "checkpoint_interval = 2\n"
        "spatial_average = true\n"
        f"{reset_line}"
        f"spatial_average_interval = {interval}\n"
        f"spatial_average_start_step = {start_step}\n"
        f"spatial_average_end_step = {end_step}\n"
        f"{phase_line}"
        "[species.electrons]\n"
        "charge = -1\n"
        "mass = 1\n"
        "weight = 0.125\n"
        "particles = 8\n"
        "thermal_velocity = 0.2\n"
        "loading = quiet_start\n"
        "[species.ions]\n"
        "charge = 1\n"
        "mass = 40\n"
        "weight = 0.125\n"
        "particles = 8\n"
        "thermal_velocity = 0.05\n"
        "loading = quiet_start\n"
    )


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_spatial_average_1d.py AURORAPIC_CLI", file=sys.stderr)
        return 2
    binary = Path(sys.argv[1]).resolve()
    require(binary.is_file(), f"missing AuroraPIC CLI: {binary}")
    project_tmp = ROOT / "tmp"
    project_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_spatial_average_", dir=project_tmp
    ) as temporary:
        work = Path(temporary)
        direct_output = work / "direct-output"
        direct_config = work / "direct.cfg"
        direct_config.write_text(
            deck(direct_output, 4), encoding="utf-8"
        )
        direct = run(binary, direct_config)
        require(direct.returncode == 0,
                f"direct spatial-average run failed: {direct.stderr}")

        split_output = work / "split-output"
        split_config = work / "split.cfg"
        split_config.write_text(
            deck(split_output, 2), encoding="utf-8"
        )
        split = run(binary, split_config)
        require(split.returncode == 0,
                f"split spatial-average run failed: {split.stderr}")
        checkpoint = split_output / "checkpoint_2.apc"
        require(checkpoint.is_file(),
                "split run did not write its averaging checkpoint")

        resumed_output = work / "resumed-output"
        resumed_config = work / "resumed.cfg"
        resumed_config.write_text(
            deck(resumed_output, 4, checkpoint), encoding="utf-8"
        )
        resumed = run(binary, resumed_config)
        require(resumed.returncode == 0,
                f"resumed spatial-average run failed: {resumed.stderr}")

        direct_profile = direct_output / "spatial_average.csv"
        resumed_profile = resumed_output / "spatial_average.csv"
        require(
            direct_profile.read_bytes() == resumed_profile.read_bytes(),
            "checkpoint continuation changed the spatial average",
        )
        for filename in (
            "spatial_kinetic_energy.csv",
            "spatial_field_average.csv",
        ):
            require(
                (direct_output / filename).read_bytes()
                == (resumed_output / filename).read_bytes(),
                f"checkpoint continuation changed {filename}",
            )
        metadata = json.loads(
            (resumed_output / "spatial_average_metadata.json").read_text(
                encoding="utf-8"
            )
        )
        require(
            metadata["complete"] is True
            and metadata["samples"] == 4
            and metadata["moment_samples"] == 4
            and metadata["moments_complete"] is True
            and metadata["expected_samples"] == 4
            and metadata["final_step"] == 4
            and metadata["spatial_average_version"] == 7
            and metadata["sampling_order"] == "post_collision",
            "completed spatial-average metadata is incorrect",
        )
        totals: dict[str, float] = {}
        with resumed_profile.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                totals[row["species"]] = totals.get(row["species"], 0.0) + (
                    float(row["number_density_mean_normalized"]) * 0.2
                )
        require(
            set(totals) == {"electrons", "ions"}
            and all(abs(value - 1.0) < 1e-14 for value in totals.values()),
            "spatial-average number density is not conservative",
        )
        with (resumed_output / "spatial_kinetic_energy.csv").open(
            newline="", encoding="utf-8"
        ) as stream:
            energy_rows = list(csv.DictReader(stream))
        require(
            energy_rows
            and max(float(row["mean_kinetic_energy_normalized"])
                    for row in energy_rows) > 0.0,
            "spatial kinetic-energy profile lost nonzero particle energy",
        )
        with (resumed_output / "spatial_field_average.csv").open(
            newline="", encoding="utf-8"
        ) as stream:
            field_rows = list(csv.DictReader(stream))
        require(
            len(field_rows) == 5
            and all(
                float(row["electric_field_rms_normalized"]) + 1e-14
                >= abs(float(row["electric_field_mean_normalized"]))
                for row in field_rows
            ),
            "spatial field RMS is inconsistent with its mean",
        )

        legacy_checkpoint = work / "checkpoint_2_v7.apc"
        legacy_lines = checkpoint.read_text(encoding="utf-8").splitlines()
        legacy_lines[0] = "AuroraPIC-checkpoint-v7"
        for index, line in enumerate(legacy_lines):
            if line.startswith("spatial_average "):
                fields = line.split()
                del fields[2]
                legacy_lines[index] = " ".join(fields)
        legacy_checkpoint.write_text(
            "\n".join(
                line for line in legacy_lines
                if not line.startswith((
                    "species_timestep_multipliers",
                    "subcycle_charge_deposition", "subcycle_charge_cache",
                    "collision_energy_totals",
                    "spatial_moments", "spatial_energy", "spatial_fields",
                    "spatial_phase", "phase_bin", "phase_species",
                    "phase_fields", "spatial_collision", "phase_eedf",
                    "phase_surface_flux", "wall_impact",
                ))
            ) + "\n",
            encoding="utf-8",
        )
        legacy_output = work / "legacy-output"
        legacy_config = work / "legacy.cfg"
        legacy_config.write_text(
            deck(legacy_output, 4, legacy_checkpoint), encoding="utf-8"
        )
        legacy = run(binary, legacy_config)
        require(legacy.returncode == 0,
                f"legacy v7 continuation failed: {legacy.stderr}")
        legacy_metadata = json.loads(
            (legacy_output / "spatial_average_metadata.json").read_text(
                encoding="utf-8"
            )
        )
        require(
            legacy_metadata["samples"] == 4
            and legacy_metadata["moment_samples"] == 2
            and legacy_metadata["moments_complete"] is False,
            "legacy v7 continuation did not expose partial moment history",
        )
        require(
            len((legacy_output / "spatial_kinetic_energy.csv")
                .read_text(encoding="utf-8").splitlines()) == 1
            and len((legacy_output / "spatial_field_average.csv")
                    .read_text(encoding="utf-8").splitlines()) == 1,
            "partial legacy moment history produced misleading profiles",
        )

        mismatched_config = work / "mismatched.cfg"
        mismatched_config.write_text(
            deck(work / "mismatched-output", 4, checkpoint, interval=2),
            encoding="utf-8",
        )
        mismatched = run(binary, mismatched_config)
        require(
            mismatched.returncode != 0
            and "spatial-average contract" in mismatched.stderr,
            "checkpoint accepted a changed spatial-average contract",
        )

        late_direct_output = work / "late-direct-output"
        late_direct_config = work / "late-direct.cfg"
        late_direct_config.write_text(
            deck(
                late_direct_output, 4,
                start_step=3, end_step=4,
            ),
            encoding="utf-8",
        )
        late_direct = run(binary, late_direct_config)
        require(
            late_direct.returncode == 0,
            f"late direct spatial-average run failed: "
            f"{late_direct.stderr}",
        )
        reset_output = work / "reset-output"
        reset_config = work / "reset.cfg"
        reset_config.write_text(
            deck(
                reset_output, 4, checkpoint,
                start_step=3, end_step=4, reset=True,
            ),
            encoding="utf-8",
        )
        reset_run = run(binary, reset_config)
        require(
            reset_run.returncode == 0,
            f"restart-reset spatial-average run failed: "
            f"{reset_run.stderr}",
        )
        require(
            (late_direct_output / "spatial_average.csv").read_bytes()
            == (reset_output / "spatial_average.csv").read_bytes(),
            "restart-reset spatial average differs from direct window",
        )
        for filename in (
            "spatial_kinetic_energy.csv",
            "spatial_field_average.csv",
        ):
            require(
                (late_direct_output / filename).read_bytes()
                == (reset_output / filename).read_bytes(),
                f"restart-reset {filename} differs from direct window",
            )
        reset_metadata = json.loads(
            (reset_output / "spatial_average_metadata.json").read_text(
                encoding="utf-8"
            )
        )
        require(
            reset_metadata["reset_on_restart"] is True
            and reset_metadata["samples"] == 2
            and reset_metadata["moment_samples"] == 2
            and reset_metadata["moments_complete"] is True
            and reset_metadata["expected_samples"] == 2
            and reset_metadata["complete"] is True,
            "restart-reset spatial-average metadata is incorrect",
        )

        missed_config = work / "missed-reset.cfg"
        missed_config.write_text(
            deck(
                work / "missed-reset-output", 4, checkpoint,
                start_step=2, end_step=4, reset=True,
            ),
            encoding="utf-8",
        )
        missed = run(binary, missed_config)
        require(
            missed.returncode != 0
            and "must start after the checkpoint step"
                in missed.stderr,
            "restart reset accepted a window with a missed sample",
        )

        rf_config = work / "rf.cfg"
        rf_config.write_text(
            deck(work / "rf-output", 4, phase_bins=2)
            .replace("dt = 0.01", "dt = 0.25")
            .replace(
                "spatial_average_end_step = 4\n",
                "spatial_average_end_step = 4\n"
                "spatial_average_rf_frequency = 1\n"
                "spatial_average_rf_cycles = 1\n",
            ),
            encoding="utf-8",
        )
        rf = run(binary, rf_config)
        require(rf.returncode == 0,
                f"whole-cycle RF averaging was rejected: {rf.stderr}")
        rf_output = work / "rf-output"
        rf_metadata = json.loads(
            (rf_output / "spatial_average_metadata.json").read_text(
                encoding="utf-8"))
        require(
            rf_metadata["phase_bins"] == 2
            and rf_metadata["phase_bin_samples"] == [2, 2],
            "RF phase-bin sample allocation is wrong",
        )
        with (rf_output / "spatial_phase_moments.csv").open(
            newline="", encoding="utf-8") as stream:
            phase_rows = list(csv.DictReader(stream))
        require(
            phase_rows
            and all(float(row["drift_separated_temperature_normalized"])
                    >= 0.0 for row in phase_rows)
            and all(
                float(row["drift_separated_temperature_normalized"])
                <= 2.0 * float(row["mean_kinetic_energy_normalized"])
                   + 1e-14
                for row in phase_rows
            ),
            "drift-separated phase temperature is invalid",
        )

        rf_split_output = work / "rf-split-output"
        rf_split_config = work / "rf-split.cfg"
        rf_split_config.write_text(
            deck(rf_split_output, 2, phase_bins=2)
            .replace("dt = 0.01", "dt = 0.25")
            .replace(
                "spatial_average_end_step = 4\n",
                "spatial_average_end_step = 4\n"
                "spatial_average_rf_frequency = 1\n"
                "spatial_average_rf_cycles = 1\n",
            ), encoding="utf-8")
        rf_split = run(binary, rf_split_config)
        require(rf_split.returncode == 0,
                f"phase-bin split run failed: {rf_split.stderr}")
        rf_resumed_output = work / "rf-resumed-output"
        rf_resumed_config = work / "rf-resumed.cfg"
        rf_resumed_config.write_text(
            deck(rf_resumed_output, 4,
                 rf_split_output / "checkpoint_2.apc", phase_bins=2)
            .replace("dt = 0.01", "dt = 0.25")
            .replace(
                "spatial_average_end_step = 4\n",
                "spatial_average_end_step = 4\n"
                "spatial_average_rf_frequency = 1\n"
                "spatial_average_rf_cycles = 1\n",
            ), encoding="utf-8")
        rf_resumed = run(binary, rf_resumed_config)
        require(rf_resumed.returncode == 0,
                f"phase-bin resumed run failed: {rf_resumed.stderr}")
        for filename in ("spatial_phase_moments.csv",
                         "spatial_phase_fields.csv"):
            require(
                (rf_output / filename).read_bytes()
                == (rf_resumed_output / filename).read_bytes(),
                f"checkpoint continuation changed {filename}",
            )

        invalid_rf_config = work / "invalid-rf.cfg"
        invalid_rf_config.write_text(
            rf_config.read_text(encoding="utf-8").replace(
                "spatial_average_end_step = 4",
                "spatial_average_end_step = 3",
            ).replace("steps = 4", "steps = 3"),
            encoding="utf-8",
        )
        invalid_rf = run(binary, invalid_rf_config)
        require(
            invalid_rf.returncode != 0
            and "whole RF cycles" in invalid_rf.stderr,
            "non-whole-cycle RF averaging window was accepted",
        )

        invalid_phase_config = work / "invalid-phase.cfg"
        invalid_phase_config.write_text(
            rf_config.read_text(encoding="utf-8").replace(
                "spatial_average_phase_bins = 2",
                "spatial_average_phase_bins = 3"),
            encoding="utf-8")
        invalid_phase = run(binary, invalid_phase_config)
        require(
            invalid_phase.returncode != 0
            and "phase_bins must divide" in invalid_phase.stderr,
            "non-divisor RF phase-bin count was accepted",
        )

    print("restart-safe 1D spatial averaging passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
