#!/usr/bin/env python3
"""Bounded tests for the locked Turner adaptive segment preparer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile

import prepare_turner_adaptive_segment as subject


BASE = """config_version = 1
units = si
dimension = 1
nx = 129
dt = 1e-10
steps = 800
output_interval = 400
output_dir = old
boundary = dirichlet
mode = transient
runtime_backend = serial
runtime_threads = 1
checkpoint_output = true
checkpoint_interval = 12800
spatial_average = true
spatial_average_interval = 1
spatial_average_start_step = 1
spatial_average_end_step = 800
spatial_average_rf_frequency = 13560000
spatial_average_rf_cycles = 2
restart_path = old.apc
[collisions.electron]
opportunity_sampling = single_bernoulli
[collisions.ion]
opportunity_sampling = single_bernoulli
[species.electrons]
charge = -1
mass = 1
weight = 1
particles = 8
thermal_velocity = 0
"""


def main() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        base = root / "base.cfg"
        checkpoint = root / "checkpoint_800.apc"
        solver = root / "aurorapic_cli"
        rule_path = root / "rule.json"
        lock_path = root / "lock.json"
        base.write_text(BASE, encoding="utf-8")
        checkpoint.write_text(
            "AuroraPIC-checkpoint-v24\ndimension 1\nstep 800\n",
            encoding="utf-8")
        solver.write_bytes(b"synthetic solver")
        rule = {
            "case_id": "turner-helium-ccp-2013-case-1",
            "locked_initial_states": [{
                "seed": 7,
                "checkpoint_step": 800,
                "checkpoint_sha256": subject.sha256(checkpoint),
                "base_config_sha256": subject.sha256(base),
            }],
            "rf_contract": {
                "frequency_hz": 13560000.0,
                "steps_per_cycle": 400,
                "cycles_per_block": 32,
            },
            "adaptive_equilibration": {
                "execution_segment_blocks": 8,
                "maximum_blocks_per_seed": 64,
                "minimum_nominal_blocks": 16,
                "minimum_ar1_effective_blocks_per_observable": 8.0,
                "maximum_absolute_projected_fractional_drift_per_observable": 0.01,
                "maximum_absolute_split_half_fractional_change_per_observable": 0.01,
                "maximum_relative_standard_error_per_observable": 0.01,
            },
        }
        rule_path.write_text(json.dumps(rule), encoding="utf-8")
        lock = {
            "status": "preregistered_not_launched",
            "rule": {"sha256": subject.sha256(rule_path)},
            "execution_order": {"seeds": [7]},
            "command_identity": {
                "solver_binary_sha256": subject.sha256(solver)},
        }
        lock_path.write_text(json.dumps(lock), encoding="utf-8")
        args = argparse.Namespace(
            rule=rule_path, execution_lock=lock_path, seed=7,
            base_config=base, checkpoint=checkpoint, solver=solver,
            prior_progress=None,
            output_dir=root / "run", output_config=root / "segment.cfg",
            report=root / "preflight.json",
            acknowledge_cost=subject.ACKNOWLEDGEMENT)
        report = subject.prepare(args)
        assert report["segment"]["target_step"] == 103200
        assert report["launched"] is False
        deck = args.output_config.read_text(encoding="utf-8")
        assert "mode = steady_state" in deck
        assert "max_steps = 103200" in deck
        assert "periodic_convergence_reset_on_restart = true" in deck
        assert "spatial_average_start_step = 103201" in deck
        assert "periodic_convergence_minimum_effective_blocks = 8.0" in deck

        continuation_checkpoint = root / "checkpoint_103200.apc"
        continuation_checkpoint.write_text(
            "AuroraPIC-checkpoint-v25\ndimension 1\nstep 103200\n",
            encoding="utf-8")
        progress_path = root / "progress.json"
        progress = {
            "case_id": rule["case_id"],
            "rule": {"sha256": subject.sha256(rule_path)},
            "execution_lock": {"sha256": subject.sha256(lock_path)},
            "admitted_segments": [{
                "seed": 7,
                "segment": 1,
                "blocks": 8,
                "end_step": 103200,
                "generated_deck_sha256": subject.sha256(args.output_config),
                "output_checkpoint_sha256": subject.sha256(
                    continuation_checkpoint),
                "integrity": {"admitted": True},
                "classification": "admitted_horizon_incomplete",
            }],
        }
        progress_path.write_text(json.dumps(progress), encoding="utf-8")
        continuation = argparse.Namespace(
            rule=rule_path, execution_lock=lock_path, seed=7,
            base_config=args.output_config,
            checkpoint=continuation_checkpoint, solver=solver,
            prior_progress=progress_path,
            output_dir=root / "run_2",
            output_config=root / "segment_2.cfg",
            report=root / "preflight_2.json",
            acknowledge_cost=subject.ACKNOWLEDGEMENT)
        continuation_report = subject.prepare(continuation)
        assert continuation_report["segment"]["index"] == 2
        assert continuation_report["segment"]["prior_admitted_blocks"] == 8
        assert continuation_report["segment"]["target_step"] == 205600
        assert continuation_report[
            "periodic_convergence_epoch_imports_prior_samples"] is True
        continuation_deck = continuation.output_config.read_text(
            encoding="utf-8")
        assert "max_steps = 205600" in continuation_deck
        assert "periodic_convergence_reset_on_restart = false" in (
            continuation_deck)

        tampered_checkpoint = root / "tampered_checkpoint.apc"
        tampered_checkpoint.write_text(
            "AuroraPIC-checkpoint-v25\ndimension 1\nstep 103200\nextra 1\n",
            encoding="utf-8")
        tampered = argparse.Namespace(**vars(continuation))
        tampered.checkpoint = tampered_checkpoint
        tampered.output_dir = root / "tampered_run"
        tampered.output_config = root / "tampered.cfg"
        tampered.report = root / "tampered.json"
        rejected_tamper = False
        try:
            subject.prepare(tampered)
        except subject.PreparationError:
            rejected_tamper = True
        assert rejected_tamper, "continuation must reject checkpoint tampering"

        rejected = False
        try:
            subject.prepare(args)
        except subject.PreparationError:
            rejected = True
        assert rejected, "preparer must refuse overwriting its products"

    print("Turner adaptive segment preparer tests passed")


if __name__ == "__main__":
    main()
